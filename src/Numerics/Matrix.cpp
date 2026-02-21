/******************************************************************************
 * @file Matrix.cpp
 * @brief Matrix assembly and linear system construction for equations
 *****************************************************************************/

#include <iostream>
#include <algorithm>
#include <cmath>

#include "Matrix.hpp"
#include "CellData.hpp"
#include "LinearInterpolation.hpp"


// ************************* Constructor & Destructor *************************

Matrix::Matrix
(
    const std::vector<Face>& faces,
    const std::vector<Cell>& cells,
    const BoundaryConditions& boundaryConds
) : allFaces_(faces),
    allCells_(cells),
    bcManager_(boundaryConds),
    numInternalFaces_(0),
    numBoundaryFaces_(0)
{
    for (const auto& face : allFaces_)
    {
        if (face.isBoundary()) ++numBoundaryFaces_;
        else ++numInternalFaces_;
    }
}


// ****************************** Public Methods ******************************

void Matrix::buildMatrix(const TransportEquation& equation)
{
    clear();
    reserveTripletList();

    size_t numCells = allCells_.size();

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        vectorB_(cellIdx) += equation.source[cellIdx];
    }

    for (size_t faceIdx = 0;
         faceIdx < allFaces_.size();
         ++faceIdx)
    {
        const Face& face = allFaces_[faceIdx];

        if (face.isBoundary())
        {
            assembleBoundaryFace(face, equation);
        }
        else
        {
            assembleInternalFace(face, equation);
        }
    }

    matrixA_.setFromTriplets
    (
        tripletList_.begin(),
        tripletList_.end()
    );
}

void Matrix::relax(Scalar alpha, const ScalarField& phiPrev)
{
    if (alpha <= 0.0)
    {
        throw   std::runtime_error
                (
                    "Matrix::relax: alpha must be positive"
                );
    }

    const int numCells = static_cast<int>(matrixA_.rows());

    if (phiPrev.size() != static_cast<size_t>(numCells))
    {
        throw   std::runtime_error
                (
                    "Matrix::relax: phiPrev size mismatch"
                    " with matrix size"
                );
    }

    // Cache original diagonal
    std::vector<Scalar> originalDiag(static_cast<size_t>(numCells));

    for (int idx = 0; idx < numCells; ++idx)
    {
        originalDiag[static_cast<size_t>(idx)] = matrixA_.coeff(idx, idx);
    }

    // Scale diagonal: a_P <- a_P / alpha
    for (int idx = 0; idx < numCells; ++idx)
    {
        matrixA_.coeffRef(idx, idx) = 
            originalDiag[static_cast<size_t>(idx)] / alpha;
    }

    // Update RHS: b <- b + ((1-alpha)/alpha) * a_P * phiPrev
    const Scalar factor = (S(1.0) - alpha) / alpha;

    for (int idx = 0; idx < numCells; ++idx)
    {
        vectorB_(idx) +=
            factor * originalDiag[static_cast<size_t>(idx)]
          * phiPrev[static_cast<size_t>(idx)];
    }
}


// ***************************** Private Methods *****************************

void Matrix::clear()
{
    tripletList_.clear();

    size_t numCells = allCells_.size();

    matrixA_.resize(numCells, numCells);
    matrixA_.setZero();

    vectorB_.resize(numCells);
    vectorB_.setZero();
}

void Matrix::reserveTripletList()
{
    size_t reserveSize = 4 * numInternalFaces_ + numBoundaryFaces_;

    tripletList_.reserve(reserveSize);
}

Scalar Matrix::extractBoundaryScalar
(
    const BoundaryData& bc,
    const std::string& fieldName
)
{
    if (bc.valueType() == BCValueType::VECTOR)
    {
        if (fieldName == "Ux") return bc.vectorValue().x();
        if (fieldName == "Uy") return bc.vectorValue().y();
        if (fieldName == "Uz") return bc.vectorValue().z();
    }

    return bc.scalarValue();
}


// ***************** Transport Equation Assembly Helpers ******************

void Matrix::assembleInternalFace
(
    const Face& face,
    const TransportEquation& equation
)
{
    const size_t ownerIdx = face.ownerCell();
    const size_t neighborIdx = face.neighborCell().value();

    Vector Sf = face.normal() * face.projectedArea();

    Vector dPN =
        allCells_[neighborIdx].centroid()
      - allCells_[ownerIdx].centroid();

    Scalar dPNMag = dPN.magnitude();
    Vector ePN = dPN / (dPNMag + vSmallValue);

    // Orthogonal component (over-relaxed)
    Vector Ef = (dot(Sf, Sf) / dot(Sf, ePN)) * ePN;

    // Diffusion coefficient at face
    Scalar Gammaf;

    if (equation.Gamma)
    {
        // Cell-based: harmonic interpolation
        const auto& G = equation.Gamma->get();
        Scalar d_Pf = face.dPfMag();
        Scalar d_Nf = face.dNfMag().value();

        Gammaf =
            dPNMag
          / ((d_Pf / (G[ownerIdx] + vSmallValue))
           + (d_Nf / (G[neighborIdx] + vSmallValue)));
    }
    else
    {
        // Face-based: use directly
        Gammaf = equation.GammaFace->get()[face.idx()];
    }

    Scalar a_diff =
        Gammaf * Ef.magnitude() / (dPNMag + vSmallValue);

    // Convection coefficients
    Scalar aPConv = S(0.0);
    Scalar aNConv = S(0.0);

    if (equation.flowRate)
    {
        Scalar flowRate = equation.flowRate->get()[face.idx()];

        auto [aP, aN] = 
            equation.convScheme->get().getFluxCoefficients(flowRate);

        aPConv = aP;
        aNConv = aN;
    }

    // Matrix coefficients for owner and neighbor cells
    tripletList_.emplace_back(ownerIdx, ownerIdx, a_diff + aPConv);

    tripletList_.emplace_back(ownerIdx, neighborIdx, -a_diff + aNConv);

    tripletList_.emplace_back(neighborIdx, neighborIdx, a_diff - aNConv);

    tripletList_.emplace_back(neighborIdx, ownerIdx, -a_diff - aPConv);

    // Non-orthogonal correction (explicit)
    Vector Tf = Sf - Ef;

    const Vector& gradPhiP = equation.gradPhi[ownerIdx];
    const Vector& gradPhiN = equation.gradPhi[neighborIdx];

    Vector gradPhif =
        equation.gradScheme.faceGradient
        (
            equation.fieldName,
            equation.phi,
            gradPhiP,
            gradPhiN,
            face.idx()
        );

    Scalar nonOrthogonalFlux = Gammaf * dot(gradPhif, Tf);

    vectorB_(ownerIdx)    += nonOrthogonalFlux;
    vectorB_(neighborIdx) -= nonOrthogonalFlux;

    // Deferred correction (explicit, convection only)
    if (equation.convScheme)
    {
        Scalar deferredCorrection =
            equation.convScheme->get().calculateCorrection
            (
                face,
                equation.phi,
                gradPhiP,
                gradPhiN,
                equation.flowRate->get()[face.idx()]
            );

        vectorB_(ownerIdx)    -= deferredCorrection;
        vectorB_(neighborIdx) += deferredCorrection;
    }
}

void Matrix::assembleBoundaryFace
(
    const Face& face,
    const TransportEquation& equation
)
{
    const size_t ownerIdx = face.ownerCell();

    const BoundaryData* bc =
        bcManager_.fieldBC(face.patch()->patchName(), equation.fieldName);

    if (!bc)
    {
        // No BC found: skip face (zero-gradient behaviour)
        return;
    }

    Vector Sf = face.normal() * face.projectedArea();
    const Vector ePf = face.e_Pf();
    const Scalar dPfMag = face.dPfMag();

    Vector Ef =(dot(Sf, Sf) / dot(Sf, ePf)) * ePf;

    // Diffusion coefficient at boundary face
    Scalar Gammaf;

    if (equation.Gamma)
    {
        Gammaf = equation.Gamma->get()[ownerIdx];
    }
    else
    {
        Gammaf = equation.GammaFace->get()[face.idx()];
    }

    Scalar aDiff = Gammaf * Ef.magnitude() / (dPfMag + vSmallValue);

    // Convert OptionalRef to raw pointer for GradientScheme compatibility
    const FaceData<Scalar>* bfv =
        equation.boundaryFaceValues ? 
        &equation.boundaryFaceValues->get() 
      : nullptr;

    if
    (
        bc->type() == BCType::FIXED_VALUE
     || bc->type() == BCType::NO_SLIP
     || bc->type() == BCType::OMEGA_WALL_FUNCTION
    )
    {
        // Dirichlet BC: phi_b is prescribed
        Scalar phiB = S(0.0);

        if (bc->type() == BCType::NO_SLIP)
        {
            phiB = S(0.0);
        }
        else if (bc->type() == BCType::OMEGA_WALL_FUNCTION)
        {
            bool hasDynamicValue =
                bfv
             && face.idx() < bfv->size()
             && std::isfinite((*bfv)[face.idx()]);

            if (hasDynamicValue)
            {
                phiB = (*bfv)[face.idx()];
            }
            else
            {
                phiB = equation.phi[ownerIdx];
            }
        }
        else
        {
            phiB = extractBoundaryScalar(*bc, equation.fieldName);
        }

        // Diffusion contribution
        tripletList_.emplace_back(ownerIdx, ownerIdx, aDiff);
        vectorB_(ownerIdx) += aDiff * phiB;

        // Convection contribution
        if (equation.flowRate)
        {
            Scalar a_conv = equation.flowRate->get()[face.idx()];
            vectorB_(ownerIdx) -= a_conv * phiB;
        }

        // Non-orthogonal correction
        Vector Tf = Sf - Ef;

        Vector gradPhif =
            equation.gradScheme.faceGradient
            (
                equation.fieldName,
                equation.phi,
                equation.gradPhi[ownerIdx],
                Vector(),
                face.idx(),
                bfv
            );

        Scalar nonOrthogonalFlux =
            Gammaf * dot(gradPhif, Tf);

        vectorB_(ownerIdx) += nonOrthogonalFlux;
    }
    else if
    (
        bc->type() == BCType::ZERO_GRADIENT
     || bc->type() == BCType::K_WALL_FUNCTION
     || bc->type() == BCType::NUT_WALL_FUNCTION
    )
    {
        // Zero normal gradient: only convection
        if (equation.flowRate)
        {
            Scalar a_conv = equation.flowRate->get()[face.idx()];

            tripletList_.emplace_back(ownerIdx, ownerIdx, a_conv);

            // Tangential gradient correction
            const Vector& gradPhiP = equation.gradPhi[ownerIdx];

            Vector eb = face.normal();
            Scalar normal_component = dot(gradPhiP, eb);
            Vector gradPhib = gradPhiP - normal_component * eb;

            Vector dCb = face.centroid() - allCells_[ownerIdx].centroid();

            Scalar correction = a_conv * dot(gradPhib, dCb);
            vectorB_(ownerIdx) -= correction;
        }
        // No convection + zero gradient = no contribution
    }
    else
    {
        // Unhandled BC type: default to zero gradient
        std::cerr
            << "Warning: Undefined boundary "
            << "condition type for field "
            << equation.fieldName << " on patch "
            << face.patch()->patchName()
            << ". Applying zero gradient."
            << std::endl;

        if (equation.flowRate)
        {
            Scalar a_conv = equation.flowRate->get()[face.idx()];
            tripletList_.emplace_back(ownerIdx, ownerIdx, a_conv);
        }
    }
}
