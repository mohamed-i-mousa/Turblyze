/******************************************************************************
 * @file Matrix.cpp
 * @brief Matrix assembly and linear system construction for equations
 *****************************************************************************/

#include "Matrix.hpp"
#include <iostream>
#include <algorithm>
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
        if (face.isBoundary())
            ++numBoundaryFaces_;
        else
            ++numInternalFaces_;
    }
}


// ****************************** Public Methods ******************************

void Matrix::buildMatrix
(
    const ScalarField& phi,
    const ScalarField& phiSource,
    const FaceFluxField& flowRateFace,
    const ScalarField& Gamma,
    const ConvectionScheme& convScheme,
    const GradientScheme& gradScheme,
    const VectorField& gradPhi,
    const std::string& fieldName
)
{
    clear();
    reserveTripletList();

    size_t numCells = allCells_.size();

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        vectorB_(cellIdx) += phiSource[cellIdx];
    }

    for (size_t faceIdx = 0; faceIdx < allFaces_.size(); ++faceIdx)
    {
        const Face& face = allFaces_[faceIdx];

        if (face.isBoundary())
        {
            assembleBoundaryFace
            (
                face,
                phi,
                flowRateFace,
                Gamma,
                gradScheme,
                gradPhi,
                fieldName
            );
        }
        else
        {
            assembleInternalFace
            (
                face,
                phi,
                flowRateFace,
                Gamma,
                convScheme,
                gradScheme,
                gradPhi,
                fieldName
            );
        }
    }

    matrixA_.setFromTriplets
    (
        tripletList_.begin(),
        tripletList_.end()
    );
}

void Matrix::buildPressureCorrectionMatrix
(
    const FaceFluxField& RhieChowFlowRate,
    const FaceFluxField& DUf,
    const ScalarField& pCorr,
    const GradientScheme& gradScheme,
    const VectorField& gradPCorr
)
{
    clear();
    reserveTripletList();

    // Build mass imbalance RHS
    for (size_t cellIdx = 0; cellIdx < allCells_.size(); ++cellIdx)
    {
        const Cell& cell = allCells_[cellIdx];

        Scalar massImbalance = 0.0;

        for (size_t i = 0; i < cell.faceIndices().size(); ++i)
        {
            const size_t faceIdx = cell.faceIndices()[i];
            const int sign = cell.faceSigns()[i];
            massImbalance += sign * RhieChowFlowRate[faceIdx];
        }

        vectorB_[cellIdx] = -massImbalance;
    }

    // Build matrix coefficients
    for (size_t faceIdx = 0; faceIdx < allFaces_.size(); ++faceIdx)
    {
        const Face& face = allFaces_[faceIdx];

        if (face.isBoundary())
        {
            assemblePCorrBoundaryFace(face, DUf, gradPCorr);
        }
        else
        {
            assemblePCorrInternalFace
            (
                face, 
                DUf, 
                pCorr, 
                gradScheme, 
                gradPCorr
            );
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
        originalDiag[static_cast<size_t>(idx)] =
            matrixA_.coeff(idx, idx);
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
    size_t reserveSize =
        4 * numInternalFaces_ + numBoundaryFaces_;

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
        if (fieldName == "U_x") return bc.vectorValue().x();
        if (fieldName == "U_y") return bc.vectorValue().y();
        if (fieldName == "U_z") return bc.vectorValue().z();
    }

    return bc.scalarValue();
}


// ***************** Transport Equation Assembly Helpers ******************

void Matrix::assembleInternalFace
(
    const Face& face,
    const ScalarField& phi,
    const FaceFluxField& flowRateFace,
    const ScalarField& Gamma,
    const ConvectionScheme& convScheme,
    const GradientScheme& gradScheme,
    const VectorField& gradPhi,
    const std::string& fieldName
)
{
    const size_t ownerIdx = face.ownerCell();
    const size_t neighborIdx = face.neighborCell().value();

    Vector S_f = face.normal() * face.projectedArea();

    Vector d_PN =
        allCells_[neighborIdx].centroid()
      - allCells_[ownerIdx].centroid();

    Scalar d_PN_magnitude = d_PN.magnitude();
    Vector e_PN = d_PN / (d_PN_magnitude + vSmallValue);

    // Orthogonal component (over-relaxed)
    Vector E_f = (dot(S_f, S_f) / dot(S_f, e_PN)) * e_PN;

    // Harmonic interpolation for diffusion coefficient
    Scalar d_Pf = face.dPfMag();
    Scalar d_Nf = face.dNfMag().value();

    Scalar Gamma_f =
        d_PN_magnitude
      / ((d_Pf / (Gamma[ownerIdx] + vSmallValue))
       + (d_Nf / (Gamma[neighborIdx] + vSmallValue)));

    Scalar a_diff =
        Gamma_f * E_f.magnitude()
      / (d_PN_magnitude + vSmallValue);

    // Convection coefficients (upwind)
    Scalar flowRate = flowRateFace[face.idx()];

    auto [a_P_conv, a_N_conv] =
        convScheme.getFluxCoefficients(flowRate);

    // Matrix coefficients for owner and neighbor cells
    tripletList_.emplace_back(ownerIdx, ownerIdx, a_diff + a_P_conv);

    tripletList_.emplace_back(ownerIdx, neighborIdx, -a_diff + a_N_conv);

    tripletList_.emplace_back(neighborIdx, neighborIdx, a_diff - a_N_conv);

    tripletList_.emplace_back(neighborIdx, ownerIdx, -a_diff - a_P_conv);

    // Non-orthogonal correction (explicit)
    Vector T_f = S_f - E_f;

    const Vector& gradPhi_P = gradPhi[ownerIdx];
    const Vector& gradPhi_N = gradPhi[neighborIdx];

    Vector gradPhi_f = 
        gradScheme.faceGradient
        (
            face.idx(),
            gradPhi_P,
            gradPhi_N,
            phi,
            fieldName
        );

    Scalar nonOrthogonalFlux = Gamma_f * dot(gradPhi_f, T_f);

    vectorB_(ownerIdx)    += nonOrthogonalFlux;
    vectorB_(neighborIdx) -= nonOrthogonalFlux;

    // Deferred correction (explicit)
    Scalar deferredCorrection = convScheme.calculateCorrection
        (
            face,
            phi,
            gradPhi_P,
            gradPhi_N,
            flowRate
        );

    vectorB_(ownerIdx)    -= deferredCorrection;
    vectorB_(neighborIdx) += deferredCorrection;
}

void Matrix::assembleBoundaryFace
(
    const Face& face,
    const ScalarField& phi,
    const FaceFluxField& flowRateFace,
    const ScalarField& Gamma,
    const GradientScheme& gradScheme,
    const VectorField& gradPhi,
    const std::string& fieldName
)
{
    const size_t ownerIdx = face.ownerCell();

    // Look up patch and BC for this face
    const auto& patchMap = bcManager_.faceToPatchMap();
    auto patchIt = patchMap.find(face.idx());

    if (patchIt == patchMap.end())
    {
        std::cerr
            << "ERROR: Boundary face " << face.idx()
            << " not found in patch map! Owner cell "
            << ownerIdx << " will have zero diagonal."
            << std::endl;
        return;
    }

    const BoundaryData* bc =
        bcManager_.fieldBC
        (
            patchIt->second->patchName(),
            fieldName
        );

    if (!bc)
    {
        std::cerr
            << "ERROR: No boundary condition found"
            << " for patch "
            << patchIt->second->patchName()
            << " field " << fieldName
            << " face " << face.idx()
            << " owner cell " << ownerIdx
            << std::endl;
        return;
    }

    Vector S_f = face.normal() * face.projectedArea();
    const Vector e_Pf = face.e_Pf();
    const Scalar dPfMag = face.dPfMag();

    Vector E_f = (dot(S_f, S_f) / dot(S_f, e_Pf)) * e_Pf;

    Scalar Gamma_f = Gamma[ownerIdx];

    Scalar a_diff =
        Gamma_f * E_f.magnitude()
      / (dPfMag + vSmallValue);

    if
    (
        bc->type() == BCType::FIXED_VALUE
     || bc->type() == BCType::NO_SLIP
    )
    {
        // Dirichlet BC: phi_b is prescribed (0 for NO_SLIP)
        Scalar phi_b = (bc->type() == BCType::NO_SLIP)
            ? S(0.0)
            : extractBoundaryScalar(*bc, fieldName);

        // Diffusion: a_diff * phi_P on diagonal, a_diff * phi_b to RHS
        tripletList_.emplace_back(ownerIdx, ownerIdx, a_diff);
        vectorB_(ownerIdx) += a_diff * phi_b;

        // Convection: flux F * phi_b to RHS
        Scalar a_conv = flowRateFace[face.idx()];
        vectorB_(ownerIdx) -= a_conv * phi_b;

        // Non-orthogonal correction
        Vector T_f = S_f - E_f;

        Vector gradPhi_f = gradScheme.faceGradient
            (
                face.idx(),
                gradPhi[ownerIdx],
                Vector(),
                phi,
                fieldName
            );

        Scalar nonOrthogonalFlux =
            Gamma_f * dot(gradPhi_f, T_f);

        vectorB_(ownerIdx) += nonOrthogonalFlux;
    }
    else if (bc->type() == BCType::ZERO_GRADIENT)
    {
        // Zero normal gradient: no diffusive flux, only convection
        // a_C += F_b (mass flux to diagonal)
        Scalar a_conv = flowRateFace[face.idx()];

        tripletList_.emplace_back(ownerIdx, ownerIdx, a_conv);

        // Tangential gradient extrapolation correction
        const Vector& gradPhi_P = gradPhi[ownerIdx];

        // Remove normal component: grad_b = grad_C - (grad_C . e_b) e_b
        Vector e_b = face.normal();
        Scalar normal_component = dot(gradPhi_P, e_b);
        Vector gradPhi_b =
            gradPhi_P - normal_component * e_b;

        // d_{Cb} = cell centroid to boundary face centroid
        Vector d_Cb =
            face.centroid()
          - allCells_[ownerIdx].centroid();

        // Source correction: -F_b * (grad_b . d_{Cb})
        Scalar correction = a_conv * dot(gradPhi_b, d_Cb);
        vectorB_(ownerIdx) -= correction;
    }
    else
    {
        // Unhandled BC type: default to zero gradient
        std::cerr
            << "Warning: Undefined boundary "
            << "condition type for field "
            << fieldName << " on patch "
            << patchIt->second->patchName()
            << ". Applying zero gradient."
            << std::endl;

        Scalar a_conv = flowRateFace[face.idx()];
        tripletList_.emplace_back(ownerIdx, ownerIdx, a_conv);
    }
}


// *************** Pressure Correction Assembly Helpers ******************

void Matrix::assemblePCorrInternalFace
(
    const Face& face,
    const FaceFluxField& DUf,
    const ScalarField& pCorr,
    const GradientScheme& gradScheme,
    const VectorField& gradPCorr
)
{
    const size_t ownerIdx = face.ownerCell();
    const size_t neighborIdx = face.neighborCell().value();

    Vector S_f = face.normal() * face.projectedArea();

    Vector d_PN =
        allCells_[neighborIdx].centroid()
      - allCells_[ownerIdx].centroid();

    Scalar dPNMag = d_PN.magnitude();
    Vector e_PN = d_PN / (dPNMag + vSmallValue);

    // Orthogonal component
    Vector E_f = (dot(S_f, S_f) / dot(S_f, e_PN)) * e_PN;

    Scalar a_diff =
        DUf[face.idx()] * E_f.magnitude()
      / (dPNMag + vSmallValue);

    tripletList_.emplace_back(ownerIdx, ownerIdx, a_diff);

    tripletList_.emplace_back(neighborIdx, neighborIdx, a_diff);

    tripletList_.emplace_back (ownerIdx, neighborIdx, -a_diff);

    tripletList_.emplace_back(neighborIdx, ownerIdx, -a_diff);

    // Non-orthogonal correction
    Vector T_f = S_f - E_f;

    const Vector& gradpCorr_P = gradPCorr[ownerIdx];
    const Vector& gradpCorr_N = gradPCorr[neighborIdx];

    Vector gradPCorr_f = 
        gradScheme.faceGradient
        (
            face.idx(),
            gradpCorr_P,
            gradpCorr_N,
            pCorr,
            "pCorr"
        );

    Scalar nonOrthogonalFlux =
        DUf[face.idx()] * dot(gradPCorr_f, T_f);

    vectorB_(ownerIdx)    += nonOrthogonalFlux;
    vectorB_(neighborIdx) -= nonOrthogonalFlux;
}

void Matrix::assemblePCorrBoundaryFace
(
    const Face& face,
    const FaceFluxField& DUf,
    const VectorField& gradPCorr
)
{
    const size_t ownerIdx = face.ownerCell();

    // Guarded patch lookup
    const auto& patchMap = bcManager_.faceToPatchMap();
    auto patchIt = patchMap.find(face.idx());

    if (patchIt == patchMap.end())
    {
        std::cerr
            << "ERROR: Boundary face " << face.idx()
            << " not found in patch map for"
            << " pressure correction!"
            << std::endl;
        return;
    }

    const BoundaryData* bc =
        bcManager_.fieldBC
        (
            patchIt->second->patchName(), "p"
        );

    // Fixed-value pressure BC: pCorr = 0 at boundary
    if (bc && bc->type() == BCType::FIXED_VALUE)
    {
        Vector S_f = face.normal() * face.projectedArea();
        const Vector e_Pf = face.e_Pf();
        const Scalar dPfMag = face.dPfMag();

        Vector E_f =
            (dot(S_f, S_f) / dot(S_f, e_Pf)) * e_Pf;

        Scalar a_diff =
            DUf[face.idx()] * E_f.magnitude()
          / (dPfMag + vSmallValue);

        tripletList_.emplace_back(ownerIdx, ownerIdx, a_diff);

        // Non-orthogonal correction (explicit)
        Vector T_f = S_f - E_f;

        const Vector& gradPCorr_P = gradPCorr[ownerIdx];

        Scalar nonOrthogonalFlux =
            DUf[face.idx()] * dot(gradPCorr_P, T_f);

        vectorB_(ownerIdx) += nonOrthogonalFlux;
    }
    // ZERO_GRADIENT: no contribution
}