/******************************************************************************
 * @file Matrix.cpp
 * @brief Matrix assembly and linear system construction for equations
 *****************************************************************************/

#include "Matrix.hpp"

#include <iostream>

#include "ErrorHandler.hpp"

// ************************* Constructor & Destructor *************************

Matrix::Matrix
(
    const Mesh& mesh,
    const BoundaryConditions& boundaryConds
) noexcept
  : mesh_(mesh),
    bcManager_(boundaryConds)
{
    for (const auto& face : mesh_.faces())
    {
        if (face.isBoundary()) ++numBoundaryFaces_;
        else ++numInternalFaces_;
    }

    size_t numCells = mesh_.numCells();

    matrixA_.resize
    (
        static_cast<Eigen::Index>(numCells),
        static_cast<Eigen::Index>(numCells)
    );
    vectorB_.resize(static_cast<Eigen::Index>(numCells));
}


// ****************************** Public Methods ******************************

void Matrix::buildMatrix(const TransportEquation& equation)
{
    clear();
    reserveTripletList();

    size_t numCells = mesh_.numCells();

    for (size_t cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        vectorB_(static_cast<Eigen::Index>(cellIdx)) += 
            equation.source[cellIdx];
    }

    for (size_t faceIdx = 0;
         faceIdx < mesh_.numFaces();
         ++faceIdx)
    {
        const Face& face = mesh_.faces()[faceIdx];

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
    if (alpha <= S(0.0))
    {
        FatalError("Matrix::relax: alpha must be positive");
    }

    const Eigen::Index numCells = matrixA_.rows();

    if (static_cast<Eigen::Index>(phiPrev.size()) != numCells)
    {
        FatalError
        (
            "Matrix::relax: phiPrev size mismatch "
            "with matrix size"
        );
    }

    // Store relaxation factor for setValues to recover
    // the pre-relaxation diagonal
    lastRelaxationFactor_ = alpha;

    // Relax diagonal and update RHS in a single pass
    const Scalar factor = (S(1.0) - alpha) / alpha;

    for (Eigen::Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar origDiag = matrixA_.coeff(cellIdx, cellIdx);

        // Scale diagonal: a_P <- a_P / alpha
        matrixA_.coeffRef(cellIdx, cellIdx) = origDiag / alpha;

        // Update RHS: b <- b + ((1-alpha)/alpha) * a_P * phiPrev
        vectorB_(cellIdx) +=
            factor * origDiag
          * phiPrev[static_cast<size_t>(cellIdx)];
    }
}


// ***************************** Private Methods *****************************

void Matrix::clear()
{
    tripletList_.clear();
    matrixA_.setZero();
    vectorB_.setZero();
    lastRelaxationFactor_ = S(0.0);
}

void Matrix::reserveTripletList()
{
    size_t reserveSize = 4 * numInternalFaces_ + numBoundaryFaces_;

    tripletList_.reserve(reserveSize);
}

Scalar Matrix::extractBoundaryScalar
(
    const BoundaryData& bc,
    std::optional<int> component
) noexcept
{
    if (component && bc.valueType() == BCValueType::VECTOR)
    {
        const Vector& v = bc.vectorValue();
        switch (*component)
        {
            case 0: return v.x();
            case 1: return v.y();
            case 2: return v.z();
            default:
                FatalError
                (
                    "extractBoundaryScalar: "
                    "component index out of range"
                );
        }
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
        mesh_.cells()[neighborIdx].centroid()
      - mesh_.cells()[ownerIdx].centroid();

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
        Scalar dPf = face.dPfMag();
        Scalar dNf = face.dNfMag().value();

        Gammaf =
            dPNMag
          / ((dPf / (G[ownerIdx] + vSmallValue))
           + (dNf / (G[neighborIdx] + vSmallValue)));
    }
    else
    {
        // Face-based: use directly
        if (!equation.GammaFace.has_value())
        {
            FatalError
            (
                "assembleInternalFace: equation has "
                "no Gamma or GammaFace"
            );
        }
        Gammaf = equation.GammaFace->get()[face.idx()];
    }

    Scalar aDiff = Gammaf * Ef.magnitude() / (dPNMag + vSmallValue);

    // Convection coefficients
    Scalar aPConv = S(0.0);
    Scalar aNConv = S(0.0);

    if (equation.flowRate)
    {
        Scalar flowRate = equation.flowRate->get()[face.idx()];

        auto [aP, aN] = ConvectionScheme::getFluxCoefficients(flowRate);

        aPConv = aP;
        aNConv = aN;
    }

    // Matrix coefficients for owner and neighbor cells
    tripletList_.emplace_back(ownerIdx, ownerIdx, aDiff + aPConv);

    tripletList_.emplace_back(ownerIdx, neighborIdx, -aDiff + aNConv);

    tripletList_.emplace_back(neighborIdx, neighborIdx, aDiff - aNConv);

    tripletList_.emplace_back(neighborIdx, ownerIdx, -aDiff - aPConv);

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
            face.idx(),
            equation.componentIdx
        );

    Scalar nonOrthogonalFlux = Gammaf * dot(gradPhif, Tf);

    vectorB_(static_cast<Eigen::Index>(ownerIdx))    += nonOrthogonalFlux;
    vectorB_(static_cast<Eigen::Index>(neighborIdx)) -= nonOrthogonalFlux;

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

        vectorB_(static_cast<Eigen::Index>(ownerIdx))    -= deferredCorrection;
        vectorB_(static_cast<Eigen::Index>(neighborIdx)) += deferredCorrection;
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

    const Vector ePf = face.dPf().normalized();
    const Scalar dPfMag = face.dPfMag();

    Vector Ef = (dot(Sf, Sf) / dot(Sf, ePf)) * ePf;

    // Diffusion coefficient at boundary face
    Scalar Gammaf;

    if (equation.Gamma)
    {
        Gammaf = equation.Gamma->get()[ownerIdx];
    }
    else
    {
        if (!equation.GammaFace.has_value())
        {
            FatalError
            (
                "assembleBoundaryFace: equation has "
                "no Gamma or GammaFace"
            );
        }
        Gammaf = equation.GammaFace->get()[face.idx()];
    }

    Scalar aDiff = Gammaf * Ef.magnitude() / (dPfMag + vSmallValue);

    using enum BCType;
    if
    (
        bc->type() == FIXED_VALUE
     || bc->type() == NO_SLIP
    )
    {
        // Dirichlet BC: phi_b is prescribed
        Scalar phiB = S(0.0);

        if (bc->type() == NO_SLIP)
        {
            phiB = S(0.0);
        }
        else
        {
            phiB = extractBoundaryScalar(*bc, equation.componentIdx);
        }

        // Diffusion contribution
        tripletList_.emplace_back(ownerIdx, ownerIdx, aDiff);
        vectorB_(static_cast<Eigen::Index>(ownerIdx)) += aDiff * phiB;

        // Convection contribution
        if (equation.flowRate)
        {
            Scalar aConv = equation.flowRate->get()[face.idx()];
            vectorB_(static_cast<Eigen::Index>(ownerIdx)) -= aConv * phiB;
        }

        // Non-orthogonal correction
        Vector Tf = Sf - Ef;

        Vector gradPhif =
            equation.gradScheme.faceGradient
            (
                equation.fieldName,
                equation.phi,
                equation.gradPhi[ownerIdx],
                Vector{},
                face.idx(),
                equation.componentIdx
            );

        Scalar nonOrthogonalFlux =
            Gammaf * dot(gradPhif, Tf);

        vectorB_(static_cast<Eigen::Index>(ownerIdx)) += nonOrthogonalFlux;
    }
    else if
    (
        bc->type() == ZERO_GRADIENT
     || bc->type() == K_WALL_FUNCTION
     || bc->type() == NUT_WALL_FUNCTION
     || bc->type() == OMEGA_WALL_FUNCTION
    )
    {
        // Zero normal gradient: only convection
        if (equation.flowRate)
        {
            Scalar aConv = equation.flowRate->get()[face.idx()];

            tripletList_.emplace_back(ownerIdx, ownerIdx, aConv);

            // Tangential gradient correction
            const Vector& gradPhiP = equation.gradPhi[ownerIdx];

            Vector eb = face.normal();
            Scalar normalComponent = dot(gradPhiP, eb);
            Vector gradPhib = gradPhiP - normalComponent * eb;

            Vector dCb = face.centroid() - mesh_.cells()[ownerIdx].centroid();

            Scalar correction = aConv * dot(gradPhib, dCb);
            vectorB_(static_cast<Eigen::Index>(ownerIdx)) -= correction;
        }
        // No convection + zero gradient = no contribution
    }
    else
    {
        // Unhandled BC type: default to zero gradient
        Warning
        (
            "Undefined boundary condition type for field "
          + equation.fieldName + " on patch "
          + face.patch()->patchName()
          + ". Applying zero gradient."
        );

        if (equation.flowRate)
        {
            Scalar aConv = equation.flowRate->get()[face.idx()];
            tripletList_.emplace_back(ownerIdx, ownerIdx, aConv);
        }
    }
}

void Matrix::setValues
(
    std::span<const size_t> cellIndices,
    std::span<const Scalar> values,
    std::span<const Scalar> fractions
)
{
    const bool hasFractions = !fractions.empty();

    for (size_t i = 0; i < cellIndices.size(); ++i)
    {
        const size_t cellIdx = cellIndices[i];
        const Scalar f = hasFractions ? fractions[i] : S(1.0);

        if (f < vSmallValue)
        {
            continue;
        }

        const Scalar diag =
            matrixA_.coeff
            (
                static_cast<Eigen::Index>(cellIdx),
                static_cast<Eigen::Index>(cellIdx)
            );

        if (f > S(1.0) - vSmallValue)
        {
            for 
            (
                Eigen::SparseMatrix<Scalar>::InnerIterator 
                it(matrixA_, static_cast<Eigen::Index>(cellIdx));
                it; 
                ++it
            )
            {
                const size_t row = static_cast<size_t>(it.row());
                if (row == cellIdx) continue;

                // Transfer coupling to neighbor source
                vectorB_(static_cast<Eigen::Index>(row)) -= it.value() * values[i];
                it.valueRef() = S(0.0);
            }

            for (size_t neighborIdx : mesh_.cells()[cellIdx].neighborCellIndices())
            {
                matrixA_.coeffRef
                (
                    static_cast<Eigen::Index>(cellIdx),
                    static_cast<Eigen::Index>(neighborIdx)
                ) = S(0.0);
            }

            vectorB_(static_cast<Eigen::Index>(cellIdx)) = diag * values[i];
        }
        else
        {
            Scalar diagPre = diag;
            if (lastRelaxationFactor_ > S(0.0))
            {
                diagPre = lastRelaxationFactor_ * diag;
            }

            Scalar coeff = f / (S(1.0) - f) * diagPre;

            matrixA_.coeffRef
            (
                static_cast<Eigen::Index>(cellIdx),
                static_cast<Eigen::Index>(cellIdx)
            ) += coeff;
            
            vectorB_(static_cast<Eigen::Index>(cellIdx)) += coeff * values[i];
        }
    }
}
