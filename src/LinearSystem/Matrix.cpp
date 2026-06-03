/******************************************************************************
 * @file Matrix.cpp
 * @brief Matrix assembly and linear system construction for equations
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Matrix.h"

// Standard library headers
#include <iostream>
#include <iterator>

// External library headers
#include <omp.h>

// Project headers
#include "ErrorHandler.h"
#include "LinearInterpolation.h"

// ************************* Special Member Functions *************************

Matrix::Matrix
(
    const Mesh& mesh,
    const BoundaryConditions& boundaryConds
) noexcept
:
    mesh_(mesh),
    bcManager_(boundaryConds)
{
    for (const auto& face : mesh_.faces())
    {
        if (face.isBoundary()) ++numBoundaryFaces_;
        else ++numInternalFaces_;
    }

    const Count numCells = mesh_.numCells();

    matrixA_.resize
    (
        eIdx(numCells),
        eIdx(numCells)
    );
    vectorB_.resize(eIdx(numCells));

    // Reserve reusable per-thread assembly buffers.
    const Count T = static_cast<Count>(omp_get_max_threads());
    const Count estimatedTriplets = 4 * numInternalFaces_ + numBoundaryFaces_;
    const Count reservePerThread = (estimatedTriplets / T) + 1;

    tripletList_.reserve(estimatedTriplets);

    perThreadTriplets_.resize(T);
    for (auto& triplets : perThreadTriplets_)
    {
        triplets.reserve(reservePerThread);
    }

    perThreadB_.assign(T, Vec::Zero(eIdx(numCells)));
}

// ****************************** Public Methods ******************************

void Matrix::buildMatrix(const TransportEquation& equation)
{
    clear();

    const Count numCells = mesh_.numCells();
    const Count numFaces = mesh_.numFaces();

    // Cell-source contribution: race-free per-cell write directly to vectorB_
    #pragma omp parallel for schedule(static)
    for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        vectorB_(eIdx(cellIdx)) += equation.source[cellIdx];
    }

    // Reset per-thread scratch (capacity retained from constructor)
    for (auto& triplets : perThreadTriplets_)
    {
        triplets.clear();
    }
    for (auto& localB : perThreadB_)
    {
        localB.setZero();
    }

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        auto& triplets = perThreadTriplets_[static_cast<Index>(tid)];
        auto& localB = perThreadB_[static_cast<Index>(tid)];

        #pragma omp for schedule(static)
        for (Index faceIdx = 0; faceIdx < numFaces; ++faceIdx)
        {
            const Face& face = mesh_.faces()[faceIdx];

            if (face.isBoundary())
            {
                assembleBoundaryFace(face, equation, triplets, localB);
            }
            else
            {
                assembleInternalFace(face, equation, triplets, localB);
            }
        }
    }

    // Merge per-thread triplet lists into the member tripletList_
    Count totalTriplets = 0;
    for (const auto& v : perThreadTriplets_) totalTriplets += v.size();
    tripletList_.reserve(totalTriplets);
    for (auto& v : perThreadTriplets_)
    {
        tripletList_.insert
        (
            tripletList_.end(),
            std::make_move_iterator(v.begin()),
            std::make_move_iterator(v.end())
        );
    }

    // Merge per-thread RHS contributions into vectorB_
    for (const auto& v : perThreadB_)
    {
        vectorB_ += v;
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

    if (eIdx(phiPrev.size()) != numCells)
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

    const Scalar factor = (S(1.0) - alpha) / alpha;

    #pragma omp parallel for schedule(static)
    for (Eigen::Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar origDiag = matrixA_.coeff(cellIdx, cellIdx);

        // Scale diagonal: a_P <- a_P / alpha
        matrixA_.coeffRef(cellIdx, cellIdx) = origDiag / alpha;

        // Update RHS: b <- b + ((1-alpha)/alpha) * a_P * phiPrev
        vectorB_(cellIdx) +=
            factor * origDiag
          * phiPrev[static_cast<Index>(cellIdx)];
    }
}

// ****************************** Private Methods *****************************

void Matrix::clear()
{
    tripletList_.clear();
    matrixA_.setZero();
    vectorB_.setZero();
    lastRelaxationFactor_ = S(0.0);
}


Scalar Matrix::faceDiffusionCoefficient
(
    const Face& face,
    const TransportEquation& equation
) const
{
    if (equation.Gamma)
    {
        const auto& gamma = equation.Gamma->get();

        if (face.isBoundary())
        {
            return gamma[face.ownerCell()];
        }

        return interpolateToFace(face, gamma);
    }

    if (!equation.GammaFace.has_value())
    {
        FatalError
        (
            "faceDiffusionCoefficient: equation has no Gamma or GammaFace"
        );
    }

    return equation.GammaFace->get()[face.idx()];
}

// ****************************** Private Methods *****************************

void Matrix::assembleInternalFace
(
    const Face& face,
    const TransportEquation& equation,
    std::vector<Eigen::Triplet<Scalar>>& triplets,
    Vec& localB
) const
{
    const Index ownerIdx = face.ownerCell();
    const Index neighborIdx = face.neighborCell().value();
    const Vector Sf = face.normal() * face.projectedArea();
    const Vector dPN =
        mesh_.cells()[neighborIdx].centroid()
      - mesh_.cells()[ownerIdx].centroid();
    const Scalar dPNMag = magnitude(dPN);
    const Vector ePN = dPN / (dPNMag + vSmallValue);

    // Orthogonal component (over-relaxed)
    const Vector Ef = (dot(Sf, Sf) / dot(Sf, ePN)) * ePN;
    const Scalar Gammaf = faceDiffusionCoefficient(face, equation);
    const Scalar aDiff = Gammaf * magnitude(Ef) / (dPNMag + vSmallValue);

    // Convection coefficients
    Scalar flowRate = S(0.0);
    Scalar aPConv = S(0.0);
    Scalar aNConv = S(0.0);

    if (equation.convection)
    {
        const ConvectionTerm& convection = *equation.convection;

        flowRate = convection.flowRate[face.idx()];

        const auto [aP, aN] = ConvectionSchemes::getFluxCoefficients(flowRate);

        aPConv = aP;
        aNConv = aN;
    }

    // Matrix coefficients for owner and neighbor cells
    triplets.emplace_back(ownerIdx, ownerIdx, aDiff + aPConv);
    triplets.emplace_back(ownerIdx, neighborIdx, -aDiff + aNConv);
    triplets.emplace_back(neighborIdx, neighborIdx, aDiff - aNConv);
    triplets.emplace_back(neighborIdx, ownerIdx, -aDiff - aPConv);

    // Non-orthogonal correction (explicit)
    const Vector Tf = Sf - Ef;
    const Vector& gradPhiP = equation.gradPhi[ownerIdx];
    const Vector& gradPhiN = equation.gradPhi[neighborIdx];
    const Vector gradPhif =
        equation.gradScheme.faceGradient
        (
            equation.field,
            equation.phi,
            gradPhiP,
            gradPhiN,
            face.idx()
        );
    const Scalar nonOrthogonalFlux = Gammaf * dot(gradPhif, Tf);

    localB(eIdx(ownerIdx)) += nonOrthogonalFlux;
    localB(eIdx(neighborIdx)) -= nonOrthogonalFlux;

    // Deferred correction (explicit, convection only)
    if (equation.convection)
    {
        const ConvectionTerm& convection = *equation.convection;

        const Scalar deferredCorrection =
            convection.scheme.correction
            (
                face,
                equation.phi,
                gradPhiP,
                gradPhiN,
                flowRate
            );

        localB(eIdx(ownerIdx)) -= deferredCorrection;
        localB(eIdx(neighborIdx)) += deferredCorrection;
    }
}


void Matrix::assembleBoundaryFace
(
    const Face& face,
    const TransportEquation& equation,
    std::vector<Eigen::Triplet<Scalar>>& triplets,
    Vec& localB
) const
{
    const Index ownerIdx = face.ownerCell();
    const BoundaryData& bc =
        bcManager_.fieldBC(face.patch()->get().patchName(), equation.field);
    const Vector Sf = face.normal() * face.projectedArea();
    const Vector ePf = normalized(face.dPf());
    const Scalar dPfMag = face.dPfMag();
    const Vector Ef = (dot(Sf, Sf) / dot(Sf, ePf)) * ePf;
    const Scalar Gammaf = faceDiffusionCoefficient(face, equation);
    const Scalar aDiff = Gammaf * magnitude(Ef) / (dPfMag + vSmallValue);
    const ConvectionTerm* convection =
        equation.convection ? &*equation.convection : nullptr;

    using enum BCType;
    if
    (
        bc.type() == fixedValue
     || bc.type() == noSlip
    )
    {
        // Dirichlet BC: phiB is prescribed
        Scalar phiB = S(0.0);

        if (bc.type() != noSlip)
        {
            phiB = bc.fixedScalarValue();
        }

        // Diffusion contribution
        triplets.emplace_back(ownerIdx, ownerIdx, aDiff);
        localB(eIdx(ownerIdx)) += aDiff * phiB;

        // Convection contribution
        if (convection != nullptr)
        {
            const Scalar aConv = convection->flowRate[face.idx()];

            localB(eIdx(ownerIdx)) -= aConv * phiB;
        }

    }
    else if (bc.type() == fixedGradient)
    {
        const Scalar gradient = bc.fixedScalarGradient(); 
        const Scalar dn = dot(face.dPf(), face.normal());

        localB(eIdx(ownerIdx)) +=
            Gammaf * gradient * face.projectedArea();

        // Boundary value for convection: phi_b = phi_P + grad * dn
        if (convection != nullptr)
        {
            const Scalar aConv = convection->flowRate[face.idx()];

            triplets.emplace_back(ownerIdx, ownerIdx, aConv);
            localB(eIdx(ownerIdx)) -= aConv * gradient * dn;
        }
    }
    else if
    (
        bc.type() == zeroGradient
     || bc.type() == kWallFunction
     || bc.type() == nutWallFunction
     || bc.type() == omegaWallFunction
    )
    {
        // Zero normal gradient: only convection
        if (convection != nullptr)
        {
            const Scalar aConv = convection->flowRate[face.idx()];

            triplets.emplace_back(ownerIdx, ownerIdx, aConv);
        }
        // No convection + zero gradient = no contribution
    }
    else
    {
        // Unhandled BC type: default to zero gradient
        Warning
        (
            "Undefined boundary condition type for field "
          + Name(fieldToString(equation.field)) + " on patch "
          + face.patch()->get().patchName()
          + ". Applying zero gradient."
        );

        if (convection != nullptr)
        {
            const Scalar aConv = convection->flowRate[face.idx()];

            triplets.emplace_back(ownerIdx, ownerIdx, aConv);
        }
    }
}


void Matrix::setValues
(
    IndexListRef cellIndices,
    ScalarListRef values,
    ScalarListRef fractions
)
{
    const bool hasFractions = !fractions.empty();

    for (Index i = 0; i < cellIndices.size(); ++i)
    {
        const Index cellIdx = cellIndices[i];
        const Scalar f = hasFractions ? fractions[i] : S(1.0);

        if (f < rootSmallValue_)
        {
            continue;
        }

        const Scalar diag =
            matrixA_.coeff
            (
                eIdx(cellIdx),
                eIdx(cellIdx)
            );

        if (f > S(1.0) - rootSmallValue_)
        {
            const Eigen::Index c = eIdx(cellIdx);

            for
            (
                Index neighborIdx
              :  mesh_.cells()[cellIdx].neighborCellIndices()
            )
            {
                const Eigen::Index r = eIdx(neighborIdx);

                const Scalar coupling = matrixA_.coeff(r, c);
                if (coupling != S(0.0))
                {
                    vectorB_(r) -= coupling * values[i];
                    matrixA_.coeffRef(r, c) = S(0.0);
                }

                matrixA_.coeffRef(c, r) = S(0.0);
            }

            vectorB_(c) = diag * values[i];
        }
        else
        {
            const Scalar diagPre =
                lastRelaxationFactor_ > S(0.0)
              ? lastRelaxationFactor_ * diag
              : diag;

            const Scalar coeff = f / (S(1.0) - f) * diagPre;

            matrixA_.coeffRef
            (
                eIdx(cellIdx),
                eIdx(cellIdx)
            ) += coeff;
            
            vectorB_(eIdx(cellIdx)) += coeff * values[i];
        }
    }
}
