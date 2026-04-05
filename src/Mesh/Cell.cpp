/******************************************************************************
 * @file Cell.cpp
 * @brief Implementation of cell geometric properties and operations
 *****************************************************************************/

#include "Cell.hpp"

#include <cmath>
#include <ostream>

#include "ErrorHandler.hpp"
#include <iomanip>

// *********************** Geometric Property Methods ***********************

void Cell::calculateGeometricProperties
(
    std::span<const Face> allFaces,
    std::span<const FaceIntegrals> allFaceIntegrals
)
{
    geometricPropertiesCalculated_ = false;
    volume_ = 0.0;
    centroid_ = Vector{};
    Vector centroidSum;

    for (size_t i = 0; i < faceIndices_.size(); ++i)
    {
        size_t faceIndex = faceIndices_[i];
        const Face& face = allFaces[faceIndex];

        if(!face.geometricPropertiesCalculated())
        {
            FatalError
            (
                "Cell " + std::to_string(idx_)
              + " calculation: Geometric properties for "
                "bounding Face " + std::to_string(face.idx())
              + " were not calculated."
            );
        }

        Scalar faceSign = S(faceSigns_[i]);
        const FaceIntegrals& integrals = allFaceIntegrals[faceIndex];

        volume_ += faceSign * integrals.volume;

        centroidSum += Vector
        (
            faceSign * integrals.x2,
            faceSign * integrals.y2,
            faceSign * integrals.z2
        );
    }

    volume_ /= S(3.0);

    if (std::abs(volume_) < smallValue)
    {
        FatalError
        (
            "Cell " + std::to_string(idx_) + " has zero volume"
        );
    }

    centroid_ = centroidSum / (S(2.0) * volume_);
    geometricPropertiesCalculated_ = true;
}

std::ostream& operator<<(std::ostream& os, const Cell& c)
{
    os  << "Cell(ID: " << c.idx() << ", Faces: [";

    const auto& faces = c.faceIndices();
    for (size_t i = 0; i < faces.size(); ++i)
    {
        os  << faces[i]
            << (i == faces.size() - 1 ? "" : ", ");
    }

    os  << "], Neighbors: [";

    const auto& neighbors = c.neighborCellIndices();
    for (size_t i = 0; i < neighbors.size(); ++i)
    {
        os  << neighbors[i]
            << (i == neighbors.size() - 1 ? "" : ", ");
    }

    os  << ']';

    if (c.geometricPropertiesCalculated())
    {
        auto flags = os.flags();
        auto prec = os.precision();

        os  << std::fixed
            << std::setprecision(6);

        os  << ", Volume: " << c.volume()
            << ", Centroid: " << c.centroid();

        os.flags(flags);
        os.precision(prec);
    }
    else
    {
        os  << ", Geometry: N/A";
    }

    os  << ')';

    return os;
}
