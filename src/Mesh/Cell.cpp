/******************************************************************************
 * @file Cell.cpp
 * @brief Implementation of cell geometric properties and operations
 *****************************************************************************/

#include "Cell.hpp"

#include <stdexcept>
#include <cmath>
#include <ostream>
#include <iomanip>

// *********************** Geometric Property Methods ***********************

void Cell::calculateGeometricProperties(std::span<const Face> allFaces)
{
    geometricPropertiesCalculated_ = false;
    volume_ = 0.0;
    centroid_ = Vector(0.0, 0.0, 0.0);
    Vector centroidSum(0.0, 0.0, 0.0);

    for (size_t i = 0; i < faceIndices_.size(); ++i)
    {
        size_t faceIndex = faceIndices_[i];
        const Face& face = allFaces[faceIndex];

        if(!face.geometricPropertiesCalculated())
        {
            throw
                std::runtime_error
                (
                    "Error in Cell " + std::to_string(idx_)
                  + " calculation:"
                  + " Geometric properties for bounding Face "
                  + std::to_string(face.idx())
                  + " were not calculated."
                );
        }

        Scalar faceSign = S(faceSigns_[i]);

        volume_ += faceSign * face.volumeContribution();

        centroidSum += Vector
        (
            faceSign * face.x2Integral(),
            faceSign * face.y2Integral(),
            faceSign * face.z2Integral()
        );
    }

    volume_ /= S(3.0);

    if (std::abs(volume_) > smallValue)
    {
        centroid_ = centroidSum / (S(2.0) * volume_);
    }
    else
    {
        throw
            std::runtime_error
            (
                "Cell " + std::to_string(idx_) + " has zero volume"
            );
    }

    if (volume_ < S(0.0))
    {
        throw
            std::runtime_error
            (
                "Error: Cell " + std::to_string(idx_)
              + " calculated negative volume ("
              + std::to_string(volume_)
              + "). Check face normal conventions"
              + " and mesh connectivity."
            );
    }
    else
    {
        geometricPropertiesCalculated_ = true;
    }
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

    os  << "]";

    if (c.geometricPropertiesCalculated())
    {
        std::ios_base::fmtflags flags = os.flags();
        auto prec = os.precision();

        os  << std::fixed
            << std::setprecision(6);

        os  << ", Volume: " << c.volume() << ", Centroid: " << c.centroid();

        os.flags(flags);
        os.precision(prec);
    }
    else
    {
        os  << ", Geometry: N/A";
    }

    os  << ")";

    return os;
}
