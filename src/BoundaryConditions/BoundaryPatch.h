#ifndef BOUNDARYPATCH_H
#define BOUNDARYPATCH_H

#include <string>
#include <vector>
#include <map>

#include "Face.h"

// Enum contains the common boundary condition
enum class BoundaryConditionType {
  VELOCITY_INLET,
  PRESSURE_INLET,
  PRESSURE_OUTLET,
  WALL,
  SYMMETRY,
  PERIODIC,
  MASS_FLOW_INLET,
  OUTFLOW,
  INTERFACE,
  INTERIOR,
  SOLID,
  FLUID,
  UNDEFINED
};

// Function reads the boundary type from the mesh file
inline BoundaryConditionType mapFluentBCToEnum(const std::string& fluentType) {
  if (fluentType == "velocity-inlet") return BoundaryConditionType::VELOCITY_INLET;
  if (fluentType == "pressure-inlet") return BoundaryConditionType::PRESSURE_INLET;
  if (fluentType == "pressure-outlet") return BoundaryConditionType::PRESSURE_OUTLET;
  if (fluentType == "wall") return BoundaryConditionType::WALL;
  if (fluentType == "symmetry") return BoundaryConditionType::SYMMETRY;
  if (fluentType == "periodic" || fluentType == "periodic-shadow") return BoundaryConditionType::PERIODIC;
  if (fluentType == "mass-flow-inlet") return BoundaryConditionType::MASS_FLOW_INLET;
  if (fluentType == "outflow") return BoundaryConditionType::OUTFLOW;
  if (fluentType == "interface") return BoundaryConditionType::INTERFACE;
  if (fluentType == "interior") return BoundaryConditionType::INTERIOR;
  if (fluentType == "solid") return BoundaryConditionType::SOLID;
  if (fluentType == "fluid") return BoundaryConditionType::FLUID;

  // If none of the above
  std::cerr << "Warning: Unknown Fluent boundary type encountered: " << fluentType << std::endl;
  return BoundaryConditionType::UNDEFINED;
}

struct BoundaryPatch {
  // ----- Members ----- //
  std::string patchName;
  std::string fluentType;
  BoundaryConditionType type;
  size_t zoneID;
  size_t firstFaceIndex;
  size_t lastFaceIndex;

  // ----- Constructor ----- //
  BoundaryPatch(size_t id, size_t start_id, size_t end_id)
      : zoneID(id), firstFaceIndex(start_id), lastFaceIndex(end_id) {}

  // ----- Functions ----- //

  size_t getNumberOfBoundaryFaces() const {
      size_t numFaces;
      numFaces = lastFaceIndex - firstFaceIndex + 1;
      return numFaces;
  }
};

#endif
