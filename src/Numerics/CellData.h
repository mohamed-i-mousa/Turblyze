#ifndef CELLDATA_H
#define CELLDATA_H

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>

#include "Scalar.h"
#include "Vector.h"

template<typename T>
class CellData {
public:
  std::string name;
  std::vector<T> internalField;

  // ----- Constructors ----- //
  CellData(const std::string& fieldName, size_t numCells)
    : name(fieldName), internalField(numCells) {}

  CellData(const std::string& fieldName, size_t numCells, const T& initialValue)
    : name(fieldName), internalField(numCells, initialValue) {}

  // ----- Access operators ----- //
  T& operator[](size_t cellIndex) {
    if (cellIndex >= internalField.size()) {
      throw std::out_of_range("Cell index out of range in CellData '" + name + "'");
    }
    return internalField[cellIndex];
  }

  const T& operator[](size_t cellIndex) const{
    if (cellIndex >= internalField.size()) {
      throw std::out_of_range("Cell index out of range in CellData '" + name + "'");
    }
    return internalField[cellIndex];
  }

  // ----- Methods ----- //

  // Get the size of the field
  size_t size() const {
    return internalField.size();
  }

  // Initialize values of the field
  void setAll(const T& value) {
    for (size_t i = 0; i < internalField.size(); i++) {
      internalField[i] = value;
    }
  }

  // Print a summary
  void printSummary(size_t itemsToShow) const {
    std::cout << "CellData: " << name << " (Size: " << internalField.size() << ")" << std::endl;
    for (size_t i = 0; i < std::min(internalField.size(), itemsToShow); ++i) {
        std::cout << "  Cell " << i << ": " << internalField[i] << std::endl; // Requires T to be streamable
    }
    if (internalField.size() > itemsToShow) {
        std::cout << "  ..." << std::endl;
    }
  }
};

using VelocityField = CellData<Vector>;
using PressureField = CellData<Scalar>;
using ScalarField = CellData<Scalar>;
using VectorField = CellData<Vector>;

#endif
