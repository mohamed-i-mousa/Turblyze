#ifndef FACEDATA_H
#define FACEDATA_H

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>

#include "Scalar.h"
#include "Vector.h"

template<typename T>
class FaceData {
public:
  std::string name;
  std::vector<T> allFacesValues;

  // ----- Constructors ----- //
  FaceData(const std::string& fieldName, size_t numFaces)
    : name(fieldName), allFacesValues(numFaces) {}

  FaceData(const std::string& fieldName, size_t numFaces, const T& initialValue)
    : name(fieldName), allFacesValues(numFaces, initialValue) {}

  // ----- Access operators ----- //
  T& operator[](size_t faceIndex) {
    if (faceIndex >= allFacesValues.size()) {
      throw std::out_of_range("Face index out of range in FaceField '" + name + "'");
    }
    return allFacesValues[faceIndex];
  }

  const T& operator[](size_t faceIndex) const {
    if (faceIndex >= allFacesValues.size()) {
      throw std::out_of_range("Face index out of range in FaceField '" + name + "'");
    }
    return allFacesValues[faceIndex];
  }

  // ----- Methods ----- //

  // Get the size of the FaceData
  size_t size() const {
    return allFacesValues.size();
  }

  // Initialize values
  void setAll(const T& value) {
    for (size_t i = 0; i < allFacesValues.size(); i++) {
      allFacesValues[i] = value;
    }
  }
  
  // Print a summary
  void printSummary(size_t itemsToShow = 5) const {
    std::cout << "FaceField: " << name << " (Size: " << allFacesValues.size() << ")" << std::endl;
    for (size_t i = 0; i < std::min(allFacesValues.size(), itemsToShow); ++i) {
      std::cout << "  Face " << i << ": " << allFacesValues[i] << std::endl; // Requires T to be streamable
    }
    if (allFacesValues.size() > itemsToShow) {
      std::cout << "  ..." << std::endl;
    }
  }
};

using FaceFluxField = FaceData<Scalar>; // e.g., for phi in p-U coupling
using FaceVectorField = FaceData<Vector>;

#endif
