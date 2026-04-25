cmake_minimum_required(VERSION 3.10)
project(Turblyze VERSION "1.0.0" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Find Eigen3 package (required for linear algebra)
find_package(Eigen3 REQUIRED NO_MODULE)

# Find VTK package with required components
find_package(VTK REQUIRED COMPONENTS
  CommonCore
  CommonDataModel
  IOLegacy
  IOXML
)

# --- List your source files ---
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/Primitives/Vector.cpp
    src/Primitives/Tensor.cpp
    src/Primitives/Logger.cpp
    src/Mesh/Face.cpp
    src/Mesh/Cell.cpp
    src/Mesh/MeshReader.cpp
    src/Mesh/MeshChecker.cpp
    src/Fields/CellData.cpp
    src/Fields/FaceData.cpp
    src/BoundaryConditions/BoundaryData.cpp
    src/BoundaryConditions/BoundaryConditions.cpp
    src/Schemes/GradientScheme.cpp
    src/Schemes/ConvectionSchemes.cpp
    src/LinearSystem/Matrix.cpp
    src/LinearSystem/LinearSolvers.cpp
    src/Solver/SIMPLE.cpp
    src/Solver/Constraint.cpp
    src/Models/kOmegaSST.cpp
    src/PostProcessing/VtkWriter.cpp
    src/PostProcessing/VtkCellOrdering.cpp
    src/PostProcessing/VtkBoundaryWriter.cpp
    src/PostProcessing/DerivedFields.cpp
    src/PostProcessing/PvdTimeSeries.cpp
    src/Case/CaseReader.cpp
    src/Application/CFDApplication.cpp
)

# --- Include directories ---
target_include_directories(${PROJECT_NAME} PUBLIC
    src/Mesh/
    src/Primitives/
    src/Fields/
    src/BoundaryConditions/
    src/Schemes/
    src/Models/
    src/Case/
    src/Application/
    src/PostProcessing/
    src/Solver/
    src/LinearSystem/
)

# Set output directory for executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

target_compile_options(${PROJECT_NAME} PRIVATE
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wconversion
    -Wsign-conversion
    -Wdouble-promotion
    -Wformat=2
    -Wcast-align
    -Wunused
    -Wnull-dereference
    -Wimplicit-fallthrough
    -Wold-style-cast
    -Woverloaded-virtual
    -Wnon-virtual-dtor
)

target_compile_definitions(${PROJECT_NAME} PUBLIC PROJECT_USE_DOUBLE_PRECISION)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(${PROJECT_NAME} PRIVATE -O3)
endif()

# Link your executable against the VTK and Eigen3 libraries
target_link_libraries(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES} Eigen3::Eigen)