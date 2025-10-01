cmake_minimum_required(VERSION 3.10)
project(MyCFDCode VERSION "1.0.0" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

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
    src/Mesh/Vector.cpp
    src/Mesh/Face.cpp
    src/Mesh/Cell.cpp
    src/Mesh/CellData.cpp
    src/Mesh/FaceData.cpp
    src/Mesh/MeshReader.cpp
    src/Mesh/checkMesh.cpp
    src/BoundaryConditions/BoundaryData.cpp
    src/BoundaryConditions/BoundaryPatch.cpp
    src/BoundaryConditions/BoundaryConditions.cpp
    src/Numerics/GradientScheme.cpp
    src/Numerics/ConvectionScheme.cpp
    src/Numerics/Matrix.cpp
    src/Numerics/flowRate.cpp
    src/Numerics/LinearSolvers.cpp
    src/Numerics/SIMPLE.cpp
    src/Numerics/Constraint.cpp
    src/Numerics/linearInterpolation.cpp
    src/Models/KOmegaSST.cpp
    src/PostProcessing/VtkWriter.cpp
    src/Setup/SetupReader.cpp
)

# --- Tell CMake about your include directories ---
target_include_directories(${PROJECT_NAME} PUBLIC
    include/Core/
    include/Mesh/
    include/BoundaryConditions/
    include/Numerics/
    include/Models/
    include/PostProcessing/
    include/Setup/
)

# Set output directory for executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)

target_compile_definitions(MyCFDCode PUBLIC PROJECT_USE_DOUBLE_PRECISION)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(${PROJECT_NAME} PRIVATE -O3)
endif()

# Link your executable against the VTK libraries
target_link_libraries(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES})
