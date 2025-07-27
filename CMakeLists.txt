cmake_minimum_required(VERSION 3.10)
project(MyCFDCode LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- List your source files ---
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/Mesh/Scalar.h 
    src/Mesh/Vector.h
    src/Mesh/Face.h
    src/Mesh/Cell.h
    src/Mesh/MeshReader.h
    src/Mesh/MeshReader.cpp
    src/BoundaryConditions/BoundaryPatch.h
    src/BoundaryConditions/BoundaryConditions.h
    src/BoundaryConditions/BoundaryData.h
    src/Numerics/CellData.h
    src/Numerics/FaceData.h
    src/Numerics/GradientScheme.h
    src/Numerics/GradientScheme.cpp
    src/Numerics/ConvectionScheme.h
    src/Numerics/ConvectionScheme.cpp
    src/Numerics/MatrixConstructor.h
    src/Numerics/MatrixConstructor.cpp
    src/Numerics/LinearSolvers.h
    src/Numerics/SIMPLE.h
    src/Numerics/SIMPLE.cpp
    src/PostProcessing/VtkWriter.h
    src/PostProcessing/VtkWriter.cpp
)

target_compile_definitions(MyCFDCode PUBLIC PROJECT_USE_DOUBLE_PRECISION)

# --- Tell CMake about your subdirectories for includes ---
target_include_directories(${PROJECT_NAME} PUBLIC
    src 
    src/Mesh
    src/BoundaryConditions
    src/Numerics
    src/PostProcessing
)

# Set output directory for executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
