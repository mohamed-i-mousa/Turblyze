cmake_minimum_required(VERSION 3.10)
project(MyCFDCode VERSION "1.0.0" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- List your source files ---
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/Mesh/MeshReader.cpp
    src/Numerics/GradientScheme.cpp
    src/Numerics/ConvectionScheme.cpp
    src/Numerics/MatrixConstructor.cpp
    src/Numerics/LinearSolvers.cpp
    src/Numerics/SIMPLE.cpp
    src/Numerics/KOmegaSST.cpp
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

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(${PROJECT_NAME} PRIVATE -O3)
endif()
