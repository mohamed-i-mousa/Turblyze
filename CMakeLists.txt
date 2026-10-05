# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: 2025-2026 Mohamed Mousa
# SPDX-License-Identifier: Apache-2.0
# -----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.20)
project(Turblyze VERSION "1.0.0" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set a default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Build options
option(TURBLYZE_DOUBLE_PRECISION
       "Build with double-precision Scalar" ON)
option(TURBLYZE_NATIVE_ARCH
       "Optimize for the build host CPU (-march=native)" ON)

# Find PETSc via pkg-config
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

# Find MPI
find_package(MPI REQUIRED COMPONENTS CXX)

# Find METIS (mesh partitioning)
find_path(METIS_INCLUDE_DIR metis.h REQUIRED)
find_library(METIS_LIBRARY metis REQUIRED)

# Find the HDF5 C library
set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED COMPONENTS C)

# Warning flags
set(TURBLYZE_WARNING_FLAGS
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wconversion
    -Wdouble-promotion
    -Wformat=2
    -Wcast-align
    -Wnull-dereference
    -Wimplicit-fallthrough
    -Wold-style-cast
    -Woverloaded-virtual
    -Wnon-virtual-dtor
)

# Release optimization flags
set(TURBLYZE_RELEASE_FLAGS -O3)
if(TURBLYZE_NATIVE_ARCH)
    list(APPEND TURBLYZE_RELEASE_FLAGS -march=native)
endif()

# Apply compiler warning and Release optimization flags to a target
function(turblyze_apply_build_flags target)
    target_compile_options(${target} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:${TURBLYZE_WARNING_FLAGS}>
    )
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(${target} PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:${TURBLYZE_RELEASE_FLAGS}>
        )
    endif()
endfunction()

# Source directories
add_library(turblyze_core OBJECT
    src/Primitives/Logger.cpp
    src/Mesh/Face.cpp
    src/Mesh/Cell.cpp
    src/Mesh/MeshReader.cpp
    src/Mesh/MeshCreator.cpp
    src/Mesh/MeshChecker.cpp
    src/Fields/Field.cpp
    src/BoundaryConditions/BoundaryConditions.cpp
    src/BoundaryConditions/BoundaryTypeFactory.cpp
    src/BoundaryConditions/BCLoader.cpp
    src/BoundaryConditions/BoundaryTypes/FixedValue.cpp
    src/BoundaryConditions/BoundaryTypes/FixedGradient.cpp
    src/BoundaryConditions/BoundaryTypes/ZeroGradient.cpp
    src/BoundaryConditions/BoundaryTypes/Symmetry.cpp
    src/Schemes/GradientSchemes/GradientScheme.cpp
    src/Schemes/GradientSchemes/LeastSquares.cpp
    src/Schemes/ConvectionSchemes/ConvectionScheme.cpp
    src/Schemes/ConvectionSchemes/Upwind.cpp
    src/Schemes/ConvectionSchemes/CentralDifference.cpp
    src/Schemes/ConvectionSchemes/SecondOrderUpwind.cpp
    src/Schemes/ConvectionSchemes/LUST.cpp
    src/Schemes/TimeSchemes/TimeScheme.cpp
    src/LinearSystem/Matrix.cpp
    src/LinearSystem/LinearSolvers.cpp
    src/LinearSystem/PETScRuntime.cpp
    src/Parallel/MPIRuntime.cpp
    src/Parallel/Comm.cpp
    src/Parallel/Reduce.cpp
    src/Parallel/GlobalIndex.cpp
    src/Parallel/MeshDecomposer.cpp
    src/Parallel/MeshDistributor.cpp
    src/Parallel/HaloExchange.cpp
    src/Solver/MomentumTransport.cpp
    src/Solver/Segregated.cpp
    src/Solver/SIMPLE.cpp
    src/Solver/PISO.cpp
    src/Models/Turbulence/TurbulenceModel.cpp
    src/Models/Turbulence/RANS/RANS.cpp
    src/Models/Turbulence/RANS/kOmegaSST.cpp
    src/PostProcessing/VTK/HDF5Wrapper.cpp
    src/PostProcessing/VTK/HDF5CellData.cpp
    src/PostProcessing/VTK/HDF5BoundaryData.cpp
    src/PostProcessing/PostProcess.cpp
    src/PostProcessing/DerivedFields.cpp
    src/PostProcessing/Forces.cpp
    src/Case/CaseReader.cpp
    src/Case/CaseConfiguration.cpp
    src/Case/RuntimeSelection.cpp
    src/Application/SolverSetup.cpp
    src/Application/CFDApplication.cpp
)

# Include directories
target_include_directories(turblyze_core PUBLIC
    src/Mesh/
    src/Primitives/
    src/Fields/
    src/BoundaryConditions/
    src/BoundaryConditions/BoundaryTypes/
    src/Schemes/
    src/Schemes/Interpolation/
    src/Schemes/GradientSchemes/
    src/Schemes/ConvectionSchemes/
    src/Schemes/TimeSchemes/
    src/Models/
    src/Models/Turbulence/
    src/Models/Turbulence/RANS/
    src/Models/Turbulence/Laminar/
    src/Case/
    src/Application/
    src/PostProcessing/
    src/PostProcessing/VTK/
    src/Solver/
    src/LinearSystem/
    src/Parallel/
)

# Include METIS headers as SYSTEM to suppress third-party compiler warnings
target_include_directories(turblyze_core SYSTEM PUBLIC
    ${METIS_INCLUDE_DIR}
)

turblyze_apply_build_flags(turblyze_core)

# Define TURBLYZE_DOUBLE_PRECISION macro if requested
if(TURBLYZE_DOUBLE_PRECISION)
    target_compile_definitions(turblyze_core PUBLIC TURBLYZE_DOUBLE_PRECISION)
endif()

# Link the solver core against HDF5, PETSc, MPI, and METIS
target_link_libraries(turblyze_core PUBLIC
    HDF5::HDF5
    PkgConfig::PETSC
    MPI::MPI_CXX
    ${METIS_LIBRARY}
)

# main()
add_executable(${PROJECT_NAME} src/main.cpp)

# Link the solver executable against the core library
target_link_libraries(${PROJECT_NAME} PRIVATE turblyze_core)

# Apply compiler warning and Release optimization flags to the executable
turblyze_apply_build_flags(${PROJECT_NAME})

# Unit tests (Catch2 + CTest)
include(CTest)

# Build Catch2 unit, MPI, and regression test suites when enabled
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
