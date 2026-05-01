cmake_minimum_required(VERSION 3.20)
project(Turblyze VERSION "1.0.0" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Build options
option(TURBLYZE_USE_DOUBLE_PRECISION
       "Build with double-precision Scalar (otherwise float)" ON)
option(TURBLYZE_NATIVE_ARCH
       "Optimize for the build host CPU (-march=native)" ON)

# Find Eigen3 package (required for linear algebra)
find_package(Eigen3 REQUIRED NO_MODULE)

# find OpenMP package (required for parallelization)
if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(LIBOMP_PREFIX "/opt/homebrew/opt/libomp")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
endif()

find_package(OpenMP REQUIRED)

# Find VTK package with required components
find_package(VTK 9 REQUIRED COMPONENTS
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
target_include_directories(${PROJECT_NAME} PRIVATE
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

# Warning flags — GCC/Clang/AppleClang only
target_compile_options(${PROJECT_NAME} PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:
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
    >
)

if(TURBLYZE_USE_DOUBLE_PRECISION)
    target_compile_definitions(${PROJECT_NAME} PRIVATE PROJECT_USE_DOUBLE_PRECISION)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(${PROJECT_NAME} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-O3>
    )
    if(TURBLYZE_NATIVE_ARCH)
        target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-march=native>
        )
    endif()
endif()

# Link executable against the VTK, Eigen3, and OpenMP libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::IOLegacy
    VTK::IOXML
    Eigen3::Eigen
    OpenMP::OpenMP_CXX
)

# VTK 9 needs autoinit to register module factories at link time.
vtk_module_autoinit(
    TARGETS ${PROJECT_NAME}
    MODULES VTK::CommonCore VTK::CommonDataModel VTK::IOLegacy VTK::IOXML
)