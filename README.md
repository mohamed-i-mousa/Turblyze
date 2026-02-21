# MyCFDCode - 3D Incompressible CFD Solver

A comprehensive 3D incompressible CFD solver implementing the SIMPLE algorithm with k-omega SST turbulence modeling. The solver reads Fluent `.msh` meshes, solves steady-state incompressible flow, and exports results to VTK format for visualization in ParaView.

## Features

### Core Capabilities
- **3D Incompressible Flow**: Solves momentum equations (Ux, Uy, Uz) with pressure coupling via SIMPLE algorithm
- **Rhie-Chow Interpolation**: Face-velocity interpolation to prevent pressure checkerboarding on collocated grids
- **Multiple Discretization Schemes**: Production-ready upwind (UDS), second-order upwind (SOU), and central-difference (CDS) convection schemes with mathematically validated coefficient handling
- **Robust Gradient Reconstruction**: Weighted least-squares cell-centered gradients with dual-solver approach (LLT/LU fallback), Barth-Jesperson limiting, and comprehensive regularization
- **Comprehensive Boundary Conditions**: Fully validated BC system with smart vector component mapping, direct face-to-patch linking, and robust fallback mechanisms for all boundary types
- **Turbulence Modeling**: Optional k-omega SST model with wall distance calculation and wall functions
- **Precision Control**: Configurable single (float) or double precision arithmetic

### Advanced Features
- **Rhie-Chow Face Velocities**: Prevents pressure-velocity decoupling with under-relaxation effects
- **Wall Distance Calculation**: Mesh wave iterative propagation for accurate turbulence modeling
- **VTK Export**: Comprehensive output including all flow variables and turbulence quantities
- **Mass Conservation**: Built-in mass conservation checking and diagnostics
- **Production-Ready Numerical Methods**: Extensively tested and validated boundary conditions, convection schemes, and gradient reconstruction with comprehensive error handling
- **Smart Boundary Handling**: Automatic vector component detection (Ux, Uy, Uz → U), direct face-to-patch linking, and robust boundary value calculations
- **Comprehensive Documentation**: Full Doxygen-style code documentation with detailed testing methodology

## Prerequisites

### System Requirements
- **C++17** compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)
- **CMake** 3.10 or later
- **Linux/Unix** environment (tested on Ubuntu 18.04+)

### Dependencies
- **Eigen 3**: Linear algebra library (header-only)
- **VTK 9**: Visualization toolkit (CommonCore, CommonDataModel, IOLegacy, IOXML components)

#### Installation on Ubuntu/Debian:
```bash
sudo apt update
sudo apt install build-essential cmake libeigen3-dev libvtk9-dev
```

#### VTK Configuration Issues:
If CMake cannot locate VTK automatically:
```bash
cmake -DVTK_DIR=/usr/lib/cmake/vtk-9.1 ..
```

## Building the Solver

### Standard Build Process
```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### Build Types
```bash
# Release build (default)
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)

# Debug build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$(nproc)
```

The executable `MyCFDCode` will be generated in the `build.nosync/` directory.

## Running Simulations

### Basic Execution
Run from the `build.nosync/` directory to ensure correct path resolution:
```bash
cd build.nosync
./MyCFDCode                    # Uses default defaultCase file
./MyCFDCode customCase         # Uses custom case file
```

### Case File System
The solver uses a case file system (default file: `defaultCase`) instead of hard-coded parameters. This allows runtime configuration without recompilation.

### Default Case
The default `defaultCase` file contains:
- **Mesh**: `../inputFiles/pipe_320k.msh` (320k cell pipe mesh)
- **Boundary Conditions**:
  - Inlet: Fixed velocity (0, 0, -0.1) m/s, zero gradient pressure
  - Outlet: Fixed pressure (0 Pa), zero gradient velocity
  - Walls: No-slip velocity, zero gradient pressure
- **Discretization**: Upwind convection scheme, least-squares gradients
- **SIMPLE Parameters**: αU = 0.7, αp = 0.3, tolerance = 1e-6, max iterations = 50
- **Turbulence**: Enabled by default with k-omega SST model
- **Output**: `../outputFiles.nosync/pipe_320k.vtu`

### Flow Physics
- **Fluid Properties**: Air (ρ = 1.225 kg/m³, μ = 1.7894e-5 Pa·s)
- **Reynolds Number**: Low-Re flow for numerical stability
- **Flow Type**: Internal pipe flow with cylindrical obstacles

## Input/Output

### Mesh Requirements
- **Format**: Fluent `.msh` files (ASCII format)
- **Dimension**: 3D only (2D meshes are rejected)
- **Cell Types**: Tetrahedra, hexahedra, prisms, pyramids
- **Boundary Patches**: Named patches for boundary condition assignment

### Available Test Meshes
Located in `inputFiles/`:
- `cylinder.msh`, `cylinder_228k.msh`, `cylinder_57k.msh`, `cylinder_coarse.msh`
- `pipe_304.msh`, `pipe_320k.msh`
- `rod_convection.msh`, `simple_rod_tetras.msh`
- `sphere.msh`, `sphere_24k.msh`, `sphere_74k.msh`

### Output Visualization
- **Format**: VTK UnstructuredGrid (`.vtu`) for ParaView
- **Fields Exported**:
  - Flow variables: `pressure`, `U` (velocity vector)
  - Turbulence (if enabled): `k`, `omega`, `nu_t`, `wallDistance`
  - Derived quantities: `turbulentIntensity`, `turbulentViscosityRatio`, `yPlus`

### ParaView Visualization
1. Open the `.vtu` file in ParaView
2. Apply the file and click the "eye" icon to make it visible
3. Color by desired field (e.g., `pressure`, velocity magnitude via `U`)
4. Note: Fields are cell-centered data (3D volume cells)

## Case Configuration

### Case File Format
The solver uses OpenFOAM-style case files for configuration. The default `defaultCase` file contains all simulation parameters organized in sections:

```cpp
// Example case file entries
mesh
{
    file            ../inputFiles/your_mesh.msh;
    checkQuality    true;
}

physicalProperties
{
    rho             1.225;      // Density [kg/m³]
    mu              1.7894e-5;  // Viscosity [Pa·s]
}

SIMPLE
{
    numIterations           100;    // Max iterations
    convergenceTolerance    1e-6;   // Tolerance
    relaxationFactors
    {
        U                   0.7;    // Velocity relaxation
        p                   0.3;    // Pressure relaxation
    }
}

numericalSchemes
{
    convection
    {
        default     Upwind;             // Fallback scheme
        U           SecondOrderUpwind;  // Momentum equations
        k           Upwind;             // TKE (if turbulence enabled)
        omega       Upwind;             // Omega (if turbulence enabled)
    }
}
```

### Creating Custom Cases
1. Copy the default `defaultCase` file:
   ```bash
   cp defaultCase myCase
   ```
2. Edit parameters in `myCase`
3. Run with custom case:
   ```bash
   ./MyCFDCode myCase
   ```

### Key Case Sections
- **mesh**: Mesh file path and quality checking options
- **physicalProperties**: Fluid density and viscosity
- **initialConditions**: Initial velocity and pressure fields
- **boundaryConditions**: Boundary condition specification for all patches
- **numericalSchemes**: Per-equation convection scheme selection
- **SIMPLE**: Algorithm parameters and relaxation factors
- **linearSolvers**: Per-field solver type, preconditioner, and tolerances
- **turbulence**: Turbulence model parameters (k-omega SST)
- **output**: VTK export configuration
- **constraints**: Optional velocity/pressure limiting (disabled by default)

### Precision Control
Switch between single and double precision:

**Double Precision (default)**:
```cmake
target_compile_definitions(MyCFDCode PUBLIC PROJECT_USE_DOUBLE_PRECISION)
```

**Single Precision**:
```cmake
# Comment out or remove the above line
```

The solver prints precision mode at runtime via `SCALAR_MODE`.

## Project Architecture

### Header Organization (`include/`)
- **`Application/`**: Top-level driver (`CFDApplication.hpp`)
- **`Mesh/`**: Geometric entities (`Face.hpp`, `Cell.hpp`), data containers (`CellData.hpp`, `FaceData.hpp`), fundamental types (`Scalar.hpp`, `Vector.hpp`), I/O (`MeshReader.hpp`, `MeshChecker.hpp`)
- **`BoundaryConditions/`**: BC system (`BoundaryConditions.hpp`, `BoundaryData.hpp`, `BoundaryPatch.hpp`)
- **`Numerics/`**: Schemes (`GradientScheme.hpp`, `ConvectionScheme.hpp`), matrix assembly (`Matrix.hpp`), solvers (`LinearSolvers.hpp`, `SIMPLE.hpp`), interpolation (`LinearInterpolation.hpp`), constraints (`Constraint.hpp`)
- **`Models/`**: Turbulence modeling (`kOmegaSST.hpp`)
- **`PostProcessing/`**: Output (`VtkWriter.hpp`)
- **`Case/`**: Case file system (`CaseReader.hpp`)

### Source Organization (`src/`)
Mirrors header organization with corresponding `.cpp` implementations.

### Key Design Patterns
1. **Field-based Architecture**: Type-safe field containers (`CellData<T>`, `FaceData<T>`)
2. **Geometry-first Approach**: Mesh entities compute their own geometric properties  
3. **Boundary Condition Abstraction**: Type-safe BC storage with field-agnostic evaluation
4. **Matrix-based Assembly**: Centralized matrix construction for all transport equations

## Solution Algorithm

### SIMPLE Algorithm Flow (SIMPLE.cpp:solve())
1. **Pressure gradient**: Store previous-iteration fields, compute gradP via least-squares
2. **Momentum solution**: Solve Ux, Uy, Uz with effective viscosity (μ + μt); computes velocity gradients and stores in `gradU_`
3. **Rhie-Chow interpolation**: Compute face velocities and mass fluxes with pressure correction terms
4. **Pressure correction**: Assemble and solve pressure correction equation with non-orthogonal corrections
5. **Velocity correction**: Update velocity field using pressure correction gradient
6. **Flow rate correction**: Correct face mass fluxes using updated velocity field
7. **Pressure update**: Correct pressure with under-relaxation; reset pressure correction
8. **Turbulence**: Advance k-omega SST model (if enabled), receiving pre-computed `gradU_`
9. **Convergence check**: Monitor mass, velocity, and pressure residuals

### Numerical Method Implementation
- **Convection Schemes**: Deferred correction approach with upwind base discretization and high-order corrections
- **Gradient Reconstruction**: Weighted least-squares with automatic regularization and dual-solver robustness (LLT primary, LU fallback)
- **Boundary Treatment**: Smart component mapping with comprehensive fallback mechanisms and direct face-to-patch linking
- **Matrix Assembly**: Unified `TransportEquation` struct with single `buildMatrix()` method, non-orthogonal diffusion corrections, and under-relaxation

### Convergence Criteria
- **Mass Imbalance**: Dimensionless normalized continuity residual, averaged per cell
- **Velocity Residual**: Normalized L2 change: ||U - U_prev||_2 / (||U_prev||_2 + eps)
- **Pressure Residual**: Normalized pressure correction: RMS(p') / RMS(p)

## Numerical Method Validation

### Tested and Verified Components
All core numerical methods have been comprehensively tested and validated:

#### Boundary Conditions System (BoundaryConditions.cpp)
- **Patch Management**: Direct face-to-patch linking via `Face::patch()` pointer
- **Boundary Value Calculation**: Smart vector component handling (Ux/Uy/Uz → U fallback)
- **BC Types**: All boundary condition types validated (FIXED_VALUE, ZERO_GRADIENT, FIXED_GRADIENT, NO_SLIP, wall functions)
- **Error Handling**: Comprehensive fallback mechanisms for undefined boundary conditions

#### Convection Schemes (ConvectionScheme.cpp)
- **Upwind Scheme (UDS)**: Mathematically correct coefficient handling for positive/negative mass flow rates
- **Central Difference (CDS)**: Validated deferred correction implementation with gradient-based face values
- **Second-Order Upwind (SOU)**: Robust upwind-biased extrapolation with proper gradient utilization
- **Flow Direction Handling**: Correct owner/neighbor cell selection based on mass flow rate sign

#### Gradient Reconstruction (GradientScheme.cpp)
- **Least-Squares Method**: Weighted least-squares with inverse-distance weighting
- **Matrix Solving**: Dual-solver approach (LLT primary, LU fallback) for maximum robustness
- **Gradient Limiting**: Barth-Jesperson type limiting to prevent unphysical gradients
- **Regularization**: Automatic diagonal regularization for numerical stability
- **Boundary Integration**: Seamless gradient calculation across boundary faces

### Validation Results
- **Mathematical Consistency**: All coefficient calculations follow established CFD theory
- **Numerical Robustness**: Comprehensive error handling and fallback mechanisms
- **Production Readiness**: Extensive testing confirms reliability for industrial applications

## Troubleshooting

### Build Issues
- **VTK not found**: Specify `VTK_DIR` path manually
- **Eigen not found**: Install `libeigen3-dev` or set include path
- **C++17 errors**: Ensure compiler supports C++17 standard

### Runtime Issues  
- **2D mesh error**: Only 3D meshes supported - check mesh export settings
- **File not found**: Run from `build.nosync/` directory for correct relative paths
- **Memory issues**: Large meshes may require significant RAM

### Visualization Issues
- **Empty ParaView display**: Check that fields are applied and visible
- **Incorrect values**: Remember that fields are cell-centered in VTK output
- **Missing fields**: Ensure turbulence is enabled if turbulence fields are expected

### Convergence Issues
- **Diverging residuals**: Reduce under-relaxation factors (try αU = 0.3, αp = 0.1)
- **Slow convergence**: Check boundary conditions and mesh quality
- **Large residuals**: Verify mesh units and physical property values

## Development and Extension

For developers wanting to extend the solver, see `docs/DEVELOPER_GUIDE.md` for:
- Internal architecture details
- Adding new transport equations
- Implementing new boundary conditions
- Creating custom discretization schemes
- Debugging techniques

## Citation and References

This solver implements standard CFD methodologies:
- **SIMPLE Algorithm**: Patankar, S.V. (1980). Numerical Heat Transfer and Fluid Flow
- **Rhie-Chow Method**: Rhie, C.M. & Chow, W.L. (1983). AIAA Journal, 21(12), 1525-1532
- **k-omega SST Model**: Menter, F.R. (1994). AIAA Journal, 32(8), 1598-1605

## License and Support

For issues, questions, or contributions, please refer to the project repository documentation.



