# MyCFDCode - 3D Incompressible CFD Solver

A comprehensive 3D incompressible CFD solver implementing the SIMPLE algorithm with k-omega SST turbulence modeling. The solver reads Fluent `.msh` meshes, solves steady-state incompressible flow, and exports results to VTK format for visualization in ParaView.

## Features

### Core Capabilities
- **3D Incompressible Flow**: Solves momentum equations (Ux, Uy, Uz) with pressure coupling via SIMPLE algorithm
- **Rhie-Chow Interpolation**: Face-velocity interpolation to prevent pressure checkerboarding on collocated grids
- **Multiple Discretization Schemes**: Upwind (UDS), second-order upwind (SOU), and central-difference (CDS) convection schemes
- **Gradient Reconstruction**: Least-squares cell-centered gradients with non-orthogonal diffusion treatment
- **Boundary Conditions**: Fixed value, fixed gradient, zero gradient, and no-slip wall conditions
- **Turbulence Modeling**: Optional k-omega SST model with wall distance calculation and wall functions
- **Precision Control**: Configurable single (float) or double precision arithmetic

### Advanced Features
- **Rhie-Chow Face Velocities**: Prevents pressure-velocity decoupling with under-relaxation effects
- **Wall Distance Calculation**: Poisson equation-based wall distance for accurate turbulence modeling  
- **VTK Export**: Comprehensive output including all flow variables and turbulence quantities
- **Mass Conservation**: Built-in mass conservation checking and diagnostics
- **Comprehensive Documentation**: Full Doxygen-style code documentation

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
# Debug build (default)
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Optimized release build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

The executable `MyCFDCode` will be generated in the `build/` directory.

## Running Simulations

### Basic Execution
Run from the `build/` directory to ensure correct path resolution:
```bash
cd build
./MyCFDCode
```

### Current Default Configuration (src/main.cpp:103)
- **Mesh**: `../inputFiles/sphere_24k.msh` (24k cell sphere mesh)
- **Boundary Conditions**: 
  - Inlet: Fixed velocity (-0.1 m/s in z-direction), zero gradient pressure
  - Outlet: Fixed pressure (0 Pa), zero gradient velocity
  - Sphere: No-slip velocity, zero gradient pressure
  - Walls: No-slip velocity, zero gradient pressure
- **Discretization**: Upwind convection scheme, least-squares gradients
- **SIMPLE Parameters**: αU = 0.3, αp = 0.1, tolerance = 1e-6, max iterations = 10
- **Turbulence**: Disabled by default (line 207)
- **Output**: `../outputFiles/sphere_24k.vtp`

### Flow Physics Setup
- **Fluid Properties**: Air (ρ = 1.225 kg/m³, μ = 1.7894e-5 Pa·s)
- **Reynolds Number**: Low-Re flow for numerical stability
- **Flow Type**: External flow around sphere

## Input/Output

### Mesh Requirements
- **Format**: Fluent `.msh` files (ASCII format)
- **Dimension**: 3D only (2D meshes are rejected)
- **Cell Types**: Tetrahedra, hexahedra, prisms, pyramids
- **Boundary Patches**: Named patches for boundary condition assignment

### Available Test Meshes
Located in `inputFiles/`:
- `cylinder.msh`, `cylinder_228k.msh`, `cylinder_57k.msh`, `cylinder_coarse.msh`
- `sphere.msh`, `sphere_24k.msh`, `sphere_74k.msh`

### Output Visualization
- **Format**: VTK PolyData (`.vtp`) for ParaView
- **Fields Exported**:
  - Flow variables: `pressure`, `velocityMagnitude`, `U_x`, `U_y`, `U_z`
  - Turbulence (if enabled): `k`, `omega`, `mu_t`, `wallDistance`
  - Derived quantities: `turbulentIntensity`, `turbulentViscosityRatio`, `yPlus`

### ParaView Visualization
1. Open the `.vtp` file in ParaView
2. Apply the file and click the "eye" icon to make it visible
3. Color by desired field (e.g., `velocityMagnitude`, `pressure`)
4. Note: Fields are face-centered data (mesh faces become ParaView cells)

## Configuration and Customization

### Modifying Simulation Parameters
Edit `src/main.cpp` to customize:

```cpp
// Mesh selection (line 103)
std::string meshFilePath = "../inputFiles/your_mesh.msh";

// Physical properties (lines 148-149) 
const Scalar rho = 1.225;      // Density [kg/m³]
const Scalar mu = 1.7894e-5;   // Viscosity [Pa·s]

// Boundary conditions (lines 161-180)
bcManager.setFixedValue("inlet", U_field, Vector(1.0, 0.0, 0.0));  // 1 m/s in x-direction

// SIMPLE parameters (lines 202-207)
simpleSolver.setRelaxationFactors(0.7, 0.3);   // More aggressive relaxation
simpleSolver.setMaxIterations(100);            // More iterations
simpleSolver.enableTurbulenceModeling(true);   // Enable k-omega SST

// Discretization scheme (line 196)
SIMPLE simpleSolver(allFaces, allCells, bcManager, gradScheme, CDS);  // Central difference
```

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
- **`Core/`**: `Scalar.h`, `Vector.h`, utility functions (`linearInterpolation.h`, `massFlowRate.h`)
- **`Mesh/`**: Geometric entities (`Face.h`, `Cell.h`), data containers (`CellData.h`, `FaceData.h`), I/O (`MeshReader.h`)
- **`BoundaryConditions/`**: BC system (`BoundaryConditions.h`, `BoundaryData.h`, `BoundaryPatch.h`)
- **`Numerics/`**: Schemes (`GradientScheme.h`, `ConvectionScheme.h`), matrix assembly (`Matrix.h`), solvers (`LinearSolvers.h`, `SIMPLE.h`)
- **`Models/`**: Turbulence modeling (`KOmegaSST.h`)
- **`PostProcessing/`**: Output (`VtkWriter.h`)

### Source Organization (`src/`)
Mirrors header organization with corresponding `.cpp` implementations.

### Key Design Patterns
1. **Field-based Architecture**: Type-safe field containers (`CellData<T>`, `FaceData<T>`)
2. **Geometry-first Approach**: Mesh entities compute their own geometric properties  
3. **Boundary Condition Abstraction**: Type-safe BC storage with field-agnostic evaluation
4. **Matrix-based Assembly**: Centralized matrix construction for all transport equations

## Solution Algorithm

### SIMPLE Algorithm Flow (SIMPLE.cpp:solve())
1. **Cache Refresh**: Update gradients and mass fluxes for current iteration
2. **Momentum Solution**: Solve Ux, Uy, Uz with effective viscosity (μ + μt)
3. **Rhie-Chow Interpolation**: Compute face velocities and mass fluxes
4. **Pressure Correction**: Assemble and solve pressure correction equation
5. **Velocity Correction**: Update velocity field using pressure correction
6. **Field Updates**: Correct mass fluxes and pressure with under-relaxation
7. **Turbulence**: Advance k-omega SST model (if enabled)
8. **Convergence Check**: Monitor mass, velocity, and pressure residuals

### Convergence Criteria
- **Mass Imbalance**: RMS of mass conservation residuals
- **Velocity Residual**: RMS change in velocity field
- **Pressure Residual**: RMS of pressure correction field

## Troubleshooting

### Build Issues
- **VTK not found**: Specify `VTK_DIR` path manually
- **Eigen not found**: Install `libeigen3-dev` or set include path
- **C++17 errors**: Ensure compiler supports C++17 standard

### Runtime Issues  
- **2D mesh error**: Only 3D meshes supported - check mesh export settings
- **File not found**: Run from `build/` directory for correct relative paths
- **Memory issues**: Large meshes may require significant RAM

### Visualization Issues
- **Empty ParaView display**: Check that fields are applied and visible
- **Incorrect values**: Remember that fields are face-centered in VTK output
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



