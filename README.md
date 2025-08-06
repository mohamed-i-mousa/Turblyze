# 3D CFD Solver with SIMPLE Algorithm

A comprehensive 3D Computational Fluid Dynamics (CFD) solver implementing the SIMPLE (Semi-Implicit Method for Pressure Linked Equations) algorithm with support for both laminar and turbulent flows using the k-omega SST turbulence model.

## Features

### 3D Flow Solver
- **Complete 3D Implementation**: Solves all three momentum equations (U, V, W components)
- **SIMPLE Algorithm**: Pressure-velocity coupling for incompressible flows
- **Rhie-Chow Interpolation**: Prevents checkerboard pressure fields
- **Under-relaxation**: Stable convergence for complex flows

### Turbulence Modeling
- **k-omega SST Model**: Shear Stress Transport turbulence model
- **Automatic Detection**: 3D mesh analysis and turbulence parameter calculation
- **Wall Functions**: Proper near-wall treatment

### Discretization Schemes
- **Convection**: Central difference, upwind, and hybrid schemes
- **Diffusion**: Central difference with non-orthogonal correction
- **Gradients**: Least squares method for unstructured meshes

### Boundary Conditions
- **Velocity**: Fixed value, fixed gradient, no-slip wall
- **Pressure**: Fixed value, fixed gradient, zero gradient
- **Turbulence**: Automatic wall treatment, inlet/outlet conditions

## 3D Implementation Details

### Momentum Equations
The solver now properly handles all three velocity components:

```cpp
// U-momentum (x-component)
∂(ρU)/∂t + ∇·(ρUU) = -∂p/∂x + ∇·[(μ + μₜ)∇U]

// V-momentum (y-component)  
∂(ρV)/∂t + ∇·(ρVV) = -∂p/∂y + ∇·[(μ + μₜ)∇V]

// W-momentum (z-component)
∂(ρW)/∂t + ∇·(ρWW) = -∂p/∂z + ∇·[(μ + μₜ)∇W]
```

### 3D Velocity Components
- **getVelocityX()**: Returns U_x component field
- **getVelocityY()**: Returns U_y component field  
- **getVelocityZ()**: Returns U_z component field
- **getVelocity()**: Returns complete 3D velocity vector field

### Mesh Analysis
The solver automatically detects 3D mesh characteristics:
- Analyzes z-range vs. cell size
- Warns if mesh appears 2D
- Provides 3D solution statistics

## Usage

### Basic Setup
```cpp
#include "SIMPLE.h"

// Create solver
SIMPLE simpleSolver(faces, cells, boundaryConditions, gradientScheme, convectionScheme);

// Configure parameters
simpleSolver.setRelaxationFactors(0.7, 0.3);  // U=0.7, p=0.3
simpleSolver.setConvergenceTolerance(1e-6);
simpleSolver.setMaxIterations(500);

// Enable turbulence modeling
simpleSolver.enableTurbulenceModeling(true);

// Solve
simpleSolver.solve();
```

### Extract 3D Results
```cpp
// Get complete velocity field
const VectorField& velocity = simpleSolver.getVelocity();

// Get individual velocity components
ScalarField U_x = simpleSolver.getVelocityX();
ScalarField U_y = simpleSolver.getVelocityY();
ScalarField U_z = simpleSolver.getVelocityZ();

// Get pressure and other fields
const ScalarField& pressure = simpleSolver.getPressure();
const FaceFluxField& massFlux = simpleSolver.getMassFlux();
```

## Output

### 3D Solution Statistics
The solver provides comprehensive 3D statistics:
- Maximum and average velocity magnitudes
- Individual velocity component statistics (U_x, U_y, U_z)
- Pressure range and distribution
- Turbulence parameters (k, ω, μₜ)
- Mass conservation verification

### VTK Export
All solution fields are exported to VTK format for visualization:
- Velocity components (U_x, U_y, U_z)
- Pressure and pressure coefficient
- Turbulence fields (k, ω, μₜ, wall distance)
- Derived quantities (turbulent intensity, y+)

## Building

```bash
mkdir build
cd build
cmake ..
make
```

## Running

```bash
./cfd_solver
```

The solver will:
1. Read 3D mesh from input files
2. Detect mesh dimensionality
3. Solve 3D flow with SIMPLE algorithm
4. Apply turbulence modeling if enabled
5. Export results to VTK format
6. Provide comprehensive 3D statistics

## Mesh Requirements

- **3D Unstructured Mesh**: Tetrahedral, hexahedral, or mixed elements
- **Boundary Patches**: Properly defined inlet, outlet, and wall boundaries
- **Mesh Quality**: Good aspect ratios and non-degenerate elements

## Performance

The 3D implementation includes:
- Efficient sparse matrix solvers (BiCGSTAB)
- Optimized memory usage for large 3D meshes
- Parallel-ready data structures
- Convergence monitoring and diagnostics

## Validation

The solver has been validated for:
- 3D laminar flows
- 3D turbulent flows with k-omega SST
- Mass conservation
- Energy conservation
- Boundary condition implementation

## Future Enhancements
- Advanced no-slip boundary condition 
- Parallel processing for large 3D meshes
- Additional turbulence models (LES)

