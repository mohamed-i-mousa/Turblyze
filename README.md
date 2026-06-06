# Turblyze - 3D Incompressible CFD Solver

A 3D incompressible CFD solver implementing the SIMPLE algorithm with k-omega SST turbulence modeling. The solver reads Fluent `.msh` meshes, solves steady-state incompressible flow, and exports results to VTK format for visualization in ParaView.

## Features

### Core Capabilities
- **3D Incompressible Flow**: Solves momentum equations (Ux, Uy, Uz) with the pressure correction via the SIMPLE algorithm

- **Collocated Grid**: Uses the Rhie-Chow Face-velocity  interpolation to prevent pressure checkerboarding

- **Multiple Convection Schemes**: Upwind (UDS), Second-Order Upwind (SOU), and Central-Difference (CDS) convection schemes with deferred-correction approach to improve stability

- **Gradient Reconstruction**: Weighted least-squares cell-centered gradients

- **Boundary Conditions**: flexible per-field BC system with direct face-to-patch linking and per-component velocity handling; unrecognized BC types are rejected with a fatal error listing the valid types

- **Turbulence Modeling**: k-omega SST model with wall distance calculation and wall functions

- **Wall Distance Calculation**: Mesh wave iterative propagation for accurate turbulence modeling

- **Shared-Memory Parallelism (OpenMP)**: All major compute loops (matrix assembly, gradient reconstruction, cell-update sweeps, turbulence transport) use OpenMP. Eigen's `BiCGSTAB` and `ConjugateGradient` solvers also use OpenMP via a RowMajor sparse-matrix layout. Thread count is set via `parallelism.numThreads` in the case file.

- **VTK Export**: Comprehensive output including all flow variables and turbulence quantities

- **Aerodynamic Force Reporting**: Optional post-solve integration of pressure and skin-friction loads over a named wall patch, decomposed into drag/lift along user-supplied directions, with non-dimensional `Cd`/`Cl` coefficients; printed to the console and written to a `<name>_forces.txt` file

- **Precision Control**: Configurable single (float) or double precision arithmetic

- **Documentation**: Full Doxygen-style code documentation


## Prerequisites

### System Requirements
- **C++20** compatible compiler (GCC 11+, Clang 15+, MSVC 2022+)
- **CMake** 3.20 or later
- **Linux/Unix** environment

### Dependencies
- **Eigen 3**: Linear algebra library (header-only)
- **VTK 9**: Visualization toolkit

#### Installation on Ubuntu/Debian:
```bash
sudo apt install build-essential cmake libeigen3-dev libvtk9-dev
```

#### Installation on MacOS:
```bash
brew install cmake eigen vtk libomp
```

> **Note (OpenMP setup)**: The `CMakeLists.txt` is tailored to two configurations
> out of the box:
> - **Linux with GCC/Clang** — OpenMP ships with the compiler; nothing to do.
> - **Apple Silicon macOS with AppleClang + Homebrew `libomp`** — the build
>   wires in `-Xpreprocessor -fopenmp` and the libomp path
>   `/opt/homebrew/opt/libomp` (hardcoded). Run `brew install libomp` first.
>
> If you are using a different setup (Intel Mac, Homebrew GCC, Homebrew
> LLVM/Clang, Intel oneAPI, MSVC, cross-compiler, custom libomp install, …),
> remove or adapt the `if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")` block
> in `CMakeLists.txt` so that the stock `find_package(OpenMP REQUIRED)` can
> discover your compiler's native OpenMP runtime. On Intel Macs the only
> change needed is updating `LIBOMP_PREFIX` to `/usr/local/opt/libomp`.

#### VTK Configuration Issues:
If CMake cannot locate VTK automatically:
```bash
cmake -S . -B build.nosync -DVTK_DIR=/usr/lib/cmake/vtk-9.1
```

## Building the Solver

### Standard Build Process
```bash
cmake -S . -B build.nosync -DCMAKE_BUILD_TYPE=Release
cmake --build build.nosync
```

### Build Types
```bash
# Release build
cmake -S . -B build.nosync -DCMAKE_BUILD_TYPE=Release
cmake --build build.nosync

# Debug build
cmake -S . -B build.nosync -DCMAKE_BUILD_TYPE=Debug
cmake --build build.nosync
```

### Build Options
| Option                          | Default | Effect                                              |
|---------------------------------|---------|-----------------------------------------------------|
| `TURBLYZE_USE_DOUBLE_PRECISION` | `ON`    | Defines `PROJECT_USE_DOUBLE_PRECISION` (double `Scalar`); turn `OFF` for single precision. |
| `TURBLYZE_NATIVE_ARCH`          | `ON`    | Adds `-march=native` in Release. Turn `OFF` for portable / cluster-deployable binaries. |

```bash
# Single-precision Release build
cmake -S . -B build.nosync -DCMAKE_BUILD_TYPE=Release \
    -DTURBLYZE_USE_DOUBLE_PRECISION=OFF

# Portable (no -march=native) Release build
cmake -S . -B build.nosync -DCMAKE_BUILD_TYPE=Release \
    -DTURBLYZE_NATIVE_ARCH=OFF
```

The executable `Turblyze` will be generated in the `build.nosync/` directory.

## Running Simulations

### Basic Execution
Run from the `build.nosync/` directory to ensure correct path resolution:
```bash
cd build.nosync
./Turblyze                     # Uses default defaultCase file
./Turblyze customCase          # Uses custom case file
```

### Case File System
The solver uses a case file system (default file: `defaultCase`). This allows runtime configuration without recompilation.

### Default Case
The default `defaultCase` file contains:
- **Mesh**: `../inputFiles/sphere.msh` (sphere in a channel)
- **Boundary Conditions**:
  - Inlet: Fixed velocity (0, 0, -20.0) m/s, zero gradient pressure
  - Outlet: Zero gradient velocity, fixed pressure (0 Pa)
  - Walls (`sphere`, `wall1`–`wall4`): No-slip velocity, zero gradient pressure; `kWallFunction`, `omegaWallFunction`, `nutWallFunction` for turbulence
- **Discretization**: Upwind convection scheme for momentum and turbulence equations. Least-squares for gradients computation
- **SIMPLE Parameters**: αU = 0.7, αp = 0.3, αk = 0.5, αω = 0.5, tolerance = 1e-3 (scaled residuals), max iterations = 500
- **Turbulence**: Enabled by default with k-omega SST model
- **Output**: `../outputFiles.nosync/sphere.vtu`

### Flow Physics
- **Fluid Properties**: Air (ρ = 1.225 kg/m³, μ = 1.7894e-5 Pa·s)
- **Flow Type**: Flow over a sphere at 20 m/s
- **Turbulence Inlet Conditions**: Turbulence intensity 5%, hydraulic diameter 0.01 m

## Input/Output

### Mesh Requirements
- **Format**: Fluent `.msh` files (ASCII format)
- **Dimension**: 3D only (2D meshes are not supported)
- **Cell Types**: Tetrahedra, hexahedra, prisms, pyramids
- **Boundary Patches**: Named patches for boundary condition assignment

### Output Visualization
- **Format**: VTK UnstructuredGrid (`.vtu`) and boundary PolyData
  (`_boundary.vtp`) for ParaView
- **Fields Exported**:
  - Main `.vtu`: `pressure`, `velocityMagnitude`, vector `velocity`,
    and, when turbulence is enabled, `k`, `omega`, `nut`, `wallDistance`
  - Boundary `.vtp` (e.g. `sphere_boundary.vtp`): all boundary patches with
    integer `patchID`, `patchZoneID`, `patchTypeID`, and `isWall` metadata.
    `wallShearStress` is included for all runs; `yPlus` is added only when
    turbulence is enabled
- **Cell Encoding**: volume cells are written as `VTK_POLYHEDRON` to preserve
  Turblyze's face topology. This is more robust for mixed/polyhedral meshes,
  but files can be larger and some ParaView filters may run slower than with
  native tetra/hex/wedge/pyramid cells.

### ParaView Visualization
1. Open the `.vtu` file in ParaView
2. Apply the file and click the "eye" icon to make it visible
3. Color by desired field (e.g., `pressure`, `velocityMagnitude`)
4. Open the corresponding `_boundary.vtp` file to inspect boundary patches
   or wall quantities (`wallShearStress` always, `yPlus` when turbulent)
5. Note: volume fields are cell-centered data; boundary metadata and wall
   diagnostics are boundary-face data

## Case Configuration

### Case File Format
The solver uses case files for configuration. The default `defaultCase` file contains all simulation parameters organized in sections. See `docs/CASE.md` for more details

```cpp
// Example case file entries
mesh
{
    file            ../inputFiles/your_mesh.msh;
    checkQuality    true;
}

physicalProperties
{
    rho             1.225;
    mu              1.7894e-5;
}

SIMPLE
{
    numIterations           300;
    convergenceTolerance    1e-3;       // Scaled residuals
    relaxationFactors
    {
        U                   0.7;
        p                   0.3;
        k                   0.5;
        omega               0.5;
    }
}

numericalSchemes
{
    convection
    {
        default     Upwind;
        U           SecondOrderUpwind;
        k           Upwind;
        omega       Upwind;
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
   ./Turblyze myCase
   ```


### Precision Control
Switch between single and double precision:

**Double Precision (default)**:
```bash
cmake -S . -B build.nosync -DTURBLYZE_USE_DOUBLE_PRECISION=ON
```

**Single Precision**:
```bash
cmake -S . -B build.nosync -DTURBLYZE_USE_DOUBLE_PRECISION=OFF
```

The solver prints precision mode at runtime via `SCALAR_MODE`.

## Project Structure

Headers and implementations live together under `src/`, following the
OpenFOAM convention:

- **`src/Application/`**: Top-level orchestration and solver assembly
  (`CFDApplication.h/.cpp`, `SolverSetup.h/.cpp`)
- **`src/Primitives/`**: Core scalar/vector/tensor types, logging, errors
- **`src/Mesh/`**: Geometric entities, topology, mesh I/O, mesh checking
- **`src/Fields/`**: Typed cell and face field containers
- **`src/BoundaryConditions/`**: Patch-based boundary condition storage,
  evaluation, and case loading
- **`src/Schemes/`**: Gradient, interpolation, and convection schemes
- **`src/LinearSystem/`**: Matrix assembly, transport equations, linear solvers
- **`src/Solver/`**: SIMPLE algorithm and solution constraints
- **`src/Models/`**: Physical models
  - **`src/Models/Turbulence/`**: Turbulence modeling
    (`kOmegaSST.h/.cpp`)
- **`src/PostProcessing/`**: Derived fields and output orchestration
  (`PostProcess.h/.cpp`)
  - **`src/PostProcessing/VTK/`**: VTK/PVD volume and boundary writers
- **`src/Case/`**: Case file system
  (`CaseReader.h/.cpp`, `CaseConfiguration.h/.cpp`)

## Documentation

Browsable HTML API documentation is generated from the Doxygen comments in the source via:

```bash
doxygen Doxyfile
```

Output is written to `docs/doxygen/html/`. Open `docs/doxygen/html/index.html` in a browser to navigate classes, call graphs, and collaboration diagrams. The `docs/doxygen/` tree is generated and is not tracked in git — regenerate it locally after pulling changes.

Requires `doxygen` and (for diagrams) `graphviz`:
```bash
# macOS
brew install doxygen graphviz

# Ubuntu/Debian
sudo apt install doxygen graphviz
```

## Numerical Method Validation

There is no committed automated test suite. Changes are validated by a
successful build and a representative case run. For numerics or solver-behavior
changes, output is compared against the reference OpenFOAM case under
`verification/`. See `docs/DEVELOPER_GUIDE.md` (§ "Testing and Debugging") for
the validation workflow and the component-level checks (boundary conditions,
convection coefficients, gradient conditioning).

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

No license has been specified for this project yet. For questions or bug
reports, contact the author or open an issue on the project's issue tracker.
