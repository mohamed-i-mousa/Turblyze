<!--
SPDX-FileCopyrightText: 2025-2026 Mohamed Mousa
SPDX-License-Identifier: Apache-2.0
-->

<h1 align="center">Turblyze</h1>
<p align="center">
  <a href="https://github.com/mohamed-i-mousa/Turblyze/actions/workflows/ci.yml"><img src="https://github.com/mohamed-i-mousa/Turblyze/actions/workflows/ci.yml/badge.svg" alt="CI Build & Test"></a>
  <a href="https://www.codefactor.io/repository/github/mohamed-i-mousa/turblyze"><img src="https://www.codefactor.io/repository/github/mohamed-i-mousa/turblyze/badge" alt="CodeFactor"></a>
  <img src="https://img.shields.io/badge/Platform-Linux%20%7C%20macOS%20%7C%20Windows%20(WSL2)-lightgrey.svg" alt="Platforms">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-Apache_2.0-blue.svg" alt="License: Apache 2.0"></a>
</p>

<p align="center">3D Incompressible CFD Solver parallelized with MPI and solves steady-state and transient flows. It offers several turbulence modeling options, reads unstructured Fluent `.msh` meshes, and export results in VTKHDF format (`.vtkhdf`) for visualization in ParaView. </p>



<p align="center">
  <a href="validation/cylinder/animation/cylinderVortexShedding.mp4">
    <img src="validation/cylinder/animation/cylinderVortexShedding.gif" alt="Turblyze: 2D Cylinder Vortex Shedding (Re ≈ 150)" width="100%">
  </a>
</p>
<p align="center">
  <em>Von Kármán vortex street past a circular cylinder (<a href="validation/cylinder/README.md">Re ≈ 150</a>) simulated with Turblyze</em>
</p>

## Features

### Core Capabilities
- **3D Incompressible Flow**: Solves momentum equations with pressure correction with the SIMPLE/PISO algorithm

- **Steady-state and Transient (URANS)**: Runs a steady-state simulation with the SIMPLE algorithm or a transient simulation with implicit Euler / Crank-Nicolson / second-order implicit time schemes using the PISO algorithm

- **Convection Schemes**: Upwind (UDS), Second-Order Upwind (SOU), Central-Difference (CDS), and LUST (Linear-Upwind Stabilized Transport) convection schemes with the deferred-correction approach

- **Gradient Reconstruction**: Weighted least-squares cell-centered gradients

- **Boundary Conditions**: `fixedValue`, `zeroGradient`, `fixedGradient`, `noSlip` (velocity), wall functions (`kWallFunction`, `omegaWallFunction`, `nutWallFunction`), and `symmetry` planes

- **Turbulence Modeling**: Laminar or k-omega SST model with wall functions

- **Distributed-Memory Parallelism (MPI)**: Runtime METIS domain decomposition with ghost-cell halo exchange. The same binary runs serially or under `mpirun -np N ./Turblyze case`, with no case-file changes. Linear systems are solved with PETSc Krylov solvers.

- **VTKHDF Export**: Comprehensive output including all flow variables and turbulence quantities, written directly through the HDF5 library in VTK's modern `.vtkhdf` format.

- **Aerodynamic Force Reporting**: Optional post-solve integration of pressure and skin-friction loads over a named wall patch, decomposed into drag/lift along user-supplied directions, with non-dimensional `Cd`/`Cl` coefficients; printed to the console and written to a `<name>_forces.txt` file

- **Precision Control**: Configurable single (float) or double precision arithmetic

- **Documentation**: Doxygen-style code documentation


## Prerequisites

### System Requirements
- **C++20**-compatible compiler (GCC 11+, or Clang/AppleClang 15+). The warning and optimization flags target GCC, Clang, and AppleClang
- **CMake** 3.20 or later
- **Linux**, **macOS**, or **Windows (via WSL2)** environment

### Dependencies
- **PETSc**: Krylov linear solvers, located through `pkg-config` (so
  `pkg-config` itself is required at configure time)
- **MPI**: any MPI implementation providing a C++-capable compiler wrapper
  (the `CXX` component); initialized through PETSc and called directly in
  `src/Parallel/`
- **METIS**: mesh partitioning for the runtime domain decomposition
- **Catch2 v3** (tests only): the CTest suite links Catch2; required only when
  `BUILD_TESTING` is on (the default). Not needed to build or run the solver —
  configure with `-DBUILD_TESTING=OFF` to drop the dependency. See [Testing](#testing).
- **HDF5 (parallel/MPI build)**: C library used to write the VTKHDF output.
  The VTKHDF writers do collective MPI-IO into one shared file per grid, so the
  build **strongly prefers a parallel (MPI-enabled) HDF5**. A serial-only HDF5
  satisfies the configure step, but the parallel writes then **fail or misbehave
  at runtime**, so install the MPI-enabled package, or point `HDF5_ROOT` /
  `HDF5_PREFER_PARALLEL` to a parallel installation.

#### Installation on Ubuntu/Debian:
```bash
sudo apt install build-essential cmake pkg-config libopenmpi-dev petsc-dev libmetis-dev libhdf5-openmpi-dev catch2
```

#### Installation on macOS:
```bash
brew install cmake pkg-config open-mpi petsc metis hdf5-mpi catch2
```

#### Installation on Windows 10/11 (WSL2):
Enable WSL and install Ubuntu from PowerShell:
```powershell
wsl --install -d Ubuntu
```
Then open the Ubuntu terminal and install the dependencies:
```bash
sudo apt install build-essential cmake pkg-config libopenmpi-dev petsc-dev libmetis-dev libhdf5-openmpi-dev catch2
```

(`catch2` is only needed to build the tests; omit it and configure with
`-DBUILD_TESTING=OFF` for a solver-only build. The suite requires Catch2 **v3**
via `find_package(Catch2 3)` — on distributions that still ship v2, install v3
separately or build with `-DBUILD_TESTING=OFF`.)

On macOS, the serial `hdf5` formula and `hdf5-mpi` conflict: if a serial `hdf5` is linked, unlink it (`brew unlink hdf5 && brew link hdf5-mpi`) before configuring, otherwise a serial HDF5 is picked up and the parallel writes fail at runtime.


## Building the Solver

```bash
cmake -S . -B build.nosync
cmake --build build.nosync -j
```

This produces the `Turblyze` executable in `build.nosync/`. The `-j` flag
builds in parallel across all available cores. The build type defaults to
`Release`; pass `-DCMAKE_BUILD_TYPE=Debug` at configure time for a debug
build.

### Build Options

Pass these at configure time, e.g.
`cmake -S . -B build.nosync -DTURBLYZE_DOUBLE_PRECISION=OFF`:

| Option                          | Default | Effect                                              |
|---------------------------------|---------|-----------------------------------------------------|
| `TURBLYZE_DOUBLE_PRECISION`     | `ON`    | Double-precision `Scalar`; turn `OFF` for single precision. The solver prints the active mode at runtime via `SCALAR_MODE`. |
| `TURBLYZE_NATIVE_ARCH`          | `ON`    | Adds `-march=native` in Release. Turn `OFF` for portable / cluster-deployable binaries. |
| `BUILD_TESTING`                 | `ON`    | Build the Catch2 test suite (needs Catch2 v3). Turn `OFF` for a solver-only build with no test dependency. |
| `TURBLYZE_REGRESSION_CD_PIN`    | Apple M4 value | Pinned sphere Cd for the `regression` test; machine-local under `-march=native`, re-pin per CPU. See [Testing](#testing). |
| `TURBLYZE_REGRESSION_CD_TOL`    | `1e-4`  | Relative tolerance for the pinned sphere Cd. |

## Testing

The test suite uses [Catch2 v3](https://github.com/catchorg/Catch2) and CTest;
it is built whenever `BUILD_TESTING` is on (the default). Build, then run:

```bash
cmake -S . -B build.nosync -DCMAKE_BUILD_TYPE=Release
cmake --build build.nosync -j
ctest --test-dir build.nosync --output-on-failure
```

Three binaries back four CTest labels:

- **`turblyze_unit_tests`** — serial unit tests (primitives, schemes, mesh
  geometry, boundary-condition linearizations and factory, least-squares
  gradients and the Barth-Jespersen limiter, case-file parsing, runtime
  selection), one CTest entry per test case, label `serial`. Never initializes
  MPI.
- **`turblyze_mpi_tests`** — MPI + PETSc tests launched under `mpirun` at 1, 2,
  and 4 ranks: global reductions and their rank-count invariance, global
  indexing, halo exchange on a decomposed fixture, matrix assembly /
  relaxation / constraint / Jacobi update, a small diffusion solve with both
  Krylov solvers, and a standalone k-omega SST construction and solve. Label
  `mpi`. Every case runs at all three rank counts, on a box decomposed through
  the production METIS path; only the halo-exchange cases `SKIP`, at np = 1,
  where there is no cut to exchange across.
- **`turblyze_abort_probe`** — a subprocess that deterministically triggers a
  `FatalError`; the `death`-labelled tests assert it prints `FATAL ERROR` and
  aborts (a control run exits cleanly).

A fourth label, `regression`, runs the solver itself: a fixed-iteration steady
sphere case at Re = 300 whose total drag coefficient is pinned (tolerance
`1e-4`). It takes ~35 s and is excluded from the fast gate. The pin is only
reproducible in the configuration it was measured in, so the test registers
itself only for a double-precision Release build with `TURBLYZE_NATIVE_ARCH`
on, and only when the (unshipped) `inputFiles/sphere.msh` is present. Because
`-march=native` ties the pin to the build host's CPU (the default value comes
from an Apple M4), other machines re-pin with
`-DTURBLYZE_REGRESSION_CD_PIN=<value>` (and, if needed,
`-DTURBLYZE_REGRESSION_CD_TOL=<value>`) from column 4 of the `Cd` line in
`regressionOut/reg_forces.txt`.

Run a subset by label:

```bash
ctest --test-dir build.nosync -L serial       # no MPI launcher involved
ctest --test-dir build.nosync -L mpi          # np = 1, 2, 4
ctest --test-dir build.nosync -L death        # FatalError subprocess probes
ctest --test-dir build.nosync -LE regression  # everything except the slow run
```

On macOS the `mpi` tests set the `HWLOC_SYNTHETIC` environment variable per
test (prte segfaults there without it), so `mpirun` needs no extra environment
in your shell; other platforms keep the real hardware topology.

## Running Simulations

### Basic Execution
Run from the `build.nosync/` directory to ensure correct path resolution:
```bash
cd build.nosync
./Turblyze                     # Uses the default `defaultCase` file
./Turblyze customCase          # Uses a custom case file
```

### Parallel Execution
The same binary runs in parallel under MPI, with no case-file changes. The mesh is decomposed at runtime with METIS:
```bash
mpirun -np 4 ./Turblyze              # Uses the default `defaultCase` file
mpirun -np 4 ./Turblyze customCase   # Uses a custom case file
```

### Case File System
The solver uses a case file system (default file: `defaultCase`). This allows runtime configuration without recompilation.

### Default Case
The default `defaultCase` file contains:
- **Mesh**: `../inputFiles/sphere.msh` (sphere in a channel)
- **Boundary Conditions**:
  - Inlet: Fixed velocity (0, 0, -0.043821) m/s, zero gradient pressure
  - Outlet: Zero gradient velocity, fixed pressure (0 Pa)
  - Walls (`sphere`, `wall1`–`wall4`): No-slip velocity, zero gradient pressure; `kWallFunction`, `omegaWallFunction`, `nutWallFunction` for turbulence
- **Time**: Transient simulation, implicit-Euler stepping (timeStep 0.1 s, totalTime 200 s, 1 outer corrector per step)
- **Discretization**: Second-Order Upwind convection scheme for momentum and Upwind convection scheme for turbulence equations. Least-squares for gradients computation
- **SIMPLE Parameters**: αU = 1.0, αp = 1.0, αk = 0.5, αω = 0.5, tolerance = 1e-3 (scaled residuals), max iterations = 500
- **Turbulence**: `Laminar` (k-omega SST is available via `model kOmegaSST`)
- **Output**: `../outputFiles.nosync/sphere.vtkhdf` (plus `sphere_boundary.vtkhdf`, and `sphere_forces.txt` when forces are enabled). The directory must already exist, the writers do not create it.

### Flow Physics
- **Fluid Properties**: Air (ρ = 1.225 kg/m³, μ = 1.7894e-5 Pa·s)
- **Flow Type**: Low-Reynolds laminar flow over a sphere
- **Turbulence Inlet Conditions** (when k-omega SST is enabled): Turbulence intensity 5%, hydraulic diameter 0.01 m

## Input/Output

### Mesh Requirements
- **Format**: Fluent `.msh` files (ASCII and binary)
- **Dimension**: 3D only (2D meshes are not supported)
- **Cell Types**: Tetrahedra, hexahedra, prisms, pyramids
- **Boundary Patches**: Named patches for boundary condition assignment

### Output Visualization
- **Format**: VTKHDF (format version 2.4), an HDF5-based container read
  natively by ParaView. Use **ParaView 6.1 or newer** (VTK 9.6): earlier
  versions (5.13–6.0, VTK 9.5) open the files but load the polyhedral volume
  grid pathologically slowly through a per-cell decomposition path that VTK
  9.6 replaced with a bulk read. Each run writes a volume UnstructuredGrid
  (`<name>.vtkhdf`) and a boundary PolyData (`<name>_boundary.vtkhdf`)
- **Fields Exported**:
  - Volume `.vtkhdf`: `pressure`, `velocityMagnitude`, vector `velocity`,
    and, when `model` is not `Laminar`, `k`, `omega`, `nut`, `wallDistance`
  - Boundary `_boundary.vtkhdf` (e.g. `sphere_boundary.vtkhdf`): all boundary patches with integer `patchIdx`, `patchZoneIdx`, `patchTypeIdx`, and `isWall` metadata.
    `wallShearStress` is included for all runs; `yPlus` is added only when `model` is not `Laminar`
- **Cell Encoding**: volume cells are written as `VTK_POLYHEDRON` to preserve
  Turblyze's face topology. This is more robust for mixed/polyhedral meshes,
  but some ParaView filters may run slower than with native
  tetra/hex/wedge/pyramid cells.
- **Transient runs**: every written step is appended into the same two
  `.vtkhdf` files, the mesh geometry is stored once and only the field data
  grows per step, so a long series stays a small fraction of the equivalent
  per-step file size. The write cadence is set by `time.writingIntervals`;
  the initial condition and the final step are always written. The files are
  fully closed between writes, so they can be opened (e.g. in ParaView) to
  monitor a run in progress, and an interrupted run still leaves a readable
  series.

### ParaView Visualization (6.1 or newer recommended)
1. Open the `.vtkhdf` file in ParaView
2. Apply the file and click the "eye" icon to make it visible
3. Color by desired field (e.g., `pressure`, `velocityMagnitude`); for a
   transient run, scrub the time slider to animate the series
4. Open the corresponding `_boundary.vtkhdf` file to inspect boundary patches
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
    numIterations               500;
    convergenceTolerance        1e-3;   // Scaled residuals
    nNonOrthogonalCorrectors    0;      // Extra p' corrector re-solves
    relaxationFactors
    {
        U                       0.7;
        p                       0.3;
        k                       0.5;
        omega                   0.5;
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

## Project Structure

Headers and implementations live together under `src/`, following the
OpenFOAM convention:

- **`src/Application/`**: Top-level orchestration and solver assembly
  (`CFDApplication.h/.cpp`, `SolverSetup.h/.cpp`)
- **`src/Primitives/`**: Core scalar/vector/tensor types, logging, errors
- **`src/Mesh/`**: Geometric entities, topology, mesh I/O, mesh checking
- **`src/Fields/`**: Typed cell and face field containers
- **`src/BoundaryConditions/`**: Boundary-type storage, evaluation, and
  case loading
- **`src/Schemes/`**: Gradient, interpolation, and convection schemes
- **`src/LinearSystem/`**: Matrix assembly, transport equations, linear solvers
- **`src/Parallel/`**: MPI runtime, METIS mesh decomposition, ghost-cell halo
  exchange, and collective reductions
- **`src/Solver/`**: Pressure-velocity coupling (SIMPLE and PISO algorithms)
- **`src/Models/`**: Physical models
  - **`src/Models/Turbulence/`**: Turbulence modeling
    (`kOmegaSST.h/.cpp`)
- **`src/PostProcessing/`**: Derived fields and output orchestration
  (`PostProcess.h/.cpp`)
  - **`src/PostProcessing/VTK/`**: VTKHDF volume and boundary writers
    (direct HDF5, no VTK library dependency)
- **`src/Case/`**: Case file system
  (`CaseReader.h/.cpp`, `CaseConfiguration.h/.cpp`)

## Documentation

Browsable HTML API documentation is generated from the Doxygen comments in the source via:

```bash
doxygen Doxyfile
```

Output is written to `docs/doxygen/html/`. Open `docs/doxygen/html/index.html` in a browser to navigate classes, call graphs, and collaboration diagrams. The `docs/doxygen/` tree is generated and is not tracked in git, so regenerate it locally after pulling changes.

Requires `doxygen` and (for diagrams) `graphviz`:
```bash
# macOS
brew install doxygen graphviz

# Ubuntu/Debian
sudo apt install doxygen graphviz
```

## Verification and Validation

Two committed V&V studies back the solver's numerics, each with its own
`README.md` (results) and `runGuide.md` (reproduction steps):

- **`verification/`**: code-to-code comparison against OpenFOAM
  (`simpleFoam`) on the same sphere case: same mesh, boundary conditions,
  k-omega SST constants, and force-coefficient definition.
- **`validation/`**: comparison against experiment: a sphere drag curve traced against the Morrison (2013) correlation, combining a laminar low-Re sweep with subcritical k-omega SST benchmark points.

## Development and Extension

For developers wanting to extend the solver, see `docs/DEVELOPER_GUIDE.md` for:
- Internal architecture details
- Adding new transport equations
- Implementing new boundary conditions
- Creating custom discretization schemes
- Debugging techniques

## License and Support

Turblyze is released under the **Apache License 2.0**, see the [`LICENSE`](LICENSE)
file for the full text.


For questions or bug reports, open an issue on the project's issue tracker.
