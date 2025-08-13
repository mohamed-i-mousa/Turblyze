## 3D CFD Solver for Unstructure Grids

A 3D incompressible CFD solver using the SIMPLE algorithm with k-omega SST turbulence modeling. Reads Fluent `.msh` meshes, solves the steady-state flow, and exports results to VTK (`.vtp`) for ParaView.

### Highlights
- **3D solver**: Solves U, V, W momentum with pressure coupling via SIMPLE
- **Rhie–Chow**: Face-velocity interpolation to avoid pressure checkerboarding
- **Schemes**: Upwind, second-order upwind, and central-difference convection; least-squares gradients; non-orthogonal diffusion treatment
- **BCs**: Fixed value, fixed gradient, zero gradient, no-slip wall
- **Turbulence (optional)**: k-omega SST with wall distance, μ_t, and derived quantities
- **VTK output**: Exports cell fields mapped to faces for visualization in ParaView

## Requirements

- **C++**: C++17-compatible compiler
- **CMake**: >= 3.10
- **Dependencies**:
  - Eigen 3 (headers)
  - VTK (Core, DataModel, IOLegacy, IOXML)

On Debian/Ubuntu:
```bash
sudo apt update
sudo apt install -y build-essential cmake libeigen3-dev libvtk9-dev
```

If CMake cannot find VTK automatically, set `VTK_DIR` to the CMake config folder:
```bash
cmake -DVTK_DIR=/usr/lib/cmake/vtk-9.1 ..
```

## Build

```bash
mkdir build
cd build
cmake ..
make -j
```

The executable is `./MyCFDCode` (same as the project name in `CMakeLists.txt`).

## Run

Run from the `build/` directory so relative input/output paths resolve correctly:
```bash
./MyCFDCode
```

Default behavior (as defined in `src/main.cpp`):
- Reads mesh: `../input_files/cylinder_coarse.msh`
- Applies inlet/outlet/cylinder/symmetry boundary conditions
- Uses Upwind convection, least-squares gradients
- Solves with SIMPLE (under-relaxation U=0.7, p=0.3, tol=1e-6, 30 iters)
- Exports VTK: `../output_files/cylinder_flow.vtp`

You will see iteration logs, residuals (mass, velocity, pressure), statistics, and a runtime summary.

## Input mesh

- Format: Fluent `.msh` (ANSYS Meshing export)
- Dimensionality: 3D only (the reader throws on 2D meshes)
- Boundary patches are read from section 45 and mapped to internal BC types
- Example patches used by the stock case: `inlet`, `outlet`, `cylinder`, `symmetry1..4`

## Configuring a case

Edit `src/main.cpp` to change the mesh, boundary conditions, schemes, and solver controls. Example:
```cpp
// Mesh
std::string meshFilePath = "../input_files/cylinder_coarse.msh";
readMshFile(meshFilePath, allNodes, allFaces, allCells, allBoundaryPatches);

// Boundary conditions
BoundaryConditions bcManager;
for (const auto& patch : allBoundaryPatches) bcManager.addPatch(patch);
bcManager.setFixedValue("inlet", "U", Vector(0.0, 0.0, -0.1));
bcManager.setZeroGradient("inlet", "p");
bcManager.setFixedValue("outlet", "p", 0.0);
bcManager.setZeroGradient("outlet", "U");
bcManager.setNoSlip("cylinder", "U");

// Schemes
GradientScheme gradScheme;
UpwindScheme UDS;               // or CentralDifferenceScheme, SecondOrderUpwindScheme

// SIMPLE setup
SIMPLE simpleSolver(allFaces, allCells, bcManager, gradScheme, UDS);
simpleSolver.setRelaxationFactors(0.7, 0.3);
simpleSolver.setConvergenceTolerance(1e-6);
simpleSolver.setMaxIterations(100);
simpleSolver.enableTurbulenceModeling(true); // set false for laminar
simpleSolver.solve();
```

Notes:
- Physical properties (ρ, μ) currently use defaults inside `SIMPLE` (air-like). Provide setters if you need to vary fluids globally.
- Precision is double by default. See Precision section to change.

## Outputs and visualization

- File: `../output_files/cylinder_flow.vtp` (VTK PolyData)
- Fields exported by default:
  - **pressure** (cell → mapped to faces)
  - **velocityMagnitude**, **U_x**, **U_y**, **U_z**
  - With turbulence enabled: **k**, **omega**, **mu_t**, **wallDistance**, and derived: **turbulentIntensity**, **turbulentViscosityRatio**, **yPlus**

Open the `.vtp` file in ParaView. Since PolyData stores faces as cells, cell-centered quantities are mapped to faces via owner cell. Color by any of the arrays listed above.

## Convergence and diagnostics

- Iteration log prints RMS-like residuals for mass, velocity, and pressure (based on p' RMS)
- Mass conservation is summarized post-solve as average per-cell imbalance
- If residuals grow very large, the solver warns and suggests actions (reduce relaxation, check BCs, improve initial guess)

## Precision control

By default, the build defines `PROJECT_USE_DOUBLE_PRECISION` (FP64). To switch to float (FP32):
1) Remove the compile definition in `CMakeLists.txt`:
```cmake
# target_compile_definitions(MyCFDCode PUBLIC PROJECT_USE_DOUBLE_PRECISION)
```
2) Reconfigure and rebuild. At runtime the program prints the current mode via `SCALAR_MODE`.

## Project layout

- `src/Mesh/`: mesh structures and Fluent `.msh` reader
- `src/BoundaryConditions/`: patch metadata and boundary condition handling
- `src/Numerics/`: gradient/convection schemes, matrix assembly, linear solvers, SIMPLE algorithm
- `src/Models/`: k-omega SST turbulence model
- `src/PostProcessing/`: VTK exporter (PolyData writer)
- `input_files/`: put your `.msh` meshes here
- `output_files/`: VTK results written here

## Troubleshooting

- VTK not found: set `VTK_DIR` to the correct CMake config dir and re-run CMake
- Eigen not found: install `libeigen3-dev` (or provide your Eigen include path)
- 2D mesh error: the reader only supports 3D; export a 3D `.msh`
- Missing input: run from the `build/` folder so `../input_files` resolves correctly
- Empty/odd visuals: PolyData stores faces; ensure you are coloring by a face array in ParaView

## Roadmap

- Advanced no-slip boundary condition 
- Parallel processing for large 3D meshes
- Additional turbulence models (LES)



