<!--
SPDX-FileCopyrightText: 2025-2026 Mohamed Mousa
SPDX-License-Identifier: Apache-2.0
-->

# Case File Reference Guide

This document provides a comprehensive reference for configuring the Turblyze solver using case files.

## Table of Contents

- [Overview](#overview)
- [File Syntax](#file-syntax)
- [Case Sections](#case-sections)
  - [1. mesh](#1-mesh)
  - [2. physicalProperties](#2-physicalproperties)
  - [3. initialConditions](#3-initialconditions)
  - [4. boundaryConditions](#4-boundaryconditions)
  - [5. time](#5-time)
  - [6. numericalSchemes](#6-numericalschemes)
  - [7. velocityCoupling](#7-velocitycoupling-optional)
  - [8. SIMPLE](#8-simple)
  - [9. PISO](#9-piso)
  - [10. linearSolvers](#10-linearsolvers)
  - [11. turbulence](#11-turbulence)
  - [12. forces](#12-forces-optional)
  - [13. output](#13-output)

## Overview

The solver uses case files for configuration. This allows runtime parameter changes without recompilation. The default case file is `defaultCase`, but you can specify any file:

```bash
./Turblyze                     # Uses default defaultCase
./Turblyze customCase          # Uses custom case file
```

**Parallel execution**: The same binary runs in parallel under MPI, with no case-file changes. The mesh is decomposed at runtime with METIS, and the linear systems are solved with PETSc Krylov solvers:

```bash
mpirun -np 4 ./Turblyze              # Uses default defaultCase
mpirun -np 4 ./Turblyze customCase   # Uses custom case file
```

## File Syntax

### Basic Syntax
- **Keywords**: Parameter names followed by values
- **Termination**: All entries end with semicolon `;`
- **Comments**: Single-line `//` or multi-line `/* */`
- **Vectors**: `(x y z)` format for 3D vectors
- **Nested sections**: `{ }` braces for sub-sections

### Example Syntax
```cpp
// Single parameters
keyword         value;
density         1.225;
enabled         true;

// Vectors
velocity        (1.0 0.0 0.0);

// Nested sections
section
{
    parameter1      value1;
    parameter2      value2;
}

// Comments
// This is a single-line comment
/* This is a
   multi-line comment */
```

## Case Sections

### 1. mesh
Specify the mesh file and control the quality checking.

```cpp
mesh
{
    file            ../inputFiles/sphere.msh;       // Required: Path to mesh file
    checkQuality    true;                           // Optional: Enable mesh quality check
}
```

The path is resolved relative to the working directory the solver is launched
from, and may be absolute. Interior spaces are preserved, so a path such as
`/Users/me/My Meshes/sphere.msh` works unquoted; runs of whitespace collapse to
a single space.

### 2. physicalProperties
Defines fluid properties.

```cpp
physicalProperties
{
    rho             1.225;          // Required: Density [kg/m³]
    mu              1.7894e-5;      // Required: Dynamic viscosity [Pa·s]
}
```

### 3. initialConditions
Sets initial field values for all solved fields. Three initialization modes are supported via `type`:

#### A. Uniform (default)
Uniform velocity and pressure assignment throughout the entire domain:
```cpp
initialConditions
{
    type            Uniform;        // Optional (default: Uniform)
    U               (0 0 -0.1);     // Required: Initial velocity [m/s]
    p               0;              // Required: Initial pressure [Pa]
    k               3.75e-5;        // Optional: Initial TKE [m^2/s^2]
    omega           1.6;            // Optional: Initial omega [1/s]
}
```

**Notes**:
- `type` is optional and defaults to `Uniform` for full backward compatibility.
- `k` and `omega` are only used when `model` is not `Laminar`
- When omitted, they are auto-computed from initial velocity using
  `turbulenceIntensity` and `hydraulicDiameter` from the `turbulence`
  section (see [turbulence](#11-turbulence)):
  - `l = 0.07 * hydraulicDiameter`
  - `k = 1.5 * (I * |U|)^2`
  - `omega = sqrt(k) / (C_mu^0.25 * l)`
- When specified, explicit values override the auto-computed defaults
- Turbulent viscosity (`nu_t`) is derived as `k / omega`

### 4. boundaryConditions
Sets up boundary conditions for all fields.

```cpp
boundaryConditions
{
    U   // Velocity field
    {
        inlet
        {
            type        fixedValue;         // BC type
            value       (0 0 -0.1);         // BC value (for fixedValue)
        }

        outlet
        {
            type        zeroGradient;       // No value needed
        }

        walls
        {
            type        noSlip;             // Equivalent to fixedValue (0 0 0)
        }
    }

    p   // Pressure field
    {
        inlet           { type zeroGradient; }
        outlet          { type fixedValue; value 0; }
        walls           { type zeroGradient; }
    }

    k   // Turbulent kinetic energy
    {
        inlet           { type fixedValue; value calculated; }
        outlet          { type zeroGradient; }
        walls           { type kWallFunction; }
    }

    omega   // Specific dissipation rate
    {
        inlet           { type fixedValue; value calculated; }
        outlet          { type zeroGradient; }
        walls           { type omegaWallFunction; }
    }

    nut   // Turbulent kinematic viscosity
    {
        inlet           { type zeroGradient; }
        outlet          { type zeroGradient; }
        walls           { type nutWallFunction; }
    }
}
```

**Boundary Condition Types:**
- `fixedValue`: Fixed value at boundary (requires `value`)
- `zeroGradient`: Zero normal gradient
- `noSlip`: No-slip condition for velocity (equivalent to `fixedValue (0 0 0)`)
- Symmetry-plane (mirror) condition : **mesh-derived, not a case-file type**.
  The Fluent `.msh` zone type determines which patches are symmetry planes;
  at load time every field on such a patch is set to the internal symmetry
  condition automatically. Scalars get a zero normal gradient and zero normal
  flux; velocity has its wall-normal component driven to zero while the
  tangential components are mirrored (`U_f = U_P - (U_P·n) n`), so the face
  carries no mass flux. Symmetry patches should be **omitted** from
  `boundaryConditions`: any entry for one (including an OpenFOAM-style
  `type symmetry;`) is ignored with a `[WARNING]` , the user has no say over
  a symmetry boundary. Writing `type symmetry;` on a patch whose mesh zone
  is *not* a symmetry plane is a fatal error. A symmetry patch must be
  planar. Wall functions and wall-distance seeding key off the mesh zone
  type (`wall`), not this condition, so a symmetry patch, never a wall
  zone, is automatically excluded from both.
- `kWallFunction` / `omegaWallFunction` / `nutWallFunction`: Wall-function
  conditions for `k`, `omega`, and `nut` on wall patches. They must be
  configured as a complete triplet on a given wall patch (all three) or
  omitted entirely; a partial set is a fatal configuration error.
- `fixedGradient`: Fixed normal gradient at the boundary (requires
  `gradient`). The face value is reconstructed as `phi_f = phi_P + gradient *
  dn`, and the diffusive flux is the prescribed `Gamma_f * gradient * |Sf|`.
  On a velocity patch `gradient` is a vector; on a scalar field it is a
  scalar. On `p` it additionally drives the explicit boundary flux
  correction (fixed-flux pressure behaviour).
  ```cpp
  // k { walls  { type fixedGradient; gradient 100;       } }
  // U { inlet  { type fixedGradient; gradient (0 0 5);   } }
  // p { outlet { type fixedGradient; gradient 0.5;       } }
  ```

**Note**: Using an unrecognized type in a case file is a fatal configuration
error. `BCLoader` aborts with
`Unknown boundary condition type '...' for field '...' on patch '...'. Valid types: ...`
rather than falling back silently.

**Calculated values:** For `k` and `omega` boundary conditions, `value`
can be set to `calculated` instead of a numeric value. The solver will
compute the value from `turbulenceIntensity` and `hydraulicDiameter`
(see [turbulence](#11-turbulence)):
- `k = 1.5 * (I * |U|)^2`
- `omega = sqrt(k) / (C_mu^0.25 * 0.07 * D_h)`

### 5. time
**Required in every case file.** Declares the time scheme, which is
**authoritative** for the velocity-coupling algorithm: `steadyState` selects
**SIMPLE** and any transient scheme selects **PISO** (see
[velocityCoupling](#7-velocitycoupling-optional)). For a steady run, set
`timeScheme steadyState`; every other key in this section is then ignored. A
transient `timeScheme` switches to the transient (URANS) path, which runs with a
fixed number of PISO outer correctors per step.

```cpp
time
{
    // Options: steadyState, implicitEuler, CrankNicolson, secondOrderImplicit
    timeScheme          secondOrderImplicit;

    timeStep            0.1;        // Time step size [s]
    totalTime           10;         // Total simulated time [s]
    writingIntervals    1;          // Write output every N steps (>= 1, default 1)
    nOuterCorrectors    20;         // Iterations per step (>= 1, default 1)
    CrankNicolsonCoeff  0.5;        // Crank-Nicolson coefficient [0, 1], default 1.0
}
```

**Time Scheme Options:**
- `steadyState`: No time term; the run is a steady SIMPLE solve.
- `implicitEuler`: First-order backward Euler, fully implicit.
- `CrankNicolson`: Second-order trapezoidal. Blended toward implicit Euler by
  `CrankNicolsonCoeff`: `1.0` is pure Crank-Nicolson and
  `0.0` collapses to implicit Euler; values around `0.5`-`0.9` trade accuracy
  for stability.
- `secondOrderImplicit`: Three-time-level second-order
  backward implicit scheme. Unconditionally stable, recommended for LES. Automatically uses first-order Euler on the initial startup time step.

**Notes**:
- `timeStep` and `totalTime` are **required** and must be positive whenever a
  transient scheme is selected. The number of steps is
  `round(totalTime / timeStep)` (at least one).
- `writingIntervals` (default `1`) writes output every N steps; the initial
  condition (t = 0) and the final step are always written. Must be `>= 1`.
- `nOuterCorrectors` (default `1`) is the number of PISO outer correctors
  performed each time step. Must be `>= 1`. The `SIMPLE.numIterations`/
  `convergenceTolerance` keys only govern the steady path.
- `CrankNicolsonCoeff` (default `1.0`) is consumed only by `CrankNicolson`; it
  must lie in `[0, 1]`.
- Transient output is appended into one temporal `.vtkhdf` file per grid
  (geometry stored once), see [output](#13-output).

### 6. numericalSchemes
Selects discretization schemes.

**Per-Equation Format** (recommended):
```cpp
numericalSchemes
{
    gradient        leastSquares;           // Cell-gradient reconstruction

    convection
    {
        default     Upwind;                 // Fall back for unspecified scheme

        U           SecondOrderUpwind;      // Momentum equations (U_x, U_y, U_z)
        
        k           Upwind;                 // Turbulent kinetic energy (optional)
        
        omega       Upwind;                 // Specific dissipation rate (optional)
    }
}
```

**Gradient Scheme Options:**
- `leastSquares`: Weighted least-squares with inverse-distance weighting (default)

The `gradient` entry may be omitted, in which case it defaults to `leastSquares`.
An unknown name is rejected at startup.

**Convection Scheme Options:**
- `Upwind`: First-order upwind (stable, diffusive)
- `CentralDifference`: Second-order central difference (accurate, may oscillate)
- `SecondOrderUpwind`: Second-order upwind (balance of accuracy and stability)
- `LUST`: Linear-Upwind Stabilized Transport (75% CentralDifference + 25% SecondOrderUpwind; low numerical dissipation with bounded stability, recommended for LES)

### 7. velocityCoupling (Optional)
Selects the pressure–velocity coupling algorithm. The two algorithms are tied
1:1 to the time scheme (see [time](#5-time)): **SIMPLE is steady-only** and
**PISO is transient-only**. The [time](#5-time) section is authoritative. If
`algorithm` disagrees with it, the parser **switches the algorithm to match and
prints a `[WARNING]`** rather than aborting. When this section is omitted,
`algorithm` defaults to `SIMPLE` (then corrected to `PISO` if the time scheme is
transient).

```cpp
velocityCoupling
{
    algorithm           PISO;          // SIMPLE (steady) | PISO (transient)
}
```

**Algorithm Options:**
- `SIMPLE`: Semi-implicit pressure-linked algorithm. **Steady-only.** Relies on
  under-relaxation; reads its controls from the [SIMPLE](#8-simple) section.
- `PISO`: Pressure-implicit with splitting of operators. **Transient-only.**
  Each outer iteration is one implicit momentum predictor followed by
  `nPrimeCorrectors` explicit PRIME corrector steps: the momentum equation is
  re-assembled with the current flux and advanced by a single explicit Jacobi
  sweep (no linear solve) before each pressure correction. Reads its controls
  from the [PISO](#9-piso) section.

**Notes**:
- The section is **optional**; when absent, `algorithm` is inferred from the
  time scheme.
- On a mismatch the time scheme wins: a transient scheme with `algorithm SIMPLE`
  ⇒ `[WARNING]` + switch to PISO; `steadyState` with `algorithm PISO` ⇒
  `[WARNING]` + switch to SIMPLE.
- Only the active algorithm's section is read (`SIMPLE{}` or `PISO{}`); the other
  may be absent.

### 8. SIMPLE
Steady-state SIMPLE controls. Read only when SIMPLE is the active algorithm.

```cpp
SIMPLE
{
    numIterations            100;       // Maximum iterations
    convergenceTolerance    1e-6;       // Convergence criterion
    nNonOrthogonalCorrectors  1;        // p' corrector re-solves (default 0)

    relaxationFactors
    {
        U                   0.7;        // Velocity under-relaxation [0-1]
        p                   0.3;        // Pressure under-relaxation [0-1]
        k                   0.5;        // Turbulent kinetic energy relaxation
        omega               0.5;        // Specific dissipation rate relaxation
    }
}
```

**Non-orthogonal correctors** (`nNonOrthogonalCorrectors`, optional, default `0`): number of pressure-correction re-solves per SIMPLE iteration, matching simpleFoam's non-orthogonal pressure corrector loop. Because p' restarts from zero every iteration, its first solve carries no non-orthogonal correction; each corrector recomputes grad(p') and re-solves with the explicit correction term. Use `0` for orthogonal (hexahedral) meshes and `1`–`2` for tetrahedral or polyhedral meshes. Each corrector adds one pressure solve per iteration.

### 9. PISO
Transient PISO controls. Read only when PISO is the active algorithm. The
time-marching cadence (`timeStep`, `totalTime`, `nOuterCorrectors`) lives in the
[time](#5-time) section.

```cpp
PISO
{
    convergenceTolerance        1e-3;   // Per-time-step convergence criterion
    nNonOrthogonalCorrectors    0;      // p' corrector re-solves (default 0)
    nPrimeCorrectors            2;      // Explicit PRIME correctors (>= 1, default 1)

    relaxationFactors
    {
        U                   1.0;        // Velocity under-relaxation
        p                   1.0;        // Pressure under-relaxation
        k                   1.0;        // Turbulent kinetic energy relaxation
        omega               1.0;        // Specific dissipation rate relaxation
    }
}
```

**Notes**:
- `nPrimeCorrectors` (default `1`, must be `>= 1`) is the number of explicit
  PRIME corrector steps per outer iteration.
- PISO is stable because the transient term makes the momentum diagonal
  dominant; very high Courant numbers can destabilize the explicit sweep. Keep
  `U` and `p` relaxation at `1.0`.

### 11. linearSolvers
Linear solver settings for each field.

```cpp
linearSolvers
{
    U   // Momentum equations
    {
        solver              BiCGSTAB;       // Solver type
        preconditioner      Jacobi;         // Supported preconditioner
        tolerance           1e-5;           // Relative residual tolerance
        maxIter             500;            // Maximum iterations
    }

    p   // Pressure equation
    {
        solver              PCG;            // Solver type
        preconditioner      Jacobi;         // Supported preconditioner
        tolerance           1e-6;           // Relative residual tolerance
        maxIter             1000;           // Maximum iterations
    }

    k   // Turbulent kinetic energy (if turbulence enabled)
    {
        solver              BiCGSTAB;       // Solver type
        preconditioner      Jacobi;         // Supported preconditioner
        tolerance           1e-5;           // Relative residual tolerance
        maxIter             500;            // Maximum iterations
    }

    omega   // Specific dissipation rate (if turbulence enabled)
    {
        solver              BiCGSTAB;       // Solver type
        preconditioner      Jacobi;         // Supported preconditioner
        tolerance           1e-5;           // Relative residual tolerance
        maxIter             500;            // Maximum iterations
    }

    // Optional: raw PETSc options, scoped per solver by equation prefix
    petscOptions        -pressure_pc_type icc -momentum_ksp_view;
}
```

Recognized keys per section:

- `solver`: linear solver type. Defaults to `BiCGSTAB` for U/k/omega
  and `PCG` for pressure. Available options:
  - `BiCGSTAB`: Bi-Conjugate Gradient Stabilized (recommended for asymmetric
    momentum and turbulence equations)
  - `PCG`: Preconditioned Conjugate Gradient (recommended for
    symmetric positive-definite pressure Poisson systems)
  - `GMRES`: Generalized Minimal Residual
  - `FGMRES`: Flexible GMRES (supports variable/nonlinear preconditioners)
  - `TFQMR`: Transpose-Free Quasi-Minimal Residual
  - `CGS`: Conjugate Gradient Squared
  - `MINRES`: Minimal Residual method (for symmetric systems)
  - `Richardson`: Stationary Richardson iteration
  - `Chebyshev`: Chebyshev polynomial iteration
  - `PreOnly`: Applies only the preconditioner (e.g. for direct LU/Cholesky
    or pure multigrid)

- `preconditioner`: system preconditioner. Defaults to `Jacobi`.
  Available options:
  - `Jacobi`: Point-Jacobi / diagonal scaling (cheap, highly parallel scalable)
  - `None`: No preconditioning (Identity)
  - `ILU`: Incomplete LU factorization (level-0 ILU; in parallel runs as
    block-Jacobi with local ILU on each rank)
  - `ICC`: Incomplete Cholesky factorization (for symmetric positive-definite
    pressure systems)
  - `SOR`: Successive Over-Relaxation / Gauss-Seidel (alias `GaussSeidel`)
  - `AMG`: Geometric-Algebraic Multigrid (alias `GAMG`, PETSc native algebraic
    multigrid; recommended for fast $O(N)$ convergence on large pressure Poisson
    grids in serial and parallel)
  - `BlockJacobi`: Block Jacobi across MPI processes
  - `ASM`: Additive Schwarz Method (overlapping domain decomposition)
  - `LU`: Direct sparse LU factorization
  - `Cholesky`: Direct sparse Cholesky factorization (for symmetric systems)

- `tolerance`: relative residual tolerance used by the iterative
  solvers (`|r| / |b|`, true unpreconditioned residual).
- `maxIter`: iteration cap before the solver gives up.
- `petscOptions`: optional string forwarded verbatim into the PETSc
  options database at startup (read up to the terminating `;`, spaces
  preserved). Each solver reads the database through its own equation
  prefix — `momentum_`, `pressure_`, `k_`, `omega_` — so an entry targets
  one solver: `-pressure_pc_type icc` changes only the pressure
  preconditioner, `-momentum_ksp_view` prints only the momentum solver's
  configuration. Unprefixed options (`-pc_type`) match no solver.

Algorithm/equation pairing is not validated. Picking `PCG` for a
non-symmetric equation will compile and run but will not converge to the
correct solution.

### 12. turbulence
Turbulence model configuration.

```cpp
turbulence
{
    model               kOmegaSST;  // Required: Laminar or kOmegaSST
    turbulenceIntensity 0.05;       // Optional: default 0.05 (5%)
    hydraulicDiameter   0.01;       // Optional: default 0.01 [m]
    roughWall           false;      // Optional: default false
}
```

**Notes**:
- `model` is required and selects the turbulence treatment: `Laminar`
  (no turbulence transport) or `kOmegaSST`. An omitted or unknown `model`
  aborts at load time with a fatal error listing the valid options.
- `turbulenceIntensity` and `hydraulicDiameter` are used to auto-compute
  initial values for `k` and `omega` when they are not explicitly
  specified in `initialConditions` (see [initialConditions](#3-initialconditions)).
  They are ignored when `model` is `Laminar`.
- `roughWall` enables the SST $F_3$ rough-wall blending function. When
  `true`, the eddy-viscosity limiter ($F_2$) is damped inside the
  roughness sublayer to prevent artificial suppression of $\nu_t$ on
  aerodynamically rough surfaces. Leave `false` (default) for smooth-wall simulations. This skips the $F_3$ computation entirely.
- k-omega SST model constants are hardcoded in
  `src/Models/Turbulence/kOmegaSST.h`
  and cannot be changed via case file
- Wall distance is computed using meshWave iterative propagation method
  (not configurable)

### 12. forces (Optional)
Aerodynamic force and coefficient reporting on one wall patch.

```cpp
forces
{
    enabled             true;           // Enable force integration
    patch               sphere;         // Wall patch to integrate
    dragDirection       (0 0 -1);       // Direction for Cd projection
    liftDirection       (0 1 0);        // Direction for Cl projection
    referenceVelocity   (0 0 -20.0);    // Reference velocity [m/s]
    referenceArea       0.00785375;     // Reference area [m^2]
}
```

**Notes**:
- When `enabled` is `false`, the remaining entries are ignored.
- `patch`, `dragDirection`, `liftDirection`, `referenceVelocity`, and
  `referenceArea` are required when `enabled` is `true`.
- `dragDirection` and `liftDirection` are normalized by the parser before use.
- Force coefficients use `0.5 * rho * |referenceVelocity|^2 *
  referenceArea`.
- Pressure loads are integrated with the face projected area vector. Friction
  loads use the model-provided wall shear stress and face contact area.
- **Steady runs** write a single summary next to the configured volume output
  file as `<name>_forces.txt`.
- **Transient runs** additionally write a per-time-step history as
  `<name>_forces.csv` with columns
  `time,pressureDrag,frictionDrag,totalDrag,pressureLift,frictionLift,totalLift,Cd,Cl`
  (one row per step, including t = 0), and still print the final-step summary.

### 13. output
Output configuration.

```cpp
output
{
    filename        ../outputFiles/result.vtkhdf;        // Output file path
    debug           false;      // Optional: verbose console output
}
```

**Notes**:
- Output format is always VTKHDF (format version 2.4; read with ParaView 6.1
  or newer). The configured filename writes the volume
  `<name>.vtkhdf` (UnstructuredGrid); a sibling `<name>_boundary.vtkhdf`
  (PolyData) is also written for all boundary patches. A known legacy
  extension on `filename` (`.vtu`, `.vtp`, `.pvd`, `.vtkhdf`) is stripped
  before the `.vtkhdf` names are formed, so older case files keep working.
- All computed volume fields are written to the volume file. Boundary patch
  metadata (`patchIdx`, `patchZoneIdx`, `patchTypeIdx`, `isWall`) is written to
  `_boundary.vtkhdf`; `wallShearStress` is written for all runs, while
  turbulence-only wall diagnostics such as `yPlus` are added only when
  turbulence is enabled.
- Volume cells are encoded as `VTK_POLYHEDRON` to preserve the mesh's
  face-based topology. Some downstream ParaView filters may run slower than
  with native tetra/hex/wedge/pyramid output.
- **Steady runs** write a single-step VTKHDF series (one time value at t = 0).
- **Transient runs** append every written step into the same two `.vtkhdf`
  files: the mesh geometry is stored exactly once and only the per-step field
  data grows, which keeps long series far smaller than one file per step.
  The t = 0 initial condition and the final step are always written;
  intermediate steps follow `time.writingIntervals`. Open the `.vtkhdf` in
  ParaView and scrub the time slider to animate the series. The files are
  fully closed between writes, so they can be opened mid-run to monitor
  progress, and an interrupted run still leaves a readable series.
- `debug` (default: `false`): When `true`, enables verbose
  console output including mesh geometry details, boundary
  condition summaries, solver configuration, per-equation
  solver convergence, turbulence field diagnostics, and
  export statistics. When
  `false`, only essential output is shown (phase headers, iteration residuals,
  convergence status, flow statistics, and error/warning messages).
