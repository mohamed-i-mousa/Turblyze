# Case File Reference Guide

This document provides a comprehensive reference for configuring the MyCFDCode solver using case files.

## Table of Contents

- [Overview](#overview)
- [File Syntax](#file-syntax)
- [Case Sections](#case-sections)

## Overview

The solver uses case files for configuration. This allows runtime parameter changes without recompilation. The default case file is `defaultCase`, but you can specify any file:

```bash
./MyCFDCode                    # Uses default defaultCase
./MyCFDCode customCase         # Uses custom case file
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
    file            ../inputFiles/pipe_320k.msh;    // Required: Path to mesh file
    checkQuality    true;                           // Optional: Enable mesh quality check
}
```

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
Sets initial field values for all solved fields.

```cpp
initialConditions
{
    U               (0 0 -0.1);     // Required: Initial velocity [m/s]
    p               0;              // Required: Initial pressure [Pa]
    k               3.75e-5;        // Optional: Initial TKE [m^2/s^2]
    omega           1.6;            // Optional: Initial omega [1/s]
}
```

**Notes**:
- `k` and `omega` are only used when turbulence modeling is enabled
- When omitted, they are auto-computed from initial velocity using
  `turbulenceIntensity` and `hydraulicDiameter` from the `turbulence`
  section (see [turbulence](#8-turbulence)):
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

    phi_wall    // Wall distance field (for turbulence)
    {
        walls           { type fixedValue; value 0; }
    }

    k   // Turbulent kinetic energy
    {
        inlet           { type fixedValue; value calculated; }
        outlet          { type zeroGradient; }
        walls           { type fixedValue; value 0; }
    }

    omega   // Specific dissipation rate
    {
        inlet           { type fixedValue; value calculated; }
        outlet          { type zeroGradient; }
        walls           { type fixedValue; value 1000; }
    }
}
```

**Boundary Condition Types:**
- `fixedValue`: Fixed value at boundary (requires `value`)
- `zeroGradient`: Zero normal gradient
- `noSlip`: No-slip condition for velocity (equivalent to `fixedValue (0 0 0)`)
- `fixedGradient`: **⚠️ Not yet implemented** - infrastructure exists but case file parsing not available
  ```cpp
  // Example syntax (NOT FUNCTIONAL - for future reference only):
  // T { walls { type fixedGradient; gradient 100; } }
  // U { walls { type fixedGradient; gradient (0 0 0.5); } }
  ```

**Note**: Attempting to use `fixedGradient` in case files will result in the solver defaulting to `zeroGradient` behavior without warning. This feature is planned for future implementation.

**Calculated values:** For `k` and `omega` boundary conditions, `value`
can be set to `calculated` instead of a numeric value. The solver will
compute the value from `turbulenceIntensity` and `hydraulicDiameter`
(see [turbulence](#8-turbulence)):
- `k = 1.5 * (I * |U|)^2`
- `omega = sqrt(k) / (C_mu^0.25 * 0.07 * D_h)`

### 5. numericalSchemes
Selects discretization schemes.

**Per-Equation Format** (recommended):
```cpp
numericalSchemes
{
    convection
    {
        default     Upwind;                 // Fall back for unspecified scheme

        U           SecondOrderUpwind;      // Momentum equations (U_x, U_y, U_z)
        
        k           Upwind;                 // Turbulent kinetic energy (optional)
        
        omega       Upwind;                 // Specific dissipation rate (optional)
    }
}
```

**Convection Scheme Options:**
- `Upwind`: First-order upwind (stable, diffusive)
- `CentralDifference`: Second-order central difference (accurate, may oscillate)
- `SecondOrderUpwind`: Second-order upwind (balance of accuracy and stability)

**Note**: Gradient scheme is hardcoded to `leastSquares` and not configurable via case file.

### 6. SIMPLE
SIMPLE algorithm parameters.

```cpp
SIMPLE
{
    numIterations            100;       // Maximum iterations
    convergenceTolerance    1e-6;       // Convergence criterion

    relaxationFactors
    {
        U                   0.7;        // Velocity under-relaxation [0-1]
        p                   0.3;        // Pressure under-relaxation [0-1]
        k                   0.5;        // Turbulent kinetic energy relaxation
        omega               0.5;        // Specific dissipation rate relaxation
    }
}
```

### 7. linearSolvers
Linear solver settings for each field.

```cpp
linearSolvers
{
    U   // Momentum equations
    {
        solver              BiCGSTAB;       // Solver type
        preconditioner      ILUT;           // Preconditioner
        tolerance           1e-5;           // Absolute tolerance
        maxIter             500;            // Maximum iterations
        computeResiduals    false;          // Optional: compute residuals
    }

    p   // Pressure equation
    {
        solver              PCG;
        preconditioner      IncompleteCholesky;
        tolerance           1e-6;
        relTol              0.05;
        maxIter             1000;
        initialShift        1.0;
        computeResiduals    false;
    }

    k   // Turbulent kinetic energy (if turbulence enabled)
    {
        solver              BiCGSTAB;
        preconditioner      ILUT;
        tolerance           1e-5;
        maxIter             500;
        computeResiduals    false;
    }

    omega   // Specific dissipation rate (if turbulence enabled)
    {
        solver              BiCGSTAB;
        preconditioner      ILUT;
        tolerance           1e-5;
        maxIter             500;
        computeResiduals    false;
    }
}
```

**Note**: `relTol` (relative tolerance) is only supported for the pressure solver, not momentum or turbulence solvers.

### 8. turbulence
Turbulence model configuration.

```cpp
turbulence
{
    model               kOmegaSST;  // Options: none, kOmegaSST
    enabled             true;       // Enable/disable turbulence
    turbulenceIntensity 0.05;       // Optional: default 0.05 (5%)
    hydraulicDiameter   0.01;       // Optional: default 0.01 [m]
}
```

**Notes**:
- `turbulenceIntensity` and `hydraulicDiameter` are used to auto-compute
  initial values for `k` and `omega` when they are not explicitly
  specified in `initialConditions` (see [initialConditions](#3-initialconditions))
- k-omega SST model constants are hardcoded in `KOmegaSST.hpp`
  and cannot be changed via case file
- Wall distance is computed using meshWave iterative propagation method
  (not configurable)

### 9. constraints (Optional)
Solution field constraints.

```cpp
constraints
{
    velocity                        // Velocity limiting
    {
        enabled         false;      // Enable velocity constraint
        maxVelocity     0.3;        // Maximum velocity magnitude [m/s]
    }

    pressure                        // Pressure limiting
    {
        enabled         false;      // Enable pressure constraint
        minPressure     -0.05;      // Minimum pressure [Pa]
        maxPressure     0.05;       // Maximum pressure [Pa]
    }
}
```

### 10. output
Output configuration.

```cpp
output
{
    filename        ../outputFiles.nosync/result.vtu;    // Output file path
    debug           false;      // Optional: verbose console output
}
```

**Notes**:
- Output format is always VTK
- All computed fields are written to the output file
- Output is written at the end of the simulation only
- `debug` (default: `false`): When `true`, enables verbose
  console output including mesh geometry details, boundary
  condition summaries, solver configuration, per-equation
  solver convergence, turbulence field diagnostics, and VTK
  export statistics. When `false`, only essential output is
  shown (phase headers, iteration residuals, convergence
  status, flow statistics, and error/warning messages).