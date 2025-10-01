# Setup File Reference Guide

This document provides a comprehensive reference for setting up the MyCFDCode solver using dictionary files.

## Table of Contents

- [Overview](#overview)
- [Dictionary File Syntax](#dictionary-file-syntax)
- [Setup Sections](#setup-sections)
- [Complete Parameter Reference](#complete-parameter-reference)
- [Examples](#examples)
- [Common Use Cases](#common-use-cases)
- [Troubleshooting](#troubleshooting)

## Overview

The MyCFDCode solver uses OpenFOAM-style dictionary files for setup. This allows runtime parameter changes without recompilation. The default setup file is `defaultSetup`, but you can specify any file:

```bash
./MyCFDCode                    # Uses default defaultSetup
./MyCFDCode customSetup        # Uses custom setup
```

## Dictionary File Syntax

### Basic Syntax
- **Keywords**: Parameter names followed by values
- **Termination**: All entries end with semicolon `;`
- **Comments**: Single-line `//` or multi-line `/* */`
- **Vectors**: `(x y z)` format for 3D vectors
- **Nested dictionaries**: `{ }` braces for sub-sections

### Example Syntax
```cpp
// Single parameters
keyword         value;
density         1.225;
enabled         true;

// Vectors
velocity        (1.0 0.0 0.0);

// Nested dictionaries
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

## Setup Sections

### 1. mesh
Controls mesh reading and quality checking.

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
Sets initial field values.

```cpp
initialConditions
{
    U               (0 0 -0.1);     // Required: Initial velocity [m/s]
    p               0;              // Required: Initial pressure [Pa]
}
```

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
            value       (0 0 -0.1);        // BC value (for fixedValue)
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
}
```

**Boundary Condition Types:**
- `fixedValue`: Fixed value at boundary (requires `value`)
- `zeroGradient`: Zero normal gradient
- `noSlip`: No-slip condition for velocity (equivalent to `fixedValue (0 0 0)`)
- `fixedGradient`: Fixed gradient (not currently implemented)

### 5. numericalSchemes
Selects discretization schemes.

```cpp
numericalSchemes
{
    convection      Upwind;         // Options: Upwind, CentralDifference, SecondOrderUpwind
    gradient        leastSquares;   // Currently only leastSquares available
}
```

**Convection Scheme Options:**
- `Upwind`: First-order upwind (stable, diffusive)
- `CentralDifference`: Second-order central difference (accurate, may oscillate)
- `SecondOrderUpwind`: Second-order upwind (balance of accuracy and stability)

### 6. SIMPLE
SIMPLE algorithm parameters.

```cpp
SIMPLE
{
    nIterations             100;        // Maximum iterations
    convergenceTolerance    1e-6;       // Convergence criterion

    relaxationFactors
    {
        U                   0.3;        // Velocity under-relaxation [0-1]
        p                   0.1;        // Pressure under-relaxation [0-1]
        k                   0.5;        // Turbulent kinetic energy relaxation
        omega               0.5;        // Specific dissipation rate relaxation
    }

    pressureCorrector                   // Optional section
    {
        nCorrectors         1;          // Number of pressure corrections
        nNonOrthCorrections 0;          // Non-orthogonal corrections
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
        solver          BiCGSTAB;       // Solver type (currently only BiCGSTAB)
        preconditioner  ILUT;           // Preconditioner (currently only ILUT)
        tolerance       1e-5;           // Absolute tolerance
        relTol          0.1;            // Relative tolerance
        maxIter         500;            // Maximum iterations
    }

    p   // Pressure equation
    {
        solver          BiCGSTAB;
        preconditioner  ILUT;
        tolerance       1e-6;
        relTol          0.05;
        maxIter         1000;
    }

    k   // Turbulent kinetic energy (if turbulence enabled)
    {
        solver          BiCGSTAB;
        preconditioner  ILUT;
        tolerance       1e-5;
        relTol          0.1;
        maxIter         500;
    }

    omega   // Specific dissipation rate (if turbulence enabled)
    {
        solver          BiCGSTAB;
        preconditioner  ILUT;
        tolerance       1e-5;
        relTol          0.1;
        maxIter         500;
    }
}
```

### 8. turbulence
Turbulence model setup.

```cpp
turbulence
{
    model           kOmegaSST;      // Model type (currently only kOmegaSST)
    enabled         true;           // Enable/disable turbulence

    kOmegaSST                       // Model constants (optional)
    {
        a1              0.31;       // Model constant
        betaStar        0.09;       // Model constant
        sigma_k1        0.85;       // Model constant
        sigma_k2        1.0;        // Model constant
        sigma_omega1    0.5;        // Model constant
        sigma_omega2    0.856;      // Model constant
    }

    wallDistance                    // Wall distance calculation (optional)
    {
        method          Poisson;    // Method (currently only Poisson)
        maxIterations   100;        // Max iterations for wall distance
    }
}
```

### 9. constraints (Optional)
Solution field constraints.

```cpp
constraints
{
    velocity                        // Velocity limiting
    {
        enabled         false;      // Enable velocity constraint
        maxVelocity     0.3;       // Maximum velocity magnitude [m/s]
    }

    pressure                        // Pressure limiting
    {
        enabled         false;      // Enable pressure constraint
        minPressure     -0.05;     // Minimum pressure [Pa]
        maxPressure     0.05;      // Maximum pressure [Pa]
    }
}
```

### 10. output
Output setup.

```cpp
output
{
    format          VTK;            // Output format (currently only VTK)
    filename        ../outputFiles/result.vtp;    // Output file path
    writeInterval   1;              // Write every N iterations
    fields          (U p k omega nu_t wallDistance velocityMagnitude);  // Fields to export

    monitoring                      // Optional convergence monitoring
    {
        enabled     true;           // Enable monitoring
        logFile     convergence.log;    // Log file name
        residuals   (U p k omega);      // Fields to monitor
    }
}
```

### 11. runControl (Optional)
Simulation control parameters.

```cpp
runControl
{
    startFrom       latestTime;     // Options: firstTime, latestTime
    stopAt          endIteration;   // Options: endIteration, convergence
    endIteration    1000;           // Maximum iterations
    printInterval   10;             // Print status every N iterations
}
```

## Complete Parameter Reference

### Required Parameters
These parameters must be present in every setup file:

| Section | Parameter | Type | Description |
|---------|-----------|------|-------------|
| mesh | file | string | Path to mesh file |
| physicalProperties | rho | scalar | Fluid density [kg/m³] |
| physicalProperties | mu | scalar | Dynamic viscosity [Pa·s] |
| initialConditions | U | vector | Initial velocity [m/s] |
| initialConditions | p | scalar | Initial pressure [Pa] |
| numericalSchemes | convection | string | Convection scheme |
| numericalSchemes | gradient | string | Gradient scheme |
| SIMPLE | nIterations | int | Maximum iterations |
| SIMPLE | convergenceTolerance | scalar | Convergence criterion |
| SIMPLE | relaxationFactors.U | scalar | Velocity relaxation factor |
| SIMPLE | relaxationFactors.p | scalar | Pressure relaxation factor |
| turbulence | model | string | Turbulence model |
| turbulence | enabled | bool | Enable turbulence |
| output | format | string | Output format |
| output | filename | string | Output file path |

### Optional Parameters
These parameters have default values if not specified:

| Section | Parameter | Default | Description |
|---------|-----------|---------|-------------|
| mesh | checkQuality | true | Enable mesh quality check |
| SIMPLE | relaxationFactors.k | 0.5 | TKE relaxation factor |
| SIMPLE | relaxationFactors.omega | 0.5 | Omega relaxation factor |
| constraints | velocity.enabled | false | Enable velocity limiting |
| constraints | pressure.enabled | false | Enable pressure limiting |

## Examples

### Example 1: Basic Pipe Flow
```cpp
mesh { file ../inputFiles/pipe.msh; }
physicalProperties { rho 1.225; mu 1.7894e-5; }
initialConditions { U (0 0 -0.1); p 0; }
boundaryConditions
{
    U { inlet { type fixedValue; value (0 0 -0.1); } outlet { type zeroGradient; } }
    p { inlet { type zeroGradient; } outlet { type fixedValue; value 0; } }
}
numericalSchemes { convection Upwind; gradient leastSquares; }
SIMPLE { nIterations 100; convergenceTolerance 1e-6; relaxationFactors { U 0.3; p 0.1; } }
turbulence { model kOmegaSST; enabled false; }
output { format VTK; filename ../outputFiles/pipe.vtp; }
```

### Example 2: Turbulent Flow Setup
```cpp
// Copy defaultSetup to new file: cp defaultSetup turbulentSetup
// Then modify these sections:
physicalProperties { rho 998.2; mu 1.003e-3; }  // Water properties
initialConditions { U (0 0 -2.0); p 0; }        // Higher velocity
boundaryConditions
{
    U { inlet { type fixedValue; value (0 0 -2.0); } /* ... */ }
    // Add phi_wall BCs for walls
    phi_wall { walls { type fixedValue; value 0; } }
}
SIMPLE
{
    nIterations 200;
    relaxationFactors { U 0.5; p 0.2; k 0.3; omega 0.3; }  // More conservative
}
turbulence { model kOmegaSST; enabled true; }   // Enable turbulence
```

## Common Use Cases

### 1. Changing Mesh
```cpp
mesh { file ../inputFiles/new_mesh.msh; }
```

### 2. Different Fluid Properties
```cpp
physicalProperties
{
    rho             998.2;      // Water density
    mu              1.003e-3;   // Water viscosity
}
```

### 3. Higher Reynolds Number Flow
```cpp
initialConditions { U (2.0 0 0); p 0; }        // Higher velocity
boundaryConditions { U { inlet { type fixedValue; value (2.0 0 0); } /* ... */ } }
SIMPLE { relaxationFactors { U 0.5; p 0.2; } } // More conservative relaxation
turbulence { enabled true; }                   // Enable turbulence
```

### 4. More Accurate Discretization
```cpp
numericalSchemes { convection SecondOrderUpwind; gradient leastSquares; }
```

### 5. Tighter Convergence
```cpp
SIMPLE
{
    nIterations 500;
    convergenceTolerance 1e-7;
}
linearSolvers
{
    p { tolerance 1e-8; relTol 0.01; }
    U { tolerance 1e-7; relTol 0.05; }
}
```

## Troubleshooting

### Common Errors

1. **"Cannot open dictionary file"**
   - Check file path and name
   - Ensure file exists and is readable

2. **"Keyword 'X' not found"**
   - Check spelling of parameter names
   - Ensure required parameters are present

3. **"Cannot convert 'X' to Type"**
   - Check parameter value format
   - For vectors, use `(x y z)` format
   - For booleans, use `true/false`, `on/off`, or `yes/no`

4. **"Sub-dictionary 'X' not found"**
   - Check section names and nesting
   - Ensure proper brace matching `{ }`

5. **"Unknown convection scheme"**
   - Valid options: `Upwind`, `CentralDifference`, `SecondOrderUpwind`

6. **Convergence issues**
   - Try more conservative relaxation factors (lower values)
   - Increase maximum iterations
   - Check boundary conditions are physically reasonable
   - For turbulent flow, ensure `phi_wall` BCs are set

### Best Practices

1. **Start conservative**: Use low relaxation factors and simple schemes initially
2. **Validate incrementally**: Test with coarse mesh and simple setup first
3. **Copy working setups**: Use example files as starting points
4. **Check boundary conditions**: Ensure all mesh patches have appropriate BCs
5. **Monitor convergence**: Enable output monitoring to track solution progress
6. **Use appropriate precision**: Match solver setup to physics requirements

### Validation

Before running large simulations:
1. Check mesh quality with `checkQuality true`
2. Verify boundary conditions are physically reasonable
3. Start with few iterations to check stability
4. Compare results with analytical solutions when possible
5. Monitor mass conservation in the output