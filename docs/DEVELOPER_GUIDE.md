## Developer Guide

This document explains the internal architecture and implementation details of the CFD solver. It is intended for contributors and users who want to understand, extend, or debug the code.

### Table of contents
- Architecture overview
- Core data structures
- Mesh I/O and topology building
- Boundary conditions system
- Numerical schemes (gradients, convection, diffusion)
- Linear system assembly (`Matrix`)
- SIMPLE algorithm (pressure–velocity coupling)
- Rhie–Chow face-velocity interpolation
- Turbulence model (k–omega SST)
- Post-processing and VTK export
- Linear solvers
- Precision and numerical tolerances
- Extending the codebase (recipes)
- Debugging and tips

### Architecture overview (Current Structure)

**Headers (`include/`):**
- **`Application/`**: top-level driver
  - `CFDApplication.hpp`
- **`Mesh/`**: geometry, fields, fundamental types, mesh I/O
  - `Face.hpp`, `Cell.hpp`, `CellData.hpp`, `FaceData.hpp`, `Scalar.hpp`, `Vector.hpp`, `MeshReader.hpp`, `MeshChecker.hpp`
- **`BoundaryConditions/`**: patch metadata and physical BC configuration
  - `BoundaryPatch.hpp`, `BoundaryData.hpp`, `BoundaryConditions.hpp`
- **`Numerics/`**: discretization and algebraic system
  - `GradientScheme.hpp`, `ConvectionScheme.hpp`, `Matrix.hpp`, `LinearSolvers.hpp`, `SIMPLE.hpp`, `LinearInterpolation.hpp`, `Constraint.hpp`
- **`Models/`**: turbulence
  - `kOmegaSST.hpp`
- **`PostProcessing/`**: output
  - `VtkWriter.hpp`
- **`Case/`**: case system
  - `CaseReader.hpp`: OpenFOAM-style case file parser

**Sources (`src/`):**
- Corresponding `.cpp` implementations for all headers
- `main.cpp`: complete end-to-end example case that loads configuration from case file


## Core data structures

### Scalar precision
- `Scalar` is aliased to `double` by default via `PROJECT_USE_DOUBLE_PRECISION` (set in `CMakeLists.txt`).
- Switch to float by removing that definition. The program prints the mode via `SCALAR_MODE`.
- Global tolerances in `include/Mesh/Scalar.hpp` (e.g., `DIVISION_TOLERANCE`, `EQUALITY_TOLERANCE`, `AREA_TOLERANCE`, `VOLUME_TOLERANCE`, `GRADIENT_TOLERANCE`).

### Vector
- Simple 3D vector with arithmetic, `dot`, `cross`, `magnitude`, normalization, and stream IO.
- Used throughout for geometry (centroids, normals) and vector fields.

### Fields
- `CellData<T>`: typed cell-centered fields with bounds-checked access.
- `FaceData<T>`: typed face-centered fields.
- Type aliases:
  - `VectorField`, `ScalarField`, `VelocityField`, `PressureField`
  - `FaceFluxField` (Scalar), `FaceVectorField` (Vector)

### Mesh entities
- `Face`
  - Topology: `nodeIndices`, `ownerCell`, optional `neighbourCell` (boundary if empty).
  - Geometry computed in `calculateGeometricProperties(allNodes)`:
    - Triangles via cross product; polygons triangulated about the face center.
    - Fields: `centroid`, `normal` (unit), `area`, and integrals (`x2_integral`, ...).
  - Metric distances `calculateDistanceProperties(allCells)`:
    - `d_Pf`, `d_Nf` vectors and magnitudes; `e_Pf`, `e_Nf` unit vectors.
- `Cell`
  - Topology: lists of `faceIndices`, `neighbourCellIndices`, and `faceSigns` (outward normal convention).
  - `calculateGeometricProperties(allFaces)`:
    - Volume via divergence theorem: `V = (1/3) Σ (r_f · S_f)` using face integrals.
    - Centroid via second-moment accumulation.


## Mesh I/O and topology building

`MeshReader` reads Fluent `.msh` files (3D only):
- Sections: comments `(0)`, dimension `(2)`, nodes `(10)`, cells `(12)`, faces `(13)`, boundaries `(45)`.
- Fluent uses hexadecimal indices for declarations; helpers convert hex→dec robustly.
- Faces section returns owner and optional neighbor cell; neighbor absent implies boundary.
- Boundaries section maps `zoneID` to `BoundaryPatch` name/type via `mapFluentBCToEnum`.
- After reading:
  - Builds `Cell.faceIndices`, `Cell.faceSigns` (+1 owner, -1 neighbor), and unique `neighbourCellIndices`.
  - Validates: min faces per cell, min nodes per face; prints a summary.

Notes:
- 2D meshes are rejected early (`dimension == 2`).
- Errors throw `std::runtime_error` with descriptive messages.


## Boundary conditions system

### Architecture
**Classes**:
- `BoundaryPatch`: Mesh patch metadata (name, Fluent type, `zoneID`, first/last face indices)
- `BoundaryData`: Type-safe storage with robust value/gradient handling
- `BoundaryConditions`: Manager class with comprehensive BC operations

### BoundaryData Implementation
**Supported BC Types**:
- `FIXED_VALUE`: Dirichlet boundary conditions
- `FIXED_GRADIENT`: Neumann boundary conditions
- `ZERO_GRADIENT`: Natural boundary conditions
- `NO_SLIP`: Special case for velocity (vector value = 0)
- `K_WALL_FUNCTION`: Wall function for turbulent kinetic energy
- `OMEGA_WALL_FUNCTION`: Wall function for specific dissipation rate
- `NUT_WALL_FUNCTION`: Wall function for turbulent viscosity

**Value Storage**:
- `BCValueType::SCALAR` or `BCValueType::VECTOR`
- Type-safe getters: `fixedScalarValue()`, `fixedScalarGradient()`
- Vector access: `vectorValue()`, `vectorGradient()`

### BoundaryConditions Manager
**Data Structure**: `patchBoundaryData[patchName][fieldName] = BoundaryData`

**Key Features**:
1. **Direct Patch Lookup**: `Face::patch()` returns `const BoundaryPatch*` directly, linked at startup via `BoundaryConditions::linkFaces()`
2. **Smart Field Mapping**: Automatic fallback from `Ux`/`Uy`/`Uz` to parent field `U`
3. **Robust Retrieval**: `fieldBC()` with comprehensive error handling
4. **Boundary Value Calculation**: `calculateBoundaryFaceValue()` for scalars, `calculateBoundaryVectorFaceValue()` for vectors

**Vector Component Handling**:
```cpp
// Smart component extraction
if (fieldName == "Ux") val = bc->vectorValue().x();
else if (fieldName == "Uy") val = bc->vectorValue().y();
else if (fieldName == "Uz") val = bc->vectorValue().z();
```

### BC Evaluation Logic
**Scalar Boundary Values**:
- **FIXED_VALUE**: `φ_f = φ_boundary`
- **ZERO_GRADIENT**: `φ_f = φ_owner`  
- **FIXED_GRADIENT**: `φ_f = φ_owner + gradient × d_n`
  where `d_n = dot(d_Pf, face_normal)`

**Vector Boundary Values**:
- **FIXED_VALUE**: `U_f = U_boundary`
- **NO_SLIP**: `U_f = (0, 0, 0)`
- **ZERO_GRADIENT**: `U_f = U_owner`
- **FIXED_GRADIENT**: `U_f = U_owner + gradient × d_n`

**Graceful Fallbacks**:
- Missing BC specifications default to zero-gradient
- Unknown field names default to cell gradient
- Invalid patches use cell values


## Numerical schemes

### Gradient reconstruction (`GradientScheme`)

#### Cell Gradient Computation (`CellGradient`)
**Method**: Weighted least-squares gradient reconstruction

**Algorithm**:
1. **Neighbor Analysis**: Validate neighbor cells and compute distance vectors
2. **Weight Calculation**: `w = 1/r²` for inverse-distance-squared weighting
3. **Matrix Assembly**: Form normal equations `ATA·∇φ = ATb`
   - `ATA = Σ w·(r ⊗ r)` (3×3 matrix)
   - `ATb = Σ w·Δφ·r` (3×1 vector)
4. **Regularization**: Add small diagonal term to prevent singularity
5. **Solution**: Eigen LLT decomposition with LU fallback
6. **Gradient Limiting**: Barth-Jesperson limiter prevents unphysical gradients

**Robustness Features**:
- **Dual solver**: LLT primary, LU fallback for poorly conditioned systems
- **Regularization**: `totalWeight × 1e-12` prevents singular matrices
- **Gradient limiting**: Prevents overshoots in high-gradient regions
- **Error handling**: Comprehensive validation and graceful failures

#### Face Gradient Computation (`FaceGradient`)
**Method**: Corrected interpolation of cell gradients

**Algorithm**:
1. **Boundary Check**: Use `calculateBoundaryFaceGradient()` for boundary faces
2. **Distance Calculation**: `d_PN = centroid_N - centroid_P`
3. **Average Gradient**: Distance-weighted interpolation via `averageFaceGradient()`
4. **Consistency Correction**: `correction = (φ_N - φ_P)/|d_PN| - (∇φ_avg · e_PN)`
5. **Final Result**: `∇φ_f = ∇φ_avg + correction × e_PN`

**Face Gradient Averaging (`averageFaceGradient`)**:
- **Weights**: `g_P = d_Nf/(d_Pf + d_Nf)`, `g_N = d_Pf/(d_Pf + d_Nf)`
- **Formula**: `∇φ_f = g_P × ∇φ_P + g_N × ∇φ_N`
- **Physical meaning**: Closer cell has more influence

#### Boundary Face Gradients (`calculateBoundaryFaceGradient`)
**Approach**: Normal/tangential decomposition

**FIXED_VALUE BC**:
1. Calculate normal gradient: `∂φ/∂n = (φ_boundary - φ_cell)/d_n`
2. Extract tangential components: `∇φ_tan = ∇φ_cell - (∇φ_cell·n)n`
3. Combine: `∇φ_f = ∇φ_tan + (∂φ/∂n)n`

**ZERO_GRADIENT BC**: `∇φ_f = ∇φ_cell`

**FIXED_GRADIENT BC**: 
1. Extract tangential: `∇φ_tan = ∇φ_cell - (∇φ_cell·n)n`
2. Apply normal gradient: `∇φ_f = ∇φ_tan + gradient_specified×n`

### Convection schemes (`ConvectionScheme`)

#### Upwind Differencing Scheme (UDS)
**Coefficients**: 
- `a_P_conv = max(massFlowRate, 0.0)`
- `a_N_conv = min(massFlowRate, 0.0)`

**Flow Direction Logic**:
- **Forward flow** (`mdot > 0`): Use owner cell value
- **Reverse flow** (`mdot < 0`): Use neighbor cell value
- **Sign handling**: `a_N_conv` correctly receives negative flow rates

**Properties**: First-order accurate, unconditionally stable

#### Central Difference Scheme (CDS)
**Implementation**: Deferred correction approach

**Matrix Coefficients**: Same as UDS for stability
**Correction Term**: `mdot × (φ_central - φ_upwind)`

**Face Value Calculation**:
```cpp
φ_f = φ_P × w + φ_N × (1-w) + (∇φ_f · d_Pf)
```
where `w = d_N/(d_P + d_N)` (inverse distance weighting)

**Features**:
- Second-order accurate on structured grids
- Requires face gradients for non-orthogonal correction
- Stable via deferred correction approach

#### Second-Order Upwind (SOU)
**Implementation**: Gradient-based extrapolation

**Face Value Calculation**:
```cpp
if (upwind_cell == owner)
    φ_f = φ_P + (∇φ_P · d_Pf)
else
    φ_f = φ_N + (∇φ_N · d_Nf)
```

**Correction Term**: `mdot × (φ_SOU - φ_UDS)`

**Properties**: Second-order accurate, bounded, TVD-like behavior

### Diffusion treatment
**Orthogonal Component**: Handled implicitly via `E_f = (S_f · e_PN) e_PN`
**Non-orthogonal Correction**: Explicit via `T_f = S_f - E_f` using face gradients
**Formula**: `∇φ_f · T_f` added to RHS for non-orthogonal meshes


## Linear system assembly (`Matrix`)

Uses a unified `TransportEquation` struct to bundle all data for any transport equation (momentum, pressure correction, turbulence). Gradients and mass fluxes are stored in `SIMPLE`, not in `Matrix`.

### TransportEquation struct
```cpp
struct TransportEquation
{
    const ScalarField& phi;             // Current field values
    std::string fieldName;              // "Ux", "k", "pCorr", etc.
    OptionalRef<FaceFluxField> flowRate;    // Face flow rates (nullopt = no convection)
    OptionalRef<ConvectionScheme> convScheme;
    OptionalRef<ScalarField> Gamma;         // Cell-based diffusion coefficient
    OptionalRef<FaceFluxField> GammaFace;   // Face-based diffusion coefficient
    const ScalarField& source;          // Explicit source term
    const VectorField& gradPhi;         // Pre-computed cell gradients
    const GradientScheme& gradScheme;
    OptionalRef<FaceData<Scalar>> boundaryFaceValues;
};
```

### Unified build method
`buildMatrix(const TransportEquation& eq)`:
- Single method handles all equation types:
  - **Momentum**: convection + diffusion via cell-based `Gamma` (nuEff)
  - **Pressure correction**: face-based diffusion via `GammaFace` (DUf), no convection (flowRate = nullopt)
  - **Turbulence k/omega**: convection + diffusion
- Internal faces: assembles diffusion and convection with non-orthogonal correction
- Boundary faces: handles FIXED_VALUE, ZERO_GRADIENT, NO_SLIP, and wall function types
- Deferred-correction for CDS/SOU added to RHS

### Under-relaxation
`relax(α, φ_prev)` performs Patankar-style implicit relaxation by scaling the diagonal and adjusting RHS with the previous state.


## SIMPLE algorithm

Entry point: `SIMPLE::solve()` performs the outer iteration until convergence or `maxIterations`:
1) Store previous-iteration fields (U, face velocities, flow rates), compute gradP.
2) `solveMomentumEquations()`: extracts Ux/Uy/Uz components, computes velocity gradients once into `gradU_` member, solves each component via `solveMomentumComponent()` with `buildMatrix()` + Patankar relaxation.
3) `calculateRhieChowFlowRate()`: compute face velocities and mass fluxes.
4) `solvePressureCorrection()`: pre-compute mass imbalance source, build and solve p' equation using `buildMatrix()` with face-based diffusion (DUf), no convection.
5) `correctVelocity()`: update U using `U = U* - D ∇p'`.
6) `correctFlowRate()`: update face mass fluxes.
7) `correctPressure()`: apply `p = p + α_p p'`; reset p'.
8) `solveTurbulence()`: advance k–ω SST using current fields and pre-computed `gradU_` (if enabled).
9) `checkConvergence()`: monitor mass imbalance (normalized per-cell average), velocity residual (normalized L2), and pressure correction residual (normalized RMS).

Controls:
- `setRelaxationFactors(α_U, α_p, α_k, α_omega)`, `setConvergenceTolerance(tol)`, `setMaxIterations(n)`. Turbulence is enabled by passing `enableTurbulence = true` to `initialize()`.


## Rhie–Chow face-velocity interpolation

Used in `calculateRhieChowFlowRate()` to prevent pressure checkerboarding:
- Start with linear-interpolated face velocity `U_f_lin`.
- Compute face D-like coefficient from interpolated momentum diagonals and geometry.
- Apply correction with face pressure gradient: `U_f = U_f_lin + D_f (∇p_f_lin - ∇p_f_cache)`.
- Add previous-iteration face under-relaxation term `(1-α_U)(U_f_prev - U_f_lin_prev)`.
- Boundary faces use centralized BC evaluation.


## Turbulence model (k–omega SST)

Class `kOmegaSST`:
- Initializes `k`, `ω`, `nut`, and computes `wallDistance` via mesh wave iterative propagation from wall boundaries.
- `solve(U, flowRateFace, gradU)` accepts pre-computed velocity gradients from SIMPLE.
- Solves ω and k transport with variable diffusion (`ν + σ·ν_t`), production/destruction, and cross-diffusion for SST.
- Calculates blending functions `F1`/`F2`, turbulent viscosity `ν_t = a1 k / max(a1 ω, S F2)`, and applies wall corrections.
- Provides getters used by SIMPLE to form effective viscosity and for post-processing: `k`, `ω`, `nut`, `wallDistance`, `wallShearStress`.


## Post-processing and VTK export

`writeVtkUnstructuredGrid(filename, allNodes, allCells, allFaces, scalarCellFields, vectorCellFields)`:
- Writes VTK UnstructuredGrid (`.vtu`) with 3D volume cells (tetrahedra, hexahedra, wedges, pyramids).
- Exports cell-centered scalar fields (pressure, turbulence quantities) and vector fields (velocity).
- Uses topology-aware node ordering for proper VTK cell types (hexahedron, wedge, pyramid).
- Used in `main.cpp` to export pressure, velocity vector `U`, and when available: `k`, `ω`, `μ_t`, `wallDistance`, and derived quantities (`turbulentIntensity`, `turbulentViscosityRatio`, `yPlus`).


## Linear solvers

`LinearSolver` provides both `solveWithBiCGSTAB()` (ILUT preconditioner) and `solveWithPCG()` (Incomplete Cholesky preconditioner):
- Per-field solver instances with independent convergence parameters.
- Configurable absolute tolerance, relative tolerance, and max iterations.
- Optional exact RMS residual computation with convergence diagnostics.
- Throws on non-finite errors; returns success boolean.


## Precision and numerical tolerances

- Precision selected at compile time via `PROJECT_USE_DOUBLE_PRECISION`.
- Tolerance constants adapt to `Scalar` (e.g., comparisons, divisions, gradient detection).
- Many algorithms include small epsilons to guard against degeneracy.


## Extending the codebase (recipes)

### Add a new scalar transport equation
1) Create a `ScalarField phi("phi", numCells, initial)` in your driver.
2) Build an effective diffusion field `Gamma` and a source `phi_source` per cell.
3) Pre-compute cell gradients `gradPhi` via `GradientScheme::cellGradient()`.
4) Create a `TransportEquation` struct with all required fields:
   ```cpp
   TransportEquation eq{phi, "phi", flowRate, convScheme,
       Gamma, std::nullopt, source, gradPhi, gradScheme,
       std::nullopt};
   ```
5) Call `matrix.buildMatrix(eq)`, then solve with `LinearSolver::solveWithBiCGSTAB()`.
6) Apply under-relaxation via `matrix.relax(alpha, phiPrev)` if needed.

### Add a new convection scheme
1) Derive from `ConvectionScheme` (base `getFluxCoefficients` returns `FluxCoefficients` struct).
2) Optionally add high-order face value and correction methods (see CDS/SOU) and integrate as deferred-correction in `Matrix`.

### Add a new boundary condition
1) **Extend enums**: Add new type to `BCType` enum in `BoundaryData.hpp`
2) **Update BoundaryData**: Add setters/getters for new BC type
3) **Extend evaluation**: Update `calculateBoundaryFaceValue()` and `calculateBoundaryFaceGradient()`
4) **Matrix integration**: Update boundary handling in `Matrix::buildMatrix()`
5) **Case file parsing**: Add parsing support in `CFDApplication::setupBoundaryConditions()`

**Example Implementation**:
```cpp
// 1. Add to BCType enum
PERIODIC,  // New BC type

// 2. Add BoundaryData methods
void setPeriodicValue(Scalar offset) {
    type = BCType::PERIODIC;
    scalarValue = offset;
}

// 3. Update evaluation logic
case BCType::PERIODIC:
    return calculatePeriodicValue(face, phi, bc->scalarValue);
```

### Expose fluid properties as inputs
- Add setters on `SIMPLE` for `rho` and `mu`, thread through to `Matrix` and `KOmegaSST` calls as needed.


## Testing and Debugging

### Comprehensive Testing Methodology

#### Boundary Conditions Testing
**Testing Strategy**: Add comprehensive std::cout debugging to trace:
1. **Patch Registration**: Verify patch names, zones, face ranges
2. **BC Storage**: Confirm type-safe storage of scalar/vector values
3. **Field Lookup**: Test `fieldBC()` with various field names
4. **Vector Components**: Verify `Ux`/`Uy`/`Uz` → `U` fallback
5. **Boundary Values**: Test `calculateBoundaryFaceValue()` for all BC types
6. **Patch Linking**: Verify `Face::patch()` returns correct patch after `linkFaces()`

**Key Tests**:
```cpp
// Test all BC types
setFixedValue("inlet", "U", Vector(1,0,0));
setZeroGradient("outlet", "p");
setNoSlip("wall", "U");
setFixedGradient("interface", "T", 100.0);
```

#### Convection Schemes Testing
**Testing Strategy**: Verify coefficient calculation and face values:
1. **Coefficient Logic**: Test `getFluxCoefficients()` for +/- mass flow rates
2. **Flow Direction**: Verify upwind cell selection
3. **Face Values**: Test interpolation and extrapolation methods
4. **Correction Terms**: Verify deferred correction calculations
5. **Boundary Integration**: Test BC application in schemes

**Critical Tests**:
```cpp
// Test flow direction handling
massFlowRate = +1.0: a_P_conv = 1.0, a_N_conv = 0.0  // Owner→Neighbor
massFlowRate = -1.0: a_P_conv = 0.0, a_N_conv = -1.0 // Neighbor→Owner
```

#### Gradient Schemes Testing
**Testing Strategy**: Verify mathematical correctness:
1. **Neighbor Validation**: Check distance calculations and weighting
2. **Matrix Assembly**: Verify ATA matrix conditioning and determinant
3. **Solver Robustness**: Test LLT/LU fallback mechanisms
4. **Gradient Limiting**: Verify limiter activation in high-gradient regions
5. **Face Interpolation**: Test averaging weights and corrections
6. **Boundary Gradients**: Verify normal/tangential decomposition

**Matrix Verification**:
```cpp
// Check matrix properties
ATA.determinant() > 0  // Well-conditioned system
regularization = totalWeight × 1e-12  // Prevents singularity
gradMag * maxDistance < 10.0 * phiRange  // Gradient limiting
```

### Debugging Strategies

#### Adding Debug Output
1. **Method Tracing**: Add entry/exit logging for key methods
2. **Parameter Logging**: Log input parameters and intermediate calculations
3. **Validation Checks**: Add assertions for mathematical consistency
4. **Performance Monitoring**: Track solver iterations and convergence

#### Common Issues and Solutions

**Boundary Condition Issues**:
- **Symptom**: "No BC specified" warnings
- **Solution**: Check patch names match mesh exactly
- **Debug**: Use `printSummary()` to list all patches and BCs

**Convection Scheme Issues**:
- **Symptom**: Incorrect flow direction or instability
- **Solution**: Verify mass flow rate signs and upwind logic
- **Debug**: Log `massFlowRate`, `a_P_conv`, `a_N_conv` values

**Gradient Issues**:
- **Symptom**: "Gradient computation failed" errors
- **Solution**: Check mesh quality and neighbor connectivity
- **Debug**: Log ATA matrix condition number and rank

#### Best Practices
1. **Modular Testing**: Test individual components before integration
2. **Mathematical Verification**: Verify algorithms against literature
3. **Boundary Case Testing**: Test with extreme parameter values
4. **Performance Profiling**: Monitor computational efficiency
5. **Regression Testing**: Maintain test cases for future validation

### Development Tips

- **BCs**: Use `BoundaryConditions::printSummary()` to inspect configuration
- **Gradients**: Check matrix conditioning with `ATA.determinant()`
- **Convection**: Verify upwind logic with simple 1D test cases
- **Solver logs**: High residuals indicate BC or relaxation issues
- **Mesh validation**: Reader throws early for malformed `.msh` files
- **ParaView**: UnstructuredGrid cells are 3D volume cells; color by cell arrays
- **Debugging**: Use comprehensive std::cout for method tracing


## Case System

The solver uses `CaseReader` for runtime configuration instead of hard-coded parameters.

### CaseReader Implementation
- **Location**: `include/Case/CaseReader.hpp` and `src/Case/CaseReader.cpp`
- **Parser**: OpenFOAM-style format with nested sections
- **Features**:
  - Type-safe template-based lookups: `lookup<Scalar>("keyword")`
  - Optional parameters with defaults: `lookupOrDefault<bool>("key", false)`
  - Nested sections: `section("sectionName")`
  - Vectors: `(x y z)` format automatically converted to `Vector`
  - Comments: Single-line `//` and multi-line `/* */`

### Case File Structure
The default `defaultCase` file is organized into logical sections:

```cpp
mesh { file path; checkQuality bool; }
physicalProperties { rho scalar; mu scalar; }
initialConditions { U vector; p scalar; }
boundaryConditions { U { patch { type value; } } p { ... } }
numericalSchemes { convection { default scheme; U scheme; k scheme; omega scheme; } }
SIMPLE { numIterations int; convergenceTolerance scalar; relaxationFactors { U scalar; p scalar; k scalar; omega scalar; } }
linearSolvers { U { solver type; preconditioner type; tolerance scalar; maxIter int; } p { ... } }
turbulence { model string; enabled bool; }
output { filename string; }
constraints { velocity { enabled bool; maxVelocity scalar; } pressure { enabled bool; ... } }
```

### Adding New Case Parameters
1. Add entry to appropriate section in `defaultCase`
2. Read in `main.cpp` using `caseReader.lookup<Type>("parameter")`
3. Apply to solver/model as needed

Example:
```cpp
// In defaultCase
SIMPLE
{
    newParameter    0.5;    // New parameter
}

// In main.cpp
auto simpleDict = caseReader.section("SIMPLE");
Scalar newParam = simpleDict.lookup<Scalar>("newParameter");
simpleSolver.setNewParameter(newParam);
```

### Error Handling
- **File not found**: Throws `std::runtime_error` with clear message
- **Parse errors**: Reports file name and line number
- **Type conversion**: Fails gracefully with conversion error messages
- **Missing parameters**: `lookup()` throws, `lookupOrDefault()` uses fallback


## Call flow

```mermaid
flowchart TD
  A[main.cpp / CFDApplication] --> B[readMshFile]
  B --> C[compute face geometry]
  C --> D[compute face distances]
  D --> E[compute cell geometry]
  E --> F[BoundaryConditions + linkFaces]
  F --> G[SIMPLE constructor]
  G --> H[SIMPLE.solve]
  H --> I[compute gradP]
  I --> J[solveMomentumEquations]
  J --> K[calculateRhieChowFlowRate]
  K --> L[solvePressureCorrection]
  L --> M[correctVelocity]
  M --> N[correctFlowRate]
  N --> O[correctPressure]
  O -.-> T[solveTurbulence]
  T --> P[checkConvergence]
  P -->|loop| I
  P --> Q[postProcess + VTK]
```


