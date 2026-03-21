# C++ Code Style Guide

Style conventions for Turblyze (C++20). These rules apply to all `.hpp`
and `.cpp` files in `include/` and `src/`.

## File Headers

**Header files (.hpp)** use the extended form with `@class` and description:
```cpp
/******************************************************************************
 * @file BoundaryPatch.hpp
 * @brief Boundary patch representation and mesh connectivity management
 *
 * This header defines the BoundaryPatch class, which represents a set of
 * faces on the domain boundary. A patch is identified by a name
 * (e.g., "inlet", "wall") and a geometric zone ID from the mesh file.
 *
 * @class BoundaryPatch
 *
 * The BoundaryPatch class provides:
 * - Identification of boundary zones (name, ID, type)
 * - Topological range definitions (start face index, end face index)
 * - Mapping between mesh file types (e.g., Fluent strings) and internal enums
 * - Helper methods for querying patch size and face validity
 *****************************************************************************/
```

**Source files (.cpp)** use the shorter form:
```cpp
/******************************************************************************
 * @file FileName.cpp
 * @brief Brief description of file purpose
 *****************************************************************************/
```

## Section Separators

**Header files (.hpp)** use plain comment separators inside classes:
```cpp
// Setter methods

// Accessor methods

// Private members

// Private methods
```

**Source files (.cpp)** use asterisk-decorated separators for major sections:
```cpp
// ****************************** Section Name ******************************
```

## Documentation Style

All documentation lives in **headers only**. Source files have no Doxygen on method implementations.

**In headers (.hpp)** — full Doxygen `/** */` blocks on every method declaration:
```cpp
/**
 * @brief Calculate boundary face value for scalar field
 * @param face Boundary face
 * @param phi Scalar field
 * @param fieldName Name of the field
 * @return Boundary value based on boundary condition
 * @throws std::runtime_error if face not found in boundary patches
 */
Scalar calculateBoundaryFaceValue
(
    const Face& face,
    const ScalarField& phi,
    const std::string& fieldName
) const;
```

Use `///` for trivial one-liners (simple constructors, member variables, and methods with no return value and no parameters):
```cpp
/// Default constructor
BoundaryConditions() = default;

/// Print summary of all boundary conditions
void printSummary() const;

/// Nested map: patch name → field name → boundary data
std::map<std::string, std::map<std::string, BoundaryData>> patchBoundaryData_;
```

Use `@details` for extended method explanations, starting on a new line:
```cpp
/**
 * @brief Calculate mass fluxes using Rhie-Chow interpolation
 * @param flowRate Face volume flow rates
 *
 * @details
 * Implements the Rhie-Chow interpolation method to prevent
 * checkerboard pressure oscillations in collocated grids.
 */
```

**In source files (.cpp)** — no Doxygen on implementations; use inline `//` comments only where logic needs explanation:
```cpp
// CASE 1: Face is "Triangle" (numNodes == 3)
if (numNodes == 3)
{
    // For planar triangles, contact = projected
    contactArea_ = projectedArea_;
}
```

## Function Signature Formatting
Multi-parameter functions use Allman-style with each parameter on its own line:
```cpp
bool setFixedValue
(
    const std::string& patchName,
    const std::string& fieldName,
    Scalar value
);
```

Single-parameter or short-signature functions can stay on one line:
```cpp
void addPatch(const BoundaryPatch& patch);
```

## Brace Style
Opening brace on a new line (Allman style) for classes, functions, and control flow:
```cpp
class BoundaryConditions
{
public:
    // ...
};

if (condition)
{
    // ...
}
else
{
    // ...
}
```

## Error Throwing
`throw` on its own line, exception type and parenthesized message indented below, Allman-style:
```cpp
throw
    std::runtime_error
    (
        "Error in Cell " + std::to_string(idx_)
      + " calculation: something went wrong."
    );
```

## Return Statement Formatting
For long return statements that exceed 80 characters or span multiple lines, place the expression on a new line after `return`, indented:
```cpp
return
    S_11*S_11 + S_22*S_22 + S_33*S_33
  + 2.0 * (S_12*S_12 + S_13*S_13 + S_23*S_23);
```

Short return statements stay on one line:
```cpp
return phi[face.ownerCell()];
```

## Line Length
Lines should not exceed 80 characters. Break long lines using the Allman-style parameter formatting, string concatenation with `+`, or continuation on the next line.

## Naming Conventions
- **Classes**: PascalCase (e.g., `LinearSolver`, `KOmegaSST`)
- **Methods**: camelCase (e.g., `solveMomentum`, `computeGradient`)
- **Member variables**: camelCase with trailing underscore (e.g., `tolerance_`, `fieldName_`)
- **Local variables**: camelCase (e.g., `cellVolume`, `faceArea`)
- **Type aliases**: PascalCase (e.g., `VectorField`, `ScalarField`)

## Stream Output Alignment

All stream insertion operators (`<<`) must align at **column 4** (position 4 from line start), following OpenFOAM conventions. This creates visual alignment of the `<<` symbols themselves.

**Core Rule:** The `<<` operator always appears at column 4, regardless of the stream object name length.

**Short stream objects (≤4 chars):**
```cpp
Info<< "Message" << value << std::endl;
os  << "Data: " << x << ", " << y << std::endl;
log << "Status" << std::endl;
```

**Longer stream objects:**
When the stream object name plus spacing to column 4 would exceed the line, place `<<` on the next line at column 4:
```cpp
std::cout
    << "Message" << value << std::endl;

std::cerr
    << "Error: " << errorMsg << std::endl;

WarningInFunction
    << "Warning message"
    << std::endl;
```

**Multi-line continuations:**
Continuation lines indent to column 4 for the `<<` operator:
```cpp
std::cout
    << "Long message part 1: " << value1
    << " part 2: " << value2
    << std::endl;
```

**In operator<< overloads:**
```cpp
std::ostream& operator<<(std::ostream& os, const Vector& v)
{
    os  << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return os;
}
```

**With format manipulators:**
```cpp
std::cout
    << std::scientific << std::setprecision(6)
    << value << " m²"
    << std::endl;
```

**Rationale:** This OpenFOAM-style alignment provides visual consistency where all `<<` operators align vertically at column 4, making stream operations easy to identify while keeping code reasonably compact.
