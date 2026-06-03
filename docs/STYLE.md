# C++ Code Style Guide

Style conventions for Turblyze (C++20). These rules apply to all `.h`
and `.cpp` files under `src/`.

## File Headers

**Header files (.h)** use the extended form with `@class` and description:
```cpp
/******************************************************************************
 * @file BoundaryPatch.h
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

Both `.h` and `.cpp` files divide their contents with **79-character asterisk
banners**. The label is centred in a field of `*`:
```cpp
// ****************************** Setter Methods ******************************
```

**Exact geometry** (so every banner lines up identically):
- The line is `// ` + left `*`s + ` ` + label + ` ` + right `*`s, **exactly 79
  characters** total (never exceed 80).
- After the `// ` prefix, 76 columns remain. With
  `stars = 76 - len(label) - 2`, split them `left = (stars + 1) / 2` and
  `right = stars - left` — so on an odd split the **extra `*` goes before the
  label**.
- A blank line precedes and follows every banner.

### Headers banner

The first banner in every file (after `#pragma once` in a header) is `Headers`.
It groups the includes with `//` sub-comments — only the groups that apply, in
this order:
```cpp
// ********************************** Headers *********************************

// Implementation header        (.cpp only — the matching .h)
#include "MeshReader.h"

// Standard library headers
#include <vector>

// External library headers     (Eigen, OpenMP, VTK)
#include <eigen3/Eigen/SparseCore>

// Project headers
#include "ErrorHandler.h"
```
A file whose includes form a single group needs no sub-comments. A header with
no includes at all (e.g. one holding only forward declarations) has no `Headers`
banner.

### Type and section banners

After `Headers`, banner each top-level definition and each member group:
```cpp
// ******************************* class Vector *******************************
// ************************* struct TransportEquation *************************
// ***************************** enum class Field *****************************
// *************************** concept CellFieldType **************************
// ******************************* namespace VTK ******************************
```
Each definition gets its own banner — a `concept` that constrains a class is
banner-separated from the `class` it precedes.

Inside a class, the canonical member-group labels are:

| Label | Covers |
|---|---|
| `Special Member Functions` | constructors, the rule-of-five block, destructor |
| `Setter Methods` / `Accessor Methods` | mutators / getters |
| `Operator Methods` | in-class `operator…` overloads |
| `Public Methods` | other public members |
| `Protected Methods` | the `protected:` section |
| `Private Methods` | private member functions |
| `Private Members` | private data |
| `Aliases` | trailing `using` type aliases |
| `Non-Member Functions` | free functions declared alongside the type |

`Private Members` / `Private Methods` (and `Protected Members` / `Protected
Methods`) banners sit **before** the access specifier, not after it:
```cpp
// ****************************** Private Members *****************************

private:

    Scalar x_ = S(0.0);
```
When a class has both a private-data group and a private-method group, **repeat
the `private:` specifier** so each banner still precedes one:
```cpp
// ****************************** Private Members *****************************

private:

    Scalar x_ = S(0.0);

// ****************************** Private Methods *****************************

private:

    void helper() const;
```
This differs from the `public:` section, which opens immediately after the class
brace and may carry several mid-section banners (`Special Member Functions`,
`Setter Methods`, `Accessor Methods`, …) under a single `public:`.

A `.cpp` reuses the **same** group names as its header, but only for the
sections it actually defines — typically `Special Member Functions` plus
`Private Methods`; inline accessors and data members never reappear in the
`.cpp`.

### "Methods" vs "Functions"

`Methods` is reserved for members **inside** a class or struct. Free functions
are never "Methods", and never "Public" (there is no access specifier for them
to be public *to*):
- A free function tied to a type or enum → `Non-Member Functions`.
- A header of standalone free functions → name it for what the functions do,
  suffixed with `Functions`: e.g. `Error Handling Functions`,
  `Interpolation Functions`.

## Documentation Style

All documentation lives in **headers only**. Source files have no Doxygen on method implementations.

**In headers (.h)** — use `///` for all method and function declarations:
```cpp
/// Calculate boundary face value for a scalar field
[[nodiscard]] Scalar boundaryFaceValue
(
    Field field,
    const ScalarField& phi,
    const Face& face
) const;

/// Default constructor
BoundaryConditions() = default;

/// Print summary of all boundary conditions
void printSummary() const;

/// Nested map: patch name → field → boundary data
std::map<Name, std::map<Field, BoundaryData>> patchBoundaryData_;
```

Use multiple `///` lines only when a non-obvious invariant or formula genuinely cannot fit on one line:
```cpp
/// Symmetric second-moment polynomial for triangle integration
/// Evaluates a² + b² + c² + ab + ac + bc
/// ∫∫_triangle x² dA = (area / 6) × secondMoment(x₁, x₂, x₃)
[[nodiscard]] static Scalar secondMoment(Scalar a, Scalar b, Scalar c);
```

Do not add trailing periods to `///` method or function doc comments — they are labels, not
sentences. Full prose sentences on member-variable comments may end with a period when
explaining a non-obvious constraint. Do not use `@param`, `@return`, or `@note` tags on
individual method declarations — the function signature and name already carry that information.
The full `/** @brief … @param … @return */` Doxygen form is reserved for the file-level
`@file`/`@class` block at the top of each header.

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
void setFixedValue
(
    const Name& patchName,
    Field field,
    Scalar value
);
```

Single-parameter or short-signature functions can stay on one line:
```cpp
void addPatch(BoundaryPatch patch);
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

## Error Handling
Two macros from `ErrorHandler.h` for all error reporting:

**Fatal errors** (unrecoverable — prints message with file/line and aborts):
```cpp
FatalError
(
    "Cell " + std::to_string(idx_)
  + " calculation: something went wrong."
);
```

**Warnings** (non-fatal — prints message with file/line and continues):
```cpp
Warning
(
    "Degenerate least-squares matrix in "
  + std::to_string(count) + " cells"
);
```

Do not use `throw`, `assert()`, or raw `std::cerr` for error reporting.

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

## Character Literals
Use `'\n'`, `':'`, `' '`, etc. for single characters — not `"\n"`, `":"`, `" "`.
Single-character string literals carry an unnecessary null terminator and are
less semantically precise:
```cpp
// Correct
std::cerr << '\n' << "Error: " << msg << '\n';

// Wrong
std::cerr << "\n" << "Error: " << msg << "\n";
```

## Naming Conventions
- **Classes**: PascalCase (e.g., `LinearSolver`, `GradientScheme`)
- **Methods**: camelCase (e.g., `solveMomentum`, `computeGradient`)
- **Member variables**: camelCase with trailing underscore (e.g., `tolerance_`, `fieldName_`)
- **Local variables**: camelCase (e.g., `cellVolume`, `faceArea`)
- **Type aliases**: PascalCase (e.g., `VectorField`, `ScalarField`)
- **Enum class names**: PascalCase (e.g., `Field`, `BCType`, `PatchType`)
- **Enumerators**: lowerCamelCase (e.g., `Field::Ux`, `BCType::fixedValue`,
  `PatchType::wall`). Avoid `ALL_CAPS` — it collides with preprocessor macros
  (C++ Core Guidelines Enum.5). Prefer `enum class` over plain `enum` (Enum.3).

### Intent-revealing aliases

Prefer the foundation aliases over bare standard-library types when the *role*
of the value is meaningful — they signal intent to the reader (the compiler
still sees the underlying type, so they are documentation, not type safety):

- **`Integer.h`** — `Index` (addresses an element) and `Count` (a size or
  quantity), both `std::size_t`; plus `IndexList`/`CountList` and the
  `IndexListRef` (`std::span<const Index>`) view.
- **`StringTypes.h`** — owned text `Name` / `Token` / `FilePath` / `Message`
  (all `std::string`) and their borrowed `std::string_view` views with the
  `*Ref` suffix (`NameRef`, `TokenRef`, `FilePathRef`, `MessageRef`).
- **`MeshContainers.h`** — owning `NodeList`/`FaceList`/`CellList`/`PatchList`
  and the borrowed `*Ref` span views (`FaceListRef`, `MutableFaceListRef`, …).

The `*Ref` suffix marks a non-owning view in the name. Domain-narrow aliases
(`FaceIndex`, `PatchName`) are intentionally **not** used — a single `Index` /
`Name` keeps the vocabulary small. Local aliases (e.g. `Face::OptionalIndex`,
`CaseReader::EntryMap`) live next to the class that needs them.

## Special Member Functions

### Declaration order

Always declare special members in this order inside the class `public` section:

```cpp
/// Copy constructor and assignment — <reason>
ClassName(const ClassName&) = delete;
ClassName& operator=(const ClassName&) = delete;

/// Move constructor and assignment — <reason>
ClassName(ClassName&&) = delete;
ClassName& operator=(ClassName&&) = delete;

/// Destructor
~ClassName() noexcept = default;
```

The brief reason on the `///` comment documents *why* the operation is restricted:
- `Not copyable (const T& members)` — reference members cannot be rebound
- `Not movable (const T& members)` — same
- `Not movable (self-referential Eigen solver)` — unsafe default move

### Choosing the right rule

| Member type | Copy | Move | Destructor |
|---|---|---|---|
| `const T&` or `T&` reference | `= delete` | `= delete` | `= default` (noexcept) |
| `std::unique_ptr<T>` | `= delete` | `= default` | `= default` (noexcept) |
| Eigen iterative solver | `= delete` | `= delete` | `= default` (noexcept) |
| No non-trivial members | *(declare none — rule of zero)* | | |

**Rule of zero**: If the compiler-generated defaults are correct, declare nothing.
Applies to value-only classes (`Face`, `Cell`, `BoundaryData`).

**Rule of five**: Whenever one special member is declared, declare all five.
Declaring only some leaves the class in a partially specified state that
surprises readers and may silently break under future refactors.

### `{}` vs `()` in initializer lists

Prefer `{}` over `()` for member initialization in constructor initializer lists.
`{}` prevents narrowing conversions and makes intent explicit:

```cpp
// Preferred
ClassName(int n, Scalar tol)
:
    count_{n},
    tolerance_{tol}
{}

// Acceptable only when initializer_list constructor would be selected instead
// (e.g. std::vector fill constructor: vector(n, value) ≠ vector{n, value})
internalField_(Mesh::cellCount(), T{})
```

Initializer list order must match member declaration order exactly.
Members absent from the list are initialized in declaration order before
the listed ones — omitting a member is not an order violation, but all
*listed* members must appear in the same relative sequence as their declarations.

## Stream Output Alignment

All stream insertion operators (`<<`) must align at **column 4** (position 4 from line start), following OpenFOAM conventions. This creates visual alignment of the `<<` symbols themselves.

**Core Rule:** The `<<` operator always appears at column 4, regardless of the stream object name length.

**Short stream objects (≤4 chars):**
```cpp
Info<< "Message" << value << '\n';
os  << "Data: " << x << ", " << y << '\n';
log << "Status" << '\n';
```

**Longer stream objects:**
When the stream object name plus spacing to column 4 would exceed the line, place `<<` on the next line at column 4:
```cpp
std::cout
    << "Message" << value << '\n';

std::cerr
    << "Error: " << errorMsg << std::endl; // Explicit flush for diagnostics

WarningInFunction
    << "Warning message"
    << std::endl; // Explicit flush for diagnostics
```

**Multi-line continuations:**
Continuation lines indent to column 4 for the `<<` operator:
```cpp
std::cout
    << "Long message part 1: " << value1
    << " part 2: " << value2 << '\n';
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
    << value << " m²" << '\n';
```

**Newlines:** Use `'\n'` for ordinary line endings. Use `std::endl`
only when the flush is intentional, such as fatal errors or rare diagnostic
messages that must be visible immediately before termination.

```cpp
// Preferred for one ordinary line
std::cout
    << "SIMPLE Loop" << '\n';
```

**Rationale:** This OpenFOAM-style alignment provides visual consistency where all `<<` operators align vertically at column 4, making stream operations easy to identify while keeping code reasonably compact.
