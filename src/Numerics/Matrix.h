#ifndef MATRIXCONSTRUCTOR_H
#define MATRIXCONSTRUCTOR_H

#include <vector>
#include <string>
#include <map>

#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Dense>

#include "Cell.h"
#include "Face.h"
#include "BoundaryConditions.h"
#include "CellData.h"
#include "ConvectionScheme.h"
#include "GradientScheme.h"

class KOmegaSST;

// Enum to select the desired time scheme
enum class TimeScheme {
    Steady,         // For steady-state problems (no time term)
    Transient       // For transient problems
};

class Matrix {
public:
    // References to mesh data and schemes
    const std::vector<Face>& allFaces;
    const std::vector<Cell>& allCells;
    const BoundaryConditions& bcManager;
    const GradientScheme& gradScheme;

    // Gradient fields
    VectorField gradP;
    VectorField gradUx;
    VectorField gradUy;
    VectorField gradUz;
    VectorField gradk;
    VectorField gradOmega;

    // Face-level data
    FaceFluxField mdotFaces;
    FaceVectorField gradP_f;
    FaceVectorField gradUx_f;
    FaceVectorField gradUy_f;
    FaceVectorField gradUz_f;
    FaceVectorField gradk_f;
    FaceVectorField gradOmega_f;

    // Flag to indicate caches are up-to-date for current iteration
    bool cachesValid = false;

    // The matrix (A) and vector (b) that represent the linear system Ax=b
    Eigen::SparseMatrix<Scalar> A_matrix;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> b_vector;

    // Constructor
    Matrix(
        const std::vector<Face>& faces,
        const std::vector<Cell>& cells,
        const BoundaryConditions& boundaryConds,
        const GradientScheme& gradientScheme
    );

    void refreshIterationCaches(const PressureField& p, const VelocityField& U, Scalar rho, const KOmegaSST* turbulenceModel);

    void constructScalarTransportMatrix(
        const std::string& fieldName,
        const ScalarField& phi,
        const ScalarField& phi_old,
        const VectorField& U_field,
        const ScalarField& phi_source,
        Scalar rho,
        Scalar Gamma,
        TimeScheme timeScheme,
        Scalar dt,
        Scalar theta,
        const ConvectionDiscretization& convScheme
    );

    const Eigen::SparseMatrix<Scalar>& getMatrixA() const { return A_matrix; }
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& getVectorB() const { return b_vector; }

private:
    std::vector<Eigen::Triplet<Scalar>> tripletList;
    std::map<size_t, const BoundaryPatch*> faceToPatchMap;

    void clear();
};

#endif