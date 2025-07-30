#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <vector>
#include <iostream>
#include <cmath>
#include "Scalar.h"

namespace LinearSolvers {

bool BiCGSTAB(Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x,
              const Eigen::SparseMatrix<Scalar>& A,
              const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& B,
              Scalar tolerance,
              int max_iterations,
              const std::string& fieldName);

}
#endif