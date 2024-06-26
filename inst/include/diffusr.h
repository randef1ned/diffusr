#ifdef INTEL_MKL_VERSION 
#define EIGEN_USE_MKL_ALL
#endif
#define EIGEN_VECTORIZE_SSE4_2

// include standard C++ headers
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <algorithm>
#include <omp.h>
#include <queue>
#include <vector>
#include <set>

// headers in this file are loaded in RcppExports.cpp

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

// other headers are loaded when C++ functions in src/ are being compiled.
// using namespace RcppSparse;
using namespace Rcpp;
using namespace Eigen;
using namespace std;

using Eigen::Map;           // 'maps' rather than copies
using Eigen::MatrixXd;      // variable size matrix, double precision
using Eigen::VectorXd;      // variable size vector, double precision
using Eigen::ArrayXd;
using Eigen::SparseMatrix;

typedef Eigen::MappedSparseMatrix<double> MSpMat;
// typedef Eigen::Map<MatrixXd> MMatrixXd;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<VectorXd> MVectorXd;

#ifndef RCPP_diffusr_H_GEN_
#define RCPP_diffusr_H_GEN_

#include "diffusr_RcppExports.h"

#endif // RCPP_diffusr_H_GEN_
