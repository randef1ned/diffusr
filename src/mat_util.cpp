/** diffusr: network diffusion algorithms in R
 *
 * Copyright (C) 2016 Simon Dirmeier
 * @author Simon Dirmeier
 * @email simon.dirmeier@bsse.ethz.ch
 *
 * This file is part of diffusr.
 *
 * diffusr is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * diffusr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with diffusr. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../inst/include/diffusr.h"

template <typename T>
T stoch_col_norm_t(const T& W) {
    size_t n = W.rows();
    T res(n, n);
    VectorXd colsums      = W.transpose() * VectorXd::Ones(n);
    //res.reserve((colsums * n).sum());
    const double    empt_col_val = 1.0 / n;
    const double    zero_col     = 0.00001;
    
    for (unsigned int i = 0; i < n; ++i) {
        if (colsums[i] <= zero_col) {
            for (size_t j = 0; j < n; ++j) {
                res.coeffRef(j, i) = empt_col_val;
            }
        } else {
            for (size_t j = 0; j < n; ++j) {
                res.coeffRef(j, i) = W.coeff(j, i) / colsums(i);
            }
        }
    }

    return res;
}

//' Column normalize a matrix, so that it is stochastic.
//'
//' @noRd
//' @param W  the adjacency matrix to be normalized
//' @return  returns the normalized matrix
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
MatrixXd stoch_col_norm_(const MatrixXd &W) {
    return stoch_col_norm_t(W);
}

//' Column normalize a matrix, so that it is stochastic.
//'
//' @noRd
//' @param W  the adjacency matrix to be normalized
//' @return  returns the normalized matrix
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
SpMat stoch_col_norm_s(const SpMat &W) {
    SpMat res = stoch_col_norm_t(W);
    res.makeCompressed();
    return res;
}

template <typename T> MatrixXd compute_laplacian(const T &W, const ArrayXd rowsums_m, const size_t P);

template <> MatrixXd compute_laplacian(const MatrixXd &W, const ArrayXd rowsums, const size_t P) {
    ArrayXd rowsums_m = rowsums.replicate(1, P).array().sqrt().inverse();
    // remove 1/0 = Inf
    rowsums_m = (rowsums_m.isInf()).select(0, rowsums_m);
    MatrixXd res = - W.array() * (rowsums_m.transpose() * rowsums_m) * W.array().cast<bool>().cast<double>();
    return res;
}

template <> MatrixXd compute_laplacian(const SpMat &W, const ArrayXd rowsums, const size_t P) {
    VectorXd rowsums_m = rowsums.replicate(1, P).cwiseSqrt().cwiseInverse();
    MatrixXd res = rowsums_m.transpose().cwiseProduct(rowsums_m).cwiseProduct(-W).cwiseProduct(W.cast<bool>().cast<double>());
    return res;
}

template <typename T> MatrixXd laplacian_t(const T& W) {
    const size_t       P = W.rows();
    //MatrixXd res(P, P);
    // Equivalent with: VectorXd rowsums = W.rowwise().sum();
    ArrayXd rowsums = (W * VectorXd::Ones(P)).array();
    
    // for diagnoals
    ArrayXd W_diags = 1 - W.diagonal().array() / (rowsums == 0).select(1, rowsums);
    W_diags *= (rowsums != 0.0).cast<double>();
    
    // for others
    MatrixXd res = compute_laplacian(W, rowsums, P);
    res.diagonal().array() = W_diags;
    return res;
}

//' Calculate the Laplacian of a weighted matrix.
//'
//' @noRd
//' @param W  the adjacency matrix for which the Laplacian is calculated
//' @return  returns the Laplacian of a matrix
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
MatrixXd laplacian_(const MatrixXd &W) {
    return laplacian_t(W);
}

//' Calculate the Laplacian of a weighted matrix.
//'
//' @noRd
//' @param W  the adjacency matrix for which the Laplacian is calculated
//' @return  returns the Laplacian of a matrix
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
MatrixXd laplacian_s(const SpMat &W) {
    return laplacian_t(W);
}

template <typename T>
VectorXd node_degrees_t(const T &W) {
    MatrixXd Wi = W.template cast<bool>().template cast<double>();
    VectorXd node_degrees = Wi.rowwise().sum();

    return node_degrees;
}

//' Get the unweighted node degrees of a adjacency matrix.
//'
//' @noRd
//' @param W  the adjacency matrix to be normalized
//' @return  returns the node degrees as vectors.
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
VectorXd node_degrees_(const MatrixXd &W) {
    return node_degrees_t(W);
}

//' Get the unweighted node degrees of a adjacency matrix.
//'
//' @noRd
//' @param W  the adjacency matrix to be normalized
//' @return  returns the node degrees as vectors.
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
VectorXd node_degrees_s(const MSpMat &W) {
    return node_degrees_t(W);
}

template <typename temp>
MatrixXd hub_normalize_t(const temp &W) {
    // auto &fc(identity_);
    size_t n = W.rows();
    // Equivalent with: MatrixXd res = MatrixXd::Constant(W.rows(), W.cols(), 0.0);
    temp res = temp(n, n);
    
    VectorXd node_degrees = node_degrees_t(W);
    
    #pragma omp parallel for
    for (unsigned int i = 0; i < n; ++i) {
        double mh;
        for (unsigned int j = 0; j < n; ++j) {
            if (W.coeff(i, j) != 0) {
                // Equivlant with: double mh = fc(node_degrees[i] / node_degrees[j]);
                mh = node_degrees[i] / node_degrees[j];
                res.coeffRef(i, j) = min(1.0, mh) / node_degrees[i];
            }
        }
    }
    return res;
}

//' Normalize the hub bias in a matrix.
//'
//' @noRd
//' @param W  the adjacency matrix to be normalized
//' @return  returns the hub corrected matrix
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
MatrixXd hub_normalize_(const MatrixXd &W) {
    return hub_normalize_t(W);
}

//' Normalize the hub bias in a matrix.
//'
//' @noRd
//' @param W  the adjacency matrix to be normalized
//' @return  returns the hub corrected matrix
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
MatrixXd hub_normalize_s(const SpMat &W) {
    return hub_normalize_t(W);
}
