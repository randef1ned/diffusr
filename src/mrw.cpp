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

template<typename type>
type inverse_matrix(const type& mat, const type &I);

template<typename type> VectorXd mrwr_o(const MatrixXd& p0, const type &W, const double r, const double thresh, const int niter, const bool do_analytical);

template<> SpMat inverse_matrix(const SpMat &mat, const SpMat &I) {
    SparseLU<SpMat> luDecomposition(mat);
    SpMat ret = luDecomposition.solve(I);
    ret.makeCompressed();
    return ret;
}
template<> MatrixXd inverse_matrix(const MatrixXd &mat, const MatrixXd &I) {
    FullPivLU<MatrixXd> luDecomposition(mat);
    return luDecomposition.inverse();
}

// [[Rcpp::plugins("cpp17")]]
template <typename temp>
VectorXd mrwr_t(const MatrixXd& p0,
                const temp& W,
                const double    r,
                const double    thresh,
                const int       niter,
                const bool      do_analytical) {
    MatrixXd pt;
    if (do_analytical) {
        size_t n = W.rows();
        // Equivalence with: MatrixXd I  = MatrixXd::Identity(W.rows(), W.cols());
        temp I(n, n);
        I.setIdentity();
        // Equivalence with: MatrixXd T  = r * (I - (1 - r) * W).inverse();
        temp T = I - (1 - r) * W;
        // then inverse the matrix T
        temp T_inv = inverse_matrix(T, I);
        I.resize(0, 0);

        pt = T_inv * r * p0;
    } else {
        MatrixXd pold;

        int iter = 0;
        pt = p0;
        do {
            if (iter % 25 == 0)
                Rcpp::checkUserInterrupt();
            pold = pt;
            pt   = (1 - r) * W * pold + r * p0;
        }
        while ((pt - pold).norm() > thresh && iter++ < niter);
    }
    VectorXd pt_v(MVectorXd(pt.data(), pt.rows()));
    return pt_v;
}

//' Do a Markon random walk (with restart) on an column-normalised adjacency
//' matrix.
//'
//' @noRd
//' @param p0  matrix of starting distribution
//' @param W  the column normalized adjacency matrix
//' @param r  restart probability
//' @param thresh  threshold to break as soon as new stationary distribution
//'   converges to the stationary distribution of the previous timepoint
//' @param niter  maximum number of iterations for the chain
//' @param do_analytical  boolean if the stationary distribution shall be
//'  computed solving the analytical solution or iteratively
//' @return  returns the matrix of stationary distributions p_inf
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
VectorXd mrwr_(const MatrixXd& p0, const MatrixXd &W, const double r, const double thresh, const int niter, const bool do_analytical) {
    return mrwr_t(p0, W, r, thresh, niter, do_analytical);
}

//' Do a Markon random walk (with restart) on an column-normalised adjacency
//' matrix.
//'
//' @noRd
//' @param p0  matrix of starting distribution
//' @param W  the column normalized adjacency matrix
//' @param r  restart probability
//' @param thresh  threshold to break as soon as new stationary distribution
//'   converges to the stationary distribution of the previous timepoint
//' @param niter  maximum number of iterations for the chain
//' @param do_analytical  boolean if the stationary distribution shall be
//'  computed solving the analytical solution or iteratively
//' @return  returns the matrix of stationary distributions p_inf
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
VectorXd mrwr_s(const MatrixXd &p0, const SpMat &W, const double r, const double thresh, const int niter, const bool do_analytical) {
    return mrwr_t(p0, W, r, thresh, niter, do_analytical);
}
