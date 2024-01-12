# diffusr: network diffusion algorithms in R
#
# Copyright (C) 2016 Simon Dirmeier
#
# This file is part of diffusr.
#
# diffusr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# diffusr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with diffusr. If not, see <http://www.gnu.org/licenses/>.


#' Create a stochastically normalized matrix/vector
#'
#' @export
#'
#' @param obj  \code{\link[base]{matrix}} (or
#'   \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}, \link[base]{vector}) that
#'   is stochstically normalized
#' @param ...  additional params
#' @return  returns the normalized \code{\link[base]{matrix}} (or
#'   \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}, \link[base]{vector})
#'
#' @importFrom checkmate assert check_matrix test_numeric test_atomic_vector
#'
#' @examples
#' W <- matrix(abs(rnorm(10000)), 100, 100)
#' stoch.W <- normalize.stochastic(W)
normalize.stochastic <- function(obj, ...) {
  is_matrix <- FALSE
  if (test_numeric(obj, lower = 0, finite = TRUE, any.missing = FALSE,
                   all.missing = FALSE, null.ok = FALSE) &&
      test_atomic_vector(obj)) {
    if (!.equals.double(sum(obj), 1, .001)) {
      message("normalizing vector!")
    }
  } else if (is.dgCMatrix(obj)) {
    assert_dgCMatrix(obj)
    is_matrix <- TRUE
  } else {
    assert(
      check_matrix(obj, mode = "numeric", any.missing = FALSE,
                   all.missing = FALSE, null.ok = FALSE),
      any(obj >= 0),
      combine = "and"
    )
    is_matrix <- TRUE
  }
  if (is_matrix) {
    sums <- colSums3(obj, !is_matrix)
    if (!all(.equals.double(sums, 1, .001))) {
      message("normalizing column vectors!")
      empt_col_val <- 1.0 / ncol(obj)

      obj <- obj / sums[col(obj)]
      # check if need wipe zeros
      zeros <- which(sums < empt_col_val)
      if (length(zeros)) {
        obj[, zeros] <- 0.00001
      }
    }
  } else {
    obj <- obj / sum(obj)
  }
  return(obj)
}

#' Calculate the Laplacian of a matrix
#'
#' @export
#'
#' @param obj  \code{\link[base]{matrix}} (or
#'   \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}) for which the Laplacian is
#'   calculated
#' @param ...  additional params
#' @return  returns the Laplacian
#'
#' @importFrom checkmate assert check_matrix
#' @importFrom Rcpp sourceCpp
#'
#' @examples
#' W <- matrix(abs(rnorm(10000)), 100, 100)
#' lapl.W <- normalize.laplacian(W)
normalize.laplacian <- function(obj, ...) {
  if (is.dgCMatrix(obj)) {
    assert_dgCMatrix(obj)
    return(laplacian_s(obj))
  } else {
    assert(
      check_matrix(obj, mode = "numeric", nrows = ncol(obj), ncols = nrow(obj),
                   min.rows = 3, any.missing = FALSE, all.missing = FALSE,
                   null.ok = FALSE),
      any(obj >= 0),
      combine = "and"
    )
    return(laplacian_(obj))
  }
}

#' Correct for hubs in an adjacency matrix
#'
#' @export
#'
#' @param obj  \code{\link[base]{matrix}} (or
#'   \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}) for which hubs
#'   are corrected
#'
#' @return  returns the \code{\link[base]{matrix}} (or
#'   \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}) with hub correction
#'
#' @useDynLib diffusr
#'
#' @importFrom checkmate assert check_matrix
#' @importFrom Rcpp sourceCpp
#'
#' @examples
#' W <- matrix(abs(rnorm(10000)), 100, 100)
#' cor.hub <- hub.correction(W)
hub.correction <- function(obj) {
  n_elements <- nrow(obj)
  if (is.dgCMatrix(obj)) {
    assert_dgCMatrix(obj)
    return(hub_normalize_s(obj))
  } else {
    assert(
      check_matrix(obj, mode = "numeric", min.rows = 3, nrows = n_elements,
                   ncols = n_elements, any.missing = FALSE, all.missing = FALSE,
                   null.ok = FALSE),
      any(obj >= 0),
      combine = "and"
    )
    return(hub_normalize_(obj))
  }
}

#' @noRd
#' @importFrom igraph graph_from_adjacency_matrix components
.is.ergodic <- function(obj) {
  adj   <- graph_from_adjacency_matrix(obj, mode = "directed", weighted = TRUE)
  comps <- components(adj)
  ifelse(length(comps$csize) == 1, TRUE, FALSE)
}

#' @noRd
#' @importFrom matrixStats colSums2
colSums3 <- function(mat, is.sparse = NULL) {
  if (is.null(is.sparse)) {
    is.sparse <- is.dgCMatrix(mat)
  }
  if (is.sparse) {
    sums <- sparseMatrixStats::colSums2(mat)
  } else {
    sums <- colSums2(mat)
  }
}
