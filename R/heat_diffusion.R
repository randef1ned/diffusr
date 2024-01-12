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


#' Graph diffusion using a heat diffusion process on a Laplacian matrix.
#'
#' @description An amount of starting heat gets distribution using the
#' Laplacian matrix of a graph. Every iteration (or time interval) \code{t}
#'heat streams from the starting nodes into surrounding nodes.
#'
#' @export
#'
#' @param h0   an \code{n x p}-dimensional numeric non-negative vector/matrix
#'  of starting temperatures
#' @param graph  an (\code{n x n})-dimensional numeric non-negative adjacence
#'  matrix representing the graph
#' @param t  time point when heat is measured
#' @param ...  additional parameters
#' @return  returns the heat on every node as numeric vector
#'
#' @useDynLib diffusr
#'
#' @importFrom checkmate assert_int assert_integer assert test_matrix
#'                       test_numeric test_atomic_vector
#' @importFrom Rcpp sourceCpp
#'
#' @references
#' \url{https://en.wikipedia.org/wiki/Laplacian_matrix} \cr
#' \url{https://en.wikipedia.org/wiki/Heat_equation}
#'
#' @examples
#' # count of nodes
#' n <- 5
#' # starting distribution (has to sum to one)
#' h0 <- as.vector(rmultinom(1, 1, prob=rep(.2, n)))
#' # adjacency matrix (either normalized or not)
#' graph <- matrix(abs(rnorm(n*n)), n, n)
#' # computation of stationary distribution
#' heat <- heat.diffusion(h0, graph)
heat.diffusion <- function(h0, graph, t = 0.5, ...) {
  assert_number(t, lower = 0, na.ok = FALSE, finite = TRUE, null.ok = FALSE)
  n_elements <- nrow(graph)

  # convert h0 if h0 is vector
  if (test_numeric(h0, lower = 0, len = n_elements, finite = TRUE,
                   any.missing = FALSE, all.missing = FALSE, null.ok = FALSE) &&
      test_atomic_vector(h0, len = n_elements)) {
    h0 <- as.matrix(h0)
  } else {
    assert(
      test_matrix(h0, mode = "numeric", nrows = n_elements, any.missing = FALSE,
                  all.missing = FALSE, null.ok = FALSE),
      any(h0 >= 0),
      combine = "and"
    )
  }

  # graph must be either matrix or dgCMatrix
  diag(graph) <- 0
  if (is.dgCMatrix(graph)) {
    assert_dgCMatrix(graph)
  } else {
    assert(
      test_matrix(graph, mode = "numeric", nrows = n_elements,
                  ncols = n_elements, min.rows = 3, any.missing = FALSE,
                  all.missing = FALSE, null.ok = FALSE),
      any(graph >= 0),
      combine = "and"
    )
  }
  heat <- heat_diffusion_(h0, laplacian_(graph), t)
  return(heat)
}
