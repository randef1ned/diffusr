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


#' Graph diffusion using nearest neighbors
#'
#' @description
#' For every node in a set of nodes the graph gets traversed along the node's
#' shortest paths to its neighbors. Nearest neighbors are added until a maximum
#' depth of \code{k} is reached. For settings where there are more than \code{k}
#' neighbors having the same distance, all neighbors are returned.
#'
#' @export
#'
#' @param nodes a \eqn{n}-dimensional integer vector of node indexes (1-based)
#'   for which the algorithm is applied iteratively
#'
#' @param graph an (\eqn{n \times n})-dimensional numeric non-negative adjacence
#'   \code{\link[base]{matrix}} (or
#'   \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}, \link[base]{vector})
#'   representing the graph
#'
#' @param k the depth of the nearest neighbor search, e.g. the depth of the
#'   graph traversal
#'
#' @param ... additional parameters
#'
#' @return returns the kNN nodes as list of integer vectors of node indexes
#'
#' @useDynLib diffusr
#'
#' @importFrom checkmate assert_int assert_integer assert test_matrix
#' @importFrom Rcpp sourceCpp
#'
#' @examples
#' # count of nodes
#' n <- 10
#' # indexes (integer) of nodes for which neighbors should be searched
#' node.idxs <- c(1L, 5L)
#' # the adjaceny matrix (does not need to be symmetric)
#' graph <- rbind(cbind(0, diag(n-1)), 0)
#' # compute the neighbors until depth 3
#' neighs <- nearest.neighbors(node.idxs, graph, 3)
nearest.neighbors <- function(nodes, graph, k = 1L, ...) {
  ## Check the fucking inputs
  n_elements <- nrow(graph)
  assert_int(k, lower = 1, upper = n_elements, na.ok = FALSE, coerce = TRUE,
             null.ok = FALSE)
  assert_integer(nodes, lower = 1, upper = n_elements, max.len = n_elements,
                 any.missing = FALSE, all.missing = FALSE, null.ok = FALSE)
  nodes <- unique(nodes)

  # graph must be either matrix or dgCMatrix
  diag(graph) <- 0
  if (is.dgCMatrix(graph)) {
    assert_dgCMatrix(graph)
    neighbors <- neighbors_s(nodes, graph, k)
  } else {
    assert(
      test_matrix(graph, mode = "numeric", min.rows = 3, nrows = n_elements,
                  ncols = n_elements, any.missing = FALSE, all.missing = FALSE,
                  null.ok = FALSE),
      any(graph >= 0),
      combine = "and"
    )
    neighbors <- neighbors_(nodes, graph, k)
  }
  names(neighbors) <- as.character(nodes)
  return(neighbors)
}
