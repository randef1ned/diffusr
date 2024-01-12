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


#' Graph diffusion using a Markov random walk
#'
#' @description
#' A Markov Random Walk takes an inital distribution \eqn{p_0}
#' and calculates the stationary distribution of that.
#' The diffusion process is regulated by a restart probability \eqn{r} which
#' controls how often the MRW jumps back to the initial values.
#'
#' @param p0  an \eqn{n \times p}-dimensional numeric non-negative vector/matrix
#'  representing the starting distribution of the Markov chain
#'  (does not need to sum to one).
#'
#' @param graph  an \eqn{n \times p}-dimensional numeric non-negative adjacence
#'   \code{\link[base]{matrix}} (or
#'   \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}) representing the graph
#'
#' @param r  a scalar between \eqn{(0, 1)}. restart probability if a Markov
#'   random walk with restart is desired
#'
#' @param thresh  threshold for breaking the iterative computation of the
#'  stationary distribution. If the absolute difference of the distribution at
#'  time point \eqn{t-1} and \eqn{t} is less than \code{thresh}, then the
#'  algorithm stops. If \code{thresh} is not reached before \code{niter}, then
#'  the algorithm stops as well.
#'
#' @param niter  maximal number of iterations for computation of the
#'  Markov chain. If \code{thresh} is not reached, then \code{niter} is used as
#'  stop criterion.
#'
#' @param do.analytical  boolean if the stationary distribution shall be
#'  computed solving the analytical solution or rather iteratively
#'
#' @param correct.for.hubs if \code{TRUE} multiplies a correction factor to the
#'  nodes, such that the random walk gets not biased to nodes with high
#'  degree. In that case the original input matrix will be normalized as:
#'  \deqn{ P(j | i) = \dfrac{1}{\text{degree}(i)} \times
#'    \min \left(1, \dfrac{\text{degree}(j)}{\text{degree}(i)}\right)}
#'  \emph{Note that this will not consider edge weights.}
#'
#' @param allow.ergodic Allow multiple components in a graph.
#'
#' @param return.pt.only Return pt only.
#'
#' @return  returns a list with the following elements
#'  \itemize{
#'   \item \code{p.inf}  the stationary distribution as numeric vector
#'   \item \code{transition.matrix} the column normalized transition matrix used
#'         for the random walk
#'  }
#'
#' @references
#' Tong, H., Faloutsos, C., & Pan, J. Y. (2006),
#' Fast random walk with restart and its applications.
#'
#' Koehler, S., Bauer, S., Horn, D., & Robinson, P. N. (2008),
#' Walking the interactome for prioritization of candidate disease genes.
#' \emph{The American Journal of Human Genetics}
#'
#' @export
#'
#' @useDynLib diffusr
#'
#' @importFrom checkmate assert_number assert_int assert_logical test_numeric
#'                       assert test_matrix check_numeric
#' @importFrom Rcpp sourceCpp
#'
#' @examples
#' # count of nodes
#' n <- 5
#' # starting distribution (has to sum to one)
#' p0    <- as.vector(rmultinom(1, 1, prob=rep(.2, n)))
#' # adjacency matrix (either normalized or not)
#' graph <- matrix(abs(rnorm(n*n)), n, n)
#' # computation of stationary distribution
#' pt    <- random.walk(p0, graph)
#'
random.walk <- function(p0, graph, r = 0.5, niter = 1e4, thresh = 1e-4,
                        do.analytical = FALSE, correct.for.hubs = FALSE,
                        allow.ergodic = FALSE, return.pt.only = FALSE) {
  ## Check the fucking inputs
  assert_number(r, lower = 0, upper = 1, na.ok = FALSE, finite = TRUE,
                null.ok = FALSE)
  assert_int(niter, lower = 2, na.ok = FALSE, coerce = TRUE, null.ok = FALSE)
  assert_number(thresh, lower = 0, na.ok = FALSE, finite = TRUE,
                null.ok = FALSE)
  assert_logical(do.analytical, len = 1, any.missing = FALSE,
                 all.missing = FALSE, null.ok = FALSE)
  assert_logical(correct.for.hubs, len = 1, any.missing = FALSE,
                 all.missing = FALSE, null.ok = FALSE)

  # graph must be either matrix or dgCMatrix
  n_elements <- nrow(graph)
  if (is.dgCMatrix(graph)) {
    assert_dgCMatrix(graph)
    sparse <- TRUE
  } else {
    assert(
      test_matrix(graph, mode = "numeric", min.rows = 3, nrows = n_elements,
                  ncols = n_elements, any.missing = FALSE, all.missing = FALSE,
                  null.ok = FALSE),
      any(graph >= 0),
      combine = "and"
    )
    sparse <- FALSE
  }

  # convert p0 if p0 is vector
  if (test_numeric(p0, lower = 0, len = n_elements, finite = TRUE,
                   any.missing = FALSE, all.missing = FALSE, null.ok = FALSE)) {
    p0 <- as.matrix(p0)
  } else {
    assert(
      test_matrix(p0, mode = "numeric", nrows = n_elements, any.missing = FALSE,
                  all.missing = FALSE, null.ok = FALSE),
      any(p0 >= 0),
      combine = "and"
    )
  }

  # begin program
  diag(graph) <- 0
  if (correct.for.hubs) {
    graph <- hub.correction(graph)
  }
  stoch.graph <- normalize.stochastic(graph)
  if ((!allow.ergodic) && (!.is.ergodic(stoch.graph))) {
    stop(paste("the provided graph has more than one component.",
               "It is likely not ergodic."))
  }

  if (sparse) {
    # sparse matrix
    l <- mrwr_s(normalize.stochastic(p0),
                stoch.graph, r, thresh, niter, do.analytical)
  } else {
    # dense matrix
    l <- mrwr_(normalize.stochastic(p0),
               stoch.graph, r, thresh, niter, do.analytical)
  }
  if (!return.pt.only) {
    l <- list(p.inf = l, transition.matrix = stoch.graph)
  }

  return(l)
}
