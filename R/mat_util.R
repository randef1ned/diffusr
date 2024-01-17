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
#'   \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}},
#'   \code{\link[base]{vector}}) that is stochstically normalized
#' @param no_r Do not use R for normalization
#' @param ...  additional params
#' @return  returns the normalized matrix/vector)
#'
#' @useDynLib diffusr
#'
#' @importFrom checkmate assert check_matrix test_numeric test_atomic_vector
#'                       test_logical
#' @importFrom memuse Sys.meminfo Sys.swapinfo howbig
#' @importFrom pryr object_size
#' @importFrom Rcpp sourceCpp
#'
#' @examples
#' W <- matrix(abs(rnorm(10000)), 100, 100)
#' stoch.W <- normalize.stochastic(W)
normalize.stochastic <- function(obj, no_r = NULL, ...) {
  is_matrix <- FALSE
  is_sparse <- FALSE
  if (!test_logical(no_r, len = 1, any.missing = FALSE, all.missing = FALSE,
                    null.ok = FALSE)) {
    no_r <- FALSE
  }
  if (test_numeric(obj, lower = 0, finite = TRUE, any.missing = FALSE,
                   all.missing = FALSE, null.ok = FALSE) &&
      test_atomic_vector(obj)) {
    if (!.equals.double(sum(obj), 1, .001)) {
      message("normalizing vector!")
    }
  } else if (is.dgCMatrix(obj)) {
    assert_dgCMatrix(obj)
    is_matrix <- is_sparse <- TRUE
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
    if (no_r) {
      if (is_sparse) {
        obj <- as(stoch_col_norm_s(obj), "dgCMatrix")
      } else {
        obj <- stoch_col_norm_(obj)
      }
    } else {
      # check memory usage;
      # if there is a memory shortage, then call C function directly
      n <- as.numeric(ncol(obj))
      memory_usage <- Sys.meminfo()
      swap_usage <- Sys.swapinfo()
      free_ram <- memory_usage$freeram@size
      free_ram <- free_ram * switch(substring(memory_usage$freeram@unit, 1, 1),
                                    "B" = 1 / 1048576, "K" = 1 / 1024, "M" = 1,
                                    "G" = 1024, "T" = 1048576,
                                    .default = 1073741824)
      swap_ram <- swap_usage$freeswap@size
      swap_ram <- swap_ram * switch(substring(swap_usage$freeswap@unit, 1, 1),
                                    "B" = 1 / 1048576, "K" = 1 / 1024, "M" = 1,
                                    "G" = 1024, "T" = 1048576,
                                    .default = 1073741824)
      free_ram <- free_ram + swap_ram
      object_ram_p <- howbig(n, n, unit = "MiB")@size     # size in practice
      object_ram_t <- as.numeric(object_size(obj)) / 1e6  # size in theory (MiB)

      # if memory is bigger than the temporary variables, then use R
      if ((free_ram > object_ram_t * 4)) {
        sums <- colSums3(obj, is_sparse)
        if (!all(.equals.double(sums, 1, .001))) {
          message("normalizing column vectors!")
          empt_col_val <- 1.0 / n

          obj <- obj / sums[col(obj)]
          # check if need wipe zeros
          zeros <- which(sums < 0.00001)
          if (length(zeros)) {
            obj[, zeros] <- empt_col_val
          }
        }
      } else if (free_ram < object_ram_p) {
        stop("You don't have sufficient memory to normalize. Required: ",
             round(object_ram_p / 1024, digits = 3), " GiB, but ",
             round(free_ram / 1024, digits = 3), " available.")
      } else {
        warning("You have just enough memory to normalize; consider ",
                "increasing your physical memory capacity in the future!")
        obj <- stoch_col_norm_s(obj)
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
#' @useDynLib diffusr
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
    sums <- matrixStats::colSums2(mat)
  }
  return(sums)
}
