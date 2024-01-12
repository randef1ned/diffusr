#' Check if a matrix is \code{dgCMatrix} class
#'
#' @param mat The matrix to check
#'
#' @return A boolean variable
#' @export
#'
#' @importFrom methods is
#'
is.dgCMatrix <- function(mat) {
  return(is(mat, "dgCMatrix"))
}

#' Check if the sparse matrix valid
#'
#' @param adj_matrix The sparse matrix to check
#' @param non_negative Check if any non negative vlues
#'
#' @export
#'
assert_dgCMatrix <- function(adj_matrix, non_negative = TRUE) {
  if (adj_matrix@Dim[1] != adj_matrix@Dim[2]) {
    stop(paste("Error: Assertion on 'adj_matrix' failed: Must have exactly",
               adj_matrix@Dim[2], "rows, but has",
               adj_matrix@Dim[1], "rows."))
  } else if (adj_matrix@Dim[1] < 4) {
    stop(paste("Error: Assertion on 'adj_matrix' failed: Must have at least",
               "3 rows, but has", adj_matrix@Dim[1], "rows."))
  } else if (non_negative && any(adj_matrix@x < 0)) {
    stop(paste("Error: Assertion on 'adj_matrix' failed: Element",
               seq_along(adj_matrix@x)[adj_matrix@x < 0][1], "is not >= 0."))
  }
}
