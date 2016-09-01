#' Find the closest neighbors of a group of nodes in a graph.
#'
#' @export
#' @author Simon Dirmeier
#' @import igraph
#' @import Matrix
#'
#' @param node.idxs  vector of node indexes (1-based) for which the algorithm is applied iteratively
#' @param graph  an <code>igraph</code> object
#' @param k  the depth of the nearest neighbor search
#' @param use.edge.weights  boolean flags if the edge weights should be considered when doing nearest neighbor lookup
#' @param ...  additional params
#' @return  returns the kNN graph as an an <code>igraph</code> object
#' @examples
#' \dontrun{
#'  TODO
#' }
neighbors <- function(node.idxs, graph, k=1L, use.edge.weights=F, ...) UseMethod("neighbors")

#' @noRd
#' @export
#' @import igraph
#' @import Matrix
neighbors.integer <- function(node.idxs, graph, k=1L, use.edge.weights=F, ...)
{
  if (any(node.idxs <= 0)) stop("Node idxs have to be 1-indexed!")
  if (!igraph::is_igraph(graph) & !is.matrix(graph) & !.is.Matrix(graph))
    stop("'graph' is not a graph object or matrix!")
  if (!is.numeric(k)) stop("k is not numeric!")
  if (!is.integer(k)) {
    k <- as.integer(k)
    message("Casting k to int!")
  }
  if (k < 1) stop("k must be greater than 0!")
  if (!is.logical(use.edge.weights)) stop("Use edge weughts should be boolean!")
  if (use.edge.weights) stop("Not yet implemented!")
  .neighbors(node.idxs, graph, k, use.edge.weights)
}

#' @noRd
#' @import igraph
#' @import Matrix
.neighbors <-  function(node.idxs, graph, k, use.edge.weights)
{
  mat <- .as.matrix(graph)
  if (any(mat < 0)) stop("graph has to be non-negative")
  if (dim(mat)[1] != dim(mat)[2]) stop("graph has to be of dimension (n x n)!")
  invisible(.neighbors_cpp(as.integer(node.idxs), graph, as.integer(k), use.edge.weights))
}