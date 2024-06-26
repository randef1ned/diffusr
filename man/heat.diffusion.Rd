% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heat_diffusion.R
\name{heat.diffusion}
\alias{heat.diffusion}
\title{Graph diffusion using a heat diffusion process on a Laplacian matrix.}
\usage{
heat.diffusion(h0, graph, t = 0.5, ...)
}
\arguments{
\item{h0}{an \eqn{n \times p}-dimensional numeric non-negative
\code{\link[base]{matrix}} (or
\code{\link[Matrix:dgCMatrix-class]{dgCMatrix}},
\code{\link[base]{vector}}) of starting temperatures}

\item{graph}{an (\eqn{n \times n})-dimensional numeric non-negative
adjacence matrix representing the graph}

\item{t}{time point when heat is measured}

\item{...}{additional parameters}
}
\value{
returns the heat on every node as numeric vector
}
\description{
An amount of starting heat gets distribution using the
  Laplacian matrix of a graph. Every iteration (or time interval) \eqn{t}
  heat streams from the starting nodes into surrounding nodes.
}
\examples{
# count of nodes
n <- 5
# starting distribution (has to sum to one)
h0 <- as.vector(rmultinom(1, 1, prob=rep(.2, n)))
# adjacency matrix (either normalized or not)
graph <- matrix(abs(rnorm(n*n)), n, n)
# computation of stationary distribution
heat <- heat.diffusion(h0, graph)
}
\references{
\url{https://en.wikipedia.org/wiki/Laplacian_matrix}

\url{https://en.wikipedia.org/wiki/Heat_equation}
}
