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

p0 <- c(1, rep(0, 4))
adja <- matrix(1, 5, 5)
graph <-  igraph::graph_from_adjacency_matrix(adja)

test_that("random walk if w restarts", {
  s <- suppressMessages(random.walk(p0, as.matrix(igraph::get.adjacency(graph)), 1))
  expect_equal(as.vector(s$p.inf), p0, tolerance = 0.001)
})

test_that("random walk with vectors is same as for matrices", {
  p0 <- matrix(runif(10*20), nrow=5)
  suppressMessages({
    mat.walk <- random.walk(p0, as.matrix(igraph::get.adjacency(graph)), .5)
    vec.walk <- sapply(seq(ncol(p0))[1], function(e) {
      p <- random.walk(p0[, e], as.matrix(igraph::get.adjacency(graph)), .5)
      p$p.inf
    })
  })
  expect_equal(mat.walk$p.inf, as.vector(vec.walk), tolerance = 0.001)
})

test_that("random walk analytical is same as iterative", {
  suppressMessages({
    ana.walk <- random.walk(p0, as.matrix(igraph::get.adjacency(graph)),
                            .5, do.analytical=TRUE)
    it.walk  <- random.walk(p0, as.matrix(igraph::get.adjacency(graph)),
                            .5, do.analytical=FALSE)
  })
  expect_equal(ana.walk$p.inf, it.walk$p.inf, tolerance = 0.001)
})

test_that("random walk if w/o restarts", {
  s <- suppressMessages(random.walk(p0, as.matrix(igraph::get.adjacency(graph)), 0))
  s$p.inf <- s$p.inf / sum(s$p.inf)
   expect_equal(as.vector(s$p.inf), rep(.2, 5), tolerance = 0.001)
})

test_that("random walk with mat for graph", {
  s <- suppressMessages(random.walk(p0, adja, 1))
  expect_equal(as.vector(s$p.inf), p0, tolerance = 0.001)
})

test_that("random walk if false p0", {
  expect_error(random.walk(p0 + .1, graph, 1))
})

test_that("random walk if false mat", {
  expect_error(random.walk(p0, 2, 1))
})

test_that("random walk if p0 not numeric", {
  expect_error(random.walk(1L, graph, 1))
})

test_that("random walk if r not in", {
  expect_error(random.walk(p0, graph, 2))
})

test_that("random walk if r not numeric", {
  expect_error(random.walk(p0, graph, "s"))
})

test_that("random walk p0 smaller zero", {
  expect_error(random.walk(p0 - 5, graph, 1))
})

test_that("hub correction workse", {
  set.seed(23)
  g <- igraph::barabasi.game(10, directed=F)
  ad <- as.matrix(igraph::get.adjacency(g))
  pt <- runif(10)
  suppressMessages({
    r1 <- random.walk(pt, ad, .5, correct.for.hubs=FALSE)
    r2 <- random.walk(pt, ad, .5, correct.for.hubs=TRUE)
  })
  expect_true(all(r1$p.inf != r2$p.inf))
  expect_true(any(r1$transition.matrix != r2$transition.matrix))
})
