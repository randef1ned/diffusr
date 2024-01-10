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
graph <-  igraph::make_graph(c(1, 2, 2, 3, 3, 4, 4, 5), directed = FALSE)
r <- .5

test_that("stochastic col norm", {
  ad <- igraph::get.adjacency(graph)
  stoch.norm <- suppressMessages(normalize.stochastic(as.matrix(ad)))
  cor.mat <- scale(ad, center=F, scale=Matrix::colSums(ad))[,]
  expect_equal(stoch.norm, cor.mat)
})

test_that("laplacian", {
  expected <- matrix(-.1, 10, 10)
  diag(expected) <- .9
  lapl <- normalize.laplacian(matrix(1, 10, 10))
  expect_equal(lapl, expected)
})

test_that("laplacian two random elements", {
  m <- matrix(abs(rnorm(100)), 10, 10)
  lapl <- normalize.laplacian(m)
  rs <- rowSums(m)
  expect_equal(lapl[1,1], 1 - (m[1,1]/rs[1]), tolerance = .001)
  expect_equal(lapl[5,3], -m[5,3]/sqrt(rs[3] * rs[5]), tolerance = .001)
})

test_that("equals double is true", {
  expect_true(.equals.double(0.9, 1, .1))
})

test_that("equals double is false", {
  expect_false(.equals.double(0.9, 1, .05))
})

test_that("normalize vector", {
  vec <- suppressMessages(sum(normalize.stochastic(1:10)))
  expect_equal(vec, 1)
})

test_that("matrix is ergodic", {
  m <- matrix(1, 10, 10)
  expect_true(.is.ergodic(m))
})

test_that("matrix is not ergodic raises error", {
  m <- matrix(0, 10, 10)
  m[1:3, 1:3]    <- 1
  m[7:10, 7:10] <- 1
  expect_false(.is.ergodic(m))
})


test_that("hub correction has no influence", {
  ad <- as.matrix(igraph::get.adjacency(graph))
  ad.hub <- hub.correction(ad)
  ad1 <- suppressMessages(normalize.stochastic(ad))
  ad2 <- suppressMessages(normalize.stochastic(ad.hub))
  expect_equal(ad1, ad2, tolerance = 0.0001)
})


test_that("hub correction has influence", {
  set.seed(23)
  g <- igraph::barabasi.game(10, directed=F)
  ad <- as.matrix(igraph::get.adjacency(g))
  ad.hub   <- hub.correction(ad)
  norm.ad  <- suppressMessages(normalize.stochastic(ad))
  norm.cor <- suppressMessages(normalize.stochastic(ad.hub))
  expect_true(norm.ad[2,1] > norm.cor[2,1])
  for (c in colSums(norm.cor)) {
    expect_equal(c, 1, tolerance = 0.0001)
  }
})
