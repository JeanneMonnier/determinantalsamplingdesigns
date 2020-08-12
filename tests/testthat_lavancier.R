library(testthat)
library(determinantalsamplingdesigns)

test_check("determinantalsamplingdesigns")

test_that("Dimensions", {
  v <- loonismary(c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8))
  n = runif(1, 1, 10)
  nrow(lavancier(v, n)) == ncol(v)
  length(lavancier(v)) == ncol(v)
  ncol(lavancier(v, n)) == n
})


test_that("Samples correctly", {
  V <- loonismary(c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8))
  Vbis <- periodicdsd(4, 1, 7)
  m = 1000
  indices <- c(1:nrow(V))
  PIk <- Re(diag(V%*%t(Conj(V))))
  PIkbis <- Re(diag(Vbis%*%t(Conj(Vbis))))
  L <- lavancier(V, m, F, rep(F, ncol(V)))
  Lbis <- lavancier(Vbis, m, F, rep(F, ncol(V)))
  compte <- c()
  comptebis <- c()
  for (l in 1:nrow(V)) {
    nombre <- L == l
    nombrebis <- Lbis == l
    compte[l] <- length(L[nombre])
    comptebis[l] <- length(Lbis[nombrebis])
  }

  effectifs <- compte/m
  effectifsbis <- comptebis/m
  expect_equal(abs(round(effectifs, 1) - round(PIk, 1)) <= 0.2, rep(TRUE, 7))
  expect_equal(abs(round(effectifsbis, 1) - round(PIkbis, 1)) <= 0.2, rep(TRUE, 7))
})
