library(testthat)
library(determinantalsamplingdesigns)

test_check("determinantalsamplingdesigns")

test_that("Probabilities respected", {
  n = 4
  r = 1
  N = 7
  round(Re(diag(periodicdsd(n, r, N, T))), 10) == round(rep(n/N, N), 10)
})

test_that("Matrice type", {
  n = 4
  r = 1
  N = 7
  #hermitian matrix
  expect_equal(round(periodicdsd(n, r, N, T) %*% t(Conj(periodicdsd(n, r, N, T))), 10) , round(periodicdsd(n, r, N, T), 10))
  #projection matrix
  expect_equal(round(eigen(periodicdsd(n, r, N, T) %*% t(periodicdsd(n, r, N, T)), symmetric = TRUE)$values, 10) == 1 | round(eigen(periodicdsd(n, r, N, T) %*% t(periodicdsd(n, r, N, T)), symmetric = TRUE)$values, 10) == 0, rep(TRUE, N))
})
