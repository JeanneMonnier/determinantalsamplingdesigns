library(testthat)
library(determinantalsamplingdesigns)

test_check("determinantalsamplingdesigns")

test_that("Probabilities respected", {
  Pi = c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8)
  expect_equal(round(diag(loonismary(Pi, T)), 2), Pi)
})

test_that("Matrice type", {
  Pi = c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8)
  #hermitian matrix
  expect_equal(round(loonismary(Pi, T) %*% t(loonismary(Pi, T)), 10), round(loonismary(Pi, T), 10))
  #projection matrix
  expect_equal(round(eigen(loonismary(Pi, T) %*% t(loonismary(Pi, T)), symmetric = TRUE)$values, 10) == 1 | round(eigen(loonismary(Pi, T) %*% t(loonismary(Pi, T)), symmetric = TRUE)$values, 10) == 0, rep(TRUE, length(Pi)))
})

test_that("Probabilities vector needed", {
  Pi1 = c(0.5, 0.75, 1, 0.2, 0.4, 0.6, 0.8)
  Pi2 = c(0.5, 0.75, 0.75, 0.2, -0.4, 0.6, 0.8)
  Pi3 = c(0.6)
  expect_error(loonismary(Pi1, T))
  expect_error(loonismary(Pi2, T))
  expect_error(loonismary(Pi3, T))
})

