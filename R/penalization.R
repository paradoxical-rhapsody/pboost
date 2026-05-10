#' @name penalization
#' 
#' @title Penalized Methods for Feature Selection
#' 
#' @description Penalized methods for feature selection.
#' 
#' @param x Feature matrix
#' @param y Response vector.
#' @param family "gaussian", "binomial", or "cox".
#' @param penalty "lasso", "SCAD" or "MCP".
#' @param tau Quantiles to be modeled.
#' @param rho Spatial autoregressive parameter. If missing or NULL, it will be estimated.
#' @param w Weight matrix (row-sum scaled being one).
#' 
#' @return Set of selected features.
#' 
#' @examples
#' set.seed(2026)
#' n <- 50
#' p <- 20
#' 
#' x <- replicate(p, rnorm(n))
#' b0 <- runif(3, 1.5, 2.0)
#' eta <- drop(x[, 1:3] %*% b0)
#' 
#' y <- rnorm(n, eta)
#' 
NULL
#> NULL
