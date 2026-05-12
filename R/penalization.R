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
#' library(survival)
#' set.seed(2026)
#' n <- 10
#' p <- 20
#' 
#' x <- replicate(p, rnorm(n))
#' b0 <- runif(3, 1.5, 2.0)
#' eta <- drop(x[, 1:3] %*% b0)
#' 
#' ## ---------- linear ----------
#' y <- rnorm(n, eta)
#' 
#' pen_glmnet(x=x, y=y, family="gaussian") |> coef()
#' pen_ncvreg(x=x, y=y, family="gaussian", penalty="MCP") |> coef()
#' pen_ncvreg(x=x, y=y, family="gaussian", penalty="SCAD") |> coef()
#' 
#' 
#' ## ---------- logistic ----------
#' y <- rbinom(n, 1, 1.0 / (1.0 + exp(-eta)))
#' 
#' pen_glmnet(x=x, y=y, family="binomial") |> coef()
#' pen_ncvreg(x=x, y=y, family="binomial", penalty="MCP") |> coef()
#' pen_ncvreg(x=x, y=y, family="binomial", penalty="SCAD") |> coef()
#' 
#' 
#' ## ---------- cox ----------
#' h0 <- 0.01
#' censoringRate <- 0.3
#' survivalTime <- -log(runif(n)) / (h0 * exp(eta))
#' censoringTime <- rexp(n, rate = -log(1 - censoringRate)/median(survivalTime))
#' y <- cbind(
#'     time = pmin(survivalTime, censoringTime),
#'     status = as.numeric(survivalTime <= censoringTime)
#' )
#' 
#' pen_glmnet(x=x, y=y, family="cox") |> coef()
#' pen_ncvreg(x=x, y=y, family="cox", penalty="MCP") |> coef()
#' pen_ncvreg(x=x, y=y, family="cox", penalty="SCAD") |> coef()
#' 
#' 
#' ## ---------- quantile ----------
#' y <- eta + rt(n, 2)
#' 
#' pen_rq(x=x, y=y, tau=0.5, penalty="lasso") |> coef()
#' pen_rq(x=x, y=y, tau=0.5, penalty="MCP") |> coef()
#' pen_rq(x=x, y=y, tau=0.5, penalty="SCAD") |> coef()
#' 
#' 
#' ## ---------- sar ----------
#' w0 <- set_rook_matrix(5, n/5)
#' rho0 <- 0.5
#' y <- solve(diag(n) - rho0 * w0, rnorm(n, eta))
#' 
#' pen_sar(x=x, y=y, rho=rho0, w=w0, penalty="lasso") |> coef()
#' pen_sar(x=x, y=y, rho=rho0, w=w0, penalty="MCP") |> coef()
#' pen_sar(x=x, y=y, rho=rho0, w=w0, penalty="SCAD") |> coef()
#' 
#' @noRd
#' 
NULL
#> NULL
