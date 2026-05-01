#' @name sar.model
#' @title Functions for Spatial Auto-regressive Model
#' 
#' @description
#' - `get_rho()`: MLE of \eqn{\rho} in SAR model. See [pboost::psar].
#' - `sar.fit()`: Fit SAR model with given \eqn{\rho}.
#' - `logLik()`, `BIC()`, `residuals()`.
#' 
#' @param x Numeric feature matrix.
#' @param y Response vector.
#' @param w Weight matrix (row-sum scaled).
#' @param object Object.
#' @param ... Arguments.
#' 
#' @return Value or list.
#' 
NULL
#> NULL


#' @name sar.model
#' @export
get_rho <- function(x, y, w) {
    stopifnot( is.matrix(x) )
    stopifnot( NROW(x) == length(y) )
    stopifnot( all( abs(rowSums(w) - 1.0) <= 1e-12 ) )
    stopifnot( all(w >= 0.0) )

    if (NCOL(x) > NROW(x))
        warning("Singularity Possible!")

    minusloglik <- function(rho) {
        A.rho <- diag(NROW(x)) - rho * w
        sig2.hat <- mean( residuals(lm.fit(x, A.rho %*% y))^2 ) # singular.ok=FALSE
        0.5*NROW(x) * log(sig2.hat) - log(det(A.rho))
    }

    # return(optimize(minusloglik, c(-1, 1))$minimum)
    rho.min <- optimize(minusloglik, c(-1, 1))
    egg <- rho.min$minimum
    attr(egg, "minusloglik") <- rho.min$objective
    return(egg)
}

#' @name sar.model
#' @export
sar.fit <- function(x, y, w) {
    stopifnot( NROW(x) == length(y) )
    stopifnot( is.matrix(x) )

    rho.hat  <- get_rho(x, y, w)
    y.tilde  <- y - rho.hat * (w %*% y)
    beta.hat <- coef(lm.fit(x, y.tilde, singular.ok = FALSE))
    err      <- drop(y.tilde - x %*% beta.hat)
    sig2.hat <- mean( err^2 )

    egg <- list(
        beta = beta.hat,
        sig2 = sig2.hat,
        rho = rho.hat,
        residuals = err
    )
    attr(egg, "dof") <- NCOL(x)
    attr(egg, "nobs") <- NROW(x)
    attr(egg, "minusloglik") <- attr(rho.hat, "minusloglik")
    attr(egg, "bic") <- 2 * attr(rho.hat, "minusloglik") + log(NROW(x)) * NCOL(x)

    class(egg) <- "sarpboost"
    return(egg)
}


#' @name sar.model
#' @export
residuals.sarpboost <- function(object, ...)
    return(object[["residuals"]])


#' @name sar.model
#' @export
logLik.sarpboost <- function(object, ...)
    return(-1.0*attr(object, "minusloglik"))


#' @name sar.model
#' @export
BIC.sarpboost <- function(object, ...)
    return(attr(object, "bic"))
