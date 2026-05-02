#' @title Lasso Penalized Quantile Regression
#' 
#' @description Fit a quantile regression model with Lasso penalty.
#' 
#' @param x Feature matrix.
#' @param y Response vector.
#' @param tau Quantile level (between 0 and 1).
#' @param lambda.vec Vector of lambda values.
#' 
#' @return A list containing the fitted model parameters and statistics.
#' 
#' @export
lasso_rq <- function(x, y, tau, lambda.vec=2.0^seq(-5, 5, by=0.01)) {
    n <- length(y)
    p <- NCOL(x)

    ebic.r <- max(0.0, 1.0 - log(n) / (2.0*log(p)))
    stopifnot( ebic.r >= 0.0 )

    ## taken from `rq()`
    Rho <- function(u, tau) u * (tau - (u < 0))

    ebic <- rep(NA_real_, length(lambda.vec))
    for (i in seq_along(lambda.vec)) {
        egg <- rq.fit.lasso(cbind(intercept=1, x), y, tau=tau, lambda=lambda.vec[i])

        dof <- sum(abs(coef(egg)[-1]) > 1e-4)
        stopifnot( dof >= 0 && dof <= p )
        fid <- sum(Rho(residuals(egg), tau))
        minusloglik <- -n * (log(tau * (1 - tau)) - 1 - log(fid/n))
        ebic.penalty <- 2.0 * ebic.r * lchoose(p, dof)
        ebic[i] <- 2*minusloglik + log(n) * dof + ebic.penalty
    }

    idx <- which.min(ebic)
    lambda.opt <- lambda.vec[idx]
    egg <- rq.fit.lasso(cbind(intercept=1, x), y, tau=tau, lambda=lambda.opt)
    egg[["beta"]] <- coef(egg)[-1][abs(coef(egg)[-1]) >= 1e-4]
    class(egg) <- c("lassorq", class(egg))
    return(egg)
}



#' @title Extract Coefficients from Lasso Penalized Quantile Regression
#'
#' @param object A fitted lasso_rq model object.
#' @param ... Additional arguments (not used).
#' @return A named vector of coefficients.
#' @export
coef.lassorq <- function(object, ...) {
    if (!inherits(object, "lassorq"))
        stop("Object must be of class 'lassorq'")
    return(object[["beta"]])
}
