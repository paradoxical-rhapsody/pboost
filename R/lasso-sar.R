#' @title Lasso Penalized Spatial Autoregressive Model
#' 
#' @description Fit a spatial autoregressive model with Lasso penalty.
#' 
#' @param x Feature matrix.
#' @param y Response vector.
#' @param w Spatial weight matrix.
#' @param rho Spatial autoregressive parameter. If NULL, it will be estimated.
#' @param lambda.vec Optional vector of lambda values for glmnet. If NULL, glmnet will choose its own sequence.
#' 
#' @return A list containing the fitted model parameters and statistics.
#' 
#' @export
lasso_sar <- function(x, y, w, rho, lambda.vec=NULL) {
    n <- length(y)
    p <- NCOL(x)

    ebic.r <- max(0.0, 1.0 - log(n) / (2.0*log(p)))
    stopifnot( ebic.r >= 0.0 )

    rho <- ifelse(missing(rho) || is.null(rho), get_rho(x, y, w), rho)
    stopifnot( length(rho) == 1 & is.numeric(rho) )

    A.rho <- diag(n) - rho * w
    y <- drop(A.rho %*% y)
    stopifnot( length(y) == n )

    standardize <- FALSE
    intercept <- FALSE
    if (missing(lambda.vec) || is.null(lambda.vec))
        models <- glmnet(x, y, "gaussian", standardize = standardize, intercept = intercept)
    else
        models <- glmnet(x, y, "gaussian", standardize = standardize, intercept = intercept, lambda = lambda.vec)

    lambda <- models[["lambda"]]
    beta.hat <- models[["beta"]]

    sig2.hat <- (as.numeric(y) - as.matrix(x %*% beta.hat))^2 |> colMeans() |> as.numeric()
    stopifnot( length(sig2.hat) == length(lambda) )

    dof <- models[["df"]]
    stopifnot( length(dof) == length(lambda) )

    BIC <- n * log(sig2.hat) + log(n) * dof # - 2*log(det(A.rho))
    EBIC <- BIC + 2.0 * ebic.r * lchoose(p, dof)

    opt.ebic <- which.min(EBIC)
    ebic.flag <- abs(as.numeric(beta.hat[, opt.ebic])) > 1e-6
    stopifnot( length(ebic.flag) == p )

    egg <- list(
        beta = as.numeric(beta.hat[, opt.ebic]),
        sig2 = sig2.hat[opt.ebic],
        rho = rho,
        flag = ebic.flag,
        lambda = lambda[opt.ebic]
    )

    class(egg) <- "sarlasso"
    return(egg)
}


#' @title Extract Coefficients from Lasso SAR Model
#' @description Extract the estimated coefficients from a fitted Lasso SAR model.
#' @param object A fitted Lasso SAR model object.
#' @param ... Additional arguments (not used).
#' @return A numeric vector of estimated coefficients.
#' @export
coef.sarlasso <- function(object, ...) {
    if (!inherits(object, "sarlasso"))
        stop("Object must be of class 'sarlasso'")
    idx <- which(object$flag)
    return(object$beta[idx])
}
