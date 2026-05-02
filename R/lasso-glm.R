#' @title Lasso Penalized Generalized Linear Models
#' @description Fit a generalized linear model with Lasso penalty.
#' @param x Feature matrix.
#' @param y Response vector.
#' @param family Family of the GLM (e.g., "gaussian", "binomial").
#' @return A list containing the fitted model parameters and statistics.
#' @export
lasso_glm <- function(x, y, family) {
    n <- length(y)
    p <- NCOL(x)

    ebic.r <- max(0.0, 1.0 - log(n) / (2.0*log(p)))
    stopifnot( ebic.r >= 0.0 )

    standardize <- FALSE
    intercept <- FALSE
    models <- glmnet(x, y, family, standardize = standardize, intercept = intercept)

    lambda <- models[["lambda"]]
    beta.hat <- models[["beta"]]
    dof <- models[["df"]]
    stopifnot( NCOL(beta.hat) == length(lambda) )
    stopifnot( length(dof) == length(lambda) )

    eta <- as.matrix(x %*% beta.hat)
    loglik <- apply(eta, 2, function(eta.col) 
        logLik(glm(y ~ 0 + offset(eta.col), family=family))
    )
    stopifnot( length(loglik) == length(lambda) )

    EBIC <- -2*loglik + dof*log(n) + 2.0 * ebic.r * lchoose(p, dof)
    idx <- which( abs(as.numeric(beta.hat[, which.min(EBIC)])) > 1e-4 )
    egg <- glm(y ~ 0 + ., family=family, data=data.frame(y, x[, idx, drop=FALSE]))
    attr(egg, "lambda") <- lambda[which.min(EBIC)]
    class(egg) <- c("lassoglm", class(egg))
    return(egg)
}
