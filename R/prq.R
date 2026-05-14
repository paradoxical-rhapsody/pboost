#' @name prq
#' @title Profile Boosting for Quantile Regression Models
#' 
#' @description
#' [prq] inherits the usage of the function [quantreg::rq].
#' 
#' @param formula Parameter passed to [quantreg::rq].
#' @param tau Parameter passed to [quantreg::rq].
#' @param data Parameter passed to [quantreg::rq].
#' @param x Parameter passed to [quantreg::rq.fit].
#' @param y Parameter passed to [quantreg::rq.fit].
#' @param subset Parameter passed to [quantreg::rq].
#' @param weights Parameter passed to [quantreg::rq].
#' @param na.action Parameter passed to [quantreg::rq].
#' @param method Parameter passed to [quantreg::rq] or [quantreg::rq.fit].
#' @param model Parameter passed to [quantreg::rq].
#' @param contrasts Parameter passed to [quantreg::rq].
#' @param ... Parameters passed to [quantreg::rq] or [quantreg::rq.fit].
#' 
#' @param stopFun Parameter passed to [pboost::pboost].
#' @param keep Parameter passed to [pboost::pboost].
#' @param maxK Parameter passed to [pboost::pboost].
#' @param verbose Parameter passed to [pboost::pboost].
#' 
#' @return A `rq` model object fitted on the selected features.
#' 
#' @examples
#' library(quantreg)
#' set.seed(2026)
#' n <- 300
#' p <- 20
#' x <- matrix(rnorm(n*p), n)
#' 
#' eta <- drop( x[, 1:3] %*% runif(3, 1.0, 1.5) )
#' y <- eta + (1.0 + x[, 3]) * rnorm(n)
#' DF <- data.frame(y, x)
#' 
#' tau <- 0.5
#' prq(y ~ ., tau, DF, verbose=TRUE)
#' 
#' BIC <- function(obj) AIC(obj, k=-1)
#' prq(y ~ ., tau, DF, stopFun=BIC, verbose=TRUE)
#' frq(y ~ ., tau, DF, stopFun=BIC, verbose=TRUE)
#' 
#' scoreFun <- function(object)
#'    return(ifelse(object[["y"]] < fitted(object), tau - 1, tau))
#'
#' pboost(x, y, rq, scoreFun, BIC, tau=tau, verbose=TRUE)
#' 
#' prq.fit(x, y, verbose=TRUE)
#' frq.fit(x, y, verbose=TRUE)
#' 
NULL
#> NULL



#' @rdname prq
#' @order 1
#' @export
prq <- function(formula, tau = 0.5, data, subset, weights, na.action,
                method = "br", model = TRUE, contrasts = NULL, ...,
                stopFun = "EBIC",
                keep = NULL, maxK = NULL, verbose = FALSE) {

    stopifnot( length(tau) == 1 )

    fml <- Formula(formula)
    mf <- model.frame(fml, data=data)

    yvec <- model.part(fml, data=mf, lhs=1, drop=TRUE)
    xmat <- model.part(fml, data=mf, rhs=1) |> as.matrix()
    use.intercept <- attr(terms(mf), "intercept") == 1L


    scoreFun <- function(object)
        return(ifelse(object[["y"]] < fitted(object), tau - 1, tau))

    n <- NROW(xmat)
    p <- NCOL(xmat)
    if (is.character(stopFun) && stopFun == "EBIC")
        stopFun <- function(object) {
            ebic.r <- max( 0.0, 1.0 - log(n) / (2.0 * log(p)) )
            stopifnot( ebic.r >= 0 )

            dof <- attr(logLik(object), "df")
            ebic.penalty <- 2.0 * ebic.r * lchoose(p - length(keep), dof - length(keep))
            stopifnot( is.finite(ebic.penalty) )

            # `BIC` for class `rq` is equivalent to `AIC` with negative `k`
            # return( -2*as.numeric(obj.loglik) + dof * log(n) + ebic.penalty )
            return( AIC(object, k=-1) + ebic.penalty )
        }

    mc <- match.call(expand.dots = TRUE)
    provided_args <- as.list(mc)[-1]
    provided_args <- provided_args[!(names(provided_args) %in% c("formula", "data", "scoreFun", "stopFun"))]

    args <- list(
        fitFun = rq,
        yvec = yvec,
        xmat = xmat,
        use.intercept = use.intercept,
        scoreFun = scoreFun,
        stopFun = stopFun
    )
    args <- c(args, provided_args)

    return(do.call(pboost, args))
}



#' @rdname prq
#' @order 2
#' @export
prq.fit <- function(x, y, tau = 0.5, method = "br", ...,
                    stopFun = "EBIC",
                    keep = NULL, maxK = NULL, verbose = TRUE) {

    stopifnot( length(tau) == 1 )

    scoreFun <- function(object)
        return(ifelse(object[["y"]] < fitted(object), tau - 1, tau))

    n <- NROW(x)
    p <- NCOL(x)
    if (is.character(stopFun) && stopFun == "EBIC")
        stopFun <- function(object) {
            ebic.r <- max( 0.0, 1.0 - log(n) / (2.0 * log(p)) )
            stopifnot( ebic.r >= 0 )

            dof <- length(object$coefficients)
            ebic.penalty <- 2.0 * ebic.r * lchoose(p - length(keep), dof - length(keep))
            stopifnot( is.finite(ebic.penalty) )

            Rho <- function(u, tau) u * (tau - (u < 0))
            fid <- sum(Rho(object$residuals, tau))
            loglik <- n * (log(tau * (1 - tau)) - 1 - log(fid/n))

            return( -2*loglik + dof * log(n) + ebic.penalty )
        }

    mc <- match.call(expand.dots = TRUE)
    provided_args <- as.list(mc)[-1]
    provided_args <- provided_args[!(names(provided_args) %in% c("x", "y", "scoreFun", "stopFun"))]

    weights <- provided_args$weights
    fitFun <- ifelse(!is.null(weights) && length(weights), rq.wfit, rq.fit)

    args <- list(
        fitFun = fitFun,
        yvec = y,
        xmat = x,
        scoreFun = scoreFun,
        stopFun = stopFun,
        use.formula = FALSE
    )
    args <- c(args, provided_args)

    return(do.call(pboost, args))
}
