#' @name frq
#' @title Forward Regression Selection for Quantile Regression Models
#' 
#' @description
#' `frq()` inherits the usage of the function [quantreg::rq], and performs forward regression selection for quantile regression models.
#' 
#' @param formula Parameter passed to [quantreg::rq].
#' @param tau Parameter passed to [quantreg::rq].
#' @param data Parameter passed to [quantreg::rq].
#' @param x Parameter passed to [quantreg::rq.fit].
#' @param y Parameter passed to [quantreg::rq.fit].
#' @param subset Parameter passed to [quantreg::rq].
#' @param weights Parameter passed to [quantreg::rq].
#' @param na.action Parameter passed to [quantreg::rq].
#' @param method Parameter passed to [quantreg::rq].
#' @param model Parameter passed to [quantreg::rq].
#' @param contrasts Parameter passed to [quantreg::rq].
#' @param ... Parameters passed to [quantreg::rq].
#' 
#' @param selectFun Parameter passed to [pboost::frs].
#' @param stopFun Parameter passed to [pboost::frs].
#' @param keep Parameter passed to [pboost::frs].
#' @param maxK Parameter passed to [pboost::frs].
#' @param verbose Parameter passed to [pboost::frs].
#' 
#' @return A `rq` model object fitted on the selected features.
#' 
#' NULL
NULL
#> NULL



#' @rdname frq
#' @order 1
#' @export 
frq <- function(formula, tau = 0.5, data, subset, weights, na.action,
                method = "br", model = TRUE, contrasts = NULL, ...,
                selectFun = logLik, stopFun = "EBIC",
                keep = NULL, maxK = NULL, verbose = FALSE) {

    stopifnot( length(tau) == 1 )

    fml <- Formula(formula)
    mf <- model.frame(fml, data=data)

    yvec <- model.part(fml, data=mf, lhs=1, drop=TRUE)
    xmat <- model.part(fml, data=mf, rhs=1) |> as.matrix()
    use.intercept <- attr(terms(mf), "intercept") == 1L


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
    provided_args <- provided_args[!(names(provided_args) %in% c("formula", "data", "stopFun"))]

    args <- list(
        fitFun = rq,
        yvec = yvec,
        xmat = xmat,
        use.intercept = use.intercept,
        stopFun = stopFun
    )
    args <- c(args, provided_args)

    return(do.call(frs, args))
}




#' @rdname frq
#' @order 2
#' @export 
frq.fit <- function(x, y, tau = 0.5, method = "br", ...,
                    selectFun = "logLik", stopFun = "EBIC",
                    keep = NULL, maxK = NULL, verbose = FALSE) {

    n <- NROW(x)
    p <- NCOL(x)

    if (is.character(selectFun) && selectFun == "logLik")
        selectFun <- function(object) {
            Rho <- function(u, tau) u * (tau - (u < 0))
            fid <- sum(Rho(object$residuals, tau))
            loglik <- n * (log(tau * (1 - tau)) - 1 - log(fid/n))
        }

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
    provided_args <- provided_args[!(names(provided_args) %in% c("x", "y", "selectFun", "stopFun"))]

    weights <- provided_args$weights
    fitFun <- ifelse(!is.null(weights) && length(weights), rq.wfit, rq.fit)

    args <- list(
        fitFun = fitFun,
        yvec = y,
        xmat = x,
        selectFun = selectFun,
        stopFun = stopFun,
        use.formula = FALSE
    )
    args <- c(args, provided_args)

    return(do.call(frs, args))
}
