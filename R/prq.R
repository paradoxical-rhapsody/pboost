#' @name prq
#' @title Profile Boosting for Quantile Regression
#' 
#' @description
#' [prq] inherits the usage of the function [quantreg::rq].
#' 
#' @param formula See [pboost].
#' @param data See [pboost].
#' @param tau Parameters passed to [quantreg::rq].
#' @param subset Parameters passed to [quantreg::rq].
#' @param weights Parameters passed to [quantreg::rq].
#' @param na.action Parameters passed to [quantreg::rq].
#' @param method Parameters passed to [quantreg::rq].
#' @param model Parameters passed to [quantreg::rq].
#' @param contrasts Parameters passed to [quantreg::rq].
#' @param ... Parameters passed to [quantreg::rq].
#' @param stopFun Parameters passed to [pboost].
#' @param keep Parameters passed to [pboost].
#' @param maxK Parameters passed to [pboost].
#' @param verbose Parameters passed to [pboost].
#' 
#' @return An `rq` model object fitted on the selected features.
#' 
#' @examples
#' library(quantreg)
#' set.seed(2025)
#' n <- 300
#' p <- 200
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
#' 
#' scorerq <- function(object) {
#'  return(ifelse(object[["y"]] < fitted(object), tau - 1, tau))
#' }
#' pboost(y ~ ., DF, rq, scorerq, EBIC, tau=tau, verbose=TRUE)
#' 
NULL
#> NULL



#' @rdname prq
#' @order 1
#' @export
prq <- function(
    formula, tau = 0.5, data, subset, weights, na.action,
    method = "br", model = TRUE, contrasts = NULL, ...,
    stopFun = EBIC, keep = NULL, maxK = NULL, verbose = FALSE) {
    stopifnot( !missing(formula) )
    stopifnot( !missing(data) )

    cl <- match.call()

    rq_template <- cl
    rq_template$stopFun <- NULL
    rq_template$keep <- NULL
    rq_template$maxK <- NULL
    rq_template$verbose <- NULL
    rq_template[[1L]] <- quote(rq)

    required_paras <- c("tau", "data", "subset", "weights", "na.action")
    for (ipara in required_paras)
        if (!is.null(cl[[ipara]]))
            rq_template[[ipara]] <- eval(cl[[ipara]], envir = parent.frame())

    fitFun <- function(formula, data) {
        call <- rq_template
        call$formula <- formula
        call$data <- data
        return( eval(call, parent.frame()) )
    }

    scoreFun <- function(object) 
        return(ifelse(object[["y"]] < fitted(object), tau - 1, tau))

    return(pboost(formula, data, fitFun, scoreFun, stopFun,
                  keep = keep, maxK = maxK, verbose = verbose))
}


#' @rdname EBIC
#' @export
EBIC.rq <- function(object, p, p.keep, ...) {
    stopifnot( inherits(object, c("rq", "rqs")) )

    if (missing(p))
        p <- get("p", envir=parent.frame())
    if (missing(p.keep))
        p.keep <- get("p.keep", envir=parent.frame())

    # `logLik(rq or rqs)` has attr `n` rather than `nobs`
    obj.loglik <- logLik(object)
    n0 <- attr(obj.loglik, "n")
    dof <- attr(obj.loglik, "df")
    ebic.r <- max( 0.0, 1.0 - log(n0) / (2.0*log(p)) )
    ebic.penalty <- ifelse(
        ebic.r <= 0.0,
        0.0,
        2.0 * ebic.r * lchoose(p - p.keep, dof - p.keep)
    )

    stopifnot( is.finite(ebic.penalty) )

    # `BIC` for class `rq` is equivalent to `AIC` with negative `k`
    # return( -2 * as.numeric(obj.loglik) + dof * log(n0) + ebic.penalty )
    return( AIC(object, k=-1) + ebic.penalty )
}