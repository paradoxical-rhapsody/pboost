#' @name prq
#' @title Profile Boosting for Quantile Regression Models
#' 
#' @description
#' [prq] inherits the usage of the function [quantreg::rq].
#' 
#' @param formula Parameter passed to [quantreg::rq].
#' @param tau Parameter passed to [quantreg::rq].
#' @param data Parameter passed to [quantreg::rq].
#' @param subset Parameter passed to [quantreg::rq].
#' @param weights Parameter passed to [quantreg::rq].
#' @param na.action Parameter passed to [quantreg::rq].
#' @param method Parameter passed to [quantreg::rq].
#' @param model Parameter passed to [quantreg::rq].
#' @param contrasts Parameter passed to [quantreg::rq].
#' @param ... Parameters passed to [quantreg::rq].
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
#' 
#' scoreFun <- function(object)
#'    return(ifelse(object[["y"]] < fitted(object), tau - 1, tau))
#'
#' pboost(y, x, rq, scoreFun, BIC, tau=tau, verbose=TRUE)
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

    mf_args <- list(formula = formula, data = data)
    mf <- do.call(model.frame, mf_args)

    yvec <- model.response(mf)

    terms <- terms(formula, data = mf)
    use.intercept <- attr(terms, "intercept") == 1L

    attr(terms, "intercept") <- 0L
    xmat <- model.matrix(terms, mf)

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
    provided_args <- provided_args[!(names(provided_args) %in% c("scoreFun", "stopFun"))]
    provided_args$formula <- NULL
    provided_args$data <- NULL

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
