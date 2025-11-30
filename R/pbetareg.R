#' @name pbetareg
#' @title Profile Boosting for Beta Regression
#' 
#' @description
#' [pbetareg] has the similar usage to the function [betareg::betareg].
#' 
#' @param formula See [pboost].
#' @param data See [pboost].
#' @param subset Parameters passed to [betareg::betareg].
#' @param na.action Parameters passed to [betareg::betareg].
#' @param weights Parameters passed to [betareg::betareg].
#' @param offset Parameters passed to [betareg::betareg].
#' @param link Parameters passed to [betareg::betareg].
#' @param link.phi Parameters passed to [betareg::betareg].
#' @param type Parameters passed to [betareg::betareg].
#' @param dist Parameters passed to [betareg::betareg].
#' @param nu Parameters passed to [betareg::betareg].
#' @param control Parameters passed to [betareg::betareg].
#' @param model Parameters passed to [betareg::betareg].
#' @param y Parameters passed to [betareg::betareg].
#' @param x Parameters passed to [betareg::betareg].
#' @param ... Parameters passed to [betareg::betareg].
#' @param stopFun Parameters passed to [pboost].
#' @param keep Parameters passed to [pboost].
#' @param maxK Parameters passed to [pboost].
#' @param verbose Parameters passed to [pboost].
#' 
#' @return An `betareg` model object fitted on the selected features.
#' 
#' @examples 
#' library(betareg)
#' set.seed(2025)
#' n <- 300
#' p <- 100
#' x <- matrix(runif(n*p), n)
#' mu <- runif(n)
#' phi <- 1.0
#' 
#' shape1 <- mu * phi
#' shape2 <- (1-mu) * phi
#' y <- rbeta(n, shape1, shape2)
#' DF <- data.frame(y, x)
#' 
#' ## The function `pbetareg` has similar usage to the built-in `glm`
#' # print( pbetareg(formula=y ~ ., data=DF, verbose=TRUE) )
#' 
#' ## use `BIC` as the stopping rule, which yields too many spurious features
#' # print( pbetareg(formula=y ~ ., data=DF, stopFun=BIC, verbose=TRUE) )
#' 
#' 
NULL
#> NULL


#' @rdname pbetareg
#' @order 1
#' @export
pbetareg <- function(
    formula, data, subset, na.action, weights, offset,
    link = c("logit", "probit", "cloglog", "cauchit", "log",
        "loglog"), link.phi = NULL, type = c("ML", "BC", "BR"),
    dist = NULL, nu = NULL, control = betareg.control(...), model = TRUE,
    y = TRUE, x = FALSE,
    ...,
    stopFun=EBIC, keep=NULL, maxK=NULL, verbose=FALSE){

    cl <- match.call(expand.dots = TRUE)

    betareg_template <- cl
    betareg_template$stopFun <- NULL
    betareg_template$keep <- NULL
    betareg_template$maxK <- NULL
    betareg_template$verbose <- NULL
    betareg_template[[1L]] <- quote(betareg)
    fitFun <- function(formula, data){
        call <- betareg_template
        call$formula <- formula
        return( eval(call, parent.frame()) )
    }

    scoreFun <- function(object) {
        phi <- predict(object, type='precision')
        mu <- predict(object, type='response')
        eta <- predict(object, type='link')
        mu.eta <- object$link$mu$mu.eta
        y <- pmin(pmax(object[["y"]], .Machine$double.eps), 1 - .Machine$double.eps)

        # weights <- object[["weights"]]
        return( mu.eta(eta) * phi * ( digamma((1-mu)*phi) - digamma(mu*phi) + qlogis(y) ) )
    }


    return(pboost(formula, data, fitFun, scoreFun, stopFun,
                  keep=keep, maxK=maxK, verbose=verbose))
}




#' @rdname EBIC
#' @export
EBIC.betareg <- function(object, p, p.keep, ...){
    stopifnot( inherits(object, "betareg") )

    if (missing(p))
        p <- get("p", envir=parent.frame())
    if (missing(p.keep))
        p.keep <- get("p.keep", envir=parent.frame())

    dof <- attr(logLik(object), "df")
    ebic.r <- max( 0.0, 1.0 - log(nobs(object)) / (2.0*log(p)) )
    ebic.penalty <- ifelse(
        ebic.r <= 0.0,
        0.0,
        2.0 * ebic.r * lchoose(p - p.keep, dof - p.keep)
    )

    # stopifnot( !is.nan(ebic.penalty) )
    stopifnot( is.finite(ebic.penalty) )
    return(BIC(object) + ebic.penalty)
}
