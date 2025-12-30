#' @name plm
#' @title Profile Boosting for Linear Models.
#' 
#' @description
#' [plm] inherits the usage of the built-in function [lm].
#' 
#' @param formula See [pboost].
#' @param data See [pboost].
#' @param subset Parameters passed to [lm].
#' @param weights Parameters passed to [lm].
#' @param na.action Parameters passed to [lm].
#' @param method Parameters passed to [lm].
#' @param model Parameters passed to [lm].
#' @param x Parameters passed to [lm].
#' @param y Parameters passed to [lm].
#' @param qr Parameters passed to [lm].
#' @param singular.ok Parameters passed to [lm].
#' @param contrasts Parameters passed to [lm].
#' @param offset Parameters passed to [lm].
#' @param ... Parameters passed to [lm].
#' @param stopFun Parameters passed to [pboost].
#' @param keep Parameters passed to [pboost].
#' @param maxK Parameters passed to [pboost].
#' @param verbose Parameters passed to [pboost].
#' 
#' @return An `lm` model object fitted on the selected features.
#' 
#' @details `plm` is an equivalent implementation to the sequential lasso method
#' proposed by Luo and Chen(2014, \doi{10.1080/01621459.2013.877275}).
#' 
#' @references
#' * Zengchao Xu, Shan Luo and Zehua Chen (2022). Partial profile score feature
#' selection in high-dimensional generalized linear interaction models.
#' Statistics and Its Interface. \doi{10.4310/21-SII706}
#' 
#' * Shan Luo and Zehua Chen (2014). A Sequential Lasso Method for Feature Selection
#' with Ultra-High Dimensional Feature Space. Journal of the American Statistical
#' Association, 109(507):223–232. \doi{10.1080/01621459.2013.877275}
#' 
#' @examples
#' set.seed(2025)
#' n <- 300
#' p <- 200
#' x <- matrix(rnorm(n*p), n)
#' 
#' eta <- drop( x[, 1:3] %*% runif(3, 1.0, 1.5) )
#' y <- eta + rnorm(n, sd=sd(eta))
#' DF <- data.frame(y, x)
#' 
#' plm(y ~ ., DF, verbose=TRUE)
#' plm(y ~ ., DF, stopFun=BIC, verbose=TRUE)
#' pboost(y ~ ., DF, lm, residuals, EBIC, verbose=TRUE)
#' 
NULL
#> NULL


#' @rdname plm
#' @order 1
#' @export
plm <- function(
    formula, data, subset, weights, na.action,
    method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
    singular.ok = TRUE, contrasts = NULL, offset, ...,
    stopFun = EBIC, keep = NULL, maxK = NULL, verbose = FALSE) {
    stopifnot( !missing(formula) )
    stopifnot( !missing(data) )

    cl <- match.call()

    lm_template <- cl
    lm_template$stopFun <- NULL
    lm_template$keep <- NULL
    lm_template$maxK <- NULL
    lm_template$verbose <- NULL
    lm_template[[1L]] <- quote(lm)

    required_paras <- c("data", "subset", "weights", "na.action", "offset")
    for (ipara in required_paras)
        if (!is.null(cl[[ipara]]))
            lm_template[[ipara]] <- eval(cl[[ipara]], envir = parent.frame())

    fitFun <- function(formula, data) {
        call <- lm_template
        call$formula <- formula
        call$data <- data
        return( eval(call, parent.frame()) )
    }

    return(pboost(formula, data, fitFun, residuals, stopFun,
                  keep = keep, maxK = maxK, verbose = verbose))
}


#' @rdname EBIC
#' @export
EBIC.lm <- function(object, p, p.keep, ...) {
    stopifnot( inherits(object, "lm") )

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

    stopifnot( is.finite(ebic.penalty) )
    return(BIC(object) + ebic.penalty)
}
