#' @name plm
#' @title Profile Boosting for Linear Models.
#' 
#' @description
#' [plm] inherits the usage of the built-in function [lm].
#' 
#' @param formula Parameter passed to [lm].
#' @param data Parameter passed to [lm].
#' @param subset Parameter passed to [lm].
#' @param weights Parameter passed to [lm].
#' @param na.action Parameter passed to [lm].
#' @param method Parameter passed to [lm].
#' @param model Parameter passed to [lm].
#' @param x Parameter passed to [lm].
#' @param y Parameter passed to [lm].
#' @param qr Parameter passed to [lm].
#' @param singular.ok Parameter passed to [lm].
#' @param contrasts Parameter passed to [lm].
#' @param offset Parameter passed to [lm] or [lm.fit].
#' @param ... Parameters passed to [lm] or [lm.fit].
#' @param tol Parameter passed to [lm.fit].
#' 
#' @param stopFun Parameter passed to [pboost::pboost].
#' @param keep Parameter passed to [pboost::pboost].
#' @param maxK Parameter passed to [pboost::pboost].
#' @param verbose Parameter passed to [pboost::pboost].
#' 
#' @return A `lm` model object fitted on the selected features.
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
#' set.seed(2026)
#' n <- 300
#' p <- 200
#' x <- matrix(rnorm(n*p), n)
#' 
#' eta <- drop( x[, 1:3] %*% runif(3, 1.0, 1.5) )
#' y <- rnorm(n, eta, sd=sd(eta))
#' DF <- data.frame(y, x)
#' 
#' plm(y ~ ., DF, verbose=TRUE)
#' plm(y ~ ., DF, stopFun=BIC, verbose=TRUE)
#' pboost(y, x, lm, residuals, verbose=TRUE)
#' 
NULL
#> NULL


#' @rdname plm
#' @order 1
#' @export
plm <- function(formula, data, subset, weights, na.action,
                method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
                singular.ok = TRUE, contrasts = NULL, offset, ...,
                stopFun = "EBIC",
                keep = NULL, maxK = NULL, verbose = FALSE) {

    mf_args <- list(formula = formula, data = data)
    mf <- do.call(model.frame, mf_args)

    yvec <- model.response(mf)

    terms <- terms(formula, data = mf)
    use.intercept <- attr(terms, "intercept") == 1L

    attr(terms, "intercept") <- 0L
    xmat <- model.matrix(terms, mf)

    return(
        pboost(yvec = yvec, xmat = xmat,
            fitFun = lm,
            scoreFun = residuals,
            stopFun = stopFun,
            subset = if (!missing(subset)) subset else NULL,
            weights = if (!missing(weights)) weights else NULL,
            na.action = if (!missing(na.action)) na.action else NULL,
            method = method,
            model = model,
            x = x,
            y = y,
            qr = qr,
            singular.ok = singular.ok,
            contrasts = contrasts,
            offset = if (!missing(offset)) offset else NULL,
            ...,
            use.intercept = use.intercept,
            keep = keep,
            maxK = maxK,
            verbose = verbose
        )
    )
}


#' @rdname plm
#' @order 2
#' @export
plm.fit <- function(x, y, offset = NULL, method = "qr",
                    tol = 1e-07, singular.ok = TRUE, ...,
                    stopFun = "EBIC",
                    keep = NULL, maxK = NULL, verbose = FALSE) {

    n <- NROW(x)
    p <- NCOL(x)
    
    if (is.character(stopFun) && stopFun == "EBIC")
        stopFun <- function(object) {
            class(object) <- "lm"

            ebic.r <- max( 0.0, 1.0 - log(n) / (2.0 * log(p)) )
            stopifnot( ebic.r >= 0 )

            dof <- length(coef(object))
            ebic.penalty <- 2.0 * ebic.r * lchoose(p - length(keep), dof - length(keep))
            stopifnot( is.finite(ebic.penalty) )

            return(BIC(object) + ebic.penalty)
        }

    return(
        pboost(yvec = y, xmat = x,
            fitFun = lm.fit,
            scoreFun = residuals,
            stopFun = stopFun,
            offset = if (!missing(offset)) offset else NULL,
            method = method,
            singular.ok = singular.ok,
            contrasts = contrasts,
            ...,
            use.formula = FALSE,
            keep = keep,
            maxK = maxK,
            verbose = verbose
        )
    )
}
