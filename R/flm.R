#' @name flm
#' @title Forward Regression Selection for Linear Models.
#' 
#' @description
#' `flm()` inherits the usage of the function [lm], and performs forward regression selection for linear models.
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
#' @param selectFun Parameter passed to [pboost::frs].
#' @param stopFun Parameter passed to [pboost::frs].
#' @param keep Parameter passed to [pboost::frs].
#' @param maxK Parameter passed to [pboost::frs].
#' @param verbose Parameter passed to [pboost::frs].
#'
#' @return A `lm` model object fitted on the selected features.
#' 
#' @examples
#' set.seed(2026)
#' n <- 300
#' p <- 200
#' x <- matrix(rnorm(n*p), n)
#' 
#' eta <- drop( x[, 1:3] %*% runif(3, 1.0, 1.5) )
#' y <- rnorm(n, eta)
#' DF <- data.frame(y, x)
#' 
#' flm(y ~ ., DF, verbose=TRUE)
#' flm(y ~ ., DF, stopFun=BIC, verbose=TRUE)
#' 
#' frs(y, x, lm, verbose=TRUE)
#' 
NULL
#> NULL


#' @rdname flm
#' @order 1
#' @export
flm <- function(formula, data, subset, weights, na.action,
                method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
                singular.ok = TRUE, contrasts = NULL, offset, ...,
                selectFun = logLik, stopFun = "EBIC",
                keep = NULL, maxK = NULL, verbose = FALSE) {

    mf_args <- list(formula = formula, data = data)
    mf <- do.call(model.frame, mf_args)

    yvec <- model.response(mf)

    terms <- terms(formula, data = mf)
    use.intercept <- attr(terms, "intercept") == 1L

    attr(terms, "intercept") <- 0L
    xmat <- model.matrix(terms, mf)

    return(
        frs(yvec = yvec, xmat = xmat, fitFun = lm,
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
            selectFun = selectFun,
            stopFun = stopFun,
            keep = keep,
            maxK = maxK,
            verbose = verbose
        )
    )
}



#' @rdname flm
#' @order 2
#' @export
flm.fit <- function(x, y, offset = NULL, method = "qr",
                    tol = 1e-07, singular.ok = TRUE, ...,
                    selectFun = "logLik", stopFun = "EBIC",
                    keep = NULL, maxK = NULL, verbose = FALSE) {

    n <- NROW(x)
    p <- NCOL(x)

    if (is.character(selectFun) && selectFun == "logLik")
        selectFun <- function(object) {
            class(object) <- "lm"
            return(logLik(object))
        }
    
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
        frs(yvec = y, xmat = x, fitFun = lm.fit,
            offset = if (!missing(offset)) offset else NULL,
            method = method,
            singular.ok = singular.ok,
            contrasts = contrasts,
            ...,
            use.formula = FALSE,
            selectFun = selectFun,
            stopFun = stopFun,
            keep = keep,
            maxK = maxK,
            verbose = verbose
        )
    )
}