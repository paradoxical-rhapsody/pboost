#' @name fglm
#' @title Forward Regression Selection for Generalized Linear Models
#' 
#' @description `fglm()` inherits the usage of [glm], and performs forward regression selection for generalized linear models.
#' 
#' @param formula Parameter passed to [glm].
#' @param data Parameter passed to [glm].
#' @param family Parameter passed to [glm].
#' @param weights Parameter passed to [glm].
#' @param subset Parameter passed to [glm].
#' @param na.action Parameter passed to [glm].
#' @param start Parameter passed to [glm].
#' @param etastart Parameter passed to [glm].
#' @param mustart Parameter passed to [glm].
#' @param offset Parameter passed to [glm].
#' @param control Parameter passed to [glm].
#' @param model Parameter passed to [glm].
#' @param method Parameter passed to [glm].
#' @param x Parameter passed to [glm].
#' @param y Parameter passed to [glm].
#' @param singular.ok Parameter passed to [glm].
#' @param contrasts Parameter passed to [glm].
#' @param ... Parameters passed to [glm].
#' @param intercept Parameter passed to [glm.fit].
#' 
#' @param selectFun Parameter passed to [pboost::frs].
#' @param stopFun Parameter passed to [pboost::frs].
#' @param keep Parameter passed to [pboost::frs].
#' @param maxK Parameter passed to [pboost::frs].
#' @param verbose Parameter passed to [pboost::frs].
#' 
#' @return A `glm` model object fitted on the selected features.
#' 
#' @examples
#' set.seed(2026)
#' n <- 300
#' p <- 200
#' x <- matrix(rnorm(n*p), n)
#' 
#' eta <- drop( x[, 1:3] %*% runif(3, 1.0, 1.5) )
#' y <- rbinom(n, 1, 1/(1+exp(-eta)))
#' DF <- data.frame(y, x)
#' 
#' fglm(y ~ ., "binomial", DF, verbose=TRUE)
#' fglm(y ~ ., "binomial", DF, stopFun=BIC, verbose=TRUE)
#' 
#' frs(y, x, glm, family="binomial", verbose=TRUE)
#' 
NULL
#> NULL



#' @rdname fglm
#' @order 1
#' @export
fglm <- function(formula, family = gaussian, data, weights, subset,
                na.action, start = NULL, etastart, mustart, offset,
                control = list(...), model = TRUE, method = "glm.fit",
                x = FALSE, y = TRUE, singular.ok = TRUE, contrasts = NULL, ...,
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
        frs(yvec = yvec, xmat = xmat, fitFun = glm,
            family = family,
            weights = if (!missing(weights)) weights else NULL,
            subset = if (!missing(subset)) subset else NULL,
            na.action = if (!missing(na.action)) na.action else NULL,
            start = start,
            etastart = if (!missing(etastart)) etastart else NULL,
            mustart = if (!missing(mustart)) mustart else NULL,
            offset = if (!missing(offset)) offset else NULL,
            control = control,
            model = model,
            method = method,
            x = x,
            y = y,
            singular.ok = singular.ok,
            contrasts = contrasts,
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



#' @rdname fglm
#' @order 2
#' @export
fglm.fit <- function(x, y, weights = rep.int(1, NROW(y)), start = NULL, etastart = NULL,
                      mustart = NULL, offset = rep.int(0, NROW(y)), family = gaussian(),
                      control = list(), intercept = TRUE, singular.ok = TRUE,
                      selectFun = "logLik", stopFun = "EBIC",
                      keep = NULL, maxK = NULL, verbose = FALSE) {

    n <- NROW(x)
    p <- NCOL(x)

    if (is.character(selectFun) && selectFun == "logLik")
        selectFun <- function(object) {
            class(object) <- "glm"
            return(logLik(object))
        }
    
    if (is.character(stopFun) && stopFun == "EBIC")
        stopFun <- function(object) {
            class(object) <- "glm"

            ebic.r <- max( 0.0, 1.0 - log(n) / (2.0 * log(p)) )
            stopifnot( ebic.r >= 0 )

            dof <- length(coef(object))
            ebic.penalty <- 2.0 * ebic.r * lchoose(p - length(keep), dof - length(keep))
            stopifnot( is.finite(ebic.penalty) )

            return(BIC(object) + ebic.penalty)
        }

    return(
        frs(yvec = y, xmat = x, fitFun = glm.fit,
            weights = weights,
            start = start,
            etastart = etastart,
            mustart = mustart,
            offset = offset,
            family = family,
            control = control,
            intercept = intercept,
            singular.ok = singular.ok,
            use.formula = FALSE,
            selectFun = selectFun,
            stopFun = stopFun,
            keep = keep,
            maxK = maxK,
            verbose = verbose
        )
    )
}
