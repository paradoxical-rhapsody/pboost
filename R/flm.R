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
#' p <- 50
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

    fml <- Formula(formula)
    mf <- model.frame(fml, data=data)

    yvec <- model.part(fml, data=mf, lhs=1, drop=TRUE)
    xmat <- model.part(fml, data=mf, rhs=1) |> as.matrix()
    use.intercept <- attr(terms(mf), "intercept") == 1L


    mc <- match.call(expand.dots = TRUE)
    provided_args <- as.list(mc)[-1]
    provided_args$formula <- NULL
    provided_args$data <- NULL

    args <- list(
        fitFun = lm,
        yvec = yvec,
        xmat = xmat,
        use.intercept = use.intercept
    )
    args <- c(args, provided_args)

    return(do.call(frs, args))
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


    mc <- match.call(expand.dots = TRUE)
    provided_args <- as.list(mc)[-1]
    provided_args <- provided_args[!(names(provided_args) %in% c("selectFun", "stopFun"))]
    provided_args$x <- NULL
    provided_args$y <- NULL

    args <- list(
        fitFun = lm.fit,
        yvec = y,
        xmat = x,
        selectFun = selectFun,
        stopFun = stopFun,
        use.formula = FALSE
    )
    args <- c(args, provided_args)

    return(do.call(frs, args))
}
