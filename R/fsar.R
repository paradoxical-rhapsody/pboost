#' @name fsar
#' @title Forward Stepwise Spatial Auto-Regressive Model
#' 
#' @description Forward stepwise strategy for spatial auto-regressive model.
#' 
#' @param x Matrix of covariates.
#' @param y Vector of response.
#' @param w Weight matrix (row-sum scaled being one).
#' @param formula Parameters passed to [spatialreg::lagsarlm].
#' @param data Parameters passed to [spatialreg::lagsarlm].
#' @param listw Parameters passed to [spatialreg::lagsarlm].
#' @param na.action Parameters passed to [spatialreg::lagsarlm].
#' @param Durbin Parameters passed to [spatialreg::lagsarlm].
#' @param type Parameters passed to [spatialreg::lagsarlm].
#' @param method Parameters passed to [spatialreg::lagsarlm].
#' @param quiet Parameters passed to [spatialreg::lagsarlm].
#' @param zero.policy Parameters passed to [spatialreg::lagsarlm].
#' @param interval Parameters passed to [spatialreg::lagsarlm].
#' @param tol.solve Parameters passed to [spatialreg::lagsarlm].
#' @param trs Parameters passed to [spatialreg::lagsarlm].
#' @param control Parameters passed to [spatialreg::lagsarlm].
#' 
#' @param selectFun Parameter passed to [pboost::frs].
#' @param stopFun Parameter passed to [pboost::frs].
#' @param keep Parameter passed to [pboost::frs].
#' @param maxK Parameter passed to [pboost::frs].
#' @param verbose Parameter passed to [pboost::frs].
#' 
#' @return Model object fitted on the selected features.
#' 
NULL
#> NULL



#' @rdname fsar
#' @order 1
#' @export
flagsarlm <- function(formula, data = list(), listw, na.action, Durbin, type,
                    method = "eigen", quiet = NULL, zero.policy = NULL, interval = NULL,
                    tol.solve = .Machine$double.eps, trs = NULL, control = list(),
                    selectFun = logLik, stopFun = "EBIC",
                    keep = NULL, maxK = NULL, verbose = FALSE) {

    fml <- Formula(formula)
    mf <- model.frame(fml, data=data)

    yvec <- model.part(fml, data=mf, lhs=1, drop=TRUE)
    xmat <- model.part(fml, data=mf, rhs=1) |> as.matrix()
    use.intercept <- attr(terms(mf), "intercept") == 1L


    mc <- match.call(expand.dots = TRUE)
    provided_args <- as.list(mc)[-1]
    provided_args <- provided_args[!(names(provided_args) %in% c("formula", "data"))]

    fitFun <- function(...) {
        dots <- list(...)
        return(do.call(lagsarlm, dots))
    }

    args <- list(
        fitFun = fitFun,
        yvec = yvec,
        xmat = xmat,
        use.intercept = use.intercept
    )
    args <- c(args, provided_args)

    return(do.call(frs, args))
}



#' @rdname fsar
#' @order 2
#' @export
fsar.fit <- function(x, y, w,
                selectFun = logLik,
                stopFun = "EBIC",
                keep = NULL, maxK = NULL, verbose = FALSE) {
    stopifnot( is.matrix(x) )
    stopifnot( all( abs(rowSums(w) - 1.0) <= 1e-12 ) )
    stopifnot( all(w >= 0.0) )
    stopifnot( !anyNA(x) )
    stopifnot( !anyNA(y) )
    stopifnot( !anyNA(w) )

    n <- NROW(x)
    p <- NCOL(x)

    if (is.character(stopFun) && stopFun == "EBIC")
        stopFun <- function(object) {
            ebic.r <- max( 0.0, 1.0 - log(n) / (2.0 * log(p)) )
            stopifnot( ebic.r >= 0 )

            dof <- attr(object, "dof")
            ebic.penalty <- 2.0 * ebic.r * lchoose(p - length(keep), dof - length(keep))
            stopifnot( is.finite(ebic.penalty) )
            stopifnot( length(ebic.penalty) == 1 )

            return(BIC(object) + ebic.penalty)
        }

    mc <- match.call(expand.dots = TRUE)
    provided_args <- as.list(mc)[-1]
    provided_args <- provided_args[!(names(provided_args) %in% c("x", "y", "stopFun"))]

    args <- list(
        fitFun = sar.fit,
        yvec = y,
        xmat = x,
        stopFun = stopFun,
        use.formula = FALSE
    )
    args <- c(args, provided_args)

    return(do.call(frs, args))
}