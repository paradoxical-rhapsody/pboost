#' @title Forward Stepwise Spatial Auto-Regressive Model
#' 
#' @description Forward stepwise strategy for spatial auto-regressive model.
#' 
#' @param x Matrix of covariates.
#' @param y Vector of response.
#' @param w Weight matrix (row-sum scaled being one).
#' 
#' @param selectFun Parameter passed to [pboost::frs].
#' @param stopFun Parameter passed to [pboost::frs].
#' @param keep Parameter passed to [pboost::frs].
#' @param maxK Parameter passed to [pboost::frs].
#' @param verbose Parameter passed to [pboost::frs].
#' 
#' @return Model object fitted on the selected features.
#' 
#' @export
fsar <- function(x, y, w,
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

    stopFun <- function(object) {
        ebic.r <- max( 0.0, 1.0 - log(n) / (2.0 * log(p)) )
        stopifnot( ebic.r >= 0 )

        dof <- attr(object, "dof")
        ebic.penalty <- 2.0 * ebic.r * lchoose(p - length(keep), dof - length(keep))
        stopifnot( is.finite(ebic.penalty) )
        stopifnot( length(ebic.penalty) == 1 )

        return(BIC(object) + ebic.penalty)
    }

    egg <- frs(
        yvec = y,
        xmat = x,
        fitFun = sar.fit,
        w = w,
        use.formula = FALSE,
        selectFun = selectFun,
        stopFun = stopFun,
        keep = keep,
        maxK = maxK,
        verbose = verbose
    )
    class(egg) <- c("psar", class(egg))
    return(egg)
}