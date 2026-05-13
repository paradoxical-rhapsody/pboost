#' @name psar
#' 
#' @title Profile Boosting for Spatial Auto-regressive Model
#' 
#' @description Profile boosting variable selection for spatial auto-regressive (SAR) model
#' \deqn{\bm{y} = \rho W \bm{y} + X \beta + \bm{\varepsilon},}
#' where \eqn{\bm{y}, \bm{\varepsilon} \in \mathbb{R}^n}, \eqn{X \in \mathbb{R}^{n \times p}},
#' \eqn{\rho \in (-1, 1)}.
#' 
#' @param x Numeric feature matrix.
#' @param y Response vector.
#' @param w Weight matrix (row-sum scaled being one).
#' 
#' @param stopFun Parameter passed to [pboost::pboost].
#' @param keep Parameter passed to [pboost::pboost].
#' @param maxK Parameter passed to [pboost::pboost].
#' @param verbose Parameter passed to [pboost::pboost].
#' 
#' @return Model object fitted on the selected features.
#' 
#' @examples
#' set.seed(2026)
#' 
#' w <- set_rook_matrix(9, 9)
#' 
#' n <- NROW(w)
#' p <- 300
#' x <- matrix(rnorm(n*p), n) %*% chol(0.7^abs(outer(1:p, 1:p, "-")))
#' eta <- drop(x[, 1:3] %*% runif(3, 1.5, 2.0))
#' 
#' rho0 <- 0.2
#' sig0 <- 1.0
#' y <- solve(diag(n) - rho0 * w, rnorm(n, eta, sd=sig0)) |> drop()
#' 
#' ## ---------- pboost ----------
#' system.time( egg <- psar(x, y, w, verbose=TRUE) )
#' y.tilde <- (diag(NROW(x)) - egg[["rho"]] * w) %*% y
#' 
#' beta.hat <- egg[["beta"]]
#' idx <- as.integer(sub("[[:alpha:]]", "", names(beta.hat)))
#' sig2.hat <- mean( (y.tilde - drop(x[, idx, drop=FALSE] %*% beta.hat))^2 )
#' print( egg[["sig2"]] - sig2.hat )
#' 
#' ## ---------- frs ----------
#' fsar(x, y, w, verbose=TRUE)
#' 
NULL
#> NULL



#' @rdname psar
#' @export
psar <- function(x, y, w,
                stopFun = "EBIC",
                keep = NULL, maxK = NULL, verbose = FALSE) {
    stopifnot( is.matrix(x) )
    stopifnot( all( abs(rowSums(w) - 1.0) <= 1e-6 ) )
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

    egg <- pboost(
        yvec = y,
        xmat = x,
        fitFun = sar.fit,
        scoreFun = residuals,
        stopFun = stopFun,
        w = w,
        use.formula = FALSE,
        keep = keep,
        maxK = maxK,
        verbose = verbose
    )
    class(egg) <- c("psar", class(egg))
    return(egg)
}
