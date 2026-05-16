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
#' w <- set_rook_matrix(15, 15)
#' 
#' n <- NROW(w)
#' p <- 30
#' x <- matrix(rnorm(n*p), n) %*% chol(0.7^abs(outer(1:p, 1:p, "-")))
#' eta <- drop(x[, 1:3] %*% runif(3, 1.5, 2.0))
#' 
#' rho0 <- 0.2
#' sig0 <- 1.0
#' y <- solve(diag(n) - rho0 * w, rnorm(n, eta, sd=sig0)) |> drop()
#' 
#' ## ---------- pboost ----------
#' system.time( egg <- psar.fit(x, y, w, verbose=TRUE) )
#' y.tilde <- (diag(NROW(x)) - egg[["rho"]] * w) %*% y
#' 
#' beta.hat <- egg[["beta"]]
#' idx <- as.integer(sub("[[:alpha:]]", "", names(beta.hat)))
#' sig2.hat <- mean( (y.tilde - drop(x[, idx, drop=FALSE] %*% beta.hat))^2 )
#' print( egg[["sig2"]] - sig2.hat )
#' 
#' \dontrun{
#' system.time(
#' plagsarlm(y ~ ., data=data.frame(y, x), listw=spdep::mat2listw(w, style="W"), verbose=TRUE)
#' )
#' }
#' 
#' ## ---------- frs ----------
#' fsar.fit(x, y, w, verbose=TRUE)
#' 
#' \dontrun{
#' system.time(
#' flagsarlm(y ~ ., data=data.frame(y, x), listw=spdep::mat2listw(w, style="W"), verbose=TRUE)
#' )
#' }
#' 
NULL
#> NULL



#' @rdname psar
#' @order 1
#' @export
plagsarlm <- function(formula, data = list(), listw, na.action, Durbin, type,
                    method = "eigen", quiet = NULL, zero.policy = NULL, interval = NULL,
                    tol.solve = .Machine$double.eps, trs = NULL, control = list(),
                    stopFun = "EBIC",
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
        use.intercept = use.intercept,
        scoreFun = residuals
    )
    args <- c(args, provided_args)

    return(do.call(pboost, args))
}



#' @rdname psar
#' @order 2
#' @export
psar.fit <- function(x, y, w,
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

    
    mc <- match.call(expand.dots = TRUE)
    provided_args <- as.list(mc)[-1]
    provided_args <- provided_args[!(names(provided_args) %in% c("x", "y", "scoreFun", "stopFun"))]

    args <- list(
        fitFun = sar.fit,
        yvec = y,
        xmat = x,
        scoreFun = residuals,
        stopFun = stopFun,
        use.formula = FALSE
    )
    args <- c(args, provided_args)

    return(do.call(pboost, args))
}
