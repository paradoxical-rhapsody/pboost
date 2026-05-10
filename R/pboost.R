#' @name pboost
#' @title Profile Boosting Framework
#' 
#' @description
#' `pboost()` is the generic workhorse function of profile boosting
#' framework for parametric regression.
#' 
#' @param yvec Response vector.
#' @param xmat Numeric feature matrix.
#' @param fitFun Function to fit the empirical risk function in
#'    the form `fitFun(formula, data, ...)`.
#' @param scoreFun Function to compute the derivative, denoted by \eqn{\frac{\partial \ell(y, \eta)}{\partial \eta}}, of empirical risk function in the form `scoreFun(object)`, where `object` is returned by `fitFun`.
#' `scoreFun()` should return a vector with the same length of `y` in `data`.
#' @param stopFun Stopping rule for profile boosting, which has the form
#'    `stopFun(object)` to evaluate the performance of model `object` returned
#'    by `fitFun`, such as [EBIC] or [BIC].
#' @param ... Additional arguments to be passed to `fitFun`.
#' @param use.formula Whether to use formula interface for model fitting. Default to `TRUE`.
#' When `use.formula=TRUE`, the the model fitting function has the form `fitFun(formula, data, ...)`; otherwise, `fitFun(x, y, ...)`.
#' @param use.intercept Include intercept in the model fitting? Valid only when `use.formula=TRUE`.
#' @param maxK Maximal number of identified features.
#'    If `maxK` is specified, it will suppress `stopFun`, saying that the
#'    profile boosting continues until the procedure identifies `maxK` features.
#'    The pre-specified features in `keep` are counted toward `maxK`.
#' @param keep Vector of indices or feature names initial features to include.
#' @param verbose Print the procedure path?
#' 
#' @return Model object fitted on the selected features.
#' 
#' @seealso [pboost::pbetareg], [pboost::pglm], [pboost::plm], [pboost::prq], [pboost::psar].
#' 
#' @examples
#' set.seed(2026)
#' n <- 200
#' p <- 300
#' x <- matrix(rnorm(n*p), n)
#' eta <- drop(x[, 1:3] %*% runif(3, 1.0, 1.5))
#' y <- rbinom(n, 1, 1/(1+exp(-eta)))
#' 
#' scoreLogistic <- function(object) {
#'     eta.hat <- object[["linear.predictors"]]
#'     return(object[["y"]] - 1/(1+exp(-eta.hat)))
#' }
#' 
#' ( result <- pboost(y, x, glm, scoreLogistic, family="binomial") )
#' 
#' all.vars(formula(result)[[3]])
#' 
NULL




#' @rdname pboost
#' @order 1
#' @export
pboost <- function(yvec, xmat, fitFun, scoreFun, stopFun = "EBIC", ...,
                    use.formula = TRUE,
                    use.intercept = TRUE,
                    keep = NULL, maxK = NULL, verbose = FALSE) {
    stopifnot( NROW(yvec) == NROW(xmat) )
    stopifnot( is.matrix(xmat), is.numeric(xmat) )

    n <- NROW(xmat)
    p <- NCOL(xmat)

    xnames <- colnames(xmat)
    if (is.null(xnames)) {
        xnames <- paste0("x", seq_len(p))
        colnames(xmat) <- xnames
    }
    stopifnot( length(xnames) == p )

    if ( is.character(stopFun) && use.formula && stopFun == "EBIC")
        stopFun <- function(object){
            ebic.r <- max( 0.0, 1.0 - log(n) / (2.0 * log(p)) )
            stopifnot( ebic.r >= 0 )

            dof <- length(all.vars(formula(object)[[3]]))
            ebic.penalty <- 2.0 * ebic.r * lchoose(p - length(keep), dof - length(keep))
            stopifnot( is.finite(ebic.penalty) )

            # return(-2*logLik(object) + dof * log(n) + ebic.penalty)
            return(BIC(object) + ebic.penalty)
        }

    stopifnot( is.character(keep) || is.numeric(keep) || is.null(keep) )
    if (!is.null(keep)) {
        if (is.character(keep)) {
            stopifnot( all(keep %in% xnames) )
        } else {
            stopifnot( all(keep %in% seq_len(p)) )
            keep <- xnames[keep]
        }
    }
    stopifnot( is.character(keep) || is.null(keep) )

    if (!is.null(maxK))
        maxK <- min( maxK, p, n-1 )

    showiter <- function(verbose, x.star, level) {
        if (verbose)
            if (missing(x.star))
                message(sprintf("Initial model: level=%.3f", level))
            else 
                message(sprintf("Adding %s: level=%.3f", x.star, level))
    }

    fit_formula <- function(fml) {
        if (use.formula) {
            return(fitFun(
                formula = fml,
                data = data.frame(yvec, xmat[, all.vars(fml[[3]]), drop=FALSE]),
                ...
            ))
        } else {
            return(fitFun(x = xmat[, all.vars(fml[[3]]), drop=FALSE], y = yvec, ...))
        }
    }

    lhs <- "yvec"
    rhs <- paste(c(as.integer(use.intercept), keep), collapse=" + ")
    fml <- as.formula(paste(c(lhs, rhs), collapse=" ~ "))
    egg <- fit_formula(fml)
    level <- stopFun(egg)
    showiter(verbose, level=level)

    if ( !is.null(maxK) && !is.null(keep) && (length(keep) >= maxK) ) {
        warning("'length(keep)' is not less than 'maxK', thus fitting on 'keep'.")
        return(egg)
    }

    while (TRUE) {
        candidates <- setdiff(xnames, all.vars(fml[[3]]))
        if (length(candidates) == 0)
            break

        x.star <- crossprod(xmat[, candidates, drop=FALSE], scoreFun(egg)) |>
            drop() |> abs() |> which.max() |> names()
        stopifnot( length(x.star) == 1 )
        stopifnot( !(x.star %in% all.vars(fml[[3]])) )

        tmp.fml <- update(fml, sprintf(". ~ . + %s", x.star))
        tmp.egg <- fit_formula(tmp.fml)
        tmp.level <- stopFun(tmp.egg)
        showiter(verbose, x.star, tmp.level)

        if (is.null(maxK) && (tmp.level >= level) )
            break

        fml <- tmp.fml
        egg <- tmp.egg
        level <- tmp.level

        if ( !is.null(maxK) && (length(all.vars(fml[[3]])) >= maxK) )
            break
    }

    egg <- fit_formula(fml)
    class(egg) <- c("frs", class(egg))
    return(egg)
}
