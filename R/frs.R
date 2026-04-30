#' @name frs
#' @title Forward Regression Selection Framework
#' 
#' @description `frs()` is a generic workhorse function of forward regression selection for parametric regression.
#' 
#' @param yvec Response vector. See [pboost].
#' @param xmat Numeric feature matrix. See [pboost].
#' @param fitFun See [pboost].
#' @param ... See [pboost].
#' @param use.formula Whether to use formula interface for model fitting. Default to `TRUE`.
#' When `use.formula=TRUE`, the the model fitting function has the form `fitFun(formula, data, ...)`; otherwise, `fitFun(x, y, ...)`.
#' @param use.intercept Include intercept in the model fitting? Valid only when `use.formula=TRUE`.
#' @param selectFun A function to evaluate the importance of an unselected feature when it is added to current model.
#' The default is `logLik`, meaning that the feature with the largest post-added log-likelihood is identified as the next one to be added to the model.
#' Note that `selectFun` is only used for selecting features, and it does not affect the stopping rule of forward regression selection, which is determined by `stopFun`.
#' @param stopFun See [pboost].
#' @param maxK See [pboost].
#' @param keep Vector of indices or feature names initial features to include.
#' @param verbose See [pboost].
#' 
#' @return Model object fitted on the selected features.
#' 
#' @seealso [pboost::fbetareg], [pboost::fglm], [pboost::flm], [pboost::frq].
#' 
#' @examples
#' set.seed(2026)
#' n <- 300
#' p <- 100
#' x <- matrix(rnorm(n*p), n)
#' eta <- drop( x[, 1:3] %*% runif(3, 1.0, 1.5) )
#' 
#' ## ---------- linear model ----------
#' y <- rnorm(n, eta)
#' frs(y, x, lm, verbose=TRUE)
#' 
#' ## ---------- logistic model ----------
#' y <- rbinom(n, 1, 1/(1+exp(-eta)))
#' frs(y, x, glm, family="binomial", verbose=TRUE)
#' 
NULL
#> NULL


#' @rdname frs
#' @order 1
#' @export
frs <- function(yvec, xmat, fitFun, ...,
                use.formula = TRUE,
                use.intercept = TRUE,
                selectFun = logLik,
                stopFun = "EBIC",
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

    lhs <- deparse(substitute(yvec))
    rhs <- paste(c(as.integer(use.intercept), keep), collapse=" + ")
    fml <- as.formula(paste(c(lhs, rhs), collapse=" ~ "))
    if (use.formula) {
        egg <- fitFun(
            formula = fml,
            data = data.frame(yvec, xmat[, all.vars(fml[[3]]), drop=FALSE]),
            ...
        )
    } else {
        egg <- fitFun(x = xmat[, all.vars(fml[[3]]), drop=FALSE], y = yvec, ...)
    }
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

        importantanceCand <- rep(-Inf, length(candidates))
        levelCand <- rep(NA_real_, length(candidates))
        for (iC in seq_along(candidates)) {
            tmp.fml <- update(fml, sprintf(". ~ . + %s", candidates[iC]))
            if (use.formula) {
                tmp.egg <- fitFun(
                    formula = tmp.fml,
                    data = data.frame(yvec, xmat[, all.vars(tmp.fml[[3]]), drop=FALSE]),
                    ...
                )
            } else {
                tmp.egg <- fitFun(x = xmat[, all.vars(tmp.fml[[3]]), drop=FALSE], y = yvec, ...)
            }
            importantanceCand[iC] <- selectFun(tmp.egg)
            levelCand[iC] <- stopFun(tmp.egg)
        }

        j.star <- which.max(importantanceCand)
        stopifnot( length(j.star) == 1 )
        stopifnot( !(candidates[j.star] %in% all.vars(fml[[3]])) )
        showiter(verbose, candidates[j.star], levelCand[j.star])

        if (is.null(maxK) && (levelCand[j.star] >= level) )
            break

        level <- levelCand[j.star]
        fml <- update(fml, sprintf(". ~ . + %s", candidates[j.star]))

        if ( !is.null(maxK) && (length(all.vars(fml[[3]])) >= maxK) )
            break
    }

    if (use.formula) {
        egg <- fitFun(
            formula = fml,
            data = data.frame(yvec, xmat[, all.vars(fml[[3]]), drop=FALSE]),
            ...
        )
    } else {
        egg <- fitFun(x = xmat[, all.vars(fml[[3]]), drop=FALSE], y = yvec, ...)
    }
    class(egg) <- c("frs", class(egg))
    return(egg)
}
