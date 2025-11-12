#' @name pboost02
#' @title Profile Boosting Framework
#' 
#' @description
#' `pboost` is the generic framework of profile boosting 
#' for parametric regression.
#' 
#' @param formula An object of class ``[formula]'' in the form `y ~ x`, where the
#' RHS contains all of candidate features. 
#' @param data An data frame containing the variables in the model.
#' @param fitFun Function to fit the empirical risk function in 
#' the form `fitFun(formula, data, ...)`.
#' @param scoreFun Function to compute the derivative of empirical 
#' risk function in the form `scoreFun(object)`, where `object` is 
#' returned by `fitFun`.
#' `scoreFun` should return a vector with the same length of `y` in `data`.
#' @param stopFun Stopping rule for profile boosting, which has the form 
#' `stopFun(object)` to evaluate the performance of model `object` returned
#' by `fitFun`, such as [BIC].
#' @param ... Additional arguments to be passed to `fitFun`.
#' @param maxK Maximal number of identified features. 
#' If `maxK` is specified, it will supress `stopFun`, saying that the 
#' profile boosting continues until the procedure identifies `maxK` features.
#' The pre-specified features in `keep` are counted toward `maxK`.
#' @param keep Initial set of features that are included in model fitting.
#' **If `keep` is specified, it should also be fully included in the RHS of `formula`.**
#' @param verbose Print the procedure path?
#' 
#' @return Model object fitted on the selected features.
#' 
#' @examples
#' # scoreLogistic <- function(object) {
#' #     D0 <- object$family$mu.eta(object$linear.predictor)
#' #     S0 <- object$family$variance(object$fitted.values)
#' #     score <- drop(D0 / S0 * (y-fitted(object)))
#' #     return(score)
#' # }
#' 
#' scoreLogistic <- function(object) {
#'     eta.hat <- object[["linear.predictors"]]
#'     return(object[["y"]] - 1/(1+exp(-eta.hat)))
#' }
#' 
#' stopLogistic <- function(object) {
#'     nobs <- attr(logLik(object), "nobs")
#'     dof <- attr(logLik(object), "df")
#'     ebic.r <- max(0.0, 1.0-log(nobs)/(2.0*log(p)))
#'     ebic.penalty <- 2 * ebic.r * lchoose(p, dof)
#'     return(BIC(object) + ebic.penalty)
#' }
#' 
#' 
#' set.seed(2025)
#' n <- 200
#' p <- 500
#' x <- matrix(rnorm(n*p), n)
#' eta <- drop(x[, 1:3] %*% runif(3, 1.0, 1.5))
#' y <- rbinom(n, 1, 1/(1+exp(-eta)))
#' DF <- data.frame(y, x)
#' 
#' ( result <- pboost(y~., DF, glm, scoreLogistic, stopLogistic, family="binomial") )
#' 
#' ## extract the selected features
#' attr(terms(formula(result), data=DF), "term.labels")
#' 
#' @noRd
NULL




#' @rdname pboost02
#' @order 1
#' @noRd
pboost02 <- function(formula, data, fitFun, scoreFun, stopFun, ..., 
                   keep=NULL, maxK=NULL, verbose=FALSE){

    xnames <- attr(terms(formula, data=data), "term.labels")
    p <- length(xnames)
    p.keep <- length(keep)

    if (!is.null(keep))
        stopifnot( all(keep %in% xnames) )

    if (!is.null(maxK))
        maxK <- min( maxK, length(xnames), NROW(data)-1, NCOL(data)-1 )

    showiter <- function(verbose, x.star, level=obj.level) {
        if (verbose) 
            if (missing(x.star))
                message(sprintf("Initial model with level=%.3f", level))
            else 
                message(sprintf("Adding %s: stopping level=%.3f", x.star, level))
    }

    # Initialzation: preserve the original formula's LHS and its intercept setting
    #   then add any `keep` variables to the RHS.
    # Extract LHS as text (preserve expressions like log(y), cbind(y1,y2), etc.)
    lhs <- paste(deparse(formula[[2]]), collapse=" ")
    intercept <- attr(terms(formula, data=data), "intercept")
    rhs <- paste(c(intercept, keep), collapse=" + ")
    fml <- as.formula(paste(c(lhs, "~", rhs), collapse=" "))
    stopifnot( intercept %in% 0:1 )
    stopifnot( !is.null(fml) )
    if (!is.null(keep))
        stopifnot( setequal(attr(terms(fml, data=data), "term.labels"), keep) )

    object <- fitFun(formula=fml, data=data, ...)
    obj.level <- stopFun(object)
    showiter(verbose, level=obj.level)

    while (TRUE) {

        stopifnot( all( attr(terms(fml, data=data), "term.labels") %in% xnames ) )
        candidates <- setdiff(xnames, attr(terms(fml, data=data), "term.labels"))
        fml.cand <- as.formula(paste(c("~ 0", candidates), collapse=" + "))

        # profilescore: crossprod( model.matrix(fml.cand, data), scoreFun(object) )
        x.star <- drop(crossprod( model.matrix(fml.cand, data), scoreFun(object) )) |> 
                    abs() |> which.max() |> names()
        fml.tmp <- update(fml, paste(c(". ~ .", x.star), collapse=" + "))
        stopifnot( !identical(fml, fml.tmp) )

        object <- fitFun(formula=fml.tmp, data=data, ...)
        # if (!object$converged) break
        obj.level.tmp <- stopFun(object)
        showiter(verbose, x.star, obj.level.tmp)

        dof <- length(attr(terms(fml.tmp, data=data), "term.labels"))
        if (is.null(maxK)) {
            if (obj.level.tmp > obj.level) break
        } else {
            if (dof > maxK) break
        }

        obj.level <- obj.level.tmp
        fml <- fml.tmp
    }

    # return(fml)
    return(fitFun(formula=fml, data=data, ...))
}
