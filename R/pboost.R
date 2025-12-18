#' @name pboost
#' @title Profile Boosting Framework
#' 
#' @description
#' `pboost` is the generic workhorse function of profile boosting 
#' framework for parametric regression.
#' 
#' @param formula An object of class [formula] of the form `LHS ~ RHS`, 
#'   where the right-hand side (RHS) specifies the candidate features 
#'   for the linear predictor \eqn{\eta = \sum_j \beta_j x_j}.
#'   
#'   The following restrictions and recommendations apply:
#'   \itemize{
#'     \item All variables appearing on the RHS must be numeric in the supplied `data`.
#'     \item For computational efficiency, each term on the RHS must correspond to a
#'           single column in the resulting model matrix. Supported expressions include 
#'           main effects (`x1`), interactions (`x1:x2`), and simple transformations 
#'           (\code{log(x1)}, \code{I(x1^2)}, etc.). 
#'           Complex terms that expand into multiple columns—such as \code{poly(x, degree)},
#'           \code{bs(x)}, or \code{ns(x)}—are **not supported**.
#'     \item Offset terms should not be included in the formula. Instead, provide them 
#'           via the dedicated \code{offset} argument of \code{fitFun}.
#'   }
#' @param data An data frame containing the variables in the model.
#' @param fitFun Function to fit the empirical risk function in
#'    the form `fitFun(formula, data, ...)`.
#' @param scoreFun Function to compute the derivative of empirical
#'    risk function in the form `scoreFun(object)`, where `object` is
#'    returned by `fitFun`.
#'    `scoreFun` should return a vector with the same length of `y` in `data`.
#' @param stopFun Stopping rule for profile boosting, which has the form
#'    `stopFun(object)` to evaluate the performance of model `object` returned
#'    by `fitFun`, such as [EBIC] or [BIC].
#' @param ... Additional arguments to be passed to `fitFun`.
#' @param maxK Maximal number of identified features.
#'    If `maxK` is specified, it will suppress `stopFun`, saying that the
#'    profile boosting continues until the procedure identifies `maxK` features.
#'    The pre-specified features in `keep` are counted toward `maxK`.
#' @param keep Initial set of features that are included in model fitting.
#'    **If `keep` is specified, it should also be fully included in the RHS
#'    of `formula`.**
#' @param verbose Print the procedure path?
#' 
#' @return Model object fitted on the selected features.
#' 
#' @examples
#' \donttest{
#' set.seed(2025)
#' n <- 200
#' p <- 300
#' x <- matrix(rnorm(n*p), n)
#' eta <- drop(x[, 1:3] %*% runif(3, 1.0, 1.5))
#' y <- rbinom(n, 1, 1/(1+exp(-eta)))
#' DF <- data.frame(y, x)
#' 
#' scoreLogistic <- function(object) {
#'     eta.hat <- object[["linear.predictors"]]
#'     return(object[["y"]] - 1/(1+exp(-eta.hat)))
#' }
#' 
#' ( result <- pboost(y~., DF, glm, scoreLogistic, EBIC, family="binomial") )
#' 
#' ## Extract the selected features
#' attr(terms(formula(result), data=DF), "term.labels")
#' }
#' 
NULL




#' @rdname pboost
#' @order 1
#' @export
pboost <- function(formula, data, fitFun, scoreFun, stopFun, ...,
                   keep = NULL, maxK = NULL, verbose = FALSE) {

    formula <- as.Formula(formula)

    ## --- `all.vars`: original variables ---
    if (!all(sapply(
        data[1, all.vars(delete.response(terms(formula(formula, rhs=1L), data=data)))],
        is.numeric
    )))
        stop("'formula' contains non-numeric feature(s).")
    
    ## --- `attr(terms, "term.labels")`: features, such as `log(x)`, `x1:x2` ---
    xnames <- attr(terms(formula(formula, rhs=1L), data=data), "term.labels") |> # features
        gsub(pattern=":", replacement="*", fixed=TRUE)
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
    intercept <- attr(terms(formula(formula, rhs=1L), data=data), "intercept")
    rhs <- paste(c(intercept, keep), collapse=" + ")
    fml <- update(formula, as.Formula(paste(c(lhs, "~", rhs), collapse=" ")))
    stopifnot( intercept %in% 0:1 )
    stopifnot( !is.null(fml) )
    if (!is.null(keep))
        stopifnot( setequal(attr(terms(fml, data=data), "term.labels"), keep) )

    object <- fitFun(formula=fml, data=data, ...)
    obj.level <- stopFun(object)
    showiter(verbose, level=obj.level)

    while (TRUE) {

        stopifnot( all( attr(terms(formula(fml, rhs=1L), data=data), "term.labels") %in% xnames ) )

        x.star <- setdiff(xnames, attr(terms(formula(fml, rhs=1L), data=data), "term.labels")) |>
            vapply(function(expr) with(data, eval(parse(text=expr))),
                   FUN.VALUE=numeric(NROW(data))) |>
            crossprod(scoreFun(object)) |>
            drop() |> abs() |> which.max() |> names()
        fml.tmp <- update(fml, as.Formula(paste(c(". ~ .", x.star), collapse=" + ")))
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
