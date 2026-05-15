#' @name pglm
#' @title Profile Boosting for Generalized Linear Models.
#' 
#' @description
#' [pglm] inherits the usage of the built-in function [glm].
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
#' @param stopFun Parameter passed to [pboost::pboost].
#' @param keep Parameter passed to [pboost::pboost].
#' @param maxK Parameter passed to [pboost::pboost].
#' @param verbose Parameter passed to [pboost::pboost].
#' 
#' @return A `glm` model object fitted on the selected features.
#' 
#' @references
#' Zengchao Xu, Shan Luo and Zehua Chen (2022). Partial profile score feature selection
#' in high-dimensional generalized linear interaction models. Statistics and Its Interface.
#' \doi{10.4310/21-SII706}
#' 
#' @examples
#' set.seed(2026)
#' n <- 200
#' p <- 100
#' x <- matrix(rnorm(n*p), n)
#' 
#' eta <- drop( x[, 1:3] %*% runif(3, 1.0, 1.5) )
#' y <- rbinom(n, 1, 1/(1+exp(-eta)))
#' DF <- data.frame(y, x)
#' 
#' ## ---------- pboost ----------
#' pglm(y ~ ., "binomial", DF, verbose=TRUE)
#' pglm(y ~ ., "binomial", DF, stopFun=BIC, verbose=TRUE)
#' 
#' scoreLogistic <- function(object) {
#'    eta.hat <- object[["linear.predictors"]]
#'    return(object[["y"]] - 1/(1+exp(-eta.hat)))
#' }
#' (result <- pboost(x, y, glm, scoreLogistic, family="binomial", verbose=TRUE))
#' all.vars(formula(result)[[3]])
#' 
#' ## ---------- frs ----------
#' fglm(y ~ ., "binomial", DF, verbose=TRUE)
#' fglm(y ~ ., "binomial", DF, stopFun=BIC, verbose=TRUE)
#' 
#' frs(x, y, glm, family="binomial", verbose=TRUE)
#' 
NULL
#> NULL



#' @rdname pglm
#' @order 1
#' @export
pglm <- function(formula, family = gaussian, data, weights, subset,
                na.action, start = NULL, etastart, mustart, offset,
                control = list(...), model = TRUE, method = "glm.fit",
                x = FALSE, y = TRUE, singular.ok = TRUE, contrasts = NULL, ...,
                stopFun = "EBIC",
                keep = NULL, maxK = NULL, verbose = FALSE) {

    fml <- Formula(formula)
    mf <- model.frame(fml, data=data)

    yvec <- model.part(fml, data=mf, lhs=1, drop=TRUE)
    xmat <- model.part(fml, data=mf, rhs=1) |> as.matrix()
    use.intercept <- attr(terms(mf), "intercept") == 1L
    # use.intercept <- attr(terms(fml, data=data), "intercept") == 1L


    scoreFun <- function(object) {
        # score <- D0/S0*(y-fitted(obj))
        # profilescore <- drop(crossprod(x, score))
        # profilescore.sd <- drop(sqrt(crossprod(x*x, D0^2/S0)))
        # return( profilescore / profilescore.sd )
        D0 <- object$family$mu.eta(object$linear.predictors)
        S0 <- object$family$variance(object$fitted.values)
        # weights <- object[["prior.weights"]]
        # return( weights*D0/S0 * (y-fitted(obj)) )
        # return( weights * D0 / S0 * residuals(object, type="response") )
        return( D0 / S0 * residuals(object, type = "response") )
    }

    mc <- match.call(expand.dots = TRUE)
    provided_args <- as.list(mc)[-1]
    provided_args <- provided_args[!(names(provided_args) %in% c("formula", "data", "scoreFun"))]

    args <- list(
        fitFun = glm,
        yvec = yvec,
        xmat = xmat,
        use.intercept = use.intercept,
        scoreFun = scoreFun
    )
    args <- c(args, provided_args)

    return(do.call(pboost, args))
}



#' @rdname pglm
#' @order 2
#' @export
pglm.fit <- function(x, y, weights = rep.int(1, NROW(y)), start = NULL, etastart = NULL,
                      mustart = NULL, offset = rep.int(0, NROW(y)), family = gaussian(),
                      control = list(), intercept = TRUE, singular.ok = TRUE,
                      stopFun = "EBIC",
                      keep = NULL, maxK = NULL, verbose = FALSE) {

    n <- NROW(x)
    p <- NCOL(x)

    scoreFun <- function(object) {
        class(object) <- "glm"
        D0 <- object$family$mu.eta(object$linear.predictors)
        S0 <- object$family$variance(object$fitted.values)
        return( D0 / S0 * residuals(object, type = "response") )
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


    mc <- match.call(expand.dots = TRUE)
    provided_args <- as.list(mc)[-1]
    provided_args <- provided_args[!(names(provided_args) %in% c("x", "y", "scoreFun", "stopFun"))]

    args <- list(
        fitFun = glm.fit,
        yvec = y,
        xmat = x,
        scoreFun = scoreFun,
        stopFun = stopFun,
        use.formula = FALSE
    )
    args <- c(args, provided_args)

    return(do.call(pboost, args))
}
