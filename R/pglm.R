#' @name pglm
#' @title Profile Boosting for Generalized Linear Models.
#' 
#' @description
#' [pglm] inherits the usage of the built-in function [glm].
#' 
#' @param formula See [pboost].
#' @param data See [pboost].
#' @param family Parameters passed to [glm].
#' @param weights Parameters passed to [glm].
#' @param subset Parameters passed to [glm].
#' @param na.action Parameters passed to [glm].
#' @param start Parameters passed to [glm].
#' @param etastart Parameters passed to [glm].
#' @param mustart Parameters passed to [glm].
#' @param offset Parameters passed to [glm].
#' @param control Parameters passed to [glm].
#' @param model Parameters passed to [glm].
#' @param method Parameters passed to [glm].
#' @param x Parameters passed to [glm].
#' @param y Parameters passed to [glm].
#' @param singular.ok Parameters passed to [glm].
#' @param contrasts Parameters passed to [glm].
#' @param ... Parameters passed to [glm].
#' @param stopFun Parameters passed to [pboost].
#' @param keep Parameters passed to [pboost].
#' @param maxK Parameters passed to [pboost].
#' @param verbose Parameters passed to [pboost].
#' 
#' @return An `glm` model object fitted on the selected features.
#' 
#' @references
#' Zengchao Xu, Shan Luo and Zehua Chen (2022). Partial profile score feature selection
#' in high-dimensional generalized linear interaction models. Statistics and Its Interface.
#' \doi{10.4310/21-SII706}
#' 
#' @examples
#' set.seed(2025)
#' n <- 300
#' p <- 200
#' x <- matrix(rnorm(n*p), n)
#' 
#' eta <- drop( x[, 1:3] %*% runif(3, 1.0, 1.5) )
#' y <- rbinom(n, 1, 1/(1+exp(-eta)))
#' DF <- data.frame(y, x)
#' 
#' pglm(y ~ ., "binomial", DF, verbose=TRUE)
#' pglm(y ~ ., "binomial", DF, stopFun=BIC, verbose=TRUE)
#' 
#' scoreLogistic <- function(object) {
#'    eta.hat <- object[["linear.predictors"]]
#'    return(object[["y"]] - 1/(1+exp(-eta.hat)))
#' }
#' pboost(y ~ ., DF, glm, scoreLogistic, EBIC, family="binomial", verbose=TRUE)
#' 
NULL
#> NULL



#' @rdname pglm
#' @order 1
#' @export
pglm <- function(
    formula, family = gaussian, data, weights, subset,
    na.action, start = NULL, etastart, mustart, offset,
    control = list(...), model = TRUE, method = "glm.fit",
    x = FALSE, y = TRUE, singular.ok = TRUE, contrasts = NULL, ...,
    stopFun = EBIC, keep = NULL, maxK = NULL, verbose = FALSE) {
    stopifnot( !missing(formula) )
    stopifnot( !missing(data) )

    cl <- match.call()

    glm_template <- cl
    glm_template$stopFun <- NULL
    glm_template$keep <- NULL
    glm_template$maxK <- NULL
    glm_template$verbose <- NULL
    glm_template[[1L]] <- quote(glm)

    required_paras <- c("data", "weights", "subset", "na.action",
                        "etastart", "mustart", "offset")
    for (ipara in required_paras)
        if (!is.null(cl[[ipara]]))
            glm_template[[ipara]] <- eval(cl[[ipara]], envir = parent.frame())

    fitFun <- function(formula, data) {
        call <- glm_template
        call$formula <- formula
        call$data <- data
        return( eval(call, parent.frame()) )
    }

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


    return(pboost(formula, data, fitFun, scoreFun, stopFun,
                  keep = keep, maxK = maxK, verbose = verbose))
}



#' @rdname EBIC
#' @export
EBIC.glm <- function(object, p, p.keep, ...) {
    stopifnot( inherits(object, "glm") )

    if (missing(p))
        p <- get("p", envir=parent.frame())
    if (missing(p.keep))
        p.keep <- get("p.keep", envir=parent.frame())

    dof <- attr(logLik(object), "df")
    ebic.r <- max( 0.0, 1.0 - log(nobs(object)) / (2.0*log(p)) )
    ebic.penalty <- ifelse(
        ebic.r <= 0.0,
        0.0,
        2.0 * ebic.r * lchoose(p - p.keep, dof - p.keep)
    )

    # stopifnot( !is.nan(ebic.penalty) )
    stopifnot( is.finite(ebic.penalty) )
    return(BIC(object) + ebic.penalty)
}
