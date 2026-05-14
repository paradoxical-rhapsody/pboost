#' @name fglm
#' @title Forward Regression Selection for Generalized Linear Models
#' 
#' @description `fglm()` inherits the usage of [glm], and performs forward regression selection for generalized linear models.
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
#' @param selectFun Parameter passed to [pboost::frs].
#' @param stopFun Parameter passed to [pboost::frs].
#' @param keep Parameter passed to [pboost::frs].
#' @param maxK Parameter passed to [pboost::frs].
#' @param verbose Parameter passed to [pboost::frs].
#' 
#' @return A `glm` model object fitted on the selected features.
#' 
NULL
#> NULL



#' @rdname fglm
#' @order 1
#' @export
fglm <- function(formula, family = gaussian, data, weights, subset,
                na.action, start = NULL, etastart, mustart, offset,
                control = list(...), model = TRUE, method = "glm.fit",
                x = FALSE, y = TRUE, singular.ok = TRUE, contrasts = NULL, ...,
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

    args <- list(
        fitFun = glm,
        yvec = yvec,
        xmat = xmat,
        use.intercept = use.intercept
    )
    args <- c(args, provided_args)

    return(do.call(frs, args))
}



#' @rdname fglm
#' @order 2
#' @export
fglm.fit <- function(x, y, weights = rep.int(1, NROW(y)), start = NULL, etastart = NULL,
                      mustart = NULL, offset = rep.int(0, NROW(y)), family = gaussian(),
                      control = list(), intercept = TRUE, singular.ok = TRUE,
                      selectFun = "logLik", stopFun = "EBIC",
                      keep = NULL, maxK = NULL, verbose = FALSE) {

    n <- NROW(x)
    p <- NCOL(x)

    if (is.character(selectFun) && selectFun == "logLik")
        selectFun <- function(object) {
            class(object) <- "glm"
            return(logLik(object))
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
    provided_args <- provided_args[!(names(provided_args) %in% c("x", "y", "selectFun", "stopFun"))]

    args <- list(
        fitFun = glm.fit,
        yvec = y,
        xmat = x,
        selectFun = selectFun,
        stopFun = stopFun,
        use.formula = FALSE
    )
    args <- c(args, provided_args)

    return(do.call(frs, args))
}
