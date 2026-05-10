#' @name fbetareg
#' @title Forward Regression Selection for Beta Regression Models
#' 
#' @description
#' `fbetareg()` inherits the usage of the function [betareg::betareg], and performs forward regression selection for beta regression models.
#'
#' @param formula Parameter passed to [betareg::betareg].
#' @param data Parameter passed to [betareg::betareg].
#' @param subset Parameter passed to [betareg::betareg].
#' @param na.action Parameter passed to [betareg::betareg].
#' @param weights Parameter passed to [betareg::betareg].
#' @param offset Parameter passed to [betareg::betareg].
#' @param link Parameter passed to [betareg::betareg].
#' @param link.phi Parameter passed to [betareg::betareg].
#' @param type Parameter passed to [betareg::betareg].
#' @param dist Parameter passed to [betareg::betareg].
#' @param nu Parameter passed to [betareg::betareg].
#' @param control Parameter passed to [betareg::betareg].
#' @param model Parameter passed to [betareg::betareg].
#' @param y Parameter passed to [betareg::betareg].
#' @param x Parameter passed to [betareg::betareg].
#' @param ... Parameters passed to [betareg::betareg].
#' 
#' @param selectFun Parameter passed to [pboost::frs].
#' @param stopFun Parameter passed to [pboost::frs].
#' @param keep Parameter passed to [pboost::frs].
#' @param maxK Parameter passed to [pboost::frs].
#' @param verbose Parameter passed to [pboost::frs].
#' 
#' @return A `betareg` model object fitted on the selected features.
#' 
#' @export
fbetareg <- function(formula, data, subset, na.action, weights, offset,
                    link = c("logit", "probit", "cloglog", "cauchit", "log","loglog"),
                    link.phi = NULL, type = c("ML", "BC", "BR"),
                    dist = NULL, nu = NULL, control = betareg.control(...),
                    model = TRUE, y = TRUE, x = FALSE, ...,
                    selectFun = logLik, stopFun = "EBIC",
                    keep = NULL, maxK = NULL, verbose = FALSE) {

    fml <- Formula(formula)
    mf <- model.frame(fml, data=data)

    yvec <- model.part(fml, data=mf, lhs=1, drop=TRUE)
    xmat <- model.part(fml, data=mf, rhs=1) |> as.matrix()
    use.intercept <- attr(terms(mf), "intercept") == 1L


    mc <- match.call(expand.dots = TRUE)
    provided_args <- as.list(mc)[-1]
    provided_args$formula <- NULL
    provided_args$data <- NULL

    args <- list(
        fitFun = betareg,
        yvec = yvec,
        xmat = xmat,
        use.intercept = use.intercept
    )
    args <- c(args, provided_args)

    return(do.call(frs, args))
}