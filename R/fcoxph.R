#' @name fcoxph
#' @title Forward Regression Selection for Cox proportional hazards Model
#' 
#' @description Forward regression selection for Cox model.
#' 
#' @param formula Parameter passed to [survival::coxph].
#' @param data Parameter passed to [survival::coxph].
#' @param weights Parameter passed to [survival::coxph].
#' @param subset Parameter passed to [survival::coxph].
#' @param na.action Parameter passed to [survival::coxph].
#' @param init Parameter passed to [survival::coxph].
#' @param control Parameter passed to [survival::coxph].
#' @param ties Parameter passed to [survival::coxph].
#' @param singular.ok Parameter passed to [survival::coxph].
#' @param robust Parameter passed to [survival::coxph].
#' @param model Parameter passed to [survival::coxph].
#' @param x Parameter passed to [survival::coxph].
#' @param y Parameter passed to [survival::coxph].
#' @param tt Parameter passed to [survival::coxph].
#' @param method Parameter passed to [survival::coxph].
#' @param id Parameter passed to [survival::coxph].
#' @param cluster Parameter passed to [survival::coxph].
#' @param istate Parameter passed to [survival::coxph].
#' @param statedata Parameter passed to [survival::coxph].
#' @param nocenter Parameter passed to [survival::coxph].
#' @param ... Parameters passed to [survival::coxph].
#' 
#' @param selectFun Parameter passed to [pboost::frs].
#' @param stopFun Parameter passed to [pboost::frs].
#' @param keep Parameter passed to [pboost::frs].
#' @param maxK Parameter passed to [pboost::frs].
#' @param verbose Parameter passed to [pboost::frs].
#' 
#' @return A `coxph` model object fitted on the selected features.
#' 
NULL
#> NULL



#' @rdname fcoxph
#' @order 1
#' @export
fcoxph <- function(formula, data, weights, subset, na.action, init, control,
                    ties = c("efron", "breslow", "exact"), singular.ok = TRUE,
                    robust, model = FALSE, x = FALSE, y = TRUE, tt, method = ties,
                    id, cluster, istate, statedata, nocenter = c(-1, 0, 1), ...,
                    selectFun = logLik, stopFun = "EBIC",
                    keep = NULL, maxK = NULL, verbose = FALSE) {
    fml <- Formula(formula)
    mf <- model.frame(fml, data=data)

    yvec <- model.part(fml, data=mf, lhs=1, drop=TRUE)
    xmat <- model.part(fml, data=mf, rhs=1) |> as.matrix()
    use.intercept <- FALSE


    mc <- match.call(expand.dots = TRUE)
    provided_args <- as.list(mc)[-1]
    provided_args$formula <- NULL
    provided_args$data <- NULL

    args <- list(
        fitFun = coxph,
        yvec = yvec,
        xmat = xmat,
        use.intercept = use.intercept
    )
    args <- c(args, provided_args)

    return(do.call(frs, args))
}
