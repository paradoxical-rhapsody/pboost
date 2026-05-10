#' @rdname penalization
#' 
#' @description
#' - `pen_glmnet(x, y, family)`: LASSO penalized models.
#' 
#' @export
pen_glmnet <- function(x, y, family) {
    stopifnot( family %in% c("gaussian", "binomial", "cox") )

    suppressWarnings(
        egg <- cv.glmnet(
            x = x,
            y = y,
            family = family,
            standardize = FALSE,
            intercept = FALSE,
            type.measure = "deviance",
            nfolds = 10
        )
    )

    class(egg) <- c("penaltyglmnet", class(egg))
    return(egg)
}



#' @title Extract Coefficients from Lasso Penalized Models
#'
#' @param object Object.
#' @param ... Additional arguments (not used).
#' 
#' @return Named vector of non-zero coefficients of under cross-validation error.
#' 
#' @export
coef.penaltyglmnet <- function(object, ...) {
    idx <- object[["index"]]['min', 'Lambda']
    b0 <- object[["glmnet.fit"]][["beta"]][, idx, drop=TRUE]
    return(b0[abs(b0) >= 1e-4])
}
