#' @rdname penalization
#' 
#' @description
#' - `pen_ncvreg(x, y, family, penalty)`: Non-Convex penalized models.
#' 
#' @export
pen_ncvreg <- function(x, y, family, penalty) {
    stopifnot( family %in% c("gaussian", "binomial", "cox") )
    stopifnot( penalty %in% c("MCP", "SCAD", "lasso") )

    if (family == "cox") {
        egg <- cv.ncvsurv(
            X = x,
            y = y,
            penalty = penalty,
            nfolds = 10
        )
    } else {
        egg <- cv.ncvreg(
            X = x,
            y = y,
            family = family,
            penalty = penalty,
            nfolds = 10
        )
    }

    class(egg) <- c("penaltyncvreg", class(egg))
    return(egg)
}



#' @title Extract Coefficients from Non-Convex Penalized Models
#'
#' @param object Object.
#' @param ... Additional arguments (not used).
#' 
#' @return Named vector of non-zero coefficients of under cross-validation error.
#' 
#' @exportS3Method
coef.penaltyncvreg <- function(object, ...) {
    idx <- object[["min"]]
    b0 <- object[["fit"]][["beta"]][, idx, drop=TRUE]
    b0 <- b0[grep("(Intercept)", names(b0), value=TRUE, fixed=TRUE, invert=TRUE)]
    return(b0[abs(b0) >= 1e-4])
}
