#' @rdname penalization
#' 
#' @description
#' - `pen_sar(x, y, rho, w, penalty)`: penalized spatial auto-regressive models.
#' 
#' @noRd
pen_sar <- function(x, y, rho, w, penalty) {
    stopifnot( penalty %in% c("MCP", "SCAD", "lasso") )

    rho <- ifelse(missing(rho) || is.null(rho), get_rho(x, y, w), rho)
    stopifnot( length(rho) == 1 & is.numeric(rho) )

    n <- length(y)
    A.rho <- diag(n) - rho * w
    y <- drop(A.rho %*% y)
    stopifnot( length(y) == n )

    if (penalty == "lasso")
        egg <- pen_glmnet(x, y, "gaussian")
    else
        egg <- pen_ncvreg(x, y, "gaussian", penalty)

    class(egg) <- c("penltysar", class(egg))
    return(egg)
}
