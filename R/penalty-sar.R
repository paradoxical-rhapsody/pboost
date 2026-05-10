#' @rdname penalization
#' 
#' @description
#' - `pen_sar(x, y, rho, w, penalty)`: penalized spatial auto-regressive models.
#' 
#' @export
pen_sar <- function(x, y, rho, w, penalty) {
    stopifnot( penalty %in% c("MCP", "SCAD", "lasso") )

    n <- length(y)
    p <- NCOL(x)

    rho <- ifelse(missing(rho) || is.null(rho), get_rho(x, y, w), rho)
    stopifnot( length(rho) == 1 & is.numeric(rho) )

    A.rho <- diag(n) - rho * w
    y <- drop(A.rho %*% y)
    stopifnot( length(y) == n )

    if (penalty == "lasso")
        egg <- pen_glmnet(x, y, "gaussian")
    else
        egg <- pen_ncvreg(x, y, "gaussian", penalty)


    class(egg) <- c("pen.sar", class(egg))
    return(egg)
}
