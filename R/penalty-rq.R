#' @rdname penalization
#' 
#' @description
#' - `pen_rq(x, y, tau, penalty)`: penalized quantile regression models.
#' 
#' @export
pen_rq <- function(x, y, tau, penalty) {
    stopifnot( penalty %in% c("MCP", "SCAD", "lasso") )
    penalty <- toupper(penalty)

    # `rq.pen.cv`, `rq.pen`
    egg <- rq.pen(
        x = x,
        y = y,
        tau = tau,
        penalty = penalty,
        scalex = FALSE
    )

    class(egg) <- c("penaltyrq", class(egg))
    return(egg)
}



#' @title Extract Coefficients from Non-Convex Penalized Models
#'
#' @param object Object.
#' @param ... Additional arguments (not used).
#' 
#' @return A named vector of coefficients.
#' 
#' @exportS3Method
coef.penaltyrq <- function(object, ...) {
    class(object) <- setdiff(class(object), "penaltyrq")
    b0 <- qic.select(object, method="PBIC") |> coef() |> drop()
    b0 <- b0[grep("intercept", names(b0), value=TRUE, fixed=TRUE, invert=TRUE)]
    return(b0[abs(b0) >= 1e-4])
}
