#' @name pbetareg
#' @title Profile Boosting for Beta Regression
#' 
#' @description
#' [pbetareg] inherits the usage of [betareg::betareg].
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
#' @param stopFun Parameter passed to [pboost::pboost].
#' @param keep Parameter passed to [pboost::pboost].
#' @param maxK Parameter passed to [pboost::pboost].
#' @param verbose Parameter passed to [pboost::pboost].
#' 
#' @return A `betareg` model object fitted on the selected features.
#' 
#' @examples
#' \dontrun{
#' set.seed(2026)
#' n <- 300
#' p <- 20
#' x <- matrix(runif(n*p), n)
#' mu <- runif(n)
#' phi <- 1.0
#' 
#' shape1 <- mu * phi
#' shape2 <- (1-mu) * phi
#' y <- rbeta(n, shape1, shape2)
#' DF <- data.frame(y, x)
#' 
#' pbetareg(y ~ ., DF, verbose=TRUE)
#' fbetareg(y ~ ., DF, verbose=TRUE)
#' }
#' 
NULL
#> NULL


#' @rdname pbetareg
#' @order 1
#' @export
pbetareg <- function(formula, data, subset, na.action, weights, offset,
                    link = c("logit", "probit", "cloglog", "cauchit", "log","loglog"),
                    link.phi = NULL, type = c("ML", "BC", "BR"),
                    dist = NULL, nu = NULL, control = betareg.control(...),
                    model = TRUE, y = TRUE, x = FALSE, ...,
                    stopFun = "EBIC",
                    keep = NULL, maxK = NULL, verbose = FALSE) {

    fml <- Formula(formula)
    mf <- model.frame(fml, data=data)

    yvec <- model.part(fml, data=mf, lhs=1, drop=TRUE)
    xmat <- model.part(fml, data=mf, rhs=1) |> as.matrix()
    use.intercept <- attr(terms(mf), "intercept") == 1L


    scoreFun <- function(object) {
        phi <- predict(object, type='precision')
        mu <- predict(object, type='response')
        eta <- predict(object, type='link')

        ## object$link is a list of two elements: one for "mean" or "mu"
        mu.eta <- object$link[[grep("^m", names(object$link), value=TRUE)]]$mu.eta
        y <- pmin(pmax(object[["y"]], .Machine$double.eps), 1 - .Machine$double.eps)

        # weights <- object[["weights"]]
        return( mu.eta(eta) * phi * ( digamma((1-mu)*phi) - digamma(mu*phi) + qlogis(y) ) )
    }

    mc <- match.call(expand.dots = TRUE)
    provided_args <- as.list(mc)[-1]
    provided_args <- provided_args[!(names(provided_args) %in% c("scoreFun"))]
    provided_args$formula <- NULL
    provided_args$data <- NULL

    args <- list(
        fitFun = betareg,
        yvec = yvec,
        xmat = xmat,
        use.intercept = use.intercept,
        scoreFun = scoreFun
    )
    args <- c(args, provided_args)

    return(do.call(pboost, args))
}
