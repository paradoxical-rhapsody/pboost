#' @name pcoxph
#' @title Profile Boosting for Cox proportional hazards Model
#' 
#' @description Profile boosting for Cox model.
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
#' @param stopFun Parameter passed to [pboost::pboost].
#' @param keep Parameter passed to [pboost::pboost].
#' @param maxK Parameter passed to [pboost::pboost].
#' @param verbose Parameter passed to [pboost::pboost].
#' 
#' @return A `coxph` model object fitted on the selected features.
#' 
#' @examples
#' library(survival)
#' set.seed(2026)
#' n <- 300
#' p <- 200
#' 
#' DF <- data.frame(
#'     time = rpois(n, 5),
#'     status = rbinom(n, 1, 0.3),
#'     matrix(rnorm(n*p), n)
#' )
#' 
#' pcoxph(Surv(time, status) ~ ., DF, verbose=TRUE)
#' fcoxph(Surv(time, status) ~ ., DF, verbose=TRUE)
#' 
NULL
#> NULL



#' @rdname pcoxph
#' @order 1
#' @export
pcoxph <- function(formula, data, weights, subset, na.action, init, control,
                    ties = c("efron", "breslow", "exact"), singular.ok = TRUE,
                    robust, model = FALSE, x = FALSE, y = TRUE, tt, method = ties,
                    id, cluster, istate, statedata, nocenter = c(-1, 0, 1), ...,
                    stopFun = "EBIC",
                    keep = NULL, maxK = NULL, verbose = FALSE) {
    fml <- Formula(formula)
    mf <- model.frame(fml, data=data)

    yvec <- model.part(fml, data=mf, lhs=1, drop=TRUE)
    xmat <- model.part(fml, data=mf, rhs=1) |> as.matrix()
    use.intercept <- FALSE


    mc <- match.call(expand.dots = TRUE)
    provided_args <- as.list(mc)[-1]
    provided_args <- provided_args[!(names(provided_args) %in% c("formula", "data", "scoreFun"))]

    args <- list(
        fitFun = coxph,
        yvec = yvec,
        xmat = xmat,
        use.intercept = use.intercept,
        scoreFun = scoreCoxph
    )
    args <- c(args, provided_args)

    return(do.call(pboost, args))
}



#' @rdname pcoxph
#' @noRd
scoreCoxph <- function(object) {
    eta <- as.numeric(object$linear.predictors)
    method <- object$method
    stopifnot( method %in% c("efron", "breslow") )

    y <- object$y
    if (is.null(y)) {
        y <- model.response(model.frame(object))
    }
    if (inherits(y, 'Surv')) {
        time <- y[,1]
        status <- y[,ncol(y)]
    } else {
        stop('Cannot recognize y')
    }

    if (!(length(time) == length(status) && length(status) == length(eta)))
        stop('length mismatch')
    n <- length(time)
    o <- order(time, decreasing = TRUE)
    time_o <- time[o]
    status_o <- as.integer(status[o])
    eta_o <- eta[o]
    w <- exp(eta_o)
    rle_time <- rle(time_o)
    lengths <- rle_time$lengths
    ends <- cumsum(lengths)
    starts <- ends - lengths + 1

    grad_o <- numeric(n)

    if (method == 'breslow') {
        # Breslow: for each event block add d / S to all at-risk entries
        cw <- cumsum(w) # cumulative weight over decreasing times
        A_o <- numeric(n)
        for (k in seq_along(ends)) {
            s <- starts[k]; e <- ends[k]
            d <- sum(status_o[s:e] == 1)
            if (d == 0) next

            S <- cw[e]
            A_o[1:e] <- A_o[1:e] + (d / S)
        }
        grad_o <- status_o - w * A_o
    } else if (method == 'efron') {
        # Efron: more careful handling of ties with fractional removal
        for (k in seq_along(ends)) {
            s <- starts[k]; e <- ends[k]
            d <- sum(status_o[s:e] == 1)
            if (d == 0) next

            deaths <- which(status_o[s:e] == 1) + s - 1
            riskset <- 1:e
            S0 <- sum(w[riskset])
            sum_deaths <- sum(w[deaths])

            # for l = 0..(d-1), partially remove l/d of the deaths
            for (l in 0:(d-1)) {
                S0l <- S0 - (l / d) * sum_deaths

                # deaths contribute +1/d each in status part; then subtract w/S0l
                for (j in deaths)
                    grad_o[j] <- grad_o[j] + 1/d - w[j] / S0l

                # non-deaths in riskset only get - w / S0l contribution
                non_deaths <- setdiff(riskset, deaths)
                if (length(non_deaths))
                    grad_o[non_deaths] <- grad_o[non_deaths] - w[non_deaths] / S0l
            }
        }
    }

    grad <- numeric(n)
    grad[o] <- grad_o
    grad
}
