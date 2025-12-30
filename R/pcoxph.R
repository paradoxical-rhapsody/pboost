#' @name pcoxph
#' @title Profile Boosting for Cox proportional hazards Model
#' 
#' @description Profile boosting for Cox model.
#' 
#' @param formula See [pboost].
#' @param data See [pboost].
#' @param weights Parameters passed to [survival::coxph].
#' @param subset Parameters passed to [survival::coxph].
#' @param na.action Parameters passed to [survival::coxph].
#' @param init Parameters passed to [survival::coxph].
#' @param control Parameters passed to [survival::coxph].
#' @param ties Parameters passed to [survival::coxph].
#' @param singular.ok Parameters passed to [survival::coxph].
#' @param robust Parameters passed to [survival::coxph].
#' @param model Parameters passed to [survival::coxph].
#' @param x Parameters passed to [survival::coxph].
#' @param y Parameters passed to [survival::coxph].
#' @param tt Parameters passed to [survival::coxph].
#' @param method Parameters passed to [survival::coxph].
#' @param id Parameters passed to [survival::coxph].
#' @param cluster Parameters passed to [survival::coxph].
#' @param istate Parameters passed to [survival::coxph].
#' @param statedata Parameters passed to [survival::coxph].
#' @param nocenter Parameters passed to [survival::coxph].
#' @param ... Parameters passed to [survival::coxph].
#' @param stopFun Parameters passed to [pboost].
#' @param keep Parameters passed to [pboost].
#' @param maxK Parameters passed to [pboost].
#' @param verbose Parameters passed to [pboost].
#' 
#' 
#' @return An `coxph` model object fitted on the selected features.
#' 
#' @examples
#' library(survival)
#' set.seed(2025)
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
#' 
NULL
#> NULL



#' @rdname pcoxph
#' @order 1
#' @export
pcoxph <- function(
    formula, data, weights, subset, na.action, init, control,
    ties = c("efron", "breslow", "exact"), singular.ok = TRUE,
    robust, model = FALSE, x = FALSE, y = TRUE, tt, method = ties,
    id, cluster, istate, statedata, nocenter = c(-1, 0, 1), ...,
    stopFun = EBIC, keep = NULL, maxK = NULL, verbose = FALSE) {
    stopifnot( !missing(formula) )
    stopifnot( !missing(data) )

    cl <- match.call()

    coxph_template <- cl
    coxph_template$stopFun <- NULL
    coxph_template$keep <- NULL
    coxph_template$maxK <- NULL
    coxph_template$verbose <- NULL
    coxph_template[[1L]] <- quote(coxph)

    required_paras <- c("data", "weights", "subset", "na.action")
    for (ipara in required_paras)
        if (!is.null(cl[[ipara]]))
            coxph_template[[ipara]] <- eval(cl[[ipara]], envir = parent.frame())

    fitFun <- function(formula, data) {
        call <- coxph_template
        call$formula <- formula
        call$data <- data
        return( eval(call, parent.frame()) )
    }

    return(pboost(formula, data, fitFun, scoreCoxph, stopFun,
                  keep = keep, maxK = maxK, verbose = verbose))
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


#' @rdname EBIC
#' @export
EBIC.coxph <- function(object, p, p.keep, ...) {
    stopifnot( inherits(object, "coxph") )

    if (missing(p))
        p <- get("p", envir=parent.frame())
    if (missing(p.keep))
        p.keep <- get("p.keep", envir=parent.frame())

    dof <- attr(logLik(object), "df")
    ebic.r <- max( 0.0, 1.0 - log(nobs(object)) / (2.0*log(p)) )
    suppressWarnings(
        ebic.penalty <- ifelse(
            ebic.r <= 0.0,
            0.0,
            2.0 * ebic.r * lchoose(p - p.keep, dof - p.keep)
        )
    )

    stopifnot( is.finite(ebic.penalty) )
    return(BIC(object) + ebic.penalty)
}
