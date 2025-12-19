#' @name pggm
#' @title Profile Boosting for Gaussian Graphical Model
#' 
#' @description Profile boosting for Gaussian graphical model.
#' 
#' @param S Covariance matrix.
#' @param nObs Number of observations.
#' @param maxK Maximum number of identified edges.
#' @param digits Integer indicating the number of decimal
#'  places or significant digits to be used.
#' @param verbose Print the procedure path?
#' 
#' @return Index set of identified features.
#' 
#' @examples
#' \donttest{
#' library(MASS)
#' library(Matrix)
#' 
#' set.seed(2025)
#' n <- 1000
#' p <- 10
#' 
#' Omega <- Diagonal(p)
#' diag(Omega[1:4, 2:5]) <- diag(Omega[2:5, 1:4]) <- 0.5
#' Sigma <- chol2inv(chol(Omega))
#' X <- mvrnorm(n, rep(0, p), Sigma, empirical=TRUE)
#' S <- cov(X)
#' system.time( egg <- pggm(S, n) )
#' }
#' 
NULL
#> NULL



#' @rdname pggm
#' @export
pggm <- function(S, nObs, maxK = floor(min(nObs-1, NROW(S)-1, 50)),
                digits = 8, verbose = FALSE) {
    stopifnot( isSymmetric(S) )
    stopifnot( nObs > 0 )
    stopifnot( maxK > 0 )

    S <- Matrix(S)

    p0 <- NROW(S)
    rEBIC <- max( 0.0, 1.0 - log(nObs)/(2.0*log(p0)))

    grad <- function(Omega, nonzeroIdx) {
        egg <- S - inv.bdMat(Omega, nonzeroIdx)
        return(2*egg - Diagonal(x=diag(egg)))
    }
    ebic <- function(Omega, dof) {
        bic <- c(sum(S*Omega) - determinant(Omega)$modulus + dof*log(nObs)/nObs)
        return( bic + 2*rEBIC*lchoose(p0*(p0-1), dof)/nObs )
    }

    showiter <- function(verbose, k, idx, level)
        if (verbose) 
            message(sprintf("Step %i: add (%i, %i) with level=%.3f)",
                            k, idx[1], idx[2], level))

    nonzeroIdx <- matrix(NA, 0, 2)
    Omega <- Diagonal(x=1.0/diag(S))
    ebicVal <- ebic(Omega, NROW(nonzeroIdx))
    ebicVec <- c(ebicVal, rep(NA_real_, maxK))
    k <- 1
    while (k <= maxK) {
        pps <- abs(grad(Omega, nonzeroIdx))
        idx.tmp <- sort(which(pps == max(pps), TRUE)[1, ])
        nonzeroIdx.tmp <- rbind(nonzeroIdx, idx.tmp)

        rmle.tmp <- rmle4ggmS4(S, nonzeroIdx.tmp)
        ebicVal.tmp <- ebic(rmle.tmp$Omega, NROW(nonzeroIdx.tmp))

        # message(sprintf("%f", ebicVal.tmp))

        if (ebicVal.tmp >= ebicVal) break

        nonzeroIdx <- nonzeroIdx.tmp
        Omega <- rmle.tmp$Omega
        ebicVal <- ebicVal.tmp
        
        k <- k + 1
        ebicVec[k+1] <- ebicVal
    }

    dimnames(nonzeroIdx) <- NULL
    return(list(
        Omega=round(Omega, digits),
        nonzeroIdx=nonzeroIdx,
        ebic=c(na.omit(ebicVec))
    ))
}






#' @rdname pggm
#' @description `rmle4ggmS4.dense`: Restricted MLE for GGM.
#' @noRd
rmle4ggmS4.dense <- function(S, nonzeroIdx, tol=1e-4, maxIter=100, verbose=FALSE) {
    stopifnot( isSymmetric(S) )
    stopifnot( all(nonzeroIdx > 0) )
    stopifnot( all(nonzeroIdx <= NROW(S)) )
    stopifnot( all(c(tol, maxIter) > 0) )

    active.mat <- Matrix(0.0, NROW(S), NCOL(S))
    active.mat[nonzeroIdx] <- 1.0
    active.mat <- active.mat + t(active.mat)
    diag(active.mat) <- 1.0

    loss <- function(Omega) c(sum(S * Omega) - determinant(Omega)$modulus)
    grad <- function(Omega, nonzeroIdx) {
        egg <- active.mat * (S - inv.bdMat(Omega, nonzeroIdx))
        return(2*egg - Diagonal(x=diag(egg)))
    }

    showiter <- function(verbose)
        if (verbose) message(sprintf("iter %i: obj=%.4f", k, obj))

    Omega0 <- Diagonal(x=1.0/diag(S))
    Omega <- W <- Omega0
    obj <- loss(Omega0)

    alpha0 <- 1.0
    alpha <- (1.0 + sqrt(1.0 + 4*alpha0^2)) / 2.0
    iter <- rep(NA_real_, maxIter)
    isConvergent <- FALSE
    delta <- 1.0
    for (k in seq_len(maxIter)) {
        showiter(verbose)
        iter[k] <- obj

        ## search point
        W <- Omega + (alpha0 - 1.0)/alpha * (Omega - Omega0)

        ## approximate solution
        Omega0 <- Omega
        obj0 <- obj
        W.go <- grad(W, nonzeroIdx) / delta
        if (max(abs(W.go)) < min(1.0e-6, 0.1*tol)) {
            isConvergent <- TRUE
            break
        }

        Omega <- symmpart(nearPD(W - W.go)$mat) # `keepDiag`
        obj <- loss(Omega)
        stopifnot( isSymmetric(Omega) )

        # here maybe arise the case `obj - obj0 \approx tol > 0`
        if (obj > obj0) {
            Omega <- Omega0
            obj <- obj0
            delta <- 1.5 * delta
            next
        }

        if ( abs(obj-obj0) <= (abs(obj0) + 0.1)*tol ) {
            isConvergent <- TRUE
            break
        }

        ## update paras
        alpha0 <- alpha
        alpha <- (1.0 + sqrt(1.0 + 4*alpha0^2)) / 2.0
    }

    list(Omega=round(Omega, 8), iter=c(na.omit(iter)), isConvertent=isConvergent)
}



#' @rdname pggm
#' @description `rmle4ggmS4`: Restricted MLE for GGM
#' @noRd
rmle4ggmS4  <- function(S, nonzeroIdx, tol=1e-4, maxIter=100) {
    stopifnot( isSymmetric(S) )
    stopifnot( all(nonzeroIdx > 0) )
    stopifnot( all(nonzeroIdx <= NROW(S)) )
    stopifnot( all(c(tol, maxIter) > 0) )

    idx <- sort(unique(c(nonzeroIdx)))

    indicator.mat <- Matrix(FALSE, NROW(S), NCOL(S), doDiag=FALSE)
    indicator.mat[nonzeroIdx] <- TRUE
    nonzeroIdx.sub <- which(indicator.mat[idx, idx, drop=FALSE], TRUE)
    obj <- rmle4ggmS4.dense(S[idx, idx, drop=FALSE], nonzeroIdx.sub, tol, maxIter)

    Omega <- Diagonal(x=1.0/diag(S))
    Omega[idx, idx] <- obj$Omega

    list(Omega=Omega, iter=obj$iter, isConvergent=obj$isConvergent)
}
