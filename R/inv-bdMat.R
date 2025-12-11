#' @title Generalized Inverse
#' @description Generalized inverse of a matrix with form `DiagBlock(dense, diag)`.
#' @param x Matrix
#' @param nonzeroIdx Index of nonzero entries.
#' @noRd
inv.bdMat <- function(x, nonzeroIdx) {
    idx <- sort(unique(c(nonzeroIdx)))

    x.inv <- Diagonal(x=1.0/diag(x))
    if (NROW(nonzeroIdx) > 0)
        x.inv[idx, idx] <- ginv(nearPD(x[idx, idx, drop=FALSE], base.matrix=TRUE)$mat)

    return(x.inv)
}
