#' @name EBIC
#' @title Extended Bayesian Information Criterion
#' @description The Extended BIC possesses the selection consistency in 
#' high-dimensional model.
#'
#' It can be called by the fitted model that has standard [logLik] method
#' to access the attributes `nobs` and `df`, such as [lm], [glm].
#' 
#' @param object Fitted model object.
#' @param p Total number of candidate features, which is available in [pboost].
#' @param p.keep Number of features that are pre-specified to be kept in model.
#' @param ... Additional parameters, which is available in [pboost].
#' 
#' @return A function to obtain the EBIC value of a fitted object.
#' 
#' @details 
#' The built-in [BIC] has the definition
#' ```
#' BIC(obj) == -2*logLik(obj) + attr(logLik(obj), "df") * log(nobs(logLik(obj))).
#' ```
#' 
#' The extended BIC (EBIC) is defined as
#' ```
#' EBIC(obj) == BIC(obj) + 2 * r * log(choose(p - |keep|, df - |keep|)).
#' ```
#' 
#' @references 
#' * Jiahua Chen and Zehua Chen (2008). Extended Bayesian information criteria for model selection with large model spaces. Biometrika, 95(3):759–771. \doi{10.1093/biomet/asn034}
#' 
#' * Jiahua Chen and Zehua Chen (2012). Extended BIC for small-n-large-p sparse GLM. Statistical Sinica, 22(2):555–574. \doi{10.5705/ss.2010.216}
#' 
#' @seealso [plm], [pglm], [pcoxph], [prq], [pbetareg].
#' 
NULL
#> NULL


#' @rdname EBIC
#' @order 1
#' @export
EBIC <- function(object, p, p.keep, ...){
    UseMethod("EBIC")
}



#' @rdname EBIC
#' @export
#' @noRd
EBIC.default <- function(object, p, p.keep=NULL, ...){
    stop("No EBIC method for object of class ", paste(class(object), collapse=", "))
}
