#' @name prostate
#' @title  Prostate Tumor Gene Expression Dataset
#' @docType data
#' 
#' @description 
#' This data arises from a large study to examine EEG
#' correlates of genetic predisposition to alcoholism.
#' See <http://kdd.ics.uci.edu/databases/eeg/> for details.
#' 
#' @format A list:
#' \describe{
#'  \item{x}{Gene expression data (matrix with 102 rows and 6033 columns).}
 #' \item{y}{Class index (vector with 102 elements.}
#' }
#' 
#' @details
#' The prostate dataset consists of 52 prostate tumor and 50 normal samples.
#' Normal and tumor classes are coded in 0 and 1, respectively, in y vector.
#' Matrix x is gene expression data and arrays were normalized, log transformed,
#' and standardized to zero mean and unit variance across genes as described in
#' Dettling (2004) and Dettling and Beuhlmann (2002).
#' See Chung and Keles (2010) for more details.
#' 
#' 
#' @source
#' Singh D, Febbo P, Ross K, Jackson D, Manola J, Ladd C, Tamayo P, Renshaw A,
#' DAmico A, Richie J, Lander E, Loda M, Kantoff P, Golub T, and Sellers W (2002),
#' "Gene expression correlates of clinical prostate cancer behavior",
#' Cancer Cell, Vol. 1, pp. 203–209.
#' 
#' 
#' @examples
#' data(prostate)
#' head(prostate$x[, 1:3])
#' print(prostate$y)
#' 
#' @keywords datasets
NULL
