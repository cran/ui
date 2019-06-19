#' Loglikelohood used in sandwich estimator of average causal effect on the treated for DR
#' 
#' Loglikelohood used in sandwich estimator of average causal effect on the treated for DR, support function for \code{\link{ui.causal}}
#' @param x coefficents.
#' @param X Covariate matrix.
#' @param z Missing or not.
#' @importFrom stats pnorm
#' @export
Logl.sandACT<-function(x,X,z){
return(mean(z*log(pnorm(X%*%x))+(1-z)*log(1-pnorm(X%*%x))))
}
