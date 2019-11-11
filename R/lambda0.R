#' Inverse Mills rato
#'
#' This function allows you to calculate the inverse Mills ratio.
#' @param x Vector
#' @importFrom stats pnorm dnorm
#' @export
lambda0<-function(x){
	return(dnorm(x)/(1-pnorm(x)))
	}