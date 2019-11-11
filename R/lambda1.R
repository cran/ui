#' Inverse Mills rato
#'
#' This function allows you to calculate the inverse Mills ratio.
#' @param x Vector
#' @importFrom stats pnorm dnorm
#' @export
lambda1<-function(x){
	return(dnorm(x)/pnorm(x))
	}