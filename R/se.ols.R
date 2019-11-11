#' Calculation of se for OLS
#'
#' This function calculates the se for UI based on OLS when we have MNAR data, for \code{\link{ui.ols}}.
#' @param X Covariate matrix.
#' @param sigmaOLScor Output from sigmaOLScor1
#' @param u Fitted values from mis.model.
#' @param gridrho Values of rho.
#' @export

se.ols<-function(X,sigmaOLScor,u,gridrho){
#M<-sigmaOLScor*sqrt(1-gridrho^2*(t(u)%*%lambda1(u)+t(lambda1(u))%*%lambda1(u))/n)
M<-sigmaOLScor*sqrt(1-gridrho^2*(mean(u*lambda1(u))-mean(lambda1(u)^2)))
Di<- sqrt(diag(solve(t(X)%*%X)))
matrix(M,ncol=1)%*%matrix(Di,nrow=1)
 }
 
