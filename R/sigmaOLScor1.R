#' Correction of OLS sigma
#'
#' This function is a bias correction of the residual standard deviation under MNAR, used by \code{\link{ui.causal}} and \code{\link{ui.ols}}.
#'
#' @param X Covariate matrix outcome.
#' @param sigmaOLS Residual sd from outcome regression.
#' @param n Number of complete cases.
#' @param p Number of covariates outcome regression.
#' @param u Fitted values from propensity score regression.
#' @param gridrho Values of rho.
#' @export

sigmaOLScor1<-function(X,sigmaOLS,n,p,u,gridrho){
K<-c((-t(u)%*%lambda1(u)-t(lambda1(u))%*%X%*%solve(t(X)%*%X)%*%t(X)%*%lambda1(u))/(n-p))
sqrt(sigmaOLS^2/(1+gridrho^2*K))
 }