#' Correction of OLS sigma for causal effects
#'
#' This function is a bias correction of the residual standard deviation under MNAR, for \code{\link{ui.causal}}.
#'
#' @param X Covariate matrix outcome.
#' @param sigmaOLS Residual sd from outcome regression.
#' @param n Number of complete cases.
#' @param p Number of covariates outcome regression.
#' @param u Fitted values from propensity score regression.
#' @param gridrho Values of rho.
#' @export

sigmaOLScor0<-function(X,sigmaOLS,n,p,u,gridrho){
K<-c((t(u)%*%lambda0(u)-t(lambda0(u))%*%X%*%solve(t(X)%*%X)%*%t(X)%*%lambda0(u))/(n-p))
sqrt(sigmaOLS^2/(1+gridrho^2*K))
}