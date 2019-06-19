#' Uncertainty intervals for probit regression
#'
#' This function allows you to derive uncertainty intervals for probit regression
#'	when there is missing data in the binary outcome. The uncertainty intervals
#'	can be used as a sensitivity analysis to ignorability (missing at random), and
#'	are derived by maximum likelihood. Note that \code{rho}=0 render the same results as
#'	a complete case analysis.
#' @param out.formula Formula for outcome regression.
#' @param mis.formula Formula for missingness mechanism. If NULL the same covariates as in the outcome regression will be used. 
#' @param data data.frame containing the variables in the formula. 
#' @param rho Vector containing the values of \code{rho} for which we want to fit the
#'	likelihood.
#' @param progress If TRUE prints out process time for each maximization of the likelihood.
#' @param max.grid Maximum distance between two elements in \code{rho}, if two wide there can 
#'	difficulties with convergence of the maximum likelihood.
#' @param alpha Default 0.05 corresponding to a confidence level of 95 for CI and UI.
#' @param method Maximization method to be passed through \code{maxLik}
#'
#'@details In order to visualize the results, you can use \code{\link{plot.uiprobit}} 
#'or \code{\link{profile.uiprobit}}.
#'
#' @importFrom  Matrix rankMatrix
#' @return A list containing:
#' \item{coef}{Estimated coefficients (outcome regression) for different values of \code{rho}.}
#' \item{rho}{The values of \code{rho} for which the likelihood is maximized.}
#' \item{vcov}{Covariance matrix.}
#' \item{ci}{Confidence intervals for different values of \code{rho}.}
#' \item{ui}{Uncertainty intervals.}
#' \item{out.model}{Outcome regression model when rho=0.}
#' \item{mis.model}{Regression model for missingness mechanism (selection).}
#' \item{se}{Standard errors from outcome regression.}
#' \item{value}{Value of maximum likelihood for different values of \code{rho}.}
#' \item{y}{Outcome vector.}
#' \item{z}{Indicator variable of observed outcome.}
#' \item{X.y}{Covariate matrix for outcome regression.}
#' \item{X.z}{Covariate matrix for missingness mechanism (selection regression model).}
#' \item{max.info}{Information about the maximization procedure. Includes whether it \code{converged}, \code{message}, \code{method} and \code{number of iterations}.}
#'@author Minna Genbäck
#'@references Genbäck, M., Ng, N., Stanghellini, E., de Luna, X. (2018). Predictors of Decline in Self-reported Health: Addressing Non-ignorable Dropout in Longitudinal Studies of Aging. \emph{European journal of ageing}, 15(2), 211-220.
#' @examples
#'library(MASS)
#'n<-500
#'
#'delta<-c(0.5,0.6,0.1,-1,1)
#'beta<-c(-0.3,-0.5,0,-0.4,-0.3)
#'
#'X<-cbind(rep(1,n),rnorm(n),runif(n),rbinom(n,2,0.5),rbinom(n,1,0.5))
#'x<-X[,-1]
#'rho=0.4
#'error<-mvrnorm(n,c(0,0),matrix(c(1,rho,rho,1),2))
#'
#'zstar<-X%*%delta+error[,1]
#'z<-as.numeric(zstar>0)
#'
#'ystar<-X%*%beta+error[,2]
#'y<-as.integer(ystar>0)
#'y[z==0]<-NA
#'data=data.frame(y=y,x1=x[,1],x2=x[,2],x3=x[,3],x4=x[,4])
#'
#'m<-ui.probit(y~x1+x2+x3+x4,data=data,rho=c(0,0.5))
#'m
#'plot(m)
#'profile(m)
#' @importFrom stats qnorm
#' @export

ui.probit<-function(out.formula,mis.formula=NULL,data,rho=c(-0.3,0.3),progress=TRUE,max.grid=0.1,alpha=0.05,method='NR'){

  cll <- match.call()
  
  # Check to make sure 0 is part of Rho
  if(min(rho)>0 | max(rho)<0){
  rho<-c(0,rho)
  warning('Note that, the rho interval did not cover 0 as recommended. Hence 0 is now added into rho.')	
  }
  if(sum(rho==0)==0){
  	rho<-c(rho,0)
  }
  rho<-sort(rho) # Sort Rho from smallest to largest value
  nrho<-length(rho)


  #Make sure that the grid of rho is as fine as specified by max.grid 
  for(i in 1:(nrho-1)){
  nri<-ceiling((rho[i+1]-rho[i])/max.grid)
  if(nri>1){
  rho<-c(rho,seq(from=rho[i],to=rho[i+1],by=((rho[i+1]-rho[i])/nri))[2:nri])
  }}

  rho<-sort(rho)
  nrho<-length(rho)
	
  ML.object<-ML.probit(out.formula,mis.formula,data,rho,progress,method)
  
  p<-dim(ML.object$coef)[1]-1
  #derive standard errors from vcov matrix
  se<-lapply(ML.object$vcov,function(x) sqrt(diag(x)))
  ML.object$se<-matrix(unlist(se),ncol=nrho)[1:(p+1),]
  
  zalpha<-qnorm(1-alpha/2)
  ML.object$zalpha<-zalpha

  ML.object$ci<-array(dim=c((p+1),nrho,2),dimnames=list(rownames(ML.object$coef),rho,c('lower','upper')))
  ML.object$ci[,,1]<-ML.object$coef-zalpha*ML.object$se
  ML.object$ci[,,2]<-ML.object$coef+zalpha*ML.object$se
  ML.object$ui<-cbind(apply(ML.object$ci[,,1],1,min),apply(ML.object$ci[,,2],1,max))
  colnames(ML.object$ui)<-c('lower','upper')

  ML.object$call <- cll
  
  class(ML.object)<-"uiprobit"
  return(ML.object)
}