#' Loglikelihood used by ui.probit 
#'
#' This function derives the Loglikelihood for \code{\link{ui.probit}}. 
#' @param par Coeficient values the logliklihood should be drived at.
#'@param rho The value of the sensitivity parameter.
#'@param X.z covariate matrix for missingness mechanism
#'@param X.y covariate matrix for the outcome regression
#'@param z indicator of wether y is missing or not
#'@param y outcome vector
#'@import mvtnorm
#'@importFrom stats pnorm dnorm
#' @export
LogL.probit<-function(par,rho,X.z = X.z, X.y = X.y, y = y, z = z){
d<-dim(X.y)[2]
beta<-par[1:d]
delta<-par[(d+1):length(par)]

q <- 2*y - 1
w1<-(q*tcrossprod(beta,X.y))[z==1]
dx<-tcrossprod(delta,X.z)
dx1<-dx[z==1]
Rhos<-(q[z==1]*rho)

n1<-sum(z==1)
Phi2<-vector(length=n1)
for(i in 1:n1){
	Phi2[i]<-pmvnorm(lower=-Inf,upper=c(w1[i],dx1[i]), 
				mean=c(0,0),corr=rbind(c(1,Rhos[i]),c(Rhos[i],1)))
}
	
logl<-sum((1-z)*log(1-pnorm(dx))) + sum(log(Phi2))
return(logl)
}

