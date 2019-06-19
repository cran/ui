#' Hessian for the loglikelihood used by ui.probit 
#'
#' This function derives the hessian in order for \code{\link{ui.probit}} to run faster.
#' @param par Coefficients.
#' @param rho Rho.
#' @param X.z Covariate matrix for missingness.
#' @param X.y Covariate matrix for outcome.
#' @param z Missing or not.
#' @param y Outcome.
#' @import mvtnorm
#' @importFrom stats pnorm dnorm
#' @export

hess<- function(par,rho,X.z = X.z, X.y = X.y, y = y, z = z){
d<-dim(X.y)[2]
beta<-par[1:d]
delta<-par[(d+1):length(par)]

q <- 2*y - 1
w1<-(q*tcrossprod(beta,X.y))[z==1]
dx<-tcrossprod(delta,X.z)

q1<-q[z==1]
dx1<-dx[z==1]
Rhos1<-q1*rho
dx0<-dx[z==0]

n1<-sum(z==1)
Phi2<-vector(length=n1)
for(i in 1:n1){
		Phi2[i]<-pmvnorm(lower=-Inf,upper=c(w1[i],dx1[i]), 
						mean=c(0,0),corr=rbind(c(1,Rhos1[i]),c(Rhos1[i],1)))}
						
pnorm1 <- pnorm((dx1 - Rhos1*w1)/sqrt(1 - rho^2))
dnorm1 <- dnorm((dx1 - Rhos1*w1)/sqrt(1 - rho^2))
pnorm2 <- pnorm((w1 - Rhos1*dx1)/sqrt(1 - rho^2))
dnorm2 <- dnorm((w1 - Rhos1*dx1)/sqrt(1 - rho^2))
dnorm.w <- dnorm(w1)
dnorm.dx <- dnorm(dx1)
lambda.dx0<- dnorm(dx0)/(1-pnorm(dx0))


secder.b<- (-dnorm.w/Phi2)*(w1*pnorm1+Rhos1*dnorm1/sqrt(1-rho^2)+ dnorm.w*pnorm1^2/Phi2)
hess.b <- crossprod(c(secder.b)*X.y[z==1,], X.y[z==1,]) 

secder.cross<-(q1*dnorm.w/Phi2)*(dnorm1/sqrt(1-rho^2)-dnorm.dx*pnorm1*pnorm2/Phi2)
hess.cross<-crossprod(c(secder.cross)*X.z[z==1,], X.y[z==1,])

secder.d0<- dx0*lambda.dx0-lambda.dx0^2
secder.d1<- (-dnorm.dx/Phi2)*(dx1*pnorm2+(Rhos1*dnorm2/sqrt(1-rho^2))+dnorm.dx*pnorm2^2/Phi2)
hess.d<-crossprod(c(secder.d0)*X.z[z==0,], X.z[z==0,])+crossprod(c(secder.d1)*X.z[z==1,], X.z[z==1,])

hess<-rbind(cbind(hess.b,t(hess.cross)),cbind(hess.cross,hess.d))
return(hess)
}