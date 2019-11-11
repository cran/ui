#' Gradient for the loglikelihood used by ui.probit 
#'
#' This function derives the gradient in order for \code{\link{ui.probit}} to run faster.
#' @param par Coefficients.
#' @param rho Rho.
#' @param X.z Covariate matrix for missingness.
#' @param X.y Covariate matrix for outcome.
#' @param z Missing or not.
#' @param y Outcome.
#' @import mvtnorm
#' @importFrom stats pnorm dnorm
#' @export
grr<- function(par,rho,X.z = X.z, X.y = X.y, y = y, z = z){
d<-dim(X.y)[2]
beta<-par[1:d]
delta<-par[(d+1):length(par)]

q <- 2*y - 1
bx<-tcrossprod(beta,X.y)
bx1<-bx[z==1]
w1<-q[z==1]*bx1
dx<-tcrossprod(delta,X.z)
q1<-q[z==1]
Rhos1<-q1*rho
dx1<-dx[z==1]

n1<-sum(z==1)
Phi2<-vector(length=n1)
for(i in 1:n1){
		Phi2[i]<-pmvnorm(lower=-Inf,upper=c(w1[i],dx1[i]), 
						mean=c(0,0),corr=rbind(c(1,Rhos1[i]),c(Rhos1[i],1)))}

gr.b<-q1*dnorm(w1)*pnorm((dx1-rho*bx1)/sqrt(1-rho^2))/Phi2

gr.beta<-crossprod(gr.b,X.y[z==1,])

gr.d1 <- dnorm(dx1)*pnorm(q1*(bx1-rho*dx1)/sqrt(1-rho^2))/Phi2 
gr.d0 <- -dnorm(dx[z==0])/(1-pnorm(dx[z==0]))

gr.delta<-crossprod(gr.d1,X.z[z==1,])+crossprod(gr.d0,X.z[z==0,])

return(c(gr.beta,gr.delta))
}





