#' Calculates standard error of Average causal effect on the treated
#' 
#' This is a support function for \code{\link{ui.causal}} and calculates standard error of Average causal effect on the treated for the doubly robust estimator.
#' @param deltasigma1 Coefficients.
#' @param X Covariate matrix outcome.
#' @param Xz Covariate matrix treatment.
#' @param y Outcome vector.
#' @param z Missingness indicator.
#' @param u Fitted values from propensity score regression.
#' @param BetaOLSy0 Coefficients from non-treated regression 
#' @param phat Fitted propensity scores.
#' @param NaivEst Naiv estimates.
#' @param n1 Number of treated.
#' @param n0 Number of non-treated.
#' @param N Total number.
#' @param p Number of covariates outcome regression.
#' @param pz Number of covariates treatment regression.
#' @importFrom numDeriv hessian
#' @export


sandACT<-function(deltasigma1,X,Xz,y,z,u,BetaOLSy0,phat,NaivEst,n1,n0,N,p,pz){

a12<-apply(X, 2,function(x){mean((z-(1-z)/(1-phat))*x)})

yhat<-X%*%BetaOLSy0

a13<-apply(Xz, 2,function(x){mean((1-z)*((y-yhat)/(1-phat))*lambda0(u)*x)})

D<-matrix(nrow=(p+pz+2),ncol=(p+pz+2),data=0)

D[1:(p+1),1:(p+1)]<-t((1-z)*X)%*%((1-z)*X)/N
D[(p+2):(p+pz+2),(p+2):(p+pz+2)]<- -hessian(Logl.sandACT,x=deltasigma1,X=Xz,z=z)

Aninv<-N/n1*c(1,-c(a12,a13)%*%solve(D))

# ptm<-proc.time()
# An<-matrix(data=0,nrow=2*p+3,ncol=2*p+3)
# An[1,1]<-n1/N
# An[1,2:(2*p+3)]<-c(a12,a13)
# An[2:(2*p+3),2:(2*p+3)]<-D
# sqrt(diag(solve(An)%*%Bn%*%solve(t(An))))
# proc.time()-ptm


phi<-matrix(nrow=N,ncol=3+p+pz)
phi[,1]<-(z*(y-X%*%BetaOLSy0-NaivEst)- (1-z)*(y-X%*%BetaOLSy0)/(1-phat)) 
phi[,2:(p+2)]<-apply(X,2,function(x){(1-z)*(y-X%*%BetaOLSy0)*x}) 
phi[,(p+3):(p+pz+3)]<-apply(Xz,2,function(x){(z*lambda1(u)-(1-z)*lambda0(u))*x}) 

# ptm<-proc.time()
# nBn<-array(dim=c(3+2*p,3+2*p,n))
# for(i in 1:n){
	# nBn[,,i]<-t(t(phi[i,]))%*%t(phi[i,])
# }
# nBn<-apply(nBn,c(1,2),sum)/N^2
# proc.time()-ptm

Bn<-matrix(nrow=(3+p+pz),ncol=(3+p+pz))
 for(i in 1:(3+p+pz)){
 	for(j in 1:(3+p+pz)){
 	Bn[i,j]<-mean(phi[,i]*phi[,j])
 	}
 }
Bn<-Bn/N

return(Aninv%*%Bn%*%Aninv)
}

