#' Calculates standard error of Average causal effect on the treated
#' 
#' This is a support function for \code{\link{ui.causal}} and calculates standard error of Average causal effect on the treated for the regression imputation estimator.
#' @param X Covariate matrix.
#' @param y Outcome vector.
#' @param z missingness indicator
#' @param BetaOLSy0 Coefficients from non-treated regression 
#' @param NaivEst Naiv estimates.
#' @param n1 Number of treated.
#' @param N Total number.
#' @param p Number of covariates outcome regression.
#' @export


sandImpACT<-function(X,y,z,BetaOLSy0,NaivEst,n1,N,p){

a12<-apply(X, 2,function(x){mean(z*x)})

D<-matrix(nrow=(p+1),ncol=(p+1),data=0)

D[1:(p+1),1:(p+1)]<-t((1-z)*X)%*%((1-z)*X)/N

Aninv<-N/n1*c(1,-a12%*%solve(D))

# ptm<-proc.time()
# An<-matrix(data=0,nrow=2*p+3,ncol=2*p+3)
# An[1,1]<-n1/N
# An[1,2:(2*p+3)]<-c(a12,a13)
# An[2:(2*p+3),2:(2*p+3)]<-D
# sqrt(diag(solve(An)%*%Bn%*%solve(t(An))))
# proc.time()-ptm


phi<-matrix(nrow=N,ncol=p+2)
phi[,1]<-(z*(y-X%*%BetaOLSy0-NaivEst)) 
phi[,2:(p+2)]<-apply(X,2,function(x){(1-z)*(y-X%*%BetaOLSy0)*x}) 

# ptm<-proc.time()
# nBn<-array(dim=c(3+2*p,3+2*p,n))
# for(i in 1:n){
	# nBn[,,i]<-t(t(phi[i,]))%*%t(phi[i,])
# }
# nBn<-apply(nBn,c(1,2),sum)/N^2
# proc.time()-ptm

Bn<-matrix(nrow=(p+2),ncol=(p+2))
 for(i in 1:(p+2)){
 	for(j in 1:(p+2)){
 	Bn[i,j]<-mean(phi[,i]*phi[,j])
 	}
 }
Bn<-Bn/N
Aninv[1]*Bn[1,1]
return(Aninv%*%Bn%*%Aninv)
}

