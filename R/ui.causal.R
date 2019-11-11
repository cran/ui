#' Uncertainty intervals for Average Causal Effects
#'
#' This function allows you to derive uncertainty intervals for the average causal effect (ACE) or the average causal effect on the treated (ACT). The function uses a regression imputation estimator and a doubly robust estimator. The uncertainty intervals can be used as a sensitivity analysis to unconfoundedness. Note that \code{rho}=0 render the same results as assuming no unobserved confounding.
#' @param out.formula Formula for the outcome regression models
#' @param treat.formula Formula for the propensity score model (regression model for treatment assignment).
#' @param data data.frame containing the variables in the formula.
#' @param rho Pre-specified interval for \code{rho0} and \code{rho1}.
#' @param rho0 Pre-specified value of \code{rho0}, if an interval it has to be the same as \code{rho1}.
#' @param rho1 Pre-specified value of \code{rho1}, if an interval it has to be the same as \code{rho0}.
#' @param ACT If TRUE Average Causal effect of the Treated is calculated, if FALSE Average Causal effect is calculated. Default is FALSE.
#' @param sand Specifies which estimator of the standard errors should be used for OR, see details.
#' @param gridn Number of fixed points within the \code{rho} interval for which sigma0 and sigma1 should be estimated.
#' @param plot If TRUE the function runs slightly slower but you will be able to plot your results using \code{\link{plot.uicausal}}.
#' @param rho.plotrange an interval larger than \code{rho} for the plot using \code{\link{plot.uicausal}}.
#' @param alpha Default 0.05 corresponding to a confidence level of 95 for CI and UI.
#'
#'@details In order to visualize the results, you can use \code{\link{plot.uicausal}}. Details about estimators can be found in Genbäck and de Luna (2018)
#'
#'The standard errors are calculated with the following estimators:
#'
#'DR ACE - simplified sandwich estimator
#'
#'DR ACT - sandwich estimator
#'
#'OR ACE - if sand=TRUE sandwich estimator (default and recommended), if sand=FALSE large sample variance
#'
#'OR ACT - if sand=TRUE sandwich estimator (default and recommended), if sand=FALSE large sample variance
#'
#'
#' @return A list containing:
#' \item{call}{The matched call}
#' \item{rho0}{The rage of \code{rho0} from which the ui is calculated}
#' \item{rho1}{If ACT==FALSE,range of \code{rho1} from which the ui is calculated}
#' \item{out.model0}{Outcome regression model for non-treated.}
#' \item{out.model1}{Outcome regression model for treated.}
#' \item{treat.model}{Regression model for treatment mechanism (propensity score).}
#' \item{sigma0}{Consistent estimate of sigma0 for different values of rho0}
#' \item{sigma1}{Consistent estimate of sigma1 for different values of rho1}
#' \item{DR}{DR inference, confidence intervals for different pre-specified values of \code{rho} for the OR estimator, uncertainty interval, coefficient estimates, confounding bias, indentification interval, standard error etc.}
#' \item{OR}{OR inference, confidence intervals for different pre-specified values of \code{rho} for the OR estimator, uncertainty interval, coefficient estimates, confounding bias, indentification interval, standard error etc.}
#'
#' @importFrom  Matrix rankMatrix
#'@author Minna Genbäck
#'@references Genbäck, M., de Luna, X. (2018). Causal Inference Accounting for Unobserved Confounding after Outcome Regression and Doubly Robust Estimation. \emph{Biometrics}. DOI: 10.1111/biom.13001
#' @examples
#'library(MASS)
#'n<-500
#'delta<-c(-0.3,0.65)
#'rho<-0.3
#'X<-cbind(rep(1,n),rnorm(n))
#'x<-X[,-1]
#'s0<-2
#'s1<-3
#'error<-mvrnorm(n, c(0,0,0), matrix(c(1,0.6,0.9,0.6,4,0.54,0.9,0.54,9), ncol=3))
#'zstar<-X%*%delta+error[,1]
#'z<- zstar>0
#'y1<-ifelse(x< (-1),0.2*x-0.1*x^2, ifelse(x< 1,0.3*x, ifelse(x<3,0.4-0.1*x^2,-0.2-0.1*x)))+error[,3]
#'y0<-ifelse(x<1.5, x-0.4*x^2, ifelse(x<2, -0.15-0.25*x+0.5*x^2, 1.85-0.25*x))+error[,2]
#'y<-y0
#'y[z==1]<-y1[z==1]
#'data<-data.frame(y,z,x)
#'
#'
#'ui<-ui.causal(y~x, z~x, data=data, rho=c(0,0.3), ACT=FALSE)
#'ui
#'plot(ui)
#'profile(ui)
#'mean(y1-y0)
#'
#'ui<-ui.causal(y~x, z~x, data=data, rho=c(0,0.3), ACT=TRUE)
#'ui
#'plot(ui)
#'mean(y1[z==1]-y0[z==1])
#'
#' @importFrom stats binomial coef complete.cases cov get_all_vars glm lm model.matrix pnorm qnorm
#' @export

ui.causal<-function(out.formula,treat.formula,data,rho=c(-0.3,0.3),rho0=NULL,rho1=NULL,ACT=FALSE, sand=TRUE,gridn=21, plot=TRUE, 
rho.plotrange=c(-0.5,0.5), alpha=0.05){
 
 
 if(is.null(rho0)&is.null(rho1)){
 	rho0<-rho
 	rho1<-rho
 }else{
 	if(is.null(rho0)) rho0<-rho1
 	if(is.null(rho1)) rho1<-rho0
 }
 rho0<-sort(rho0)
 rho1<-sort(rho1)
 
 ###Warnings
 if(class(data)!="data.frame"){
	stop('Data must be a data frame')
  }
  
  y.data<-get_all_vars(out.formula,data=data)
  z.data<-get_all_vars(treat.formula,data=data)
  
#  remove individuals with partial missing data in the covariates
  if(sum(complete.cases(y.data)==FALSE)>0){
	warning(paste('Partial missing values in covariates! ',sum(complete.cases(y.data)==FALSE),'individual(s) are removed from the outcome regression.'))
	y.data<-y.data[complete.cases(y.data),]
  }
  if(sum(complete.cases(z.data)==FALSE)>0){
	warning(paste('Partial missing values in covariates. ',sum(complete.cases(z.data)==FALSE),'individual(s) are removed from the propensity score regression.'))
	z.data<-z.data[complete.cases(z.data),]
  }
    
  output<-list()
	if(plot==TRUE){
    if((length(rho0)==1&length(rho1)==1)|sum(rho0!=rho1)>0){
    	#warning('In order to plot results rho0 and rho1 must be equal intervals, plot changed to FALSE.')
    	plot=FALSE
    	}
	if(min(rho.plotrange)>=min(rho)){
	#	warning('Lower bound of plotrange is >= lower bound of UI, it needs to be less than. You will not be able to plot this object.')
		rho.plotrange[1]=(min(rho)-1)/2}
	if(max(rho.plotrange)<=max(rho)){
		#warning('Upper bound of plotrange is <= upper bound of UI, it needs to be greater than. You will not be able to plot this object.')
		rho.plotrange[2]=(max(rho)+1)/2}
	}
	
	output$plot<-list()
    output$plot$plot=plot	
  
  z<-z.data[,1]
  y<-y.data[,1]
  y0<-y[z==0]
  y1<-y[z==1]
      
  X<-model.matrix(lm(out.formula,data=y.data))
  out.model0<-lm(out.formula,data=y.data[z==0,])  
  out.model1<-lm(out.formula,data=y.data[z==1,])  
  treat.model<-glm(treat.formula,family=binomial(link="probit"),data=z.data)
  
  output$out.model0<-out.model0
  output$out.model1<-out.model1
  output$treat.model<-treat.model
  
  X1<-model.matrix(out.model1)
  X0<-model.matrix(out.model0)
  Xz<-model.matrix(treat.model)

  d0<-dim(X0)
  d1<-dim(X1)
  n0<-d0[1]
  n1<-d1[1]
  N<-n0+n1
  p<-d0[2]-1
  pz<-dim(Xz)[2]-1
	
  gamma<-coef(summary(treat.model))[,1]
	
  #Estimating BetaOLS untreated
  BetaOLSy0<-coef(summary(out.model0))[,1]
  varBetaOLSy0<-diag((coef(summary(out.model0))[,2])^2)
  sigma0hatOLS<-summary(out.model0)$sigma
  #Estimating BetaOLS treated
  BetaOLSy1<-coef(summary(out.model1))[,1]
  varBetaOLSy1<-diag((coef(summary(out.model1))[,2])^2)
  sigma1hatOLS<-summary(out.model1)$sigma

  u<-(Xz%*%gamma)
  u0<-(Xz[z==0,]%*%gamma)
  u1<-(Xz[z==1,]%*%gamma)
  OLSlambda0<-solve(t(X0)%*%X0)%*%t(X0)%*%lambda0(u0)
  l0z0<-mean(lambda0(u0))
  l0<-mean(lambda0(u))
  OLSlambda1<-solve(t(X1)%*%X1)%*%t(X1)%*%lambda1(u1)
  l1z1<-mean(lambda1(u1))
  l1<-mean(lambda1(u))

  
  t0<-gridrho.f(rho0,gridn,rho.plotrange,plot)
  gridrho0<-t0[[1]]
  nui0<-t0[[2]]
  output$plot$nui0<-nui0
  
  output$plot$sigma0<-sigmaOLScor0(X0,sigma0hatOLS,n0,p,u0,gridrho0)
  output$sigma0<-output$plot$sigma0[nui0]
  output$rho0<-rho0
  output$plot$gridrho0<-gridrho0
  output$gridrho0<-gridrho0[nui0]
  	
  if(ACT==FALSE){
  t1<-gridrho.f(rho1,gridn,rho.plotrange,plot)
  gridrho1<-t1[[1]]
  nui1<-t1[[2]]
  output$rho1<-rho1
  output$plot$gridrho1<-gridrho1
  output$gridrho1<-gridrho1[nui1]
  
  output$plot$sigma1<-sigmaOLScor1(X1,sigma1hatOLS,n1,p,u1,gridrho1)
  output$sigma1<-output$plot$sigma1[nui1]
  output$plot$nui1<-nui1
  }
  	
  output$DR$conf.bias<-vector(length=2)
  names(output$DR$conf.bias)<-c('Min', 'Max')
  output$DR$IdentInt<-output$DR$conf.bias
  output$DR$ui<-output$DR$conf.bias
  output$DR$naivci<-output$DR$conf.bias
  output$OR<-output$DR
  
  ##################### NAIV ESTIMATES #####################
  phat<-pnorm(Xz%*%gamma)
  phat0<-phat[z==0]
  if(ACT==TRUE){
  	output$OR$naiv<-mean(y1-X1%*%BetaOLSy0)
  	output$DR$naiv<-output$OR$naiv-sum((y0-X0%*%BetaOLSy0)/(1-phat0))/n1
  }else{
  	output$OR$naiv<-(sum(X1%*%BetaOLSy1-X1%*%BetaOLSy0)+sum(X0%*%BetaOLSy1-X0%*%BetaOLSy0))/N
	phat1<-phat[z==1]
	output$DR$naiv<-output$OR$naiv+(sum((y1-X1%*%BetaOLSy1)/phat1)-sum((y0-X0%*%BetaOLSy0)/(1-phat0)))/N
  }
		
	##################### CONFOUNDING BIAS #####################
	exgz1<-rep(1,n1)%*%X1/n1
	exgz0<-rep(1,n0)%*%X0/n0
	
	
	if(ACT==TRUE){
	#Counfounding bias regression imputation ACT
	f<-gridrho0*output$plot$sigma0*c(l1z1+exgz1%*%OLSlambda0)
	names(f)<-as.character(round(gridrho0,4))
	output$plot$OR$conf.bias.grid<-f
	output$OR$conf.bias.grid<-f[nui0]
	output$OR$conf.bias<-c(min(f[nui0]),max(f[nui0]))
	#Total bias DR ACT
	f<-gridrho0*output$plot$sigma0*l0*N/n1
	names(f)<-as.character(round(gridrho0,4))
	output$plot$DR$conf.bias.grid<-f
	output$DR$conf.bias.grid<-f[nui0]
	output$DR$conf.bias<-c(min(f[nui0]),max(f[nui0]))
	
	}else{
	#Counfounding bias regression imputation ACE
	f1<-gridrho1*output$plot$sigma1*c(rep(1,N)%*%X%*%OLSlambda1)/N
	f2<-gridrho0*output$plot$sigma0*c(rep(1,N)%*%X%*%OLSlambda0)/N
	output$plot$OR$conf.bias.grid<-matrix(rep(f1,length(f2)),nrow=length(f1))+t(matrix(rep(f2,length(f1)),nrow=length(f2)))
	rownames(output$plot$OR$conf.bias.grid)<-as.character(round(gridrho1,4))
	colnames(output$plot$OR$conf.bias.grid)<-as.character(round(gridrho0,4))
	output$OR$conf.bias.grid<-output$plot$OR$conf.bias.grid[nui1,nui0]
	output$OR$conf.bias<-c(min(output$OR$conf.bias.grid),max(output$OR$conf.bias.grid))
	#Total bias DR ACE
	f1<-gridrho1*output$plot$sigma1*l1
	f2<-gridrho0*output$plot$sigma0*l0
	output$plot$DR$conf.bias.grid<-matrix(rep(f1,length(f2)),nrow=length(f1))+t(matrix(rep(f2,length(f1)),nrow=length(f2)))
	dimnames(output$plot$DR$conf.bias.grid)<-dimnames(output$plot$OR$conf.bias.grid)
	output$DR$conf.bias.grid<-output$plot$DR$conf.bias.grid[nui1,nui0]
	output$DR$conf.bias<-c(min(output$DR$conf.bias.grid),max(output$DR$conf.bias.grid))
	}
	
	#Coef, NAIV ESTIMATES - BIAS (for each value of rho)
	#Want the same rownames and colnames as conf.bias.grid
	output$plot$DR$coef<-output$plot$DR$conf.bias.grid
	output$plot$DR$coef<-c(output$DR$naiv)-output$plot$DR$conf.bias.grid
	output$plot$OR$coef<-output$plot$OR$conf.bias.grid
	output$plot$OR$coef<-c(output$OR$naiv)-output$plot$OR$conf.bias.grid
	output$DR$coef<-output$DR$conf.bias.grid
	output$DR$coef<-c(output$DR$naiv)-output$DR$conf.bias.grid
	output$OR$coef<-output$OR$conf.bias.grid
	output$OR$coef<-c(output$OR$naiv)-output$OR$conf.bias.grid
	
	#Identification interval, NAIV ESTIMATES - BIAS (interval)
	output$OR$IdentInt<-rep(output$OR$naiv,2)-output$OR$conf.bias[2:1]
	output$DR$IdentInt<-rep(output$DR$naiv,2)-output$DR$conf.bias[2:1]
	
	
	#Standard errors
	if(ACT==TRUE){
		if(sand==TRUE){
		output$OR$se<-sqrt(sandImpACT(X,y,z,BetaOLSy0,output$OR$naiv,n1,N,p))
		}else{
		output$OR$se<-sqrt((2*sigma1hatOLS^2+  (BetaOLSy1-BetaOLSy0)%*%cov(X1)%*%(BetaOLSy1-BetaOLSy0))/n1)
		}
		output$DR$se<-sqrt(sandACT(gamma,X,Xz,y,z,u,BetaOLSy0,phat,output$DR$naiv,n1,n0,N,p,pz))
	}else{
		if(sand==TRUE){
		output$OR$se<-sqrt(sandImpACE(X,y,z,BetaOLSy0,BetaOLSy1,output$OR$naiv,N,p))
			}else{
		output$OR$se<-sqrt(rep(1,N)%*%X%*%(sigma1hatOLS^2*solve(t(X1)%*%X1)+sigma0hatOLS^2*solve(t(X0)%*%X0))%*%t(X)%*%rep(1,N)/N^2 + (BetaOLSy1-BetaOLSy0)%*%cov(X)%*%(BetaOLSy1-BetaOLSy0)/N)	
			}
		I<-(X%*%BetaOLSy1-X%*%BetaOLSy0+z*((y-X%*%BetaOLSy1)/phat)- (1-z)*((y-X%*%BetaOLSy0)/(1-phat)))-output$DR$naiv
		output$DR$se<-sqrt(sum(I^2))/N
		}
	zalpha<-qnorm(1-alpha/2)
	#Confidence intervals
		if(ACT){
		dim<-c(length(output$plot$OR$conf.bias.grid),2)
		dnames<-list(names(output$plot$OR$conf.bias.grid),c('lower','upper'))
	}else{
		dim<-c(dim(output$plot$OR$conf.bias.grid),2)
		dnames<-dimnames(output$plot$OR$conf.bias.grid)
		dnames[[3]]<-c('lower','upper')
		}
		
	output$plot$OR$ci<-array(c(c(output$OR$naiv-zalpha*output$OR$se)-output$plot$OR$conf.bias.grid, c(output$OR$naiv+zalpha*output$OR$se)-		
						output$plot$OR$conf.bias.grid),dim=dim,dimnames=dnames)
	output$plot$DR$ci<-array(c(c(output$DR$naiv-zalpha*output$DR$se)-output$plot$DR$conf.bias.grid, c(output$DR$naiv+zalpha*output$DR$se)-
						output$plot$DR$conf.bias.grid),dim=dim,dimnames=dnames)
						
	if(ACT){
		dim<-c(length(output$OR$conf.bias.grid),2)
		dnames<-list(names(output$OR$conf.bias.grid),c('lower','upper'))
	}else{
		dim<-c(dim(output$OR$conf.bias.grid),2)
		dnames<-dimnames(output$OR$conf.bias.grid)
		dnames[[3]]<-c('lower','upper')
		}
		
	output$OR$ci<-array(c(c(output$OR$naiv-zalpha*output$OR$se)-output$OR$conf.bias.grid, c(output$OR$naiv+zalpha*output$OR$se)-		
						output$OR$conf.bias.grid),dim=dim,dimnames=dnames)
	output$DR$ci<-array(c(c(output$DR$naiv-zalpha*output$DR$se)-output$DR$conf.bias.grid, c(output$DR$naiv+zalpha*output$DR$se)-
						output$DR$conf.bias.grid),dim=dim,dimnames=dnames)

	output$OR$ui<-c(min(output$OR$ci),max(output$OR$ci))
	output$DR$ui<-c(min(output$DR$ci),max(output$DR$ci))
	output$OR$naivci<-output$OR$naiv+c(-zalpha*output$OR$se,zalpha*output$OR$se)
	output$DR$naivci<-output$DR$naiv+c(-zalpha*output$DR$se,zalpha*output$DR$se)
	output$plot$gridn<-gridn
	
	names(output$OR$ui)<-c('lower','upper')
	names(output$OR$naivci)<-c('lower','upper')
	names(output$OR$conf.bias)<-c('lower','upper')
	names(output$OR$IdentInt)<-c('lower','upper')
	names(output$DR$ui)<-c('lower','upper')
	names(output$DR$naivci)<-c('lower','upper')
	names(output$DR$conf.bias)<-c('lower','upper')
	names(output$DR$IdentInt)<-c('lower','upper')
	
	
	output$call <- match.call()
	class(output)<-"uicausal"
return(output)
}


