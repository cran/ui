#' Uncertainty intervals for OLS regression
#'
#' This function allows you to derive uncertainty intervals for OLS regression when there is missing data in the continuous outcome. The uncertainty intervals can be used as a sensitivity analysis to ignorability (missing at random). Note that rho=0 render the same results as a complete case analysis.
#' @param out.formula Formula for outcome regression.
#' @param mis.formula Formula for missingness mechanism. If NULL the same covariates as in the outcome regression will be used. 
#' @param data data.frame containing the variables in the formula.
#' @param rho The limits of rho for which the uncertainty interval should be constructed.
#' @param alpha Default 0.05 corresponding to a confidence level of 95 for CI and UI.
#' @param gridn The number of distinct points within the interval \code{rho} at which confidence intervals should be constructed. Default is 101.
#'@details In order to visualize the results, you can use \code{\link{plot.uiols}}, 
#' or \code{\link{profile.uiols}}.
#' @importFrom  Matrix rankMatrix
#' @return A list containing:
#' \item{call}{The matched call}
#' \item{ci}{Confidence intervals for different values of \code{rho}}
#' \item{ui}{Uncertainty intervals}
#' \item{coef}{Estimated coefficients (outcome regression) for different values of \code{rho}}
#' \item{out.model}{Outcome regression model when rho=0.}
#' \item{mis.model}{Regression model for missingness mechanism (selection).}
#' \item{rho}{The range of \code{rho} for which we want to construct an uncertainty interval}
#' \item{gridrho}{The values of \code{rho} for which bias and standard errors are derived}
#' \item{sigma}{Consistant estimate of sigma}
#' \item{se}{Standard error for different values of \code{rho}}
#' \item{ciols}{Confidence intervals from a complete case analysis}
#' \item{ident.bound}{Bounds for the coefficient estimates.}
#'@author Minna Genbäck
#'@references Genbäck, M., Stanghellini, E., de Luna, X. (2015). Uncertainty Intervals for Regression Parameters with Non-ignorable Missingness in the Outcome. \emph{Statistical Papers}, 56(3), 829-847.
#' @examples
#'library(MASS)
#'n<-500
#'delta<-c(0.5,0.3,0.1)
#'beta<-c(0.8,-0.2,0.3)
#'X<-cbind(rep(1,n),rnorm(n),rbinom(n,1,0.5))
#'x<-X[,-1]
#'rho=0.4
#'error<-mvrnorm(n,c(0,0),matrix(c(1,rho*2,rho*2,4),2))
#'zstar<-X%*%delta+error[,1]
#'z<-as.numeric(zstar>0)
#'y<-X%*%beta+error[,2]
#'y[z==0]<-NA
#'data<-data.frame(y,x,z)
#' ui<-ui.ols(y~X1+X2,data=data,rho=c(-0.5,0.5))
#' ui
#' plot(ui)
#'
#' @importFrom stats binomial coef get_all_vars glm lm model.matrix qnorm as.formula
#' @export
ui.ols<-function(out.formula,mis.formula=NULL,data,rho=c(-0.3,0.3),alpha=0.05,gridn=101){
  ###Warnings
   if(class(data)!="data.frame"){
	stop('Data must be a data frame')
    }
    
    y<-get_all_vars(out.formula,data=data)[,1]
  	z <- as.numeric(!is.na(y))
  	data$z<-z

    out.model<-lm(out.formula,data=data[z==1,])
    
  if(is.null(mis.formula)){
  	mis.formula<-update.formula(out.formula, as.formula(paste('z', "~ ."))) 
    }
    
	mis.model<-glm(mis.formula,family=binomial(link="probit"),data=data)

    X<-model.matrix(out.model)
    Xz<-model.matrix(mis.model)[z==1,]

	d<-dim(X)
	n<-d[1]
	p<-d[2]-1 
	pz<-dim(Xz)[2]-1

	output=list()
	
	output$out.model<-out.model
	output$mis.model<-mis.model
    #estimating delta and u 
    deltasigma1<-coef(summary(mis.model))[,1]
	u<-(Xz%*%deltasigma1)
		
	#Estimating BetaOLS
	OLSC<-coef(summary(out.model))
	BetaOLS<-OLSC[,1]
	sigmaOLS<-summary(out.model)$sigma
	BetaOLSse<-OLSC[,2]
	
	zalpha<-qnorm(1-alpha/2)
	output$zalpha<-zalpha
	output$ciols<-cbind(BetaOLS-zalpha*BetaOLSse, BetaOLS+zalpha*BetaOLSse)
	#output$coefols<-BetaOLS
	
	if(length(rho)>1){
	output$gridrho<-rho[1]+((rho[2]-rho[1])/(gridn-1))*0:(gridn-1)
	}else{output$gridrho<-rho
		gridn=1}
	output$sigma<-sigmaOLScor1(X,sigmaOLS,n,p,u,output$gridrho)

	# Estimation of the term with lambda1 for bound
	OLSlambda<-solve(t(X)%*%X)%*%t(X)%*%lambda1(u)
	
	if(length(rho)>1){
	Mt<-matrix(output$gridrho*output$sigma,ncol=1)%*%matrix(OLSlambda,nrow=1)
	output$coef<-apply(Mt,1,function(x)BetaOLS-x)
	output$ident.bound<-cbind(apply(output$coef,1,min),apply(output$coef,1,max))
	}else{
	output$coef<-BetaOLS-output$gridrho*as.vector(output$sigma)*OLSlambda
	output$ident.bound<-cbind(output$coef,output$coef)
	}
	
	#output$seols<-BetaOLSse
	output$se<-t(se.ols(X,output$sigma,u,output$gridrho))
	

	output$ci<-array(dim=c(p+1,gridn,2) ,dimnames=list(colnames(X),rho=as.character(round(output$gridrho,4)),c('lower','upper')))
	output$ci[,,1]<-output$coef-zalpha*output$se
	output$ci[,,2]<-output$coef+zalpha*output$se
	output$ui<-cbind(apply(output$ci[,,1],1,min),apply(output$ci[,,2],1,max))
	
	rownames(output$ciols)<-colnames(X)
	#names(output$coefols)<-colnames(X)
	rownames(output$ident.bound)<-colnames(X)
	rownames(output$ui)<-colnames(X)
	rownames(output$coef)<-colnames(X)
	rownames(output$se)<-colnames(X)
	
	colnames(output$coef)<-as.character(round(output$gridrho,4))
	colnames(output$se)<-colnames(output$coef)
	colnames(output$ciols)<-c('lower','upper')
	colnames(output$ident.bound)<-c('lower','upper')	
	colnames(output$ui)<-c('lower','upper')
	output$call <- match.call()
	output$rho <- rho
	class(output)<-"uiols"
	output
	}
	