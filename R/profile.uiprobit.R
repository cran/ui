#' Plot of UI and CI
#'
#' Plot function for objects returned from \code{\link{ui.probit}}.
#' Plots confidence intervals for different values of rho and the uncertainty interval.
#' @param fitted An object of class uiprobit
#' @param plot.all If TRUE, plots all covariates.
#' @param which Specify which variables should be plotted by either sending in their names in a vector or a vector with their numbers (1 intercept, 2 for the first covariate etc.).
#' @param intercept If TRUE, also plots the intercept.
#' @param xlab Title for x-axis, default is \code{expression(rho)}.
#' @param ylab Title for y-axis, default is the variable names.
#' @param cex.lab Size of lables.
#' @param mar Margin around panels in plot.
#' @param ... Additional arguments, use is discouraged.
#'
#'
#' @importFrom graphics par plot polygon lines
#' @export


profile.uiprobit<-function(fitted, plot.all=TRUE, which=NA, intercept=FALSE, xlab=NULL,ylab=NULL, cex.lab=2,mar=c(6,6,2,2), ...){
p<-dim(fitted$coef)[1]-1
nrho<-dim(fitted$coef)[2]

if(sum(is.na(which)>0)){
	if(plot.all==FALSE) warning('Need to specify which variable in order to not plot all.')
	plot.all=TRUE
}else{
	plot.all=FALSE
	if(mode(which)=='character'){
		which<-which(NamnX %in% which)
	}
}
if(plot.all){
	if(intercept){which<- 1:(p+1)}else{which<- 2:(p+1)}
	}

if(is.null(ylab)){NamnX<-rownames(fitted$coef)
	}else{if(length(as.vector(ylab))==1){
	NamnX<-rep(ylab,length(rownames(fitted$coef)))	
	}else{
		if(length(as.vector(ylab))==length(which)){
			NamnX<-rep("",length(rownames(fitted$coef)))
	NamnX[which]=ylab	}else{
		stop('Wrong dimensions on ylab.')
	}}}
if(is.null(xlab)){xlab=expression(rho)}
	
	
if(sum(is.na(fitted$ui))>0){
	ui<-fitted$uirough
	ci<-fitted$cirough
}else{
	ui<-fitted$ui
	ci<-fitted$ci
}
n<-length(which)


dim.par<-c(floor(sqrt(n)),ceiling(sqrt(n)))

if(dim.par[1]*dim.par[2]<n){
	dim.par[1]<-dim.par[1]+1
}


par(mfrow=dim.par,mar=mar)

for(i in which){
plot(fitted$rho,fitted$coef[i,],type='l',ylim=c(min(c(0, ui[i,])),max(c(0, ui[i,]))),xlab=xlab,ylab=NamnX[i],cex.lab=cex.lab)
polygon(c(fitted$rho, rev(fitted$rho)), c(ci[i,,2],rev(ci[i,,1])), col = "grey90", border = NA)
lines(fitted$rho,fitted$coef[i,])
lines(fitted$rho,ci[i,,1],lty=2)
lines(fitted$rho,ci[i,,2],lty=2)
lines(c(-1,1),c(0,0))
}	
}
