#' Plot of UI and CI
#'
#' Plot function for objects returned from \code{\link{ui.causal}}.
#' Plots confidence intervals for different values of rho0=rho1=rho.
#' @param fitted An object of class uicausal
#' @param DR If TRUE, plots both DR if FALSE OR.
#' @param main Main title, default is no title.
#' @param xlab Title for x-axis, default is \code{expression(rho)}.
#' @param ylab Title for y-axis, default is the variable names.
#' @param ... Additional arguments, use is discouraged.
#'
#'
#' @importFrom graphics lines par plot polygon
#' @export


profile.uicausal<-function(fitted, DR=TRUE, main='', xlab=NULL,ylab='', ...){

if(DR==TRUE){
	ui.obj<-fitted$DR
	}else{
		ui.obj<-fitted$OR
	} 

if(length(fitted$rho0)==1){
	stop("Not able to plot profile since rho0 and/or rho1 is an integer.")
}	
if(!is.null(fitted$rho1)){
	if(length(fitted$rho1)==1){
		stop("Not able to plot profile since rho0 and/or rho1 is an integer.")
	}	
	if((sum(fitted$rho1!=fitted$rho0)>0)){
	stop("Not able to plot profile since rho0 and rho1 are unequal intervals.")
	}}

if(is.null(fitted$rho1)){
	coef<-ui.obj$coef
	ci<-ui.obj$ci
	}else{
	coef<-diag(ui.obj$coef)
	ci<-cbind(diag(ui.obj$ci[,,1]),diag(ui.obj$ci[,,2]))
	}

if(is.null(xlab)){xlab<-expression(rho)}

gridrho<-fitted$gridrho0
ui<-ui.obj$ui

plot(gridrho,coef,type='l',ylim=c(min(c(0, ui)),max(c(0, ui))),xlab=xlab,ylab=ylab)
polygon(c(gridrho, rev(gridrho)), c(ci[,2],rev(ci[,1])), col = "grey90", border = NA)
lines(gridrho,coef)
lines(gridrho,ci[,1],lty=2)
lines(gridrho,ci[,2],lty=2)
lines(c(-1,1),c(0,0))
	
}
