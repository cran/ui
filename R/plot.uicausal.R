#' Plot of UI and CI
#'
#' Plot function for objects returned from \code{\link{ui.causal}}.
#' Plots confidence intervals for different values of rho and the uncertainty interval.
#' @param x An object of class uicausal
#' @param DR If TRUE the doubly robust estimator is plotted, otherwise the outcome regression estimator is plotted.
#' @param main Main title, default is no title.
#' @param ylab Title for y axis, default is no title.
#' @param xlab Title for xaxis, default is \code{expression(rho)}.
#' @param ... Additional arguments, use is discouraged.
#'
#' @importFrom graphics arrows axis lines plot polygon text
#' @export


plot.uicausal<-function(x,DR=TRUE,main='',xlab=NULL,ylab='', ...){
if(DR==TRUE){
	ex<-x$plot$DR
	ui.obj<-x$DR
	}else{
		ex<-x$plot$DR
		ui.obj<-x$OR
	}

if(x$plot$plot==FALSE){
	stop("Not able to plot uicausal object, with plot==FALSE. In ui.causal() rho0 and rho1 must be equal intervals otherwise plot is changed to FALSE.")
}

rho<-x$plot$gridrho0
if(is.null(x$rho1)){
	coef<-ex$coef
	ci<-cbind(ex$ci[,1],ex$ci[,2])
	}else{
	coef<-diag(ex$coef)
	ci<-cbind(diag(ex$ci[,,1]),diag(ex$ci[,,2]))
	}

if(is.null(xlab)){xlab<-expression(rho)}

plot(rho,coef,type='l',mgp=c(2,1,0),ylim=c(min(ci[,1]),max(ci[,2])),ylab=ylab,xlab=xlab,xaxs='i',main=main)

polygon(c(rho, rev(rho)), c(ci[,2],rev(ci[,1])), col = "grey90", border = NA)

lines(rho,coef)
lines(rho,ci[,1],lty=2)
lines(rho,ci[,2],lty=2)
lines(c(-1,1),c(0,0))


Stn<-x$plot$gridn
Fin<-x$plot$gridn*2-1
n0<-which(round(rho,10)==0)
if(length(n0)>1){
	n0<-which(abs(rho)==min(abs(rho)))
}

if(n0==Stn|n0==Fin){
nv<-Stn+round((Fin-Stn)/4)*0:4
if(n0==Stn){nv<-nv[-1]}else{nv<-nv[-5]}
}else{
nv<-sort(Stn+round((Fin-Stn)/3)*0:3)
}

minR<-which(ci[Stn:Fin,1]==min(ci[Stn:Fin,1]))+(Stn-1)
maxR<-which(ci[Stn:Fin,2]==max(ci[Stn:Fin,2]))+(Stn-1)

lines(c(rho[4],rho[minR]),c(ci[minR,1],ci[minR,1]),lty=2,col='blue')
lines(c(rho[4],rho[maxR]),c(ci[maxR,2],ci[maxR,2]),lty=2,col='blue')
arrows(x0=rho[Stn],y0=min(ci[,1])-10,y1=ci[Stn,1],code=3,length=0,lty=2,col='blue')
arrows(x0=rho[Fin],y0=min(ci[,1])-10,y1=ci[Fin,1],code=3,length=0,lty=2,col='blue')



arrows(x0=0,y0=ci[n0,1],y1=ci[n0,2],code=3,length=0.1,col='red')
for(i in nv){
arrows(x0=rho[i],y0=ci[i,1],y1=ci[i,2],code=3,length=0.1,lty=3,lwd=1.5)
}

arrows(x0=rho[4],y0=min(ci[Stn:Fin,1]),y1=max(ci[Stn:Fin,2]),code=3,length=0.1,col='blue')


text(mean(rho[2:3]),mean(ci),'UI',col='blue')
axis(1,at=round(c(rho[Stn],rho[Fin]),3),col.axis='blue')
}
