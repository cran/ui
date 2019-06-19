#' Plot of UI and CI
#'
#' Plot function for objects returned from \code{\link{ui.ols}}. Plots confidence intervals,
#' coefficients and significans assuming ignorability and the uncertainty interval under non-ignorability.
#' @param x An object of class uiols
#' @param plot.all If TRUE, plots all covariates.
#' @param which Specify which variables should be plotted by either sending in their names in a vector or a vector with their numbers (1 intercept, 2 for the first covariate etc.).
#' @param intercept If TRUE, also plots the intercept.
#' @param ylab Vector of names for the y-axis, default is the variable names.
#' @param col Vector containing the color of confidence intervals (default black) and uncertainty intervals (default red).
#' @param ... Additional arguments, use is discouraged.
#'
#'
#' @importFrom graphics arrows axis lines par plot points text
#' @importFrom stats pnorm
#' @export

plot.uiols<-function(x,plot.all=TRUE, which=NA,intercept=FALSE,ylab=NULL,col=c("black","red"), ...){
p<-dim(x$coef)[1]-1


if(sum(is.na(which)>0)){
	if(plot.all==FALSE) warning('Need to specify which variable in order to not plot all.')
	plot.all=TRUE
}else{
	plot.all=FALSE
	if(mode(which)=='character'){
		which<-which(rownames(x$coef) %in% which)
	}
}
if(plot.all){
	if(intercept){which<- 1:(p+1)}else{which<- 2:(p+1)}
	}


if(is.null(ylab)){NamnX<-rownames(x$coef)[which]
	}else{
		if(length(as.vector(ylab))==length(rownames(x$coef)[which])){
	NamnX=ylab	}else{
		stop('Wrong dimensions on ylab.')
	}}
	
ci<-x$ciols[which,]
ui<-x$ui[which,]
coefols<-(ci[,1]+ci[,2])/2
pv<-coef(summary(x$out.model))[which,4]
sign<-ifelse(pv>0.1,"",ifelse(pv>0.05&pv<=0.1,".",ifelse(pv>0.01&pv<=0.05,"*",ifelse(pv>0.001&pv<=0.01,"**",'***'))))

n=length(which)
locCI<-(n:1)
locUI<-locCI-0.1
locCoef<-locCI+0.15

cexM<-0.6+4/p
mar=c(2+cexM,1+cexM*max(nchar(NamnX,type="width"))/2,2+cexM,1.5)

par(mar=mar,mgp = c(3, cexM*0.8, 0))
plot(ui[1,],c(locUI[1],locUI[1]),type="l",xlim=c(min(ui),max(ui)),ylim=c(min(locUI)-0.1,max(locCoef)+0.1),col=col[2],lwd=2*cexM,cex.main=cexM,cex.axis=cexM,xlab="",ylab="",yaxt='n',main='')

for(i in 2:n){
	lines(ui[i,],c(locUI[i],locUI[i]),col=col[2],lwd=2*cexM)
}
for(i in 1:n){
	lines(ci[i,],c(locCI[i],locCI[i]),lwd=2*cexM,col=col[1])
	points(coefols[i],locCI[i],pch=19,cex=0.6*cexM,col=col[1])
	text(coefols[i],locCoef[i],labels=paste(as.character(round(coefols[i],2)),sign[i]),cex=0.8*cexM,col=col[1])
}
arrows(0,c(min(locUI)-0.2,max(locCoef)+0.2),0,0,length=0,lwd=0.7*cexM)
axis(2,at=locCI,labels=NamnX,las=1,cex.axis=cexM)
}
