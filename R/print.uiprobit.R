#' Prints objects of class uiprobit
#'
#' @param x an objects returned from \code{\link{ui.probit}}
#' @param digits number of digits to be printed.
#' @param digitsci number of digits to be printed in the confidence interval.
#' @param digitsui number of digits to be printed in the uncertainty interval.
#' @param ... Additional arguments, use is discouraged.
#'
#' @export

print.uiprobit<-function(x,digits=3,digitsci=digits,digitsui=digits, ...){
	if(sum(is.na(x$ui))>0){
	ui<-x$uirough
	ci<-x$cirough
}else{
	ui<-x$ui
	ci<-x$ci
}
w<-which(x$rho==0)
ci<-apply(ci[,w,],1,function(x){interv.p(x,digitsci)})
ui<-apply(ui,1,function(x){interv.p(x,digitsui)})

cat("\nCall:\n", deparse(x$call), "\n\n\n", sep = "")
cat("Confidence intervals (CI) derived assuming ignorable dropout (rho=0)","\n",sep="")
cat("Uncertainty intervals (UI) derived assuming ",min(x$rho),"<=rho<=", max(x$rho),"\n","\n",sep="")
Tab<-data.frame(Est=round(x$coef[,w],digits),ci=ci,ui=ui)
print(Tab)
	if(sum(is.na(x$ui))>0){cat("\n","Note that the standard errors are not exact.","\n", "For an exact uncertainty interval choose se = TRUE instead.",sep="")}
}



