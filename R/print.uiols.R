#' Prints objects of class uiols
#'
#' @param x an objects returned from \code{\link{ui.ols}}
#' @param digits number of digits to be printed.
#' @param digitsci number of digits to be printed in the confidence interval.
#' @param digitsui number of digits to be printed in the uncertainty interval.
#' @param ... Additional arguments, use is discouraged.
#'
#' @export

print.uiols<-function(x,digits=3,digitsci=digits,digitsui=digits, ...){

w<-which(x$gridrho==0)
ci<-apply(x$ci[,w,],1,function(x){interv.p(x,digitsci)})
ui<-apply(x$ui,1,function(x){interv.p(x,digitsui)})

cat("\nCall:\n", deparse(x$call), "\n\n\n", sep = "")
cat("Confidence intervals (CI) derived assuming ignorable dropout (rho=0)","\n",sep="")
cat("Uncertainty intervals (UI) derived assuming ",min(x$rho),"<=rho<=", max(x$rho),"\n","\n",sep="")
Tab<-data.frame(Est=round(x$coef[,w],digits),ci=ci,ui=ui)
print(Tab)
	if(sum(is.na(x$ui))>0){cat("\n","Note that the standard errors are not exact.","\n", "For an exact uncertainty interval choose se = TRUE instead.",sep="")}
}



