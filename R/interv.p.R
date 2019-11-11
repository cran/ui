#' Print interval in parantesis
#'
#' This function allows you to print an interval (vector of two elements) in a parantesis single element.
#' @param v Lower and upper bounds.
#' @param digits Number of decimals.
#'
#'
#' @export


interv.p<-function(v,digits=3){
paste('(',round(v[1],digits),', ',round(v[2],digits),')',sep='')
}
