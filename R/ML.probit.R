#' Fit maximum likelihood for fixed values of rho
#'
#'This is a support function for \code{\link{ui.probit}}
#' @param out.formula Formula for outcome regression.
#' @param mis.formula Formula for regression model for the missingness mechanism.
#' @param data Data frame containing the variables in the formulas
#' @param rho Vector containing the values of rho for which we want to fit the likelihood.
#' @param progress If TRUE prints out process time for each maximazation of the likelihood.
#' @param method Maximazation method to be passed through \code{maxLik}
#' @importFrom maxLik maxLik
#' @importFrom stats get_all_vars glm binomial update.formula vcov setNames as.formula
#' @export
ML.probit<-function(out.formula,mis.formula=NULL,data,rho=c(-0.5,0.5),progress=TRUE,method='NR'){

  ###Warnings
   if(class(data)!="data.frame"){
	stop('Data must be a data frame')
    }
    
  y<-get_all_vars(out.formula,data=data)[,1]
  z <- as.numeric(!is.na(y))
  
  data$z<-z
  out.model<-glm(out.formula,family=binomial(link="probit"),data=data)  
  
  m.f<-update.formula(out.formula, as.formula(paste('z', "~ ."))) 
	#Need to keep all NA in X.y, if different covariates in mis.model and out.model 	make new model without missing in y
  if(is.null(mis.formula)){
  	mis.formula<-m.f
	mis.model<-glm(mis.formula,family=binomial(link="probit"),data=data)
	t.m<-mis.model
    }else{
	mis.model<-glm(mis.formula,family=binomial(link="probit"),data=data)
	t.m<-glm(m.f,family=binomial(link="probit"),data=data)
}

  X.y<-model.matrix(t.m)
  X.z<-model.matrix(mis.model)
  
  nrho <- length(rho)
  p1 <- length(coef(out.model))
  p2 <- length(coef(mis.model))
  d<- p1+p2

  coef <- matrix(nrow=(p1 + p2), ncol = nrho)
  rownames(coef) <- c(names(out.model$coef), names(mis.model$coef))
  colnames(coef)<-paste(rho)

  value <- vector(length = nrho)
  names(value) <- paste(rho)
  max.type <- value

  i0 <- which(rho == 0)
  coef[, i0] <- c(coef(out.model), coef(mis.model))

  max.info <- list(converged = array(dim=nrho), message = array(dim=nrho), method = array(dim=nrho),
                   iterations = array(dim=nrho))
  max.info$converged[i0] <- mis.model$converged == TRUE & out.model$converged == TRUE
  max.info$message[i0] <- " "
  max.info$method[i0] <- "glm"
  max.info$iterations[i0] <- NA

  
  glmVcov <- vcov(out.model)
  mis.glmVcov <- vcov(mis.model)

  value[i0] <- LogL.probit(coef[, i0], rho = 0, X.z = X.z, X.y = X.y, y = y, z = z)

  sigma <- list()
  sigma[[i0]] <- matrix(0, nrow = d, ncol = d)
  sigma[[i0]][1:p1, 1:p1] <- glmVcov
  sigma[[i0]][(p1 + 1):d, (p1 + 1):d] <- mis.glmVcov
  sigma[[i0]][(p1 + 1), (p1 + 1)] <- summary(out.model)$dispersion^2/(2*out.model$df.residual)
  dimnames.sigma <- c(dimnames(glmVcov)[[1]], dimnames(mis.glmVcov)[[1]])
  dimnames(sigma[[i0]]) <- list(dimnames.sigma, dimnames.sigma)

  if(sum(rho < 0) > 0){
    for(i in 1:(i0 - 1)){

      if(progress == TRUE){
        cat("Optimization for rho =", rho[i0 - i], sep=" ", "\n")
        ptm <- proc.time()
      }

      f <- function(par)
        LogL.probit(par, rho = rho[i0 - i], X.z = X.z, X.y = X.y, y = y, z=z)

      g <- function(par)
        grr(par, rho = rho[i0 - i], X.z = X.z, X.y = X.y, y = y, z=z)

      h <- function(par)
        hess(par, rho = rho[i0 - i], X.z = X.z, X.y = X.y, y = y, z=z)

      ML <- maxLik(logLik = f, grad = g, hess = h, start = coef[, i0 - i + 1], method = method)

      max.info$converged[i0 - i] <- ML$code <= 2
      max.info$message[i0 - i] <- ML$message
      max.info$method[i0 - i] <- ML$type
      max.info$iterations[i0 - i] <- ML$iterations
      coef[,i0 - i] <- ML$estimate
      value[i0 - i] <- ML$maximum

      sigma[[i0 - i]] <- -solve(ML$hessian)


      if(progress == TRUE)
        cat("   Time elapsed:", (proc.time()-ptm)[3], "s", "\n")
    }
  }

  if(sum(rho > 0) > 0){
    for(i in (i0 + 1):nrho){

      if(progress == TRUE){
        cat("Optimization for rho =", rho[i], sep=" ", "\n")
        ptm <- proc.time()
      }

      f <- function(par)
        LogL.probit(par, rho = rho[i], X.z = X.z, X.y = X.y, y = y, z = z)

      g <- function(par)
        grr(par, rho = rho[i], X.z = X.z, X.y = X.y, y = y , z = z)

      h <- function(par)
        hess(par, rho = rho[i], X.z = X.z, X.y = X.y, y = y, z = z)

      ML <- maxLik(logLik = f, grad = g, hess = h, start = coef[, i - 1], method = method)

      max.info$converged[i] <- ML$code <= 2
      max.info$message[i] <- ML$message
      max.info$method[i] <- ML$type
      max.info$iterations[i] <- ML$iterations
      coef[,i] <- ML$estimate
      value[i]  <-  ML$maximum

      sigma[[i]] <- -solve(ML$hessian)

      if(progress == TRUE)
        cat("   Time elapsed:", (proc.time()-ptm)[3], "s", "\n")
    }
  }

  mis.coef  <- as.matrix(coef[(p1 + 1):(p1 + p2), ])
  if(p2 == 1){
    mis.coef  <- t(mis.coef)
    rownames(mis.coef) <- names(mis.model$coefficients)
  }

  colnames(mis.coef) <- paste(rho)
  coef <- as.matrix(coef[1:(p1), ])
  if(p1 == 1){
    coef <- t(coef)
    rownames(coef) <- names(out.model$coefficients)
  }
  colnames(coef) <- paste(rho)

  max.info <- lapply(max.info, setNames, nm = paste(rho))
  ML.object<-list(coef = coef, rho = rho, mis.coef = mis.coef, mis.model = mis.model, out.model = out.model, X.z = X.z,
              X.y = X.y, y = y, z = z, value = value, vcov = sigma, max.info = max.info)
  return(ML.object)

}

