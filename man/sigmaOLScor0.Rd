% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigmaOLScor0.R
\name{sigmaOLScor0}
\alias{sigmaOLScor0}
\title{Correction of OLS sigma for causal effects}
\usage{
sigmaOLScor0(X, sigmaOLS, n, p, u, gridrho)
}
\arguments{
\item{X}{Covariate matrix outcome.}

\item{sigmaOLS}{Residual sd from outcome regression.}

\item{n}{Number of complete cases.}

\item{p}{Number of covariates outcome regression.}

\item{u}{Fitted values from propensity score regression.}

\item{gridrho}{Values of rho.}
}
\description{
This function is a bias correction of the residual standard deviation under MNAR, for \code{\link{ui.causal}}.
}
