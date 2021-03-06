% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/se.ols.R
\name{se.ols}
\alias{se.ols}
\title{Calculation of se for OLS}
\usage{
se.ols(X, sigmaOLScor, u, gridrho)
}
\arguments{
\item{X}{Covariate matrix.}

\item{sigmaOLScor}{Output from sigmaOLScor1}

\item{u}{Fitted values from mis.model.}

\item{gridrho}{Values of rho.}
}
\description{
This function calculates the se for UI based on OLS when we have MNAR data, for \code{\link{ui.ols}}.
}
