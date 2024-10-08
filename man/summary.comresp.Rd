% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resp.R
\name{summary.comresp}
\alias{summary.comresp}
\title{Summary of the \code{comresp} object}
\usage{
\method{summary}{comresp}(object, ...)
}
\arguments{
\item{object}{An object of class \code{comresp}}

\item{...}{Currently not used.}
}
\value{
The function returns a list with the variance components (\code{varcomp})
and a data frame informing which genotypes had the highest DGE, the highest and
the lowest IGE, and the highest TGV. This data frame also informs the number of
aggressive, homoeostatic and sensitive. If a multi-age model was fitted, the output
list will have another data frame with the same information for each age.
}
\description{
A brief summary of the \code{comresp} object.
}
\seealso{
\link{resp}
}
