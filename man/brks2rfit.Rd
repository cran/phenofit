% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roughFit.R
\name{brks2rfit}
\alias{brks2rfit}
\title{get rough fitting}
\usage{
brks2rfit(brks)
}
\arguments{
\item{brks}{returned by function \code{\link[=season_mov]{season_mov()}}}
}
\value{
\itemize{
\item \code{data}:
\itemize{
\item t
\item y
\item QC_flag
}
\item \code{tout}:
\item \code{zs}: list of iter1, ..., itern
\item \code{ws}: list of iter1, ..., itern
}
}
\description{
get rough fitting
}
\keyword{internal}
