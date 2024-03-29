% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qcFUN.R
\name{qcFUN}
\alias{qcFUN}
\alias{getBits}
\alias{qc_summary}
\alias{qc_StateQA}
\alias{qc_FparLai}
\alias{qc_5l}
\alias{qc_NDVIv4}
\alias{qc_NDVI3g}
\alias{qc_SPOT}
\title{Initial weights according to qc}
\usage{
getBits(x, start, end = start)

qc_summary(QA, wmin = 0.2, wmid = 0.5, wmax = 1)

qc_StateQA(QA, wmin = 0.2, wmid = 0.5, wmax = 1)

qc_FparLai(QA, FparLai_QC = NULL, wmin = 0.2, wmid = 0.5, wmax = 1)

qc_5l(QA, wmin = 0.2, wmid = 0.5, wmax = 1)

qc_NDVIv4(QA, wmin = 0.2, wmid = 0.5, wmax = 1)

qc_NDVI3g(QA, wmin = 0.2, wmid = 0.5, wmax = 1)

qc_SPOT(QA, wmin = 0.2, wmid = 0.5, wmax = 1)
}
\arguments{
\item{x}{Binary value}

\item{start}{Bit starting position, count from zero}

\item{end}{Bit ending position}

\item{QA}{quality control variable}

\item{wmin}{Double, minimum weigth (i.e. weight of snow, ice and cloud).}

\item{wmid}{Dougle, middle weight, i.e. marginal}

\item{wmax}{Double, maximum weight, i.e. good}

\item{FparLai_QC}{Another QC flag of \code{MCD15A3H}}
}
\value{
A list object with
\itemize{
\item \code{weigths}: Double vector, initial weights.
\item \code{QC_flag}: Factor vector, with the level of
\code{c("snow", "cloud", "shadow", "aerosol", "marginal", "good")}
}
}
\description{
\itemize{
\item \code{getBits}: Extract bitcoded QA information from bin value
\item \code{qc_summary}: Initial weigths based on Quality reliability of VI pixel,
suit for MOD13A1, MOD13A2 and MOD13Q1 (SummaryQA band).
\item \code{qc_5l}: Initial weights based on Quality control of five-level
confidence score, suit for MCD15A3H(LAI, FparLai_QC), MOD17A2H(GPP, Psn_QC)
and MOD16A2(ET, ET_QC).
\item \code{qc_StateQA}: Initial weights based on \code{StateQA}, suit for MOD09A1, MYD09A1.
\item \code{qc_FparLai}: For MODIS LAI
\item \code{qc_NDVI3g}: For AVHRR NDVI3g
\item \code{qc_NDVIv4}: For AVHRR NDVIv4
}
}
\details{
If \code{FparLai_QC} specified, \code{I_margin = SCF_QC >= 2 & SCF_QC <= 3}.
}
\note{
\code{qc_5l} and \code{qc_NDVIv4} only returns \code{weight}, without \code{QC_flag}.
}
\examples{
set.seed(100)
QA <- as.integer(runif(100, 0, 2^7))

r1 <- qc_summary(QA, wmin = 0.2, wmid = 0.5, wmax = 1)
r2 <- qc_StateQA(QA, wmin = 0.2, wmid = 0.5, wmax = 1)
r_5l <- qc_5l(QA, wmin = 0.2, wmid = 0.5, wmax = 1)
r_NDVI3g <- qc_NDVI3g(QA, wmin = 0.2, wmid = 0.5, wmax = 1)
r_NDVIv4 <- qc_NDVIv4(QA, wmin = 0.2, wmid = 0.5, wmax = 1)
}
\references{
https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13A1

https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD15A3H

Erwin Wolters, Else Swinnen, Carolien Toté, Sindy Sterckx.
SPOT-VGT COLLECTION 3 PRODUCTS USER MANUAL V1.2, 2018, P47
}
\seealso{
\code{\link[=qc_sentinel2]{qc_sentinel2()}}
}
