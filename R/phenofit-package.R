#' @title phenofit
#' @name phenofit
#' @aliases phenofit-package
#' @docType package
#' @keywords Vegetation phenology package
#' @description Vegetation phenology package
#' 
#' @import magrittr numDeriv
#' @import ggplot2 
#' @importFrom dplyr nth mutate group_by group_map group_modify select
#' ungroup across starts_with as_tibble
#' @importFrom gridExtra arrangeGrob
#' @importFrom data.table data.table as.data.table := is.data.table fwrite fread 
#' dcast
#' @importFrom zoo na.approx index zoo
#' @importFrom purrr map map_df
#' @importFrom lubridate ymd yday year month day dyears is.Date
#' @importFrom utils object.size
#' @importFrom grDevices dev.off cairo_pdf colorRampPalette
#' @import stats graphics
#' @import zeallot
#' @useDynLib phenofit, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

.onLoad <- function (libname, pkgname){
    season_input <- c("rFUN", "wFUN", "iters", "wmin", "lambda", ".lambda_vcurve", 
        "frame", "nf", "minpeakdistance", "ypeak_min", "r_max", "r_min", "rtrough_max", 
        "MaxPeaksPerYear", "MaxTroughsPerYear", "calendarYear", "adj.param", 
        "rm.closed", "is.continuous", ".check_season", "maxExtendMonth", 
        "nextend", "len_min", "len_max", "verbose")
    
    if(getRversion() >= "2.15.1") {
        utils::globalVariables(
            c(".", ".SD", ".N", "..vars", 
              "left", "len", "right", "y_peak",
              "meth", "doy", "origin", # tidy_pheno
              "DayOfYear", "SummaryQA", "site", "EVI", "w", "QC_flag", # tidy_MOD13
              "beg", "end",  # plot_curvefits
              "val", "type", "flag", "peak", # season, 
              "i", "qc", "y", "sitename",  # phenofit_TS.avhrr
              season_input, 
              "years", "nyear", "d_fit", "rfit", "info_peak"
            )
        )
    }
}

# .onUnload <- function (libpath) {
#   library.dynam.unload("phenofit", libpath)
# }

# dplyr's version now returns last row instead of last col when
# supplied a data frame
last <- function(x) {
  x[[length(x)]]
}

first <- function(x) {
  x[[1]]
}
