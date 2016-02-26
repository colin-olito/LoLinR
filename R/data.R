#' O2 consumption data for sea urchins
#'
#' A dataset containing VO2 time series for 4 indivudal sea urchins.
#' Data are provided by Colin Olito (unpublished data).
#'
#' @format A data frame with 166 rows and 5 variables:
#' \itemize{
#'		\item time: time in minutes
#'		\item A: Volume O2 (mL) data for urchin A
#'		\item B: Volume O2 (mL) data for urchin B
#'		\item C: Volume O2 (mL) data for urchin C
#'		\item D: Volume O2 (mL) data for urchin D
#' }
#' @name TestO2data
#' @docType data
NULL

#' CO2 production in a speckled cockroach
#'
#' A dataset containing CO2 production data by a speckled cockroach. This 
#' dataset represents a small example of data used in analyses described in Schimpf
#' et al. 2011. Cockroaches that exchange respiratory gases discontinuously
#' survive food and water restriction. Evolution 66: 597--604.
#'
#' @format A data frame with 392 rows and 2 variables:
#' \itemize{
#'		\item Time: time in seconds
#'		\item CO2.ppm: Rate of CO2 concentration in ppm
#' }
#'
#' @name thinned_cockroach_data
#' @docType data
NULL

#' O2 consumption data -- flow-through respirometry for cormorants
#'
#' A dataset containing O2 consumption data for a single Great cormorant. This
#' dataset represents a small example of data used for analyses described in White
#' et al. 2011. Metabolic rate throughout the annual cycle reveals the demands
#' of an Arctic existence in Great Cormorants. Ecology 92: 475--486.
#'
#' @format A data frame with 447 rows and 2 variables:
#' \itemize{
#'		\item Time: time in hours
#'		\item VO2.ml.min: Rate of oxygen consumption (mL O2/min)
#' }
#' @name thinned_cormorant_data
#' @docType data
NULL

#' O2 consumption in a cane toad
#'
#' A dataset containing O2 consumption data for a cane toad. This dataset
#' represents a small example of data used for analyses described in Halsey and
#' White. 2010. Measuring energetics and behaviour using accelerometry in cane 
#' toads Bufo marinus. PLoS One 5: e10170.
#'
#' @format A data frame with ... rows and ... variables:
#' @format A data frame with 442 rows and 2 variables:
#' \itemize{
#'		\item Time: time in seconds
#'		\item FO2: fractional O2 concentration 
#' }
#' @name thinned_toad_data
#' @docType data
NULL
