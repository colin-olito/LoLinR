#' O2 consumption data for sea urchins
#'
#' A dataset containing VO2 time series for 4 individual sea urchins.
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
#' @name UrchinData
#' @docType data
NULL

#' O2 consumption data for Bugula neritina
#'
#' A dataset containing O2 saturation time series for 4 individual 
#' Bugula neritina larvae. These data representa  small subset of
#' the data used in analyses described in Pettersen A.K., White C.R., 
#' Marshall D.J. 2015. Why does offspring size affect performance? 
#' Integrating metabolic scaling with life-history theory. Proc. Roy.
#' Soc. B 282:20151946. (doi:10.1098/rspb.2015.1946)
#'
#' @format A data frame with 107 rows and 5 variables:
#' \itemize{
#'		\item Time.s: time in seconds
#'		\item A2: O2 saturation (%) data for individual larvae in chamber A2
#'		\item A3: O2 saturation (%) data for individual larvae in chamber A3
#'		\item C2: O2 saturation (%) data for individual larvae in chamber C2
#'		\item D3: O2 saturation (%) data for individual larvae in chamber D3
#' }
#' @name BugulaData
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
#' @name CormorantData
#' @docType data
NULL

#' O2 production data for pelagic phytoplankton in ponds
#'
#' A dataset containing O2 production time series for 
#' pelagic phytoplankton in several ponds and light conditions.
#' These data represent a small example of the data used for analyses 
#' described in Yvon-Durocher et al. (2015). Five years of experimental
#' warming increases the biodiversity and productivity of phytoplankton.
#' PLoS Biology 13(12): e1002324. doi:10.1371/journal. pbio.1002324.
#'
#' @format A data frame with 242 rows and 5 variables:
#' \itemize{
#'		\item pond: pond ID
#'		\item month.num: Sampling month ID
#'		\item flux: light conditions
#'		\item Time.m: Time in minutes
#'		\item O2: O2 (mg/L) data for phytoplankton
#' }
#' @name PelagicData
#' @docType data
NULL