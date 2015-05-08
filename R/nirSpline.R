#' nirSpline
#' 
#' Fits an aggregated functional data model to near-infrared spectroscopy data.
#' 
#' @param wavelengths 
#' @param absorbance
#' @param concentration
#' @param norder
#'
#' @export
#' @return List of two items 
#'   \item{betaMatrix}{TODO}
#'   \item{fittedMatrix}{TODO}
#'
#' @examples
#' require(fda)
#' data(tecator)
#' with(tecator, nirSpline(seq(850, 1050, length = 100), training[,1:100], 
#'                         training[,123:125]/100))
#' 
#' @references
#'  Ramsay, J. O. (2006) Functional Data Analysis. New York: Springer.
#'  
#'  Dias, R., Garcia, N., Ludwig, G. and Saraiva, M. (2014) Aggregated functional 
#'  data model for Near-Infrared Spectroscopy calibration and prediction, 
#'  Journal of Applied Statistics
#'
#' @seealso \code{\link{tecator}}
#' @keywords Aggregated Functional Data Analysis
nirSpline <- function(wavelengths, absorbance, concentration, ...)
UseMethod("nirSpline")