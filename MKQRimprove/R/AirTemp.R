#' @title AirTemp
#' @description A dataset containing the air temperature extremum of 1850 - 2020 monthly.
#'  
#'
#' @format A data frame with 107 rows and 3 variables
#'\describe{
#'   the first column: year
#'   the second column: month
#'   the third column: Air temperature extremum monthly
#' }
#' @source Berkeley
#' @docType data
#' @keywords AirTemp
#' @usage data("AirTemp")
#'
#' @references
#' \itemize{
#' \item Berkeley(2011)
#'   http://berkeleyearth.lbl.gov/auto/Global
#'   }
#'
#' @examples
#' \dontrun{
#' data("AirTemp")
#' y <- AirTemp[,3]
#' x <- seq(1850,2020,0.0833) #transform into time series
#' tempmonth <- as.data.frame(cbind(x,y)) #use tempmonth to do real data analysis
#' }

