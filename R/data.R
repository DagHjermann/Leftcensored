#' Concentration of a polybrominated substance in fish    
#' 
#' A data set containing concentration of a polybrominated substance in fish captured in three 
#' locations (column 'stations') over several years. 
#' 
#' @format A data frame (in tibble format) with 699 rows and 4 variables:  
#' \describe{
#'   \item{station}{A code distinguishing the three stations}
#'   \item{year}{Year of capture}
#'   \item{concentration}{Concentration of the substance (when concentr_flag = NA), or the 
#'   detection limit (when concentr_flag = '<')}
#'   \item{LOQ_flag}{A character column indicating whether the value in concentration is the actual
#'   concentration (LOQ_flag = NA) or the upper limit for concentration (LOQ_flag = '<')}
#' }   
#' @source \url{https://ocean.ices.dk/ohat/}
"polybrom"