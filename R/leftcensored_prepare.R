#' Prepare data for lc_linear
#' 
#' This function is meant for analysis of chemical concentration data as function of year. These data are given as
#' two columns in the data: Concentration and LOQ flag. In our case, we assume that an empty LOQ flag 
#' indicates that the concentration is above the LOQ. A non-empty LOQ flag indicates that the concentration is below the LOQ,
#' and the value given in the Concentration column is actually the LOQ. This procedure organises the data so the fit for lc_linear,
#' and applies the standard column names needed. In addition, the LOQ is not given when the data is above LOQ. 
#' This procedure also sets a 'mock' LOQ for uncensored data (lowest value - 1; this does not affect the analysis).         
#' 
#' @param data Data set (a data frame)
#' @param x Variable name for the predictor variable  
#' @param y Variable name of the response variable (for chemical data, concentration data)  
#' @param censored Variable name of the variable with the LOQ flag
#' @param censored_value Value of LOQ flag variable for values under LOQ, by default "<" 
#' @param log Value (TRUE/FALSE) for whether to log-transform concentrations or not (TRUE by default)
#' @param const Constant to add before log-transformation (to avoid problems with log of zero)
#' 
#' 
#' @import dplyr
#' 
#' @examples 
#' 
#' # The 'polybrom' data set contains the variables 'year', 'concentration' (concentration of a 
#' # polybrominated substance), and 'LOQ_flag' ('<' for measurements under LOQ, the Limit Of Quantification).
#' # The measurements are from three different stations. See help("polybrom").
#' 
#' # Show the data in the original dataset:
#' library(ggplot2)
#' ggplot(polybrom, aes(year, concentration, color = LOQ_flag)) + 
#'   geom_point() + 
#'   facet_wrap(vars(station)) 
#' 
#' # Make the data ready for analysis
#' # We also choose to log-transform the data in this case 
#' polybrom_for_analysis <- lc_prepare(polybrom, 
#'                                     x = "year", 
#'                                     y = "concentration", 
#'                                     censored = "LOQ_flag",
#'                                     log = TRUE)
#'                                     
#' # x, y and LOQ variables are given standard names (other variables are left unchanged). 
#' # 'y' is set to NA for censored data:  
#' ggplot(polybrom_for_analysis, aes(x, y, color = factor(uncensored))) + geom_point() + facet_wrap(vars(station))
#' 
#' # 'threshold' is set to a random value (lower than observed values) for all uncensored data:
#' ggplot(polybrom_for_analysis, aes(x, threshold, color = factor(uncensored))) + geom_point() + facet_wrap(vars(station))
#' 
#' @export
lc_prepare <- function(data,
                       x = "x",
                       y = "y",
                       censored = "LOQ_flag",
                       censored_value = "<",
                       log = FALSE,
                       const = 0){
  # Change names in data set
  varnames_user <- c(x, y, censored)
  varnames_new <- c("x", "y", "Flag")
  for (i in 1:3){
    sel <- names(data) %in% varnames_user[i]
    if (sum(sel)==1){
      names(data)[sel] <- varnames_new[i]
    } else {
      stop(paste("No unique variable named", sQuote(varnames_user[i]), "was found in the data set"))
    } 
  }
  result <- subset(data, !is.na(y))
  result$threshold <- as.numeric(NA)
  result$uncensored <- as.numeric(NA)

  if (log){
    result$y <- log(result$y + const)
  }
  check <- sum(is.na(result$y))
  if (check > 0){
    warning("Log transformation resulted in ", check,  "NA values. These are deleted from the data.")
    result <- subset(data, !is.na(y))
  }
  
  
  #
  # Values under LOQ
  #
  sel_cens <- result$Flag %in% censored_value
  result$threshold[sel_cens] <- result$y[sel_cens]
  result$y[sel_cens] <- NA
  result$uncensored[sel_cens] <- 0
  #
  # Values over LOQ
  #
  result$threshold[!sel_cens] <- min(result$y[!sel_cens]) - 1
  result$uncensored[!sel_cens] <- 1
  
  # Returning resulting data frame
  result
}
