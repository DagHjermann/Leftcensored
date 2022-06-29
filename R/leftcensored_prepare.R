#' Prepare data for lm_linear
#' 
#' This function is meant for analysis of chemical concentration data as function of year. These data are given as
#' two columns in the data: Concentration and LOQ flag. In our case, we assume that an empty LOQ flag 
#' indicates that the concentration is above the LOQ. A non-empty LOQ flag indicates that the concentration is below the LOQ,
#' and the value given in the Concentration column is actually the LOQ. This procedure organises the data so the fit for lm_linear.
#' In addition, the LOQ is not given when the data is above LOQ. This procedure also sets the LOQ for every observation, in the following way:
#' For years with at least one observation below LOD, we use the lowest LOD. For years with no observations below LOD, LOD is set to be 
#' the lowest value - 10%.         
#' 
#' @param data Data set (a data frame)
#' @param var_year Variable name of the variable with the year data
#' @param var_concentration Variable name of the variable with the concentration data
#' @param var_LOQflag Variable name of the variable with the LOQ flag
#' @param value_under_LOQ Value of LOQ flag variable for values under LOQ 
#' @param const Constant to add before log-transformation (to avoid problems with log of zero)
#' 
#' 
#' @import dplyr
#' 
#' @export
lc_prepare <- function(data,
                       var_year = "year",
                       var_concentration = "concentration",
                       var_LOQflag = "concentr_flag",
                       value_under_LOQ = "<",
                       const = 0){
  # Change names in data set
  varnames_user <- c(var_year, var_concentration, var_LOQflag)
  varnames_new <- c("x", "y", "Flag")
  for (i in 1:3){
    sel <- names(data) %in% varnames_user[i]
    if (sum(sel)==1){
      names(data)[sel] <- varnames_new[i]
    } else {
      stop(paste("No unique variable named", sQuote(varnames_user[i]), "was found in the data set"))
    } 
  }
  result <- data %>%
    filter(!is.na(y)) %>%
    mutate(
      y = log(y + const),
      threshold = as.numeric(NA),
      uncensored = as.numeric(NA))
  #
  # Values under LOQ
  #
  sel_cens <- result$Flag %in% value_under_LOQ
  result$threshold[sel_cens] <- result$y[sel_cens]
  result$y[sel_cens] <- NA
  result$uncensored[sel_cens] <- 0
  #
  # Values over LOQ
  #
  result$threshold[!sel_cens] <- min(result$y[!sel_cens]) - 10
  result$uncensored[!sel_cens] <- 1
  result
}
