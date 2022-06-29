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
                       const = 0){
  # Change names in data set
  varnames_user <- c(var_year, var_concentration, var_LOQflag)
  varnames_new <- c("Year", "y", "Flag")
  for (i in 1:3){
    sel <- names(data) %in% varnames_user[i]
    if (sum(sel)==1){
      names(data)[sel] <- varnames_new[i]
    } else {
      stop(paste("No unique variable named", sQuote(varnames_user[i]), "was found in the data set"))
    } 
  }
  data_1 <- data %>%
    filter(!is.na(y)) %>%
    mutate(log_y = log(y + const)) %>% 
    group_by(Year) %>%
    mutate(
      n_below_loq = sum(!is.na(Flag))
    ) %>%
    ungroup()
  # Years with some data below LOQ
  data_2a <- data_1 %>%
    filter(n_below_loq > 0) %>%
    group_by(Year, n_below_loq) %>%
    mutate(
      LOQ_per_year = min(log_y[!is.na(Flag)])  # years with at least one below LOD: use lowest LOD 
    ) %>%
    ungroup()
  # Then we also must make sure that any data that are above LOD has LOQ_per_year lower than this
  sel <- with(data_2a, is.na(Flag) & LOQ_per_year > log_y)
  data_2a$LOQ_per_year[sel] <- data_2a$log_y[sel] - log(1.10)
  # Years with no data below LOQ
  data_2b <- data_1 %>%
    filter(n_below_loq == 0) %>%
    group_by(Year) %>%
    mutate(
      LOQ_per_year = min(log_y) - log(1.10)   # years with no obs below LOD: use lowest value - 10%
    ) %>%
    ungroup()
  result <- bind_rows(data_2a, data_2b) %>%
    mutate(
      x = Year,
      y = ifelse(is.na(Flag), log_y, 0),
      uncensored = ifelse(is.na(Flag), 1, 0),
      threshold = ifelse(uncensored == 0, log_y, LOQ_per_year)
    ) %>%
    select(x, y, Flag, uncensored, threshold)
  result
}
