#' Prepare data for leftcensored_lm
#' 
#' This function is meant for analysis of chemical concentration data as function of year. These data are given as
#' two columns in the data: Concentration and LOQ flag. In our case, we assume that an empty LOQ flag 
#' indicates that the concentration is above the LOQ. A non-empty LOQ flag indicates that the concentration is below the LOQ,
#' and the value given in the Concentration column is actually the LOQ. This procedure organises the data so the fit for leftcensored_lm.
#' In addition, the LOQ is not given when the data is above LOQ. This procedure also sets the LOQ for every observation, in the following way:
#' For years with at least one observation below LOD, we use the lowest LOD. For years with no observations below LOD, LOD is set to be 
#' the lowest value - 10%.         
#' 
#' @param data Data set (a data frame)
#' @param var_year Variable name of the variable with the year data
#' @param var_concentration Variable name of the variable with the concentration data
#' @param var_LOQflag Variable name of the variable with the LOQ flag
#' 
#' @import dplyr
#' 
#' @export
leftcensored_prepare <- function(data,
                                 var_year = "Year",
                                 var_concentration = "Conc",
                                 var_LOQflag = "Flag"){
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
  data %>%
    filter(!is.na(y)) %>%
    mutate(y = log(y + 0.00001)) %>% 
    group_by(Year) %>%
    mutate(
      n_below_loq = sum(!is.na(Flag)),
      LOQ_per_year = case_when(
        n_below_loq > 0 ~ min(y[!is.na(Flag)]),  # years with at least one below LOD: use lowest LOD 
        n_below_loq == 0 ~ min(y) - log(1.10)    # years with no obs below LOD: use lowest value - 10%
      )
    ) %>%
    ungroup() %>% 
    mutate(
      x = Year,
      y_aboveLOQ = ifelse(is.na(Flag), 1, 0),
      y_cens = y,   # ifelse(y_aboveLOD, y, NA)
      y_LOQ = ifelse(y_aboveLOQ == 0, y, LOQ_per_year)
    ) %>%
    select(x, y, Flag, y_aboveLOQ, y_cens, y_LOQ)
}
# prepare_data(df_test)
