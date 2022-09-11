#' Clean data before running trend analyses, using Rob's rule 1
#' 
#' Rule 1. Time series should be truncated from the left until Nplus/N >= 0.5.
#'
#' @param data 
#'
#' @return A data set with the same columns, but might have some rows deleted
#' @export
#'
#' @examples
#' 
#' # Simulate some data
#' testdata <- data.frame(
#'   x = rep(2009:2020, each = 3),
#'   threshold = NA,
#'   uncensored = 1)
#' 
#' # Set the  data for the first 7 years to be left-censored   
#' sel <- testdata$x %in% c(2010:2016)
#' testdata$threshold[sel] <- testdata$y[sel]
#' testdata$y[sel] <- NA
#' testdata$uncensored[sel] <- 0
#' 
#' testdata_cleaned <- lc_clean1(testdata)

lc_clean1 <- function(data){
  
  # Data must have variables named 'x' and 'uncensored'
  
  data_summ <- data %>%
    group_by(x) %>%
    summarise(N_over = sum(uncensored %in% 1),
              N_under = sum(uncensored %in% 0),
              N_text = paste0(N_over, "/", N_under), .groups = "drop") %>%
    as.data.frame()
  
  # 'Nplus_sum' for row i is the number of x values with uncensored data 
  # from row i to the end
  # 'N_sum' for row i is the number of x values from row i to the end 
  # 'Nplus_N_ratio' is the ratio between the two   
  data_summ$Nplus_sum <- rev(cumsum(rev(data_summ$N_over > 0)))
  data_summ$N_sum <- rev(1:nrow(data_summ))
  data_summ$Nplus_N_ratio <- with(data_summ, Nplus_sum/N_sum)
  
  N <- length(data_summ$N_over)
  Nplus <- sum(data_summ$N_over >= 1)
  
  if (Nplus/N < 0.5){
    
    x_accepted <- data_summ$x[data_summ$Nplus_N_ratio >= 0.5]
    
    if (length(x_accepted) >= 1){
      
      data_result <- subset(data, x >= x_accepted[1])
      data_result_summ <- subset(data_summ, x >= x_accepted[1])
      N <- length(data_summ$N_over)
      Nplus <- sum(data_summ$N_over > 0)
      
      message("Deleted data with x < ", x_accepted[1], " (", nrow(data) - nrow(data_result), " rows of data)\n",
              "Final data has ", length(data_result_summ$N_over), " unique x values (",
              sum(data_result_summ$N_over >= 1), " with uncensored data), and ",
              nrow(data_result), " rows of data."
      )
    } else {
      data_result <- NULL
      warning("Criterion 1 cannot be fulfilled ")
    }
  } else {
    
    data_result <- data
    message("No changes made to the data (lc_clean1)")
    
  }
  
  data_result
  
}

#' Add flag to data before running trend analyses, using Rob's rule 1
#' 
#' Rule 1. Time series should be truncated from the left until Nplus/N >= 0.5.
#'
#' @param data Must contain variables named 'x' and 'uncensored'
#'
#' @return The same data set, with added TRUE/FALSE variable "Rule1"
#' @export
#'
#' @examples
#' 
#' # Simulate some data
#' testdata <- data.frame(
#'   x = rep(2009:2020, each = 3),
#'   threshold = NA,
#'   uncensored = 1)
#' 
#' # Set the  data for the first 7 years to be left-censored   
#' sel <- testdata$x %in% c(2010:2016)
#' testdata$threshold[sel] <- testdata$y[sel]
#' testdata$y[sel] <- NA
#' testdata$uncensored[sel] <- 0
#' 
#' testdata_flagged <- lc_flag1(testdata)

lc_flag1 <- function(data, show_result = TRUE){
  
  # Data must have variables named 'x' and 'uncensored'
  
  data_summ <- data %>%
    group_by(x) %>%
    summarise(N_over = sum(uncensored %in% 1),
              N_under = sum(uncensored %in% 0),
              N_text = paste0(N_over, "/", N_under), .groups = "drop") %>%
    as.data.frame()
  
  # 'Nplus_sum' for row i is the number of x values with uncensored data 
  # from row i to the end
  # 'N_sum' for row i is the number of x values from row i to the end 
  # 'Nplus_N_ratio' is the ratio between the two   
  data_summ$Nplus_sum <- rev(cumsum(rev(data_summ$N_over > 0)))
  data_summ$N_sum <- rev(1:nrow(data_summ))
  data_summ$Nplus_N_ratio <- with(data_summ, Nplus_sum/N_sum)
  
  N <- length(data_summ$N_over)
  Nplus <- sum(data_summ$N_over >= 1)
  
  data$Rule1 <- TRUE
  data_summ$Rule1 <- TRUE
  
  if (Nplus/N < 0.5){
    
    x_accepted <- data_summ$x[data_summ$Nplus_N_ratio >= 0.5]
    
    if (length(x_accepted) >= 1){
      
      sel1 <- data$x < x_accepted[1]
      data$Rule1[sel1] <- FALSE

      if (show_result){
        sel2 <- data_summ$x < x_accepted[1]
        data_summ$Rule1[sel2] <- FALSE
        message("Flagged data (Rule1 = FALSE) with x < ", x_accepted[1], " (", sum(!data$Rule1), " rows of data)\n",
              "Unflagged data (Rule1 = TRUE) has ", sum(data_summ$Rule1), " unique x values (",
              sum(data_summ$Rule1 & data_summ$N_over >= 1), " with uncensored data), and ",
              sum(data$Rule1), " rows of data.")
      }
    } else {
      data_result <- NULL
      if (show_result)
        message("All data flagged")
    }
  } else {
    
    if (show_result)
      message("No changes made to the data (lc_clean1)")
    
  }
  
  data
  
}


#' Clean data before running trend analyses, using Rob's rule 2
#' 
#' Rule 2. If a linear/smooth trend is fitted, the first year must be non-censored 
#'
#' @param data 
#'
#' @return A data set with the same columns, but might have some rows deleted
#' @export
#'
#' @examples
#' 
#' # Simulate some data
#' testdata <- data.frame(
#'   x = rep(2009:2020, each = 3),
#'   threshold = NA,
#'   uncensored = 1)
#' 
#' # Set the  data for the first 7 years to be left-censored   
#' sel <- testdata$x %in% c(2010:2016)
#' testdata$threshold[sel] <- testdata$y[sel]
#' testdata$y[sel] <- NA
#' testdata$uncensored[sel] <- 0
#' 
#' testdata_cleaned <- lc_clean1(testdata)
#' 

lc_clean2 <- function(data){
  
  # Rule 2. If a linear/smooth trend is fitted, the first year must be non-censored 
  
  # This is typically run for linear and smooth fits (k >= 2)  
  # Data must have variables named 'x' and 'uncensored'  
  
  data_summ <- data %>%
    group_by(x) %>%
    summarise(N_over = sum(uncensored %in% 1),
              N_under = sum(uncensored %in% 0),
              N_text = paste0(N_over, "/", N_under), .groups = "drop") %>%
    as.data.frame()
  
  if (data_summ$N_over[1] == 0){
    
    x_accepted <- data_summ$x[data_summ$N_over >= 1]
    
    if (length(x_accepted) >= 1){
      data_result <- subset(data, x >= x_accepted[1])
      data_deleted <- subset(data, x < x_accepted[1])
      message("Deleted data with x < ", x_accepted[1], " (", nrow(data_deleted), " rows of data)")
    } else {
      data_result <- NULL
      warning("Criterion 2 cannot be fulfilled (no data are uncensored)")
    }
    
  } else {
    
    data_result <- data
    message("No changes made to the data (lc_clean2)")
    
  }
  
  data_result
  
}

#' Flag data before running trend analyses, using Rob's rule 2
#' 
#' Rule 2. If a linear/smooth trend is fitted, the first year must be non-censored 
#' 
#' @param data 
#'
#' @return The same data set, with added TRUE/FALSE variable "Rule1"
#' @export
#'
#' @examples
#' 
#' # Simulate some data
#' testdata <- data.frame(
#'   x = rep(2009:2020, each = 3),
#'   threshold = NA,
#'   uncensored = 1)
#' 
#' # Set the  data for the first 7 years to be left-censored   
#' sel <- testdata$x %in% c(2010:2016)
#' testdata$threshold[sel] <- testdata$y[sel]
#' testdata$y[sel] <- NA
#' testdata$uncensored[sel] <- 0
#' 
#' testdata_flagged <- lc_flag2(testdata)

lc_flag2 <- function(data, show_result = TRUE){
  
  # Rule 2. If a linear/smooth trend is fitted, the first year must be non-censored 
  
  # This is typically run for linear and smooth fits (k >= 2)  
  # Data must have variables named 'x' and 'uncensored'  
  
  data_summ <- data %>%
    group_by(x) %>%
    summarise(N_over = sum(uncensored %in% 1),
              N_under = sum(uncensored %in% 0),
              N_text = paste0(N_over, "/", N_under), .groups = "drop") %>%
    as.data.frame()
  
  data$Flag2 <- TRUE
  
  if (data_summ$N_over[1] == 0){
    
    x_accepted <- data_summ$x[data_summ$N_over >= 1]
    
    if (length(x_accepted) >= 1){
      sel <- data$x < x_accepted[1]
      data$Flag2[sel] <- FALSE
      if (show_result)
        message("Flagged data (Rule2 = FALSE() with x < ", x_accepted[1], " (", sum(!data$Flag2), " rows of data)")
    } else {
      data$Flag2 <- FALSE
      if (show_result)
        message("No uncensored values in the data")
    }
    
  } else {
    
    if (show_result)
      message("No values flagged (Rule2 set to TRUE for all rows)")
    
  }
  
  data
  
}

