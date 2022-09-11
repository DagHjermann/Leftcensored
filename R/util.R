
# Rename variables, with checking  
rename_check <- function(data, old, new){
  sel <- names(data) %in% old
  if (sum(sel) == 0)
    stop("Could not find variable ", sQuote(old))
  if (sum(sel) >= 2)
    stop("Found ", sum(sel), " variables named ", sQuote(old))
  names(data)[sel] <- new
  data
}


# Get ordered data
get_ordered_data1 <- function(data,
                             x = "x", 
                             y = "y", 
                             uncensored = "uncensored",
                             threshold = "threshold",
                             measurement_error = NULL){
  
  data <- as.data.frame(data)
  
  # browser()
  data <- rename_check(data, x, "x")
  data <- rename_check(data, y, "y")
  data <- rename_check(data, uncensored, "uncensored")
  data <- rename_check(data, threshold, "cut")
  
  if (!is.null(measurement_error))
    data <- rename_check(data, measurement_error, "meas_error")
  
  # If uncensored = FALSE/TRUE, change it to 0/1 
  data$uncensored <- as.numeric(data$uncensored)
  
  # Keep only data where uncensored is 0 or 1
  n1 <- nrow(data)
  data <- data %>% 
    filter(uncensored %in% c(0,1))
  n2 <- nrow(data)
  if (n2 < n1)
    warning("Note: ", n1-n2, " rows deleted because they had no accepted value for 'uncensored' (", sQuote(uncensored), ")")
  
  if (sum(data$uncensored == 0) > 0){
    
    # Set all censored data to NA (if not already done)
    data$y[data$uncensored == 0] <- NA
    
    # Order file with uncensored data first
    result <- bind_rows(
      data[data$uncensored == 1,],
      data[data$uncensored == 0,]
    )
    # y_comb is the combination of y and cut
    # - will be used as the response in the analysis
    # - will only affect the uncensored values
    result$y_comb <- c(data[data$uncensored == 1, "y"], 
                       data[data$uncensored == 0, "cut"])
    
  } else {
    
    # If no uncensored data
    result <- data
    result$y_comb <- data$y
    
  }
  
  result
  
}

get_ordered_data2 <- function(data_ordered, added_x_values){
  
  #
  # Make x data that will be used for making the fitted line + conf.int.
  # - to be added to the rest of the data
  # - we also make a y (could perhaps be NA?)
  #
  # Make x
  if (length(added_x_values) == 1){
    dat_for_fit <- data.frame(
      x = seq(min(data_ordered$x), max(data_ordered$x), length = added_x_values)
    )
  } else if (length(added_x_values) > 1){
    dat_for_fit <- data.frame(
      x = added_x_values
    )
  } else {
    stop("'new_x' must be supplied")
  }
  # Make y
  # Use ordinary gam (on uncensored data) to add the y
  # (Reason: perhaps the y value affects some sort of normalization)
  dat_for_fit$y_comb <- predict(
    mgcv::gam(y_comb~x, data = data_ordered),
    newdata = dat_for_fit)
  
  #
  # Add "x data for fitted line" to the "actual" data
  #
  # dat_ordered1 = data file
  # dat_ordered2 = with addition for daa for fitted line (3o points)
  
  list(
    data = bind_rows(data_ordered, dat_for_fit),
    data_for_fit = dat_for_fit
  )
  
}

get_dic <- function(runjags_object){
  
  #
  # DIC
  #
  dic.pd <- rjags::dic.samples(model = runjags::as.jags(runjags_object), n.iter=1000, type="pD")
  
  # Not used now:
  # dic.popt <- dic.samples(model=model_run, n.iter=30000, type="popt"); dic.popt
  
  # Select the observations for which we got penalties
  dic.sel.pd <- !is.nan(dic.pd$penalty )
  
  # Get penalties and deviances for those
  pd <- dic.pd$penalty[dic.sel.pd]
  deviance <- dic.pd$deviance[dic.sel.pd]
  
  # Calculate DIC
  dic <- sum(deviance) + sum(pd)
  
  list(
    deviance = deviance,
    pd = pd,
    dic = dic)

}




  