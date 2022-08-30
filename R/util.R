
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
get_dat_ordered1 <- function(data,
                             x = "x", 
                             y = "y", 
                             uncensored = "uncensored",
                             threshold = "threshold"){
  
  data <- rename_check(data, x, "x")
  data <- rename_check(data, y, "y")
  data <- rename_check(data, uncensored, "uncensored")
  data <- rename_check(data, threshold, "cut")
  
  # If uncensored = FALSE/TRUE, change it to 0/1 
  data$uncensored <- as.numeric(data$uncensored)
  
  # Keep only data where uncensored is 0 or 1
  n1 <- nrow(data)
  data <- data %>% 
    filter(uncensored %in% c(0,1))
  n2 <- nrow(data)
  if (n2 < n1)
    warning("Note: ", n1-n2, " rows deleted because they had no accepted value for 'uncensored' (", sQuote(uncensored), ")")
  
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
  
  result
  
}

get_dat_ordered2 <- function(data_ordered, fit_length){
  
  #
  # Make x data that will be used for making the fitted line + conf.int.
  # - to be added to the rest of the data
  # - we also make a y (could perhaps be NA?)
  #
  # Make x
  dat_for_fit <- data.frame(
    x = seq(min(data_ordered$x), max(data_ordered$x), length = 30)
  )
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
  bind_rows(
    data_ordered,
    dat_for_fit
  )
  
}
