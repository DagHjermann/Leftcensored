
set_x_to_censored <- function(data, x_values){
  sel <- data$x %in% x_values
  data$threshold[sel] <- data$y[sel]
  data$y[sel] <- NA
  data$uncensored[sel] <- 0
  data
}

