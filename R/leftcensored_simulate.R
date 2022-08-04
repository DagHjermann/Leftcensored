#' Simulate censored data
#' 
#' Simulates data where the real value y depends linearly on x, with random variation added.
#' The values are censored, so we only observe the actual value above some threshold.
#' The threshold has two different levels depending on x. This is a typical situation if y is chemical
#' measurements and x is time in years, as chemical methods tends to improve (decrease the limit of detection) over the years.
#' 
#' @param n Number of simulated observations
#' @param intercept Intercept of the linear relationship between y and x (default = 30)  
#' @param slope Slope of the linear relationship between y and x (default = -3)  
#' @param sigma Standard deviance around the linear relationship (the s.d. of y for a given x; default = 4)
#' @param x Values of x. If this is not set, values form 0 to 10 (uniform distribution) will be used
#' @param threshold_1 threshold when x <= threshold_change
#' @param threshold_2 threshold when x > threshold_change
#' @param threshold_change The level of x when threshold changes from threshold_1 to threshold_2
#' @param plot If TRUE, makes a plot of the data
#' @param seed If you want the simulation to be reproducible, set seed to some integer  
#' 
#' @return A data frame with 5 columns is returned invisibly. The variables are 1) x; 2) y (the real, unobserved value);
#' 3) y_uncens (the observed y values, which equals y for values above threshold, and is NA for values under threshold); 
#' 4) uncensored (equals 0 for values below threshold and 1 for values above threshold); and 5) threshold (the threshold 
#' value for each observation).
#' 
#' @seealso \code{\link{lc_linear}}
#' 
#' @examples 
#' 
#' # Default parameters
#' lc_simulate()
#'
#' # Modified parameters
#' sim <- lc_simulate(slope = -1.5, threshold_1 = 27, threshold_2 = 20)
#' 
#' # The data object of the output can be directly used in lm_leftcensored 
#' result <- lc_linear(sim$data)
#' 
#' @export
# Some R code to simulate data from the above model
lc_simulate <- function(
  n = 100,
  intercept = 30,
  slope = -3,
  sigma = 4,
  x = NULL,            # takes precedence over n (n is ignored if x is set)
  threshold_1 = 20,
  threshold_2 = 10,
  threshold_change = 4,
  plot = TRUE,
  seed = NULL){
  
  if (!is.null(seed)){
    # Set the seed so this is repeatable
    set.seed(seed)
  }

  if (is.null(x)){
    x <- sort(runif(n, 0, 10)) # Sort as it makes the plotted lines neater
  } else {
    n <- length(x)     # the given n is ignored if x is set
  } 
  y <- rnorm(n, mean = intercept + slope * x, sd = sigma)
  
  threshold <- ifelse(x <= threshold_change, threshold_1, threshold_2)
  uncensored <- ifelse(y > threshold, 1, 0)
  y_plot <- ifelse(y > threshold, y, threshold)  # if below threshold => y = threshold (for plotting)
  y_obs <- ifelse(y > threshold, y, NA)      # if below threshold => y = NA (for analysis)
  
  result_data <- data.frame(x = x,
                            y = y_obs,
                            y_plot = y_plot, 
                            y_real = y,
                            uncensored = uncensored,
                            threshold = threshold)
  
  if (plot){
    lc_plot_sim(result_data, intercept, slope, threshold_1, threshold_2, threshold_change)
  }
  result <- list(
    data = result_data,
    intercept = intercept, 
    slope = slope)
  invisible(result)
}

#' @export
#' 
#' 
lc_plot_sim <- function(data, 
                        sim_intercept, sim_slope, 
                        sim_threshold_1, sim_threshold_2, sim_threshold_change){
  # Real (partly unobserved) concentrations - NOT used in analysis:
  plot(y_real ~ x, data)                          
  # Observed concentrations - used in analysis:
  points(y_plot ~ x, data, pch = 20, col = "blue2") 
  # Actual (unknown) trend line
  abline(sim_intercept, sim_slope, col = "blue2")
  # Plot censoring limits  
  segments(x0 = c(min(data$x), sim_threshold_change), 
           x1 = c(sim_threshold_change, max(data$x)),
           y0 = c(sim_threshold_1, sim_threshold_2),
           col = "red")
}

