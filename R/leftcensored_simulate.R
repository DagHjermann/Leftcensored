#' Simulate censored data
#' 
#' Simulates data where the real value y depends linearly on x, with random variation added.
#' The values are censored, so we only observe the actual value above some threshold.
#' The threshold has two different levels depending on x. This is a typical situation if y is chemical
#' measurements and x is time in years, as chemical methods tends to improve (decrease the limit of detection) over the years.
#' 
#' @param n Number of simulated observations
#' @param intercept Intercept of the linear relationship between y and x
#' @param slope Slope of the linear relationship between y and x
#' @param sigma Standard deviance around the linear relationship (the s.d. of y for a given x)
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
#' @seealso \code{\link{lm_leftcensored}}
#' 
#' @examples 
#' 
#' # Default parameters
#' leftcensored_simulate()
#'
#' # Modified parameters
#' sim <- leftcensored_simulate(slope = -1.5, threshold_1 = 27, threshold_2 = 20)
#' 
#' # The data object of the output can be directly used in lm_leftcensored 
#' result <- leftcensored_lm(sim$data)
#' 
#' @export
# Some R code to simulate data from the above model
leftcensored_simulate <- function(
  n = 100,
  intercept = 30,
  slope = -3,
  sigma = 4,
  threshold_1 = 20,
  threshold_2 = 10,
  threshold_change = 4,
  plot = TRUE,
  seed = NULL){
  
  if (!is.null(seed)){
    # Set the seed so this is repeatable
    set.seed(seed)
  }
  x <- sort(runif(n, 0, 10)) # Sort as it makes the plotted lines neater
  y <- rnorm(n, mean = intercept + slope * x, sd = sigma)
  
  threshold <- ifelse(x <= threshold_change, threshold_1, threshold_2)
  uncensored <- ifelse(y > threshold, 1, 0)
  y_obs <- ifelse(y > threshold, y, threshold)  # if below threshold => y = threshold (for plotting)
  y_uncens <- ifelse(y > threshold, y, NA)    # if below threshold => y = NA (for analysis)
  
  if (plot){
    # Plot
    plot(x, y) # actual concentrations
    # Measured concentrations
    points(x, y_obs, pch = 19, col = "blue2")
    # Concentrations below threshold
    sel <- uncensored==0
    points(x[sel], y_obs[sel], pch = 20, col = "red")
    # Actual trend line
    lines(x, intercept + slope * x)
    # Plot censoring limits  
    segments(x0 = c(min(x), threshold_change), 
             x1 = c(threshold_change, max(x)),
             y0 = c(threshold_1, threshold_2),
             col = "red")
  }
  result <- list(
    data = data.frame(x = x,
                       y = y,
                       y_uncens = y_uncens, 
                       uncensored = uncensored,
                       threshold = threshold),
    intercept = intercept, 
    slope = slope)
  invisible(result)
}
