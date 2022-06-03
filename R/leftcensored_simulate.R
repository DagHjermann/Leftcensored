#' Simulate censored data
#' 
#' Simulates data where the real value y depends linearly on x, with random variation added.
#' The values are censored, so we only observe the actual value above some limit LOQ.
#' LOQ has two different levels depending on x. This is a typical situation if y is chemical
#' measurements and x is time in years, as chemical methods tends to improve (decrease LOQ) over the years.
#' 
#' @param n Number of simulated observations
#' @param intercept Intercept of the linear relationship between y and x
#' @param slope Slope of the linear relationship between y and x
#' @param sigma Standard deviance around the linear relationship (the s.d. of y for a given x)
#' @param LOQ_1 LOQ when x <= LOQ_change
#' @param LOQ_2 LOQ when x > LOQ_change
#' @param LOQ_change The level of x when LOQ changes from LOQ_1 to LOQ_2
#' @param plot If TRUE, makes a plot of the data
#' @param seed If you want the simulation to be reproducible, set seed to some integer  
#' 
#' @return A data frame with 5 columns is returned invisibly. The variables are 1) x; 2) y (the real, unovserved value);
#' 3) y_cens (the observed y values, which equals y for values above LOQ, and is NA for values under LOQ); 
#' 4) y_aboveLOQ (equals 0 for values below LOQ and 1 for values above LOQ); and 5) y_LOQ (the LOQ 
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
#' sim <- leftcensored_simulate(slope = -1.5, LOQ_1 = 27, LOQ_2 = 20)
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
  LOQ_1 = 20,
  LOQ_2 = 5,
  LOQ_change = 4,
  plot = TRUE,
  seed = NULL){
  
  if (!is.null(seed)){
    # Set the seed so this is repeatable
    set.seed(seed)
  }
  x <- sort(runif(n, 0, 10)) # Sort as it makes the plotted lines neater
  y <- rnorm(n, mean = intercept + slope * x, sd = sigma)
  
  y_LOQ <- ifelse(x <= LOQ_change, LOQ_1, LOQ_2)
  y_aboveLOQ <- ifelse(y > y_LOQ, 1, 0)
  y_obs <- ifelse(y > y_LOQ, y, y_LOQ)  # if below LOQ => y = LOQ (for plotting)
  y_cens <- ifelse(y > y_LOQ, y, NA)    # if below LOQ => y = NA (for analysis)
  
  if (plot){
    # Plot
    plot(x, y) # actual concentrations
    # Measured concentrations
    points(x, y_obs, pch = 19, col = "blue2")
    # Concentrations below LOQ
    sel <- y_aboveLOQ==0
    points(x[sel], y_obs[sel], pch = 20, col = "red")
    # Actual trend line
    lines(x, intercept + slope * x)
    # Plot censoring limits  
    segments(x0 = c(min(x), LOQ_change), 
             x1 = c(LOQ_change, max(x)),
             y0 = c(LOQ_1, LOQ_2),
             col = "red")
  }
  result <- list(
    data = data.frame(x = x,
                       y = y,
                       y_cens = y_cens, 
                       y_aboveLOQ = y_aboveLOQ,
                       y_LOQ = y_LOQ),
    intercept = intercept, 
    slope = slope)
  invisible(result)
}
