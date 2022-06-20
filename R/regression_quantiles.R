#' Get regression quantiles for y for a number of x values within the range of x of the given data set
#' 
#' @param data data frame containing x values in a column named 'x'
#' @param mcmc MCMC object
#' @param n Number of x values used (default is 20)
#' 
#' @return Matrix of n rows and 5 columns, containing quantiles for a given x value 
#' 
#' @examples
#' sim <- leftcensored_simulate()
#' result <- lm_linear(sim$data)
#' regression_quantiles(sim$data, result$model)
#' 
#' @export
regression_quantiles <- function(data, mcmc, n = 20){
  x <- seq(min(data$x), max(data$x), length = n)
  y <- x %>% purrr::map(predicted_y_quantiles, mcmc = mcmc,
                 quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)) %>% bind_cols() %>% t()
  colnames(y) <- c("2.5%", "25%", "50%", "75%", "97.5%") 
  list(x=x, y=y)
}
