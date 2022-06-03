#' Get data for plotting confidence interval of the regression
#'  
#' @param data data frame containing x values in a column named 'x'
#' @param mcmc MCMC object
#' @param n Number of x values used (default is 20). Increase for more smoothness of the CI
#' 
#' @return Matrix of n rows and 5 columns: x, y, y_lo and y_hi 
#' 
#' @examples
#' sim <- leftcensored_simulate()
#' result <- leftcensored_lm(sim$data)
#' plotdata <- regression_ci(sim$data, result$model)
#' plot(y ~ x, type = "l", data = plotdata, ylim = range(plotdata[,3:4]))
#' lines(y_lo ~ x, data = plotdata, lty = 2)
#' lines(y_hi ~ x, data = plotdata, lty = 2)
#'  
#' @export
regression_ci <- function(data, mcmc, n = 20){
  qs <- regression_quantiles(data, mcmc, n = n)
  data.frame(x = qs$x, y = qs$y[,"50%"], y_lo = qs$y[,"2.5%"], y_hi = qs$y[,"97.5%"])
}
