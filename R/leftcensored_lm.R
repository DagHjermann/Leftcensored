#' Linear regression for left-censored data
#' 
#' This function runs linear regression when the dependent (y) variable is left-censored. That is,
#' if the actual value is below some below some limit of quantification, LOQ, we we cannot measure it. We only know that the 
#' actual value is somewhere below LOQ. The variables need to have particular names, so you may want to use
#' prepare_data() to make your data ready for this function.
#' 
#' @param data The data must have (at least) four columns: the predictor variable, the (uncensored) response variable, 
#' a variable which is 1 for every uncensored observation, and a variable with the threhold of censoring for every censored observation
#' @param x Variable name for the predictor (independent) variable        
#' @param y Variable name for the response (dependent) variable. The censored values can be anything (they will be set to NA 
#' depending on the values of \code{uncensored}    
#' @param uncensored Variable name for a variable which is 1 for uncensored values and 0 for censored
#' values.     
#' @param threshold Variable name for a variable containing the threshold for censoring (e.g. for chemical data,
#' limit of detection or limit of quantification).It may be constant, or it may vary for each observation censored observation.   
#' @param measurement_error (optional) Variable name for a variable giving the standard value of the measurement error.
#' If given, this is taken into account, increasing the quantile range (i.e. the confidence interval) of the estimated
#' intercept and slope.   
#' @param resolution The number of points along the x axis used to describe the spline. 
#' @param n.chains The number of MCMC chains (replicates) to run. The default is 4. Using more than 1 chain enables us to say whether 
#' @param n.iter The number of iterations for each MCMC chains. The default is 5000, which is usually sufficient for this application.
#' @param n.burnin The number of iterations to remove at start of each MCMC chain, before results are collected for statistics. The default is 1000.
#' If n.burnin is too small, plots of the trace (see examples) will show whether the chains are homogeneous along the chain (there should be from 
#' no decreasing or increasing trend in the traceplot). One can also use R2jags::traceplot to assess whether the chains behave differently, 
#' depending on their different starting point. If they do, n.burnin should be increased.
#' @param n.thin The number of MCMC iterations that are kept for statistics.
#' @param type Type of Jags model formulation used. The default (type = 'Qi') uses the method of Qi et al. (2022). The alternative is 
#' type = 'dinterval', which uses the method used in the standard JAGS documentation. The advantage of the method of Qi et al. is that
#' it returns the deviance information criterion (DIC).   
#' 
#' @keywords Statistics, regression, censored
#' 
#' @references 
#' Qi X, Zhou S and Plummer M. 2022. On Bayesian modeling of censored data in JAGS. BMC Bioinformatics 23: 102
#' (doi: 10.1186/s12859-021-04496-8). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8944154/  
#' 
#' Plummer N (2017). JAGS Version 4.3.0 user manual. https://sourceforge.net/projects/mcmc-jags/files/  
#' 
#' @return The function returns a list with two parts: \code{summary} and \code{model}. \code{summary} shows a summary of the result, i.e., estimates
#' of the parameters of the linear regression: intercept, slope and sigma (the estimated standard deviation of the data around the regression line).
#' It is common to use the quantiles for parameter estimates, i.e., using the  
#' 50% quantile as the "best estimate" of the parameters, and using the 2.5% and 97.5% quantiles as endpoints of a 95% confidence interval. 
#' \code{model} is the output of the jags() command, which is what \code{lc_linear()} runs under the hood for estimation. 
#' It is an MCMC object, which has methods for functions such as plot (see ?mcmc). If you have some knowledge of the MCMC technique for     
#' state-space models, this can be used for diagnostic plots of the model. See examples.
#'   
#' @examples
#' 
#' #
#' # 1. No measurement error
#' #
#' 
#' # Simulate data
#' set.seed(11)
#' sim <- lc_simulate(n = 30)     # also plots the data
#' 
#' # Estimate regression
#' result <- lc_linear(sim$data)
#' 
#' # Get best estimates and plot its regression line on top of the plot of the data  
#' a <- result$intercept["50%"]
#' b <- result$slope["50%"]
#' abline(a, b, col = "green2")
#' lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
#' lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")
#' 
#' # Check quantiles of the parameters
#' result$summary$quantiles
#' 
#' # Make a standard MCMC plot: the trace and the density for each estimated parameter (long output)
#' # par(mar = c(2,4,3,1))
#' # plot(result$model)
#' 
#' # Plot the trace for each MCMC run (long output)  
#' # par(mfrow = c(2,2), mar = c(2,4,3,1))
#' # coda::traceplot(result$model, ask = FALSE)
#' 
#' #
#' # 2. With measurement error
#' #
#' 
#' # Simulate data
#' set.seed(11)
#' sim <- lc_simulate(n = 30)     # also plots the data
#' 
#' # Add SD of the measurement error
#' # (For simplicity, we use the same SD for all rows, bt it may vary for case to case,
#' # for instance if the error is given as a percentage of the observed value)
#' sim$data$meas_error <- 5
#' 
#' # Estimate regression
#' result_me <- lc_linear(sim$data)
#' 
#' # Get best estimates and plot the regression line  
#' a <- result_me$intercept["50%"]
#' b <- result_me$slope["50%"]
#' abline(a, b, col = "green2")
#' lines(y_lo ~ x, data = result_me$plot_data, lty = "dashed", col = "green2")
#' lines(y_hi ~ x, data = result_me$plot_data, lty = "dashed", col = "green2")
#' 
#' @export
lc_linear <- function(data,
                      x = "x", 
                      y = "y", 
                      uncensored = "uncensored",
                      threshold = "threshold",
                      resolution = 50,
                      n.chains = 4, 
                      n.iter = 5000, 
                      n.burnin = 1000, 
                      n.thin = 2,
                      mean_y = NULL,
                      sd_y = NULL,
                      plot_input = FALSE,
                      plot_norm = FALSE,
                      detailed = FALSE,
                      type = "Qi",
                      measurement_error = NULL){
  
  # Censoring vs truncation:
  # https://stats.stackexchange.com/a/144047/13380 
  # Also see here for an alternative way to set up the model:
  # https://stats.stackexchange.com/questions/6870/censoring-truncation-in-jags
  # NOTE: CHECK THIS for multi-level statistical models:
  # https://stats.stackexchange.com/questions/185254/multi-level-bayesian-hierarchical-regression-using-rjags
  
  if (type == "Qi" & is.null(measurement_error)){
    result <- lc_linear_qi(data=data, x=x, y=y, uncensored=uncensored, threshold=threshold,
                           resolution=resolution, n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin,
                           n.thin=n.thin, mean_y=mean_y, sd_y=sd_y,
                           plot_input=plot_input, plot_norm=plot_norm, detailed = FALSE)
  } else if (type == "Qi" & !is.null(measurement_error)){
    result <- leftcensored:::lc_linear_qi_measerror(
                           data=data, x=x, y=y, uncensored=uncensored, threshold=threshold,
                           resolution=resolution, n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin,
                           n.thin=n.thin, mean_y=mean_y, sd_y=sd_y,
                           plot_input=plot_input, plot_norm=plot_norm, detailed = FALSE,
                           measurement_error = measurement_error)
    
  } else {
    result <- lc_linear_dinterval(data=data, x=x, y=y, uncensored=uncensored, threshold=threshold,
                                  resolution=resolution, n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin,
                                  n.thin=n.thin, mean_y=mean_y, sd_y=sd_y,
                                  plot_input=plot_input, plot_norm=plot_norm, detailed = FALSE)
  }
  
  result
  
}

