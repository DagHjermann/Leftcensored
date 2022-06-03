#' Linear regression for left-censored data
#' 
#' This function runs linear regression when the dependent (y) variable is left-censored. That is,
#' if the actual value is below some below some limit of quantification, LOQ, we we cannot measure it. We only know that the 
#' actual value is somewhere below LOQ. The variables need to have particular names, so you may want to use
#' prepare_data() to make your data ready for this function.
#' 
#' @param data The data must have four columns, named x, y_cens, y_aboveLOQ, y_LOQ. x is the
#' predictor (independent) variable; y_cens is the observed dependent data, which typically is NA for data below LOQ;
#' y_aboveLOQ is 1 for data above LOQ; and y_LOQ is the LOQ (which may be constant, or may vary for each observation).
#' @param n.chains The number of MCMC chains (replicates) to run. The default is 4. Using more than 1 chain enables us to say whether 
#' @param n.iter The number of iterations for each MCMC chains. The default is 5000, which is usually sufficient for this application.
#' @param n.burnin The number of iterations to remove at start of each MCMC chain, before results are collected for statistics. The default is 1000.
#' If n.burnin is too small, plots of the trace (see examples) will show whether the chains are homogeneous along the chain (there should be from 
#' no decreasing or increasing trend in the traceplot). One can also use R2jags::traceplot to assess whether the chains behave differently, 
#' depending on their diferent starting point. If they do, n.burnin should be increased.
#' @param n.thin The number of MCMC iterations that are kept for statistics.
#' 
#' @keywords Statistics, regression, censored
#' 
#' @return The function returns a list with two parts: \code{summary} and \code{model}. \code{summary} shows a summary of the result, i.e., estimates
#' of the parameters of the linear regression: intercept, slope and sigma (the estimated standard deviation of the data around the regression line).
#' It is common to use the quantiles for parameter estimates, i.e., using the  
#' 50% quantile as the "best estimate" of the parameters, and using the 2.5% and 97.5% quantiles as endpoints of a 95% confidence interval. 
#' \code{model} is the output of the jags() command, which is what \code{leftcensored_lm()} runs under the hood for estimation. 
#' It is an MCMC object, which has methods for functions such as plot (see ?mcmc). If you have some knowledge of the MCMC technique for     
#' state-space models, this can be used for diagnostic plots of the model. See examples.
#'   
#' @examples
#' # Simulate data and estimate regression
#' sim <- leftcensored_simulate(n = 30)
#' result <- leftcensored_lm(sim$data)
#' 
#' # Get best estimates and plot its regression line on top of the plot  
#' a <- result$summary$quantiles["intercept", "50%"]
#' b <- result$summary$quantiles["slope", "50%"]
#' abline(a, b, col = "green2")
#' 
#' # Example with real data
#' # Prepare the data
#' data_test <- prepare_data(concentrations)
#' 
#' # Perform the analysis
#' result <- leftcensored_lm(df_test2)
#' 
#' # MCMC summary
#' result$summary
#' 
#' # Check quantiles of the parameters
#' result$summary$quantiles
#' 
#' # Make a standard MCMC plot: the trace and the density for each estimated parameter  
#' par(mar = c(2,4,3,1))
#' plot(result$model)
#' 
#' # Plot the trace for each MCMC run  
#' par(mfrow = c(2,2), mar = c(2,4,3,1))
#' traceplot(result$model, ask = FALSE)
#' 
#' @export
leftcensored_lm <- function(data,
                            x = "x", 
                            y = "y_cens", 
                            uncensored = "y_aboveLOQ",
                            threshold = "y_LOQ",
                            n.chains = 4, 
                            n.iter = 5000, 
                            n.burnin = 1000, 
                            n.thin = 2){
  
  # Censoring vs truncation:
  # https://stats.stackexchange.com/a/144047/13380 
  # Also see here for an alternative way to set up the model:
  # https://stats.stackexchange.com/questions/6870/censoring-truncation-in-jags
  # NOTE: CHECK THIS for multi-level statistical models:
  # https://stats.stackexchange.com/questions/185254/multi-level-bayesian-hierarchical-regression-using-rjags
  
  # Set all censored data to NA (if not already done)
  # Important! Otherwise all LOQ stuff is ignored
  data[[y]][!data[[uncensored]] == 1] <- NA
  
  # Jags code to fit the model to the simulated data
  # Jags code to fit the model to the simulated data
  
  model_code = '
model
{
  # Likelihood
  for (i in 1:n) {
    y_aboveLOQ[i] ~ dinterval(y_cens[i], y_LOQ[i])
    y_cens[i] ~ dnorm(intercept + slope * x[i], sigma^-2)
  }

  # Priors
  intercept ~ dnorm(0, 100^-2)
  slope ~ dnorm(0, 100^-2)
  sigma ~ dunif(0, 10)
}
'
  ### Set up data and parameters
  # Set up the data
  model_data <- list(n = nrow(data), 
                    y_cens = data[[y]], 
                    y_aboveLOQ = data[[uncensored]],
                    y_LOQ = data[[threshold]],
                    x = data[[x]])
  # Choose the parameters to watch
  model_parameters <-  c("intercept", "slope", "sigma")
  
  ### Run model
  # Run the model
  model_run <- R2jags::jags(data = model_data,
                   parameters.to.save = model_parameters,
                   model.file=textConnection(model_code),
                   n.chains=n.chains,   # Number of different starting positions
                   n.iter = n.iter,     # Number of iterations
                   n.burnin = n.burnin, # Number of iterations to remove at start
                   n.thin = n.thin)     # Amount of thinning
  
  # model_run
  model_mcmc <- coda::as.mcmc(model_run)
  # summary(model_mcmc) %>% str()
  list(summary = summary(model_mcmc),
       model = model_mcmc)
  
}
