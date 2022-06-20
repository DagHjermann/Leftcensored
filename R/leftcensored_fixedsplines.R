#' Linear regression for left-censored data
#' 
#' This function runs spline regression (with fixed number and placement of knots) when the dependent
#' (y) variable is left-censored. That is, if the actual value is below some below some limit 
#' of quantification, LOQ, we we cannot measure it. We only know that the actual value is
#' somewhere below LOQ. If the data are given as typical chemical measurements with one column
#' containing "<" for measurements below LOQ, and the value column giving the LOQ value in those cases,
#' you may want to use prepare_data() to make your data ready.
#' 
#' @param data The data must have (at least) four columns: the predictor variable, the (uncensored) response variable, 
#' a variable which is 1 for every uncensored observation, and a variable with the threshold of censoring for every censored observation
#' @param x Variable name for the predictor (independent) variable        
#' @param y Variable name for the response (dependent) variable. The censored values can be anything (they will be set to NA 
#' depending on the values of \code{uncensored}    
#' @param uncensored Variable name for a variable which is 1 for uncensored values and 0 for censored
#' values.     
#' @param threshold Variable name for a variable containing the threshold for censoring (e.g. for chemical data,
#' limit of detection or limit of quantification).It may be constant, or it may vary for each observation censored observation.   
#' @param n.chains The number of MCMC chains (replicates) to run. The default is 4. Using more than 1 chain enables us to say whether 
#' @param knots Either a single number for the number of equidistant knots, or a vector giving
#' the placement of the knots.   
#' @param resolution The number of points along the x axis used to describe the spline. 
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
leftcensored_fixedsplines <- function(data,
                            x = "x", 
                            y = "y_uncens", 
                            uncensored = "uncensored",
                            threshold = "threshold",
                            knots = 9,
                            resolution = 50,
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
  
  # All the values of 'uncensored' that are not 1, will be set to 
  data[[uncensored]][!data[[uncensored]] == 1] <- 0
  
  # For making predicted lines (with SE) 
  xmin <- min(data[[x]], na.rm = TRUE)
  xmax <- max(data[[x]], na.rm = TRUE)
  x.out <- seq(xmin, xmax, length = resolution)
  
  # Create knots and splinen basis functions  
  if (length(knots) == 1)
    knots <- seq(xmin, xmax, length = knots)

  # Generate basis splines at the observed x values, using splines::bs()  
  B <- t(splines::bs(data[[x]], knots = knots, degree=3, intercept = TRUE)) 
  B.out <- t(splines::bs(x.out, knots = knots, degree=3, intercept = TRUE)) 
  
  # Generate basis splines for describing predicted spline line   
  # x <- seq(from = xmin, to = xmax, length = resolution) # generating inputs
  # B <- t(splines::bs(x, knots = knots, degree=3, intercept = TRUE)) 
  
  # Jags code to fit the model to the simulated data

  model_code = '
model
{
  y.hat <- a0*x + a %*% B ## expected response  
  y.hat.out <- a0*x.out + a %*% B.out ## expected response  
  
  # Likelihood
  for (i in 1:n){
    uncensored[i] ~ dinterval(y[i], threshold[i]) 
    y[i] ~ dnorm(y.hat[i], tau)
  }
  sigma <- 1/tau         ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  # Linear effect  
  a0 ~ dnorm(0, lambda[1])
  
  # Specify priors for spline terms
  for (k in 1:K) {
    a[k] ~ dnorm(0, lambda[2])
  }
  
  ## smoothing parameter priors   
  for (m in 1:2) {
    lambda[m] ~ dgamma(.05,.005)
    rho[m] <- log(lambda[m])
  }
}
'
  ### Set up data and parameters
  # Set up the data
  model_data <- list('n' = nrow(data), 
                    'y' = data[[y]], 
                    'uncensored' = data[[uncensored]],
                    'threshold' = data[[threshold]],
                    'x' = data[[x]],
                    'x.out' = x.out,
                    'B' = B,
                    'B.out' = B.out,
                    'K' = dim(B)[1])
  # Choose the parameters to watch
  model_parameters <-  c('a0', 'a', 'sigma','tau','y.hat.out')
  
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
  
  summary <- summary(model_mcmc)
  
  # summary(model_mcmc) %>% str()
  
  # Get predicted line 
  quants <- summary$quantiles
  length.out <- length(x.out)
  pick_rownames <- sprintf("y.hat.out[%i]", 1:length.out)
  plot_data <- data.frame(
    x = x.out, 
    y = quants[pick_rownames,"50%"],
    y_lo = quants[pick_rownames,"2.5%"],
    y_hi = quants[pick_rownames,"97.5%"]
  )
  
  list(summary = summary,
       plot_data = plot_data,
       model = model_mcmc)
  
  
}
