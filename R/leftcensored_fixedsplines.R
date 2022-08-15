#' Linear regression for left-censored data
#' 
#' This function runs spline regression (with fixed number and placement of knots) when the dependent
#' (y) variable is left-censored. That is, if the actual value is below some below some limit 
#' of quantification, LOQ, we we cannot measure it. We only know that the actual value is
#' somewhere below LOQ. If the data are given as typical chemical measurements with one column
#' containing "<" for measurements below LOQ, and the value column giving the LOQ value in those cases,
#' you may want to use lc_prepare() to make your data ready.
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
#' @param n.iter The number of iteratioterations to remove at start of each MCMC chain, before results are collected for statistics. The default is 1000.
#' If n.burnin is too small, plons for each MCMC chains. The default is 5000, which is usually sufficient for this application.
#' @param n.burnin The number of its of the trace (see examples) will show whether the chains are homogeneous along the chain (there should be from 
#' no decreasing or increasing trend in the traceplot). One can also use R2jags::traceplot to assess whether the chains behave differently, 
#' depending on their diferent starting point. If they do, n.burnin should be increased.
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
#' # Simulate data and estimate regression
#'X <- seq(from=-1, to=1, by=.025) # generating inputs
#'B <- t(splines::bs(X, knots=seq(-1,1,1), degree=3, intercept = TRUE)) # creating the B-splines
#'num_data <- length(X); num_basis <- nrow(B)
#'a0 <- 0.2 # intercept
#'
#'set.seed(991)
#'# num_basis <- 6
#'a <- rnorm(num_basis, 0, 1) # coefficients of B-splines
#'n_param <- length(a)
#'
#'Y_true <- as.vector(a0*X + a%*%B) # generating the output
#'Y <- Y_true + rnorm(length(X),0,.1) # adding noise
#'
#'dat_sim <- data.frame(x = X, y_uncensored = Y, y_true = Y_true)
#'# ggplot(dat_sim, aes(x, y_uncensored)) +
#'#   geom_point() +
#'#   geom_line(aes(y = y_true), color = "blue")
#'
#'# Add censoring 
#'dat_sim$y <- dat_sim$y_uncensored
#'dat_sim$uncensored <- 1
#'threshold_fixed <- -0.3
#'sel <- dat_sim$y_uncensored < threshold_fixed
#'dat_sim$y[sel] <- NA  
#'dat_sim$uncensored[sel] <- 0  
#'dat_sim$threshold <- threshold_fixed
#'
#'# Plot
#'lc_plot(dat_sim)
#'
#' result_nonlin <- lc_fixedsplines(dat_sim, knots = 3)
#' 
#' # Show result
#' lc_plot(dat_sim, results = result_nonlin)
#' 
#' \dontrun{
#' # Show summary of the posterior for each variable:
#' result_nonlin$summary
#' }
#' 
#' @export
lc_fixedsplines <- function(data,
                            x = "x", 
                            y = "y", 
                            uncensored = "uncensored",
                            threshold = "threshold",
                            knots = 9,
                            resolution = 50,
                            n.chains = 4, 
                            n.iter = 5000, 
                            n.burnin = 1000, 
                            n.thin = 2,
                            type = "Qi"){
  
  # Censoring vs truncation:
  # https://stats.stackexchange.com/a/144047/13380 
  # Also see here for an alternative way to set up the model:
  # https://stats.stackexchange.com/questions/6870/censoring-truncation-in-jags
  # NOTE: CHECK THIS for multi-level statistical models:
  # https://stats.stackexchange.com/questions/185254/multi-level-bayesian-hierarchical-regression-using-rjags
  
  if (is.numeric(grep(type, c("Qi", "dbern", "Bernoulli")))){
    result <- lc_fixedsplines_qi(
      data=data, x=x, y=y, uncensored=uncensored, threshold=threshold,
      knots=knots,
      resolution=resolution, n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin,
      n.thin=n.thin)
  } else {
    result <- lc_fixedsplines_dinterval(
      data=data, x=x, y=y, uncensored=uncensored, threshold=threshold,
      knots=knots,
      resolution=resolution, n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin,
      n.thin=n.thin)
  }
  
  result
  
}
