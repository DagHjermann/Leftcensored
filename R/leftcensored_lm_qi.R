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
#' @param resolution The number of points along the x axis used to describe the spline. 
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
#' \code{model} is the output of the jags() command, which is what \code{lc_linear()} runs under the hood for estimation. 
#' It is an MCMC object, which has methods for functions such as plot (see ?mcmc). If you have some knowledge of the MCMC technique for     
#' state-space models, this can be used for diagnostic plots of the model. See examples.
#'   
#' @examples
#' # Simulate data and estimate regression
#' set.seed(11)
#' sim <- lc_simulate(n = 30)
#' result <- lc_linear(sim$data)
#' 
#' # Get best estimates and plot its regression line on top of the plot  
#' a <- result$intercept["50%"]
#' b <- result$slope["50%"]
#' abline(a, b, col = "green2")
#' lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
#' lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")
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
#' coda::traceplot(result$model, ask = FALSE)
#' 
#' @export
lc_linear_qi <- function(data,
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
                      detailed = FALSE){
  
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
  
  if (plot_input){
    plot(data$x, data$y, ylim = range(data$y, data$threshold, na.rm = TRUE))
    points(data$x[data$uncensored == 0], data$threshold[data$uncensored == 0], pch = 20, col = 2)
  }
  
  # For making predicted lines (with SE) 
  xmin <- min(data[[x]], na.rm = TRUE)
  xmax <- max(data[[x]], na.rm = TRUE)
  x.out <- seq(xmin, xmax, length = resolution)

  # Jags code to fit the model to the simulated data

  model_code = '
model
{
  # Uncensored observations 
  for (o in 1:O) {
    y.uncens[o] ~ dnorm(intercept + slope * x[o], sigma^-2)
  }
  # Censored observations 
  for (c in 1:C) {
    Z1[c] ~ dbern(p[c])
    p[c] <- pnorm(cut[c], intercept + slope * x[O+c], sigma^-2)
  }
  for (i in 1:resolution) {
    y.hat.out.norm[i] ~ dnorm(intercept + slope * x.out[i], sigma^-2)
  }
  
  for (i in 1:resolution){
    y.hat.out[i] <- y.hat.out.norm[i]*sd_y + mean_y
  }

  # Priors
  intercept ~ dnorm(0, 100^-2)
  slope ~ dnorm(0, 100^-2)
  sigma ~ dunif(0, 10)
}
'
  ### Set up data and parameters
  
  # Normalize y
  # Achieves mean = 0
  if (is.null(mean_y))
    mean_y <- mean(data[[y]], na.rm = TRUE)
  if (is.null(sd_y))
    sd_y <- sd(data[[y]], na.rm = TRUE)
  # Functions for normalization and "un-normalization" (back-transformation)
  norm_y <- function(x) (x-mean_y)/sd_y
  unnorm_y <- function(x) x*sd_y + mean_y

  # Normalize (or more correctly centralize) x
  # Achieves mean = 0
  mean_x <- mean(data[[x]], na.rm = TRUE)
  # Functions for normalization and "un-normalization" (back-transformation)
  norm_x <- function(x) x-mean_x
  
  # Set up the data for Qi's method
  # Data has normalized y values and centralized x values
  
  # Split the data into uncensored and censored parts
  data_obs <- data[data[[uncensored]] %in% 1,]
  data_cen <- data[data[[uncensored]] %in% 0,]
  data_all <- rbind(data_obs, data_cen)

  model_data <- list(x = norm_x(data[[x]]),
                     y.uncens = norm_y(data_obs[[y]]),
                     O = nrow(data_obs),
                     Z1 = rep(1, nrow(data_cen)),  # because all are left-censored, see text below 'Model 2' in Qi' et al. 2022's paper
                     cut = norm_y(data_cen[[threshold]]),
                     C = nrow(data_cen),
                    x.out = norm_x(x.out),
                    resolution = resolution,
                    mean_y = mean_y,
                    sd_y = sd_y)

  if (plot_norm){
    plot(model_data$x, model_data$y.uncens, ylim = range(model_data$y.uncens, model_data$threshold, na.rm = TRUE))
    points(model_data$x[model_data$uncensored == 0], model_data$threshold[model_data$uncensored == 0], pch = 20, col = 2)
  }
  
  # Choose the parameters to watch
  if (detailed){
    model_parameters <-  c("intercept", "slope", "sigma", 
                           "y.hat.out.norm", "y.hat.out", "y.uncens")
  } else {
    model_parameters <-  c("intercept", "slope", "sigma", "y.hat.out")
  }
  
  # Initial values  
  init_model_df <- data.frame(x = data_obs[[x]], y.uncens = data_obs[[y]])
  init_model <- lm(y.uncens ~ x, data = init_model_df)
  init_summ <- summary(init_model)$coef
  jags.inits <- function(){
    list("intercept" = rnorm(1, mean = init_summ[1,1], sd = init_summ[1,2]), 
         "slope" =  rnorm(1, mean = init_summ[2,1], sd = init_summ[2,2]),
         "sigma" = runif(1))
  }
  
  ### Run model
  # Run the model
  model_run <- R2jags::jags(
    data = model_data,
    init = jags.inits,
    parameters.to.save = model_parameters,
    model.file=textConnection(model_code),
    n.chains=n.chains,   # Number of different starting positions
    n.iter = n.iter,     # Number of iterations
    n.burnin = n.burnin, # Number of iterations to remove at start
    n.thin = n.thin)     # Amount of thinning
  
  
  # model_run
  model_mcmc <- coda::as.mcmc(model_run)
  # summary(model_mcmc) %>% str()
  
  summary <- summary(model_mcmc)

  #
  # Get predicted line 
  #
  quants <- summary$quantiles
  length.out <- length(x.out)
  pick_rownames <- sprintf("y.hat.out[%i]", 1:length.out)
  # y and lower and upper CI  values are back-transformed (un-normalized) using unnorm:
  plot_data <- data.frame(
    x = x.out, 
    y = quants[pick_rownames,"50%"],
    y_lo = quants[pick_rownames,"2.5%"],
    y_hi = quants[pick_rownames,"97.5%"]  )
  
  # Get regression coefficients, originals:
  intercept.norm <- summary$quantiles["intercept",]
  slope.norm <- summary$quantiles["slope",]
  #
  # Get regression coefficients, back-transformed
  intercept <- intercept.norm*sd_y - slope.norm*sd_y*mean_x + mean_y
  slope <- slope.norm*sd_y
  #
  # Basis for formulae above:
  #   y' = (y - mean_y)/sd_y   (1)
  #   x' = x - mean_x          (2)
  # Slope formula used for normalized data:
  #   y' = a' + b'x'           (3)
  #   - where a' and b' are the intercept and slope found for normalized data
  # To get the formulae used above, substitute (1) and (2) into (3):
  #   (y - mean_y)/sd_y = a' + b'(x - mean_x)               (4)
  # And solve for y on the left side. This results in
  #   y = [a'*sd_y - b'*sd_y*mean_x + mean_y] + [b*sd_y]*x  (5)
  # Where the two parentheses are the back-transformed intercept and slope, respectively
  #
  
  list(summary = summary(model_mcmc),
       plot_data = plot_data,
       intercept = intercept,
       slope = slope,
       model_data = model_data,     # 
       model = model_mcmc,
       mean_y = mean_y,
       sd_y = sd_y,
       norm_y = norm_y)
  
}


#' @export

lc_linear_measerror <- function(data,
                                x = "x", 
                                y = "y", 
                                uncensored = "uncensored",
                                threshold = "threshold",
                                measurement_error = 0.1,
                                resolution = 50,
                                n.chains = 4, 
                                n.iter = 5000, 
                                n.burnin = 1000, 
                                n.thin = 2,
                                mean_y = NULL,
                                sd_y = NULL,
                                # minimum_y = NULL,
                                detailed = TRUE){
  
  # Set all censored data to NA (if not already done)
  # Important! Otherwise all LOQ stuff is ignored
  data[[y]][!data[[uncensored]] == 1] <- NA
  
  # All the values of 'uncensored' that are not 1, will be set to 
  data[[uncensored]][!data[[uncensored]] == 1] <- 0
  
  # For making predicted lines (with SE) 
  xmin <- min(data[[x]], na.rm = TRUE)
  xmax <- max(data[[x]], na.rm = TRUE)
  x.out <- seq(xmin, xmax, length = resolution)
  
  # Jags code to fit the model to the simulated data
  
  model_code = '
model
{
  # Likelihood
  for (i in 1:n) {
    uncensored[i] ~ dinterval(y.uncens.error[i], threshold[i])
    y.uncens.error[i] ~ dnorm(y.uncens[i], sigma2^-2) 
    y.uncens[i] ~ dnorm(intercept + slope * x[i], sigma^-2)  
  }
  #  y.uncens.error[i] ~ dnorm(y.uncens[i], se_measurement[i]^-2)
  
  for (i in 1:resolution) {
    y.hat.out.norm[i] ~ dnorm(intercept + slope * x.out[i], sigma^-2)
  }
  
  for (i in 1:resolution){
    y.hat.out[i] <- y.hat.out.norm[i]*sd_y + mean_y
  }

  # Priors
  intercept ~ dnorm(0, 100^-2)
  slope ~ dnorm(0, 100^-2)
  sigma ~ dunif(0, 10)
  sigma2_sd <- 0.1*sigma2_mean
  sigma2 ~ dnorm(sigma2_mean, sigma2_sd^-2)
}
'
  ### Set up data and parameters
  # y_measerror[i] <- se_measurement*abs(y.uncens[i])
  # y.uncens.error[i] ~ dnorm(y.uncens[i], (0.1*abs(y.uncens[i])^2)
  
  # Normalize y
  # Achieves mean = 0
  if (is.null(mean_y))
    mean_y <- mean(data[[y]], na.rm = TRUE)
  if (is.null(sd_y))
    sd_y <- sd(data[[y]], na.rm = TRUE)
  # Functions for normalization and "un-normalization" (back-transformation)
  norm_y <- function(x) (x-mean_y)/sd_y
  unnorm_y <- function(x) x*sd_y + mean_y
  
  # Normalize (or more correctly centralize) x
  # Achieves mean = 0
  mean_x <- mean(data[[x]], na.rm = TRUE)
  # Functions for normalization and "un-normalization" (back-transformation)
  norm_x <- function(x) x-mean_x
  
  # Set up the data
  model_data <- list(n = nrow(data), 
                     y.uncens.error = norm_y(data[[y]]), 
                     uncensored = data[[uncensored]],
                     threshold = norm_y(data[[threshold]]),
                     x = norm_x(data[[x]]),
                     x.out = norm_x(x.out),
                     resolution = resolution,
                     mean_y = mean_y,
                     sd_y = sd_y,
                     sigma2_mean = measurement_error*mean(data[[y]], na.rm = TRUE))
  
  # Choose the parameters to watch
  if (detailed){
    model_parameters <-  c("intercept", "slope", "sigma", "y.hat.out", 
                           "y.uncens", "y.uncens.error", "uncensored")
  } else {
    model_parameters <-  c("intercept", "slope", "sigma", "y.hat.out")
  }
  
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
  
  summary <- summary(model_mcmc)
  
  #
  # Get predicted line 
  #
  quants <- summary$quantiles
  length.out <- length(x.out)
  pick_rownames <- sprintf("y.hat.out[%i]", 1:length.out)
  plot_data <- data.frame(
    x = x.out,
    y = quants[pick_rownames,"50%"],
    y_lo = quants[pick_rownames,"2.5%"],
    y_hi = quants[pick_rownames,"97.5%"]
  )
  
  # Get regression coefficients, originals:
  intercept.norm <- summary$quantiles["intercept",]
  slope.norm <- summary$quantiles["slope",]
  #
  # Get regression coefficients, back-transformed
  intercept <- intercept.norm*sd_y - slope.norm*sd_y*mean_x + mean_y
  slope <- slope.norm*sd_y
  
  list(summary = summary(model_mcmc),
       plot_data = plot_data,
       intercept = intercept,
       slope = slope,
       model_data = model_data,     # 
       model = model_mcmc,
       mean_y = mean_y,
       sd_y = sd_y,
       norm_y = norm_y)  
  
}


#' @export

lc_linear_measerror_min <- function(data,
                                x = "x", 
                                y = "y", 
                                uncensored = "uncensored",
                                threshold = "threshold",
                                measurement_error = 0.1,
                                resolution = 50,
                                n.chains = 4, 
                                n.iter = 5000, 
                                n.burnin = 1000, 
                                n.thin = 2,
                                mean_y = NULL,
                                sd_y = NULL,
                                minimum_y = -10,     # NEW: this sets a minimum value for y.uncens[i], using T()
                                detailed = TRUE){
  
  # See 'test_script' for a test of this. It increased the bias drastically,
  #   instead of decreasing bias
  
  
  # Set all censored data to NA (if not already done)
  # Important! Otherwise all LOQ stuff is ignored
  data[[y]][!data[[uncensored]] == 1] <- NA
  
  # All the values of 'uncensored' that are not 1, will be set to 
  data[[uncensored]][!data[[uncensored]] == 1] <- 0
  
  # For making predicted lines (with SE) 
  xmin <- min(data[[x]], na.rm = TRUE)
  xmax <- max(data[[x]], na.rm = TRUE)
  x.out <- seq(xmin, xmax, length = resolution)
  
  # Jags code to fit the model to the simulated data
  
  model_code = '
model
{
  # Likelihood
  for (i in 1:n) {
    uncensored[i] ~ dinterval(y.uncens.error[i], threshold[i])
    y.uncens.error[i] ~ dnorm(y.uncens[i], sigma2^-2) 
    y.uncens[i] ~ dnorm(intercept + slope * x[i], sigma^-2) T(y.min,)  
  }

  for (i in 1:resolution) {
    y.hat.out.norm[i] ~ dnorm(intercept + slope * x.out[i], sigma^-2)
  }
  
  for (i in 1:resolution){
    y.hat.out[i] <- y.hat.out.norm[i]*sd_y + mean_y
  }

  # Priors
  intercept ~ dnorm(0, 100^-2)
  slope ~ dnorm(0, 100^-2)
  sigma ~ dunif(0, 10)
  sigma2_sd <- 0.1*sigma2_mean
  sigma2 ~ dnorm(sigma2_mean, sigma2_sd^-2)
}
'
  ### Set up data and parameters
  # y_measerror[i] <- se_measurement*abs(y.uncens[i])
  # y.uncens.error[i] ~ dnorm(y.uncens[i], (0.1*abs(y.uncens[i])^2)
  
  # Normalize y
  # Achieves mean = 0
  if (is.null(mean_y))
    mean_y <- mean(data[[y]], na.rm = TRUE)
  if (is.null(sd_y))
    sd_y <- sd(data[[y]], na.rm = TRUE)
  # Functions for normalization and "un-normalization" (back-transformation)
  norm_y <- function(x) (x-mean_y)/sd_y
  unnorm_y <- function(x) x*sd_y + mean_y
  
  # Normalize (or more correctly centralize) x
  # Achieves mean = 0
  mean_x <- mean(data[[x]], na.rm = TRUE)
  # Functions for normalization and "un-normalization" (back-transformation)
  norm_x <- function(x) x-mean_x
  
  # Set up the data
  model_data <- list(n = nrow(data), 
                     y.uncens.error = norm_y(data[[y]]), 
                     uncensored = data[[uncensored]],
                     threshold = norm_y(data[[threshold]]),
                     x = norm_x(data[[x]]),
                     x.out = norm_x(x.out),
                     resolution = resolution,
                     mean_y = mean_y,
                     sd_y = sd_y,
                     sigma2_mean = measurement_error*mean(data[[y]], na.rm = TRUE),
                     y.min = minimum_y)
  
  # Choose the parameters to watch
  if (detailed){
    model_parameters <-  c("intercept", "slope", "sigma", "y.hat.out", 
                           "y.uncens", "y.uncens.error", "uncensored")
  } else {
    model_parameters <-  c("intercept", "slope", "sigma", "y.hat.out")
  }
  
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
  
  summary <- summary(model_mcmc)
  
  #
  # Get predicted line 
  #
  quants <- summary$quantiles
  length.out <- length(x.out)
  pick_rownames <- sprintf("y.hat.out[%i]", 1:length.out)
  plot_data <- data.frame(
    x = x.out,
    y = quants[pick_rownames,"50%"],
    y_lo = quants[pick_rownames,"2.5%"],
    y_hi = quants[pick_rownames,"97.5%"]
  )
  
  # Get regression coefficients, originals:
  intercept.norm <- summary$quantiles["intercept",]
  slope.norm <- summary$quantiles["slope",]
  #
  # Get regression coefficients, back-transformed
  intercept <- intercept.norm*sd_y - slope.norm*sd_y*mean_x + mean_y
  slope <- slope.norm*sd_y
  
  list(summary = summary(model_mcmc),
       plot_data = plot_data,
       intercept = intercept,
       slope = slope,
       model_data = model_data,     # 
       model = model_mcmc,
       mean_y = mean_y,
       sd_y = sd_y,
       norm_y = norm_y)  
  
}