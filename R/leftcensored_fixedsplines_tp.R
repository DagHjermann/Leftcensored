#' Thin plate splines for censored data  
#'
#' @param data 
#' @param x 
#' @param y 
#' @param uncensored 
#' @param threshold 
#' @param k 
#' @param resolution 
#' @param n.chains 
#' @param n.iter 
#' @param n.burnin 
#' @param n.thin 
#' @param model_parameters_for_convergence 
#'
#' @return
#' @export
#'
#' @examples
lc_fixedsplines_tp <- function(data,
                               x = "x", 
                               y = "y", 
                               uncensored = "uncensored",
                               threshold = "threshold",
                               k = 5,
                               resolution = 50,
                               n.chains = 2, 
                               n.iter = 4000, 
                               n.burnin = 4000, 
                               n.thin = 5,
                               model_parameters_for_convergence =  c("b","rho","lambda","scale"),
                               normalize = TRUE
){
  
  # Rename variables
  rename_check <- function(data, old, new){
    sel <- names(data) %in% old
    if (sum(sel) == 0)
      stop("Could not find variable ", old)
    names(data)[sel] <- new
    data
  }
  
  data <- rename_check(data, x, "x")
  data <- rename_check(data, y, "y")
  data <- rename_check(data, uncensored, "uncensored")
  data <- rename_check(data, threshold, "cut")
  
  # Set all censored data to NA (if not already done)
  # Important! Otherwise all LOQ stuff is ignored
  data$y[!data$uncensored == 1] <- NA
  
  # All the values of 'uncensored' that are not 1, will be set to 
  data$uncensored[!data$uncensored == 1] <- 0
  
  # Order file with uncensored data first
  data_obs <- data[sel_uncens,]
  data_cen <- data[!sel_uncens,]
  dat_cens_ordered1 <- bind_rows(data_obs, data_cen)

  # y_comb is the combination of y and cut
  # - will be used as the response in the analysis
  # - will only affect the uncensored values
  dat_cens_ordered1$y_comb <- c(data_obs$y, data_cen$cut)
  
  
  # Make x data that will be used for making the fittd line + conf.int.
  # - to be added to the rest of the data
  # - we also make a y (could perhaps be NA?)
  #
  # Make x
  dat_for_fit <- data.frame(
    x = seq(min(dat_cens_ordered1$x, na.rm = TRUE), 
            max(dat_cens_ordered1$x, na.rm = TRUE),
            length = resolution)
  )
  # Make y
  # Use ordinary gam (on uncensored data) to add the y
  # (Reason: perhaps the y value affects some sort of normalization)
  dat_for_fit$y_comb <- predict(
    mgcv::gam(y_comb ~ s(x), data = dat_cens_ordered1),
    newdata = dat_for_fit)
  
  #
  # Add "x data for fitted line" to the "actual" data
  #
  # dat_cens_ordered1 = data file
  # dat_cens_ordered2 = with addition for daa for fitted line (3o points)
  dat_cens_ordered2 <- bind_rows(
    dat_cens_ordered1,
    dat_for_fit
  )
  
  # Normalize data
  # Achieves mean = 0
  if (normalize){
    
  mean_y <- mean(dat_cens_ordered2$y_comb, na.rm = TRUE)
  mean_x <- mean(dat_cens_ordered2$x, na.rm = TRUE)
  sd_y <- sd(dat_cens_ordered2$y_comb, na.rm = TRUE)
  sd_x <- sd(dat_cens_ordered2$x, na.rm = TRUE)
  norm_x <- function(x) (x-mean_x)/sd_x
  unnorm_x <- function(x) x*sd_x + mean_x
  norm_y <- function(x) (x-mean_y)/sd_y
  unnorm_y <- function(x) x*sd_y + mean_y
  
  # Normalize x and y_comb   
  # norm <- normalize_lm(dat_cens_ordered2$x, c(data_obs$y_comb, data_cen$cut))
  dat_cens_ordered2$x <- norm_x(dat_cens_ordered2$x)
  dat_cens_ordered2$y_comb <- norm_y(dat_cens_ordered2$y_comb)
  dat_cens_ordered2$cut <- norm_y(dat_cens_ordered2$cut)
  
  }

  # Model file that will be made by 'jagam' 
  # - Note: this is not the model code that will be used. We will only,
  # use other parts of the output from jagam()
  jags.file <- paste(tempdir(),"/test.jags",sep="") 
  
  jags_code <- get_jags_model_code(bs = "tp", k = k, type = "leftcensored")
  
  form <- as.formula(paste0("y_comb ~ s(x, bs='tp', k=", k, ")"))
  
  # Makes (1) jags.file_tp5_orig (was used as basis for "_leftcens" file)
  # Make (2) jagam_object$jags.data (will be manpulated below)
  jagam_object <- jagam(form, 
                     data = dat_cens_ordered2,    # file no. 2 here
                     file = jags.file,   # this file will be overwritten (was used as basis for "_leftcens" file)
                     sp.prior = "gamma", 
                     diagonalize = TRUE)
  
  # Modify jags.data object
  n <- sum(dat_cens_ordered1$uncensored %in% 1)   # file no. 1 here 
  m <- sum(dat_cens_ordered1$uncensored %in% 0)   # file no. 1 here
  jagam_object$jags.data$n <- n   # - makes sure only these data are use for the likelihood
  jagam_object$jags.data$m <- m   #    - " -
  jagam_object$jags.data$Z <- c(rep(0,n), rep(1, m))
  jagam_object$jags.data$cut <- dat_cens_ordered1$cut[dat_cens_ordered1$uncensored %in% 0]  # jags.file
  
  # jagam_object$jags.data
  
  #
  # Run model - alt. 2 ---- 
  # - Using runjags::autorun.jags
  #
  
  # Choose the parameters to watch
  # Sample varaibles that have been inserted only to get the fitted line
  mu_fitted_names1 <- paste0("mu[", nrow(dat_cens_ordered1)+1, ":", nrow(dat_cens_ordered2), "]")
  mu_fitted_names2 <- paste0("mu[", seq(nrow(dat_cens_ordered1)+1, nrow(dat_cens_ordered2)), "]")
  
  ### Run model
  # Initial run  
  model_converged <- runjags::autorun.jags(
    data = jagam_object$jags.data,
    monitor = model_parameters_for_convergence,     
    inits = jagam_object$jags.ini,
    model = jags_code,
    n.chains = n.chains,    # Number of different starting positions
    startsample = n.iter,   # Number of iterations
    startburnin = n.burnin, # Number of iterations to remove at start
    thin = n.thin)          # Amount of thinning
  
  # Add all model parameters and get samples for them
  model_result <- runjags::extend.jags(model_converged, 
                                       add.monitor = mu_fitted_names1,
                                       sample = n.iter)
  
  # model_result
  model_mcmc <- coda::as.mcmc(model_result)
  
  summary <- summary(model_mcmc)
  quants <- summary$quantiles
  pick_rownames <- rownames(quants) %in% mu_fitted_names2
  # Make 'plot_data'  
  # y and lower and upper CI  values are back-transformed (un-normalized) using unnorm:
  plot_data <- data.frame(
    x = dat_for_fit$x, 
    y = quants[pick_rownames,"50%"],
    y_lo = quants[pick_rownames,"2.5%"],
    y_hi = quants[pick_rownames,"97.5%"]  )
  
  result <- list(summary = summary,
                 plot_data = plot_data,
                 mean_y = mean_y,
                 sd_y = sd_y,
                 norm_y = norm_y,
                 dic_all = dic.pd,
                 dic = dic,
                 model_data = model_data)  
  
  if (keep_model)
    result = append(result, list(model = model_mcmc))
  if (keep_model_from_jags)
    result = append(result, list(model_from_jags = model_result))
  
  result
  
}


get_jags_model_code <- function(bs = "tp",
                                k = 5,
                                type = "leftcensored"){
  
  if (bs == "tp" & k == 3 & type == "leftcensored"){
    
    code <- '
model {
  mu <- X %*% b ## expected response
  for (i in 1:n) { y[i] ~ dnorm(mu[i], tau) } ## response 
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- max(pnorm(cut[j], mu[n+j], tau), 0.01)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  K1 <- S1[1:2,1:2] * lambda[1]  + S1[1:2,3:4] * lambda[2]
  b[2:3] ~ dmnorm(zero[2:3],K1) 
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

  } else if (bs == "tp" & k == 4 & type == "leftcensored"){
    
    code <- '
model {
  mu <- X %*% b ## expected response
  for (i in 1:n) { y[i] ~ dnorm(mu[i], tau) } ## response 
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- max(pnorm(cut[j], mu[n+j], tau), 0.01)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  for (i in c(2:3)) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(4)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'
  } else if (bs == "tp" & k == 5 & type == "leftcensored"){
    
    code <- '
model {
  mu <- X %*% b ## expected response
  for (i in 1:n) { y[i] ~ dnorm(mu[i], tau) } ## response 
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- max(pnorm(cut[j], mu[n+j], tau), 0.01)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  for (i in c(2:4)) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(5)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

  } else {
    
    stop("The given combination of bs = ", sQuote(bs), ", k = ", k, ", and type = ", sQuote(type),
         " has not been implemented.")
    
  }

code

}

