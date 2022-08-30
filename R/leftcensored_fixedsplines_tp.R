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
                               normalize = TRUE,
                               make_data_only = FALSE,
                               initialize_only = FALSE,
                               measurement_error = NULL){
  
  # Rename variables, reorder data and add 'y_comb'    
  dat_ordered1 <- get_dat_ordered1(
    data = dat_test, x = "xx", y = "yx", uncensored = "uncensoredx", threshold = "cutx"
  )
  
  # Make additional data that will be used only to make the fited line,
  #  and add them to the data  
  dat_ordered2 <- get_dat_ordered2(dat_ordered1, fit_length = resolution)
  
  
  # Normalize data
  # Achieves mean = 0
  if (normalize){
    
    # Standardization to mean = 0 and sd = 1 
    # mean_x <- mean(dat_ordered2$x, na.rm = TRUE)
    # sd_x <- sd(dat_ordered2$x, na.rm = TRUE)
    # norm_x <- function(x) (x-mean_x)/sd_x
    # unnorm_x <- function(x) x*sd_x + mean_x
    
    # scale
    scale <- 10
    
    # x: standardization to min = 0 and max = scale 
    min_x <- min(dat_ordered2$x, na.rm = TRUE)
    max_x <- max(dat_ordered2$x, na.rm = TRUE)
    norm_x <- function(x) (x-min_x)/(max_x-min_x)*scale
    unnorm_x <- function(x) x*(max_x-min_x) + min_x/scale
    
    # y: standardization to max = scale 
    max_y <- max(dat_ordered2$y, na.rm = TRUE)
    norm_y <- function(y) y/max_y*scale
    unnorm_y <- function(y) y*max_y/scale
    
    # Normalize x, y and y_comb   
    # norm <- normalize_lm(dat_ordered2$x, c(data_obs$y_comb, data_cen$cut))
    dat_ordered2$x <- norm_x(dat_ordered2$x)
    dat_ordered2$y <- norm_y(dat_ordered2$y)
    dat_ordered2$y_comb <- norm_y(dat_ordered2$y_comb)
    dat_ordered2$cut <- norm_y(dat_ordered2$cut)
    
  }
  
  jd_tp5_lc <- get_jagam_object(dat_ordered2, k = k)
  
  if (make_data_only){
    
    result <- list(
      data_for_analysis = dat_ordered2,
      jagam_object= jd_tp5_lc)
    
    
  } else if (!make_data_only & initialize_only){
    
    # Runs rjags::jags.model
    # Fast check of whther the code is working  
    
    # Get JAGS code
    jagscode_txt <-  get_jags_model_code(bs = "tp", k = k, type = "leftcensored")
    
    jm <- rjags::jags.model(textConnection(jagscode_txt), 
                            data=jd_tp5_lc$jags.data, 
                            inits=jd_tp5_lc$jags.ini, 
                            n.chains=2)  # changed from n.chains=1 
    
    result <- list(
      data_for_analysis = dat_ordered2,
      jagam_object= jd_tp5_lc,
      jagscode_txt = jagscode_txt)

  } else if (!make_data_only & !initialize_only){
    
    # Get JAGS code
    jagscode_txt <-  get_jags_model_code(bs = "tp", k = 5, type = "leftcensored")
    
    form <- as.formula(paste0("y_comb ~ s(x, bs='tp', k=", k, ")"))
    
    # Makes (1) jags.file_tp5_orig (was used as basis for "_leftcens" file)
    # Make (2) jagam_object$jags.data (will be manipulated below)
    jagam_object <- mgcv::jagam(gam_formula, 
                                data = dat_ordered2,    # file no. 2 here
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
    
  }
  
  result
  
}






get_jagam_object <- function(data, k = 5){
  
  jags.file <- paste0(tempdir(), "/temporary.jags") 
  
  gam_formula <- as.formula(paste0("y_comb ~ s(x, bs='tp', k=", k, ")"))

  # Make (1) jags.file_tp5_orig (was used as basis for "_leftcens" file)
  # Make (2) jagam_object$jags.data (will be manpulated below)
  jagam_object <- jagam(gam_formula, 
                        data = data,    # file no. 2 here
                        file = jags.file,   # this file will be overwritten (was used as basis for "_leftcens" file)
                        sp.prior = "gamma", 
                        diagonalize = TRUE)
  
  # Modify jags.data object
  n <- sum(dat_ordered1$uncensored %in% 1)   # file no. 1 here 
  m <- sum(dat_ordered1$uncensored %in% 0)   # file no. 1 here
  jagam_object$jags.data$n <- n   # - makes sure only these data are use for the likelihood
  jagam_object$jags.data$m <- m   #    - " -
  jagam_object$jags.data$Z <- c(rep(0,n), rep(1, m))
  jagam_object$jags.data$cut <- dat_ordered1$cut[dat_ordered1$uncensored %in% 0]
  
  jagam_object
  
}

