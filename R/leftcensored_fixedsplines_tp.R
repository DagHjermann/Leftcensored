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
#' @return A list including \code{summary} (summary of the JAGS coda object), \code{plot_data} (x and y 
#' data for plotting the fitted model), and \code{dic} (the DIC value).  
#' 
#' @export
#'
#' @examples
#' 
#' 
#' 
#' # Simulate data ----
#' set.seed(2) ## simulate some data... 
#' n <- 50
#' dat <- mgcv::gamSim(1,n=n,dist="normal",scale=1)  
#' 
#' # we will use only x2 and y, and x2 is renamed 'x'
#' dat <- dat[c("x2", "y")]
#' names(dat)[1] <- "x"
#' 
#' # Plot original data  
#' ggplot(dat, aes(x, y)) +
#'   geom_point()
#' 
#' # Make censored data (here, using a fixed thresholdm, but that is not necessary)
#' thresh <- 4
#' dat_cens <- dat[c("x","y")]
#' dat_cens$y_orig <- dat_cens$y       # original (will not be used)
#' sel_uncens <- dat_cens$y > thresh
#' dat_cens$y[!sel_uncens] <- NA
#' dat_cens$cut <- thresh
#' dat_cens$cut[sel_uncens] <- NA
#' dat_cens$uncensored <- 0
#' dat_cens$uncensored[sel_uncens] <- 1
#' 
#' # Plot censored data
#' ggplot() +
#'   geom_point(data = dat_cens[sel_uncens,], aes(x = x, y = y)) +
#'   geom_point(data = dat_cens[!sel_uncens,], aes(x = x, y = cut), shape = 6)
#' 
#' # Quick test that JAGS runs
#' # debugonce(lc_fixedsplines_tp)
#' test <- lc_fixedsplines_tp(data = dat_cens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
#'                            normalize = TRUE, k = 3, initialize_only = TRUE)
#' 
#' # Add measurement error (here: a fixed value for all observations)
#' dat_cens$error <- 1
#' 
#' # Plot
#' ggplot() +
#'   geom_pointrange(data = dat_cens[sel_uncens,], aes(x = x, y = y, ymin = y-error, ymax = y+error)) +
#'   geom_point(data = dat_cens[!sel_uncens,], aes(x = x, y = cut), shape = 6)
#' 
#' # Quick test that JAGS runs
#' test <- lc_fixedsplines_tp(data = dat_cens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
#'                            measurement_error = "error", 
#'                            normalize = TRUE, k = 3, initialize_only = TRUE)
#' 
#'
lc_fixedsplines_tp <- function(data,
                               x = "x", 
                               y = "y", 
                               uncensored = "uncensored",
                               threshold = "threshold",
                               k = 5,
                               predict_x = 50,
                               n.chains = 2, 
                               n.iter = 4000, 
                               n.burnin = 4000, 
                               n.thin = 5,
                               model_parameters_for_convergence =  c("b","rho","lambda","scale"),
                               normalize = TRUE,
                               make_data_only = FALSE,
                               initialize_only = FALSE,
                               raftery = TRUE,
                               max.time = "2 minutes",
                               keep_jags_model = FALSE,
                               keep_mcmc_model = FALSE,
                               measurement_error = NULL){
  
  if (k <= 2){
    k_jagam <- 3
  } else {
    k_jagam <- k
  }
  
  # Rename variables, reorder data and add 'y_comb'    
  dat_ordered1 <- get_ordered_data1(
    data = data, x = x, y = y, uncensored = uncensored, threshold = threshold,
    measurement_error = measurement_error
  )
  
  # Make additional data that will be used only to make the fited line,
  #  and add them to the data  
  dat_ordered2_list <- get_ordered_data2(dat_ordered1, added_x_values = predict_x)
  dat_ordered2 <- dat_ordered2_list$data
  
  #
  # For standardiztion
  #
  
  # scale
  scale <- 10
  
  # x: standardization to min = 0 and max = scale 
  min_x <- min(dat_ordered2$x, na.rm = TRUE)
  max_x <- max(dat_ordered2$x, na.rm = TRUE)
  norm_x <- function(x) (x-min_x)/(max_x-min_x)*scale
  denorm_x <- function(x) x*(max_x-min_x) + min_x/scale
  
  # y: standardization to max = scale 
  max_y <- max(dat_ordered2$y, na.rm = TRUE)
  norm_y <- function(y) y/max_y*scale
  denorm_y <- function(y) y*max_y/scale
  
  # Normalize data
  # Achieves mean = 0
  if (normalize){
    
    # Normalize x, y and y_comb   
    # norm <- normalize_lm(dat_ordered2$x, c(data_obs$y_comb, data_cen$cut))
    dat_ordered2$x <- norm_x(dat_ordered2$x)
    dat_ordered2$y <- norm_y(dat_ordered2$y)
    dat_ordered2$y_comb <- norm_y(dat_ordered2$y_comb)
    dat_ordered2$cut <- norm_y(dat_ordered2$cut)
    if (!is.null(measurement_error))
      dat_ordered2$meas_error <- norm_y(dat_ordered2$meas_error)

  }
  
  jagam_object <- leftcensored:::get_jagam_object(dat_ordered1, dat_ordered2, 
                                   k_jagam = k_jagam, 
                                   measurement_error = measurement_error,
                                   k_orig = k)
  
  if (make_data_only){
    
    result <- list(
      data_for_analysis = dat_ordered2,
      jagam_object = jagam_object)
    
  } else {
    
    model_meas_error <- !is.null(measurement_error)
    model_leftcensored <- mean(dat_ordered1$uncensored) < 1
    
    # Get JAGS code
    if (!model_meas_error & model_leftcensored){
      jagscode_txt <-  leftcensored:::get_jags_model_code(bs = "tp", k_code = k, type = "leftcensored")
    } else if (model_meas_error & model_leftcensored){
      jagscode_txt <-  leftcensored:::get_jags_model_code(bs = "tp", k_code = k, type = "leftcensored_measerror")
    } else if (!model_meas_error & !model_leftcensored){
      jagscode_txt <-  leftcensored:::get_jags_model_code(bs = "tp", k_code = k, type = "uncensored")
    } else if (model_meas_error & !model_leftcensored){
      jagscode_txt <-  leftcensored:::get_jags_model_code(bs = "tp", k_code = k, type = "uncensored_measerror")
    }
    
  }
    
  if (!make_data_only & initialize_only){
    
    # Runs rjags::jags.model
    # Fast check of whther the code is working  
    
    jm <- rjags::jags.model(textConnection(jagscode_txt), 
                            data=jagam_object$jags.data, 
                            inits=jagam_object$jags.ini, 
                            n.chains=2)  # changed from n.chains=1   
    
    result <- list(
      data_for_analysis = dat_ordered2,
      jagam_object= jagam_object,
      jagscode_txt = jagscode_txt)

  } else if (!make_data_only & !initialize_only){
    
    # Choose the parameters to watch
    # Sample varaibles that have been inserted only to get the fitted line
    mu_fitted_names1 <- paste0("mu[", nrow(dat_ordered1)+1, ":", nrow(dat_ordered2), "]")
    mu_fitted_names2 <- paste0("mu[", seq(nrow(dat_ordered1)+1, nrow(dat_ordered2)), "]")
    
    if (raftery){
      raftery.options <- list()
    } else {
      raftery.options <- FALSE
    }
    
    ### Run model
    # Initial run  
    model_converged <- runjags::autorun.jags(
      data = jagam_object$jags.data,
      monitor = model_parameters_for_convergence,     
      inits = jagam_object$jags.ini,
      model = jagscode_txt,
      n.chains = n.chains,    # Number of different starting positions
      startsample = n.iter,   # Number of iterations
      startburnin = n.burnin, # Number of iterations to remove at start
      thin = n.thin,          # Amount of thinning
      raftery.options = raftery.options,
      max.time = max.time)
    
    # DIC
    dic_list <- get_dic(model_converged)
    
    # Add all model parameters and get samples for them
    model_result <- runjags::extend.jags(model_converged, 
                                         add.monitor = mu_fitted_names1,
                                         sample = 2000)
    
    # model_result
    model_mcmc <- coda::as.mcmc(model_result)
    
    summary <- summary(model_mcmc)
    quants <- summary$quantiles
    pick_rownames <- rownames(quants) %in% mu_fitted_names2
    
    if (normalize){
      # Denormalize predicted data
      y_med <- denorm_y(quants[pick_rownames,"50%"])
      y_lo <- denorm_y(quants[pick_rownames,"2.5%"])
      y_hi <- denorm_y(quants[pick_rownames,"97.5%"])
    } else {
      y_med <- quants[pick_rownames,"50%"]
      y_lo <- quants[pick_rownames,"2.5%"]
      y_hi <- quants[pick_rownames,"97.5%"]
    }
    
    # Make 'plot_data'  
    # y and lower and upper CI  values are back-transformed (un-normalized) using denorm:
    plot_data <- data.frame(
      x = dat_ordered2_list$data_for_fit$x, 
      y = y_med,
      y_lo = y_lo,
      y_hi = y_hi  
      )

    result <- list(summary = summary,
                   plot_data = plot_data,
                   dic = dic_list$dic,
                   deviance = dic_list$deviance,
                   pd = dic_list$pd,
                   model_data = jagam_object$jags.data,
                   dat_ordered1 = dat_ordered1,
                   dat_ordered2 = dat_ordered2)  
    
    if (keep_jags_model)
      result = append(result, list(model = model_result))
    if (keep_mcmc_model)
      result = append(result, list(model_from_jags = model_mcmc))
    
  }
  
  result
  
}




get_jagam_object <- function(data_ordered1, data_ordered2, k_jagam = 5, measurement_error, 
                             k_orig = k_orig){
  
  jags.file <- paste0(tempdir(), "/temporary.jags") 
  
  gam_formula <- as.formula(paste0("y_comb ~ s(x, bs='tp', k=", k_jagam, ")"))

  # Make (1) jags.file_tp5_orig (was used as basis for "_leftcens" file)
  # Make (2) jagam_object$jags.data (will be manpulated below)
  jagam_object <- mgcv::jagam(gam_formula, 
                        data = data_ordered2,    # file no. 2 here
                        file = jags.file,   # this file will be overwritten (was used as basis for "_leftcens" file)
                        sp.prior = "gamma", 
                        diagonalize = TRUE)
  
  # Modify jags.data object
  
  n <- sum(data_ordered1$uncensored %in% 1)   # file no. 1 here 
  m <- sum(data_ordered1$uncensored %in% 0)   # file no. 1 here
  
  jagam_object$jags.data$n <- n   # - makes sure only these data are use for the likelihood
  
  if (m > 0){
    jagam_object$jags.data$m <- m   #    - " -
    jagam_object$jags.data$Z <- c(rep(0,n), rep(1, m))
    jagam_object$jags.data$cut <- data_ordered1$cut[data_ordered1$uncensored %in% 0]
  }
  
  if (!is.null(measurement_error))
    jagam_object$jags.data$meas_error <- data_ordered1$meas_error[data_ordered1$uncensored %in% 1]
  
  if (k_orig == 1){
    jagam_object$jags.data$b2 <- 0
    jagam_object$jags.data$zero <- NULL
    jagam_object$jags.data$S1 <- NULL
    jagam_object$jags.ini$b <- jagam_object$jags.ini$b[1]
  }
  
  if (k_orig == 2){
    jagam_object$jags.data$b2 <- 0
    jagam_object$jags.ini$b <- c(jagam_object$jags.ini$b[1], jagam_object$jags.ini$b[3])
  }
  
  if (k_orig >= 4){
    jagam_object$jags.data$k <- k_orig
  }
  
  jagam_object
  
}

