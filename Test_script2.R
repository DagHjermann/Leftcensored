
#
# When starting testing
#
if (FALSE){
  library(devtools)
  load_all()
}

#
# Linear regression ----
#


# Simulate data and estimate regression 
set.seed(11)
sim <- lc_simulate(n = 30)   # also plots the data
# debugonce(lc_linear)
result <- lc_linear(sim$data)

# Get best estimate fitted line       
a <- result$intercept["50%"]
b <- result$slope["50%"]
# Add regression line to the plot  
abline(a, b, col = "green3")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green3")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green3")

# DIC
result$dic


#
# Linear regression, measurement error ----
#

#
# . Simplest case ----
#

# Simulate data and estimate regression 
set.seed(11)
sim <- lc_simulate(n = 30)   # also plots the data
sim$data$meas_error <- 5
result_me <- lc_linear(sim$data, measurement_error = "meas_error")

# Get best estimate fitted line       
a <- result_me$intercept["50%"]
b <- result_me$slope["50%"]
# Add regression line to the plot  
abline(a, b, col = "purple")
# Add confidence interval  
lines(y_lo ~ x, data = result_me$plot_data, lty = "dashed", col = "purple")
lines(y_hi ~ x, data = result_me$plot_data, lty = "dashed", col = "purple")

# DIC
result_me$dic

#
# . Additive measurement error ----
#

error_size <- c(0.1, 2, 5, 10)

# Function for setting measurement error and then estimate regression
get_result <- function(error_size){
  sim$data$meas_error <- error_size
  lc_linear(sim$data, measurement_error = "meas_error")
}

# Estimate regression for all error sizes  
result_me_list <- purrr::map(error_size, get_result)
names(result_me_list) <- paste("Error =", error_size)

# Plot results
lc_plot(sim$data, results = result_me_list)
lc_plot(sim$data, results = result_me_list, facet = "cols")


#
# . Measurement error given as percent, log-transformed data ----
#

# Data are log-transformed -> proportional error becomes additive
# Assume that error sd is proportional, e.g. 20% of the value
#   (errorfraction = 0.2)
# Log-transforming Y to Y' = log(Y) means that 
#   standard error becomes additive with sd = exp(errorfraction)-1
# See "Theory" below!


# Simulate data and estimate regression 
set.seed(11)
sim <- lc_simulate(n = 30, 
                   intercept = 3, slope = -0.2, sigma = 1, 
                   threshold_1 = 1.5, threshold_2 = 0.8, threshold_change = 4)
sim$data$y_orig <- exp(sim$data$y)
# plot(y_orig ~ x, data = sim$data)

# Assume 20% measurement error (fraction = 0.2)
sim$data <- sim$data %>%
  mutate(
    y_orig = exp(y),
    error_orig = y_orig*0.2,
    error_log = exp(0.2)-1
  )

# Plot on original scale
ggplot(sim$data, aes(x, y_orig)) +
  geom_pointrange(aes(ymin = y_orig - error_orig, 
                      ymax = y_orig + error_orig))

# Plot on log scale
ggplot(sim$data, aes(x, y)) +
  geom_pointrange(aes(ymin = y - error_log, 
                      ymax = y + error_log))

result_me <- lc_linear(sim$data, measurement_error = "error_log")

lc_plot(sim$data, results = result_me)


#
# .. same data, different percentages ----
#

get_data <- function(error_percent){
  sim$data <- sim$data %>%
    mutate(
      y_orig = exp(y),
      error_orig = y_orig*error_percent/100,
      error_log = exp(error_percent/100)-1
    )
}

library(purrr)
results_me <- data.frame(
  error_percent = c(0.1, 15, 45)) %>%
  mutate(
    data = map(error_percent, get_data),
  )

# Estimate regression for all error percentages  
result_me_list <- purrr::map(results_me$data, lc_linear, measurement_error = "error_log")
names(result_me_list) <- paste("Error =", results_me$error_percent)

# Plot results
# Narrower confidence intervals with higher error_percent??
lc_plot(sim$data, results = result_me_list, facet = "wrap")

# The data look as expected  
results_me %>%
  select(error_percent, data) %>%
  tidyr::unnest(cols = c(error_percent, data)) %>%
  ggplot(aes(x, y)) +
  geom_pointrange(aes(ymin = y - error_log, 
                      ymax = y + error_log)) +
  facet_wrap(vars(error_percent))




#
# Linear regression, actual data with proportional error (e.g. 20%) ----
#

# Get one station
data_test_orig <- subset(polybrom, station %in% "23B")

# Prepare data
# debugonce(lc_prepare)
data_test_prep <- lc_prepare(data_test_orig, 
                             x = "year",
                             y = "concentration", 
                             censored = "LOQ_flag",
                             log = TRUE,
                             keep_original_columns = TRUE)

# Plot
lc_plot(data_test_prep)



#
# . Try different percentages of measurement error ----
#

get_data <- function(error_percent, data = data_test_prep){
  data$meas_error <- exp(error_percent/100) - 1
  data
}
# get_data(15, data_test_prep)
# get_data(45, data_test_prep)

# X <- get_data(15, data_test_prep)
# ggplot(X, aes(x, y)) + 
#   geom_pointrange(aes(ymin = y-meas_error, ymax = y+meas_error))
# X <- get_data(45, data_test_prep)
# ggplot(X, aes(x, y)) + 
#   geom_pointrange(aes(ymin = y-meas_error, ymax = y+meas_error))

# Function for setting measurement error and then estimate regression
get_result <- function(data_for_analysis){
  lc_linear(data_for_analysis, measurement_error = "meas_error")
}

results_error <- data.frame(
  error_percent = c(0.1, 15, 45)) %>%
  mutate(
    data = map(error_percent, get_data),
    lc_result = map(data, lc_linear, measurement_error = "meas_error")
  )
  
  




# Estimate regression for all error sizes  
result_me_list <- purrr::map(error_percent, get_result)
names(result_me_list) <- paste("Error =", error_percent)

str(result_me_list, 1)
# Plot results
lc_plot(data_test_prep, results = result_me_list, facet = "wrap")

map_dfr(result_me_list, "slope")


result_ord <- lc_linear(data_test_prep)   
result_me <- lc_linear(data_test_prep, measurement_error = 0.2)   

lc_plot(data_test_prep)

lc_plot(data_test_prep, 
        results = list(Ordinary = result_ord,
                       Meas.error.20perc = result_me))



# Get best estimate fitted line       
a <- result_me$intercept["50%"]
b <- result_me$slope["50%"]
# Add regression line to the plot  
abline(a, b, col = "purple")
# Add confidence interval  
lines(y_lo ~ x, data = result_me$plot_data, lty = "dashed", col = "purple")
lines(y_hi ~ x, data = result_me$plot_data, lty = "dashed", col = "purple")

# DIC
result_me$dic



#
# Non-linear regression ----
#

#
# Test Qi version with simulated data ----   
#

#
# Simulate strongly non-linear data
#
set.seed(991)
X <- seq(from=-1, to=1, by=.025) # generating inputs
B <- t(splines::bs(X, knots=seq(-1,1,1), degree=3, intercept = TRUE)) # creating the B-splines
num_data <- length(X); num_basis <- nrow(B)
a0 <- 0.2 # intercept
a <- rnorm(num_basis, 0, 1) # coefficients of B-splines
n_param <- length(a)

Y_true <- as.vector(a0*X + a%*%B) # generating the output
Y <- Y_true + rnorm(length(X),0,.1) # adding noise

dat_sim <- data.frame(x = X, y_uncensored = Y, y_true = Y_true)
# ggplot(dat_sim, aes(x, y_uncensored)) +
#   geom_point() +
#   geom_line(aes(y = y_true), color = "blue")

# Add censoring 
dat_sim$y <- dat_sim$y_uncensored
dat_sim$uncensored <- 1
threshold_fixed <- -0.3
sel <- dat_sim$y_uncensored < threshold_fixed
dat_sim$y[sel] <- NA  
dat_sim$uncensored[sel] <- 0  
dat_sim$threshold <- threshold_fixed

# Plot
lc_plot(dat_sim)


#
# . fit models and test DIC values ----
#
# Does work as expected
#
# debugonce(lc_fixedsplines_qi)
result_linear <- lc_linear(dat_sim)
result_2knots <- lc_fixedsplines(dat_sim, knots = 2)
result_3knots <- lc_fixedsplines(dat_sim, knots = 3)
result_4knots <- lc_fixedsplines(dat_sim, knots = 4)
result_5knots <- lc_fixedsplines(dat_sim, knots = 5)
dic_values <- data.frame(
  Model = c("Linear", "2 knots", "3 knots", "4 knots", "5 knots"),
  DIC <- c(result_linear_qi$dic, result_2knots_qi$dic, result_3knots_qi$dic,
           result_4knots_qi$dic, result_5knots_qi$dic)
)
barplot(dic_values$DIC, names.arg = dic_values$Model)

#
# . plot data ---- 
#
#   true model (if existing), and fitted line of model(s)
#

lc_plot(dat_sim, y_true = "y_true", results = result_3knots_qi)
lc_plot(dat_sim, y_true = "y_true", results = list(result_3knots_qi))
lc_plot(dat_sim, 
        y_true = "y_true", 
        results = list(Linear = result_linear_qi,
                       Nonin_2knots = result_2knots_qi,
                       Nonlin_3knots = result_3knots_qi,
                       Nonlin_4knots = result_4knots_qi,
                       Nonlin_5knots = result_5knots_qi))


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Show that normalizing Y means that standard error is divided by sd(Y)
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Random variable with mean 100 and sd 10
n <- 5000
y <- rep(100, n)
err <- rnorm(n, 0, 10)
# Observed y (ym):
ym <- y + err 
sd(ym)

# Normalize 'ym'
ym_norm <- (ym-100)/10
mean(ym_norm)
sd(ym_norm)

# Recalculate normalized ym from y and err
ym_norm2 <- (y-100)/10 + err/10

# ym_norm and ym_norm2 seems exactly the same
mean(ym_norm2)
sd(ym_norm2)

# More 'proof':
plot(ym_norm[1:50])
points(ym_norm2[1:50], pch = 18, col = 2)

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Theory: assume that error sd is proportional, e.g. 20% of the value
#   (errorfraction = 0.2)
# Show that log-transforming Y to Y' = log(Y) means that 
#   standard error becomes additive with sd = exp(errorfraction)-1
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Random variable with mean 100 and sd 10
n <- 10000
y_mean <- rep(100, n)
sigma <- rnorm(n, 0, 10)
# Observed y (ym):
y <- y_mean + sigma 
sd(y)

logy <- log(y)
sd(logy)
log(sd(y))


# Random variable with mean 100 and sd 10
y_actual <- runif(1000, 20, 200)
errorfraction = 0.2
sd <- y_actual*errorfraction
error <- rnorm(1000, 0, sd)
y_obs <- y_actual + error

# Plot SD relative to value
plot(y_actual, sd)
plot(y_obs, sd)

# Plot some points withiut and with error
plot(head(y_actual, 50))
points(head(y_obs, 50), pch = 18, col = 2)

# Plot all 
plot(y_actual, y_obs)

# Check some large numbers
sel <- y_obs > 150
sd(y_obs[sel] - y_actual[sel])
175*0.2  # the "expected", quite close

# Check some smaller numbers
sel <- y_obs < 40
sd(y_obs[sel] - y_actual[sel])
30*0.2  # the "expected", quite close

# Log transform  
log_y_actual <- log(y_actual)
log_y_obs <- log(y_obs)

# Plot some points withiut and with error 
plot(head(log_y_actual, 50))
points(head(log_y_obs, 50), pch = 18, col = 2)

# Plot all 
plot(log_y_actual, log_y_obs)

# Expected additive approximate errors
error_logscale <- rnorm(1000, 0, exp(errorfraction) - 1)
log_y_obs2 <- log_y_actual + error_logscale

# Plot all 
points(log_y_actual, log_y_obs2, pch = 18, col = 2)

# Check sigma from variation around line - very close  
mod1 <- lm(log_y_obs ~ log_y_actual)
summary(mod1)$sigma
mod2 <- lm(log_y_obs2 ~ log_y_actual)
summary(mod2)$sigma
# [1] 0.2172967
# [1] 0.2218082

