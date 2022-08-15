
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


