
#
# Start up ----
#

library(devtools)
# library(ggplot2)
# library(mgcv)
# library(splines)
load_all()

#
# Checking and adding ..... -----
#

check()
use_package('dplyr')
use_package('ggplot2', 'suggests')
use_mit_license()
check()
use_package('ggplot2')
check()
use_package('splines')

# help(package = "leftcensored")

help(lc_linear)

#
# Simulate data and estimate regression ----
#
set.seed(10)
sim <- leftcensored_simulate(n = 30)
# debugonce(lc_linear)
result <- lc_linear(sim$data)

# Get best estimates and plot its regression line on top of the plot  
a <- result$intercept["50%"]
b <- result$slope["50%"]
abline(a, b, col = "green2")
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")


#
# With measurement error ----
#

set.seed(11)
sim <- leftcensored_simulate(n = 30)

result <- lc_linear_measerror(sim$data, measurement_error = 0.1)  

# Get best estimates and plot its regression line on top of the plot  
a <- result$intercept["50%"]
b <- result$slope["50%"]
abline(a, b, col = "green2")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")

#
# With measurement error, less censoring ----
#

set.seed(11)
sim <- leftcensored_simulate(n = 30, threshold_1 = 15, threshold_2 = 5)

# Without measurement error (green)
result <- lc_linear(sim$data)
a <- result$intercept["50%"]
b <- result$slope["50%"]
abline(a, b, col = "green2")
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")

# With measurement error (blue)
result <- lc_linear_measerror(sim$data, measurement_error = 0.1)
a <- result$intercept["50%"]
b <- result$slope["50%"]
abline(a, b, col = "blue3")
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "blue3")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "blue3")


#
# With measurement error, test: without error vs. with different error sizes ----
#

# Simulate data  
set.seed(6)
sim0 <- sim1 <- sim2 <- sim3 <- leftcensored_simulate(n = 30)

# sd_values   
sd_values <- c(0, 0.10, 0.25, 0.50)

# Make list of simulation objects
sim_list <- list(sim0, sim1, sim2, sim3)

# Estimate models
result_list <- list()
result_list[[1]] <- lc_linear(sim_list[[1]]$data)
for (i in 2:4){
  result_list[[i]] <- lc_linear_measerror(sim_list[[i]]$data, 
                                                measurement_error = sd_values[i],
                                                detailed = TRUE)
}

# Extract quantiles of intercept and slope  
quantiles <- purrr::map(result_list, c("summary", "quantiles"))
quantiles <- purrr::map(quantiles, as.data.frame)
names(quantiles) <- sprintf("SD = %.2f", sd_values)
interc <- purrr::map_dfr(quantiles, ~.["intercept",], .id = "Model")
slope <- purrr::map_dfr(quantiles, ~.["slope",], .id = "Model")

# Plot slope estimates
ggplot(slope, aes(Model, `50%`)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`))

# Plot graphs for 
for (i in 1:4){
  gg <- ggplot(sim_list[[i]]$data, aes(x, y, color = factor(uncensored))) +
    geom_point() +
    geom_abline(
      intercept = interc$`50%`[i],
      slope = slope$`50%`[i]
    ) +
    labs(title = names(quantiles)[i])
  print(gg)
}

#
# Check with detailed quantiles (including 'y_uncens' and 'uncensored') ----
#

# 1. ordinary lm  

set.seed(6)
sim <- leftcensored_simulate(n = 30)
result <- lc_linear(sim$data, detailed = TRUE)
# Check detailed quantiles
qs <- result$summary$quantiles
# rownames(qs)
sim$data[7:9,]
# plot(y~x, sim$data[7:9,])
qs[c("y_uncens[7]", "y_uncens[8]", "y_uncens[9]"),]
qs[c("uncensored[7]", "uncensored[8]", "uncensored[9]"),]

# 2. With measurement error 

set.seed(6)
sim <- leftcensored_simulate(n = 30)
sim$data$se_measurement <- 0.1*abs(sim$data$y)  
# sim$data$se_measurement[sim$data$uncensored == 0] <- 0.00001
result <- lc_linear_measerror(sim$data, sigma2 = 0.10, detailed = TRUE)

# Check detailed quantiles
qs <- result$summary$quantiles
# rownames(qs)
sim$data[7:9,]
# plot(y~x, sim$data[7:9,])
qs[c("y_uncens[7]", "y_uncens[8]", "y_uncens[9]"),]
qs[c("y_uncens_error[7]", "y_uncens_error[8]", "y_uncens_error[9]"),]
qs[c("uncensored[7]", "uncensored[8]", "uncensored[9]"),]




#
# leftcensored_prepare
#

# ?leftcensored_prepare

# Get actual data  
fn <- "../../seksjon 212/Milkys2_pc/Files_from_Jupyterhub_2020/Raw_data/109_adjusted_data_2021-09-15.rds"                         # data FROM Milkys2 on PC 
dat_all <- readRDS(fn)

# Tables for less-than's
dat_all %>%
  filter(MYEAR > 2012) %>%
  xtabs(~PARAM + addNA(FLAG1), .)

dat_all %>%
  filter(MYEAR > 2012 & substr(PARAM,1,3) == "BDE") %>%
  xtabs(~PARAM + addNA(FLAG1), .)

# Example of time series
dat1 <- dat_all %>%
  filter(PARAM == "BDE99" & TISSUE_NAME == "Lever" & STATION_CODE == "30B" & !is.na(VALUE_WW)) %>%
  mutate(Over_LOQ = !is.na(FLAG1)) %>%
  arrange(desc(Over_LOQ))

ggplot(dat1, aes(MYEAR, VALUE_WW, shape = Over_LOQ, color = Over_LOQ)) +
  geom_jitter(width = 0.1) +
  scale_y_log10()

# Prepare data
data_test <- leftcensored_prepare(dat1, var_year = "MYEAR", var_concentration = "VALUE_WW", var_LOQflag = "FLAG1")

ggplot(data_test, aes(x, y_uncens)) +
  geom_point(aes(y = threshold), shape = 1, size = 3, color = "red") +
  geom_point(aes(color = uncensored))

# Linear MCMC  
result <- lc_linear(data_test)

# MCMC summary
result$summary

# Spline MCMC  
result <- leftcensored_fixedsplines(data_test, x = "x", y = "y_uncens", uncensored = "uncensored", threshold = "threshold")

plot(y_uncens ~ x, data = data_test)
points(y_uncens ~ x, data = subset(data_test, uncensored == 0), pch = 20, col = "red")  
lines(y ~ x, data = result$plot_data, col = "red")
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "red")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "red")  

