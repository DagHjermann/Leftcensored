
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
set.seed(11)
sim <- lc_simulate(n = 30)
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
sim <- lc_simulate(n = 30)

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
sim <- lc_simulate(n = 30, threshold_1 = 15, threshold_2 = 5)

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
sim0 <- sim1 <- sim2 <- sim3 <- lc_simulate(n = 30)

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

names(result_list) <- sprintf("SD = %.2f", sd_values)

# Extract quantiles of intercept and slope  
intercepts <- purrr::map_dfr(result_list, c("intercept"), .id = "Model")
slopes <- purrr::map_dfr(result_list, c("slope"), .id = "Model")

# Plot slope estimates
ggplot(slopes, aes(Model, `50%`)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`))

# Plot graphs for 
for (i in 1:4){
  gg <- ggplot(sim_list[[i]]$data, aes(x, y_plot, color = factor(uncensored))) +
    geom_point() +
    geom_abline(
      intercept = intercepts$`50%`[i],
      slope = slopes$`50%`[i]
    ) +
    labs(title = names(result_list)[i])
  print(gg)
}


#
# lc_prepare ----
#

# ?lc_prepare

# Get actual data  
fn <- "../../seksjon 212/Milkys2_pc/Files_from_Jupyterhub_2020/Raw_data/109_adjusted_data_2021-09-15.rds"                         # data FROM Milkys2 on PC 
dat_all <- readRDS(fn)

# Tables for less-than's
dat_all %>%
  filter(MYEAR > 2012) %>%
  xtabs(~PARAM + addNA(FLAG1), .)

# as above, BDEs only 
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

# Plot  
ggplot(data_test, aes(x, y_uncens)) +
  geom_point(aes(y = threshold), shape = 1, size = 3, color = "red") +
  geom_point(aes(color = factor(uncensored)))

# Linear MCMC  
result <- lc_linear(data_test)  

# Plot with regression line  
a <- result$intercept["50%"]
b <- result$slope["50%"]
ggplot(data_test, aes(x, y_uncens)) +
  geom_point(aes(y = threshold), shape = 1, size = 3, color = "red") +
  geom_point(aes(color = factor(uncensored))) +
  geom_abline(intercept = result$intercept["50%"], slope = result$slope["50%"])

# Plot with regression line and confidence interval   
a <- result$intercept["50%"]
b <- result$slope["50%"]
ggplot(data_test, aes(x = x)) +
  geom_ribbon(data = result$plot_data, aes(ymin = y_lo, ymax = y_hi), fill = "grey80") + 
  geom_point(aes(y = threshold), shape = 1, size = 3, color = "red") +
  geom_point(aes(y = y_uncens, color = factor(uncensored))) +
  geom_abline(intercept = result$intercept["50%"], slope = result$slope["50%"])

# Slope from lc_linear
result$slope

# Slope from linear regression, ignoring "<"  
ols <- summary(lm(y_uncens ~ x, data = data_test))$coef["x",]

# Compare
df_slopes <- bind_rows(
  data.frame(
    Analysis = "lc_linear",
    slope_min = result$slope["2.5%"],
    slope = result$slope["50%"],
    slope_max = result$slope["97.5%"]),
  data.frame(
    Analysis = "lm",
    slope_min = ols["Estimate"] - 2*ols["Std. Error"],
    slope = ols["Estimate"],
    slope_max = ols["Estimate"] + 2*ols["Std. Error"])
)
ggplot(df_slopes, aes(x = Analysis)) +
  geom_pointrange(aes(y = slope, ymin = slope_min, ymax = slope_max))


#
# Spline MCMC  
#
result <- leftcensored_fixedsplines(data_test, x = "x", y = "y_uncens", uncensored = "uncensored", threshold = "threshold")

plot(y_uncens ~ x, data = data_test)
points(y_uncens ~ x, data = subset(data_test, uncensored == 0), pch = 20, col = "red")  
lines(y ~ x, data = result$plot_data, col = "red")
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "red")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "red")  

