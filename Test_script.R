
#
# Start up ----
#

library(devtools)
# library(ggplot2)
# library(mgcv)
# library(splines)
load_all()

#
# Checking, adding packages ..... -----
#

check()
use_package('dplyr')
use_package('ggplot2', 'suggests')
use_mit_license()
check()
use_package('ggplot2')
check()
use_package('splines')

# Installing package  
install()

# help(package = "leftcensored")

help(lc_linear)

#
# Simulate data and estimate regression ----
#
set.seed(11)
sim <- lc_simulate(n = 30)   # also plots the data
# debugonce(lc_linear)
result <- lc_linear(sim$data)

# Get best estimate fitted line       
a <- result$intercept["50%"]
b <- result$slope["50%"]
# Add regression line to the plot  
abline(a, b, col = "green2")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")

# The result is a list of 9 things:  
str(result, 1)
str(result, 2)
str(result$summary, 1)
head(result$summary$statistics)
head(result$summary$quantiles)
rownames(result$summary$quantiles)

pick_rownames <- sprintf("y.hat.out[%i]", 1:length.out)

#' # Make a standard MCMC plot: the trace and the density for each estimated parameter  
#' par(mar = c(2,4,3,1))
#' plot(result$model)

#
# With measurement error ----
#

set.seed(11)
sim <- lc_simulate(n = 30)

result <- lc_linear_measerror(sim$data, measurement_error = 0.1)  
# result <- lc_linear_measerror(sim$data, measurement_error = 0.1, detailed = TRUE)  

# Get best estimates and plot its regression line on top of the plot  
a <- result$intercept["50%"]
b <- result$slope["50%"]
abline(a, b, col = "green2")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")

# Truncae distribution  
y_min <- 0
(y_min - result$mean_y)/result$sd_y

# Test plot
plot(result$model_data$x, result$model_data$y.uncens.error, 
     ylim = range(result$model_data$y.uncens.error, result$model_data$threshold, na.rm = TRUE))
sel <- is.na(result$model_data$y.uncens.error)
points(result$model_data$x[sel], result$model_data$threshold[sel], pch = 18, col = "red")
  

set.seed(11)
sim <- lc_simulate(n = )
head(sim$data)

# Normalized y, centralized x - no data under threshold   
set.seed(22)
sim <- lc_simulate(x = seq(-4,4,length=50), intercept = 0, slope = -0.25, sigma = 0.75, 
                   threshold_1 = -3, threshold_2 = -3, threshold_change = 0)
mean(sim$data$y_real, na.rm = TRUE)
sd(sim$data$y_real, na.rm = TRUE)
# - analysis   
result <- lc_linear_measerror(sim$data, measurement_error = 0.1)  
a <- result$intercept["50%"]
b <- result$slope["50%"]
abline(a, b, col = "green2")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")

sim_est <- function(..., mean_y = NULL, sd_y = NULL, seed = 22){
  set.seed(seed)
  sim <- lc_simulate(...)  
  mean(sim$data$x, na.rm = TRUE) %>% cat("mean x:", ., "\n")
  mean(sim$data$y_real, na.rm = TRUE) %>% cat("mean y:", ., "\n")
  sd(sim$data$y_real, na.rm = TRUE) %>% cat("sd y:", ., "\n")
  # - analysis   
  result <- lc_linear_measerror(sim$data, measurement_error = 0.1, mean_y = mean_y, sd_y = sd_y)  
  a <- result$intercept["50%"]
  b <- result$slope["50%"]
  abline(a, b, col = "green2")
  # Add confidence interval  
  lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
  lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")
  invisible(result)
}

# 0a. Normalized y, centralized x - no data under threshold   
X <- sim_est(x = seq(-4,4,length=50), intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -1000, threshold_2 = -1000, threshold_change = 0)

# 0b. As 0a, lower N     
X <- sim_est(x = seq(-4,4,length=30), intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -1000, threshold_2 = -1000, threshold_change = 0)

# run for k samples with n = 30
if (FALSE){
  seeds <- round(runif(n = 2, min = 1, max = 10000))
  X_list <- seeds %>%
    purrr::map(~sim_est(x = seq(-4,4,length=30), intercept = 0, slope = -0.25, sigma = 0.75, 
                 threshold_1 = -1000, threshold_2 = -1000, threshold_change = 0, seed = .)
    )
}

# 1. normalized y, centralized x - some data under threshold   
X <- sim_est(x = seq(-4,4,length=50), intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -0.1, threshold_2 = -0.5, threshold_change = 0)

# 2a. as 1 but lower N  
X <- sim_est(x = seq(-4,4,length=30), intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -0.1, threshold_2 = -0.5, threshold_change = 0, seed = 40)

# 2b. as 2a (lower N) but mean_y and sd_y are given    
X <- sim_est(x = seq(-4,4,length=30), intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -0.1, threshold_2 = -0.5, threshold_change = 0, mean_y = 0, sd_y = 0)

# 3. As 1 but "un-centralize" x (adding 10)       
# - threshold parameters also adjusted
# - works fine  
X <- sim_est(x = seq(-4,4,length=50) + 10, intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -2.7, threshold_2 = -3.2, threshold_change = 10)

# 4. As 3 but also "un-centralize" y (adding 10)       
# - threshold parameters also adjusted
# - still works fine  
X <- sim_est(x = seq(-4,4,length=50) + 10, intercept = 10, slope = -0.25, sigma = 0.75, 
             threshold_1 = -2.7 + 10, threshold_2 = -3.2 + 10, threshold_change = 10)

# 5. As 4 but also "un-normalize" y (adjusting slope)       
# - threshold parameters also adjusted
# - still works fine  
X <- sim_est(x = seq(-4,4,length=50) + 10, intercept = 10, slope = -0.55, sigma = 0.75, 
             threshold_1 = 4.7, threshold_2 = 3.2, threshold_change = 10)

# 6. As 5 but lower sample size      
# - threshold parameters also adjusted
# - still works fine  
X <- sim_est(x = seq(-4,4,length=30) + 10, intercept = 10, slope = -0.55, sigma = 0.75, 
             threshold_1 = 4.7, threshold_2 = 3.2, threshold_change = 10)

# 7. As 6 but "un-centralize" y       
X <- sim_est(x = seq(-4,4,length=50), intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -0.1, threshold_2 = -0.5, threshold_change = 0)


result <- lc_linear_measerror(sim$data, measurement_error = 0.1, detailed = TRUE)  


rownames(result$summary$quantiles)

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
# Add example data ----
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

# One BDE only, stations 
param <- "BDE99"
tab <- dat_all %>%
  filter(MYEAR > 2012 & PARAM == param) %>%
  xtabs(~STATION_CODE + addNA(FLAG1), .)
# Stations with total sample size >= 80 and at least 20 less-thans: 
tab2 <- tab[apply(tab, 1, sum) >= 50 & tab[,1] >= 15,]
tab2

# One BDE only, data for plot
dat1 <- dat_all %>%
  filter(PARAM == param & TISSUE_NAME == "Lever",
         STATION_CODE %in% rownames(tab2) & !is.na(VALUE_WW)) %>%
  mutate(Over_LOQ = is.na(FLAG1)) %>%
  arrange(desc(Over_LOQ)) 

# One BDE only, plot
ggplot(dat1, aes(MYEAR, VALUE_WW, shape = Over_LOQ, color = Over_LOQ)) +
  geom_jitter(width = 0.1) +
  scale_y_log10() +
  facet_wrap(vars(STATION_CODE))+
  labs(title = param)

if (FALSE){
  
  # Make example dataset
  polybrom <- dat1 %>% 
    filter(STATION_CODE %in% c("13B", "23B", "53B")) %>%
    rename(
      station = STATION_CODE,
      year = MYEAR,
      concentration = VALUE_WW,
      concentr_flag = FLAG1,
    ) %>%
    arrange(station, year) %>%
    select(station, year, concentration, concentr_flag)
  
  # Add to package  
  usethis::use_data(polybrom)

}



#
# lc_prepare ----
#


xtabs(~station + addNA(concentr_flag), polybrom)

# Prepare data
data_test_orig <- subset(polybrom, station %in% "53B")

# debugonce(lc_prepare)
data_test_prep <- lc_prepare(data_test_orig)  
# data_test_prep <- data_test_prep %>% filter(uncensored == 1)  

# CHECK 1
check1 <- with(data_test_prep, y < threshold)
sum(check1)

# CHECK 2
check2 <- with(data_test_prep, y == threshold)
sum(check2)
data_test_prep[check2,]

# FIX CHECK 2
data_test_prep$threshold[check2] <- data_test_prep$threshold[check2] - 1

# NOW IT WORKS
result <- lc_linear(data_test_prep, plot_input = TRUE, plot_norm = TRUE)  

# Check

str(result, 1)
str(result$model_data, 1)

# Test plot

# - plot data
sel_uncens <- !data_test_orig$concentr_flag %in% "<" 
plot(log(concentration) ~ year, 
     data = data_test_orig[sel_uncens,], 
     ylim = range(log(data_test_orig$concentration), na.rm = TRUE),
     pch = 16, col = 4)
points(log(concentration) ~ year, 
     data = data_test_orig[!sel_uncens,], 
     pch = 6, col = 2)
# Get best estimate fitted line       
a <- result$intercept["50%"]
b <- result$slope["50%"]
# Add regression line to the plot  
abline(a, b, col = "green2")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")





sel_uncens <- result$model_data$uncensored == 1
plot(result$model_data$x[sel_uncens], result$model_data$x[sel_uncens], )

ggplot(data_test_orig, aes(year, concentration, color = is.na(concentr_flag))) +
  geom_point()  



ggplot(data_test_orig, aes(year, log(concentration), color = is.na(concentr_flag))) +
  geom_point()


# Plot  

ggplot() +
  geom_jitter(
    data = data_test_prep %>% filter(uncensored == 0), 
    aes(x, threshold), width = 0.3
    ) +
  geom_jitter(
    data = data_test_prep %>% filter(uncensored == 1), 
    aes(x, y), width = 0.3, color = "red"
  )
  
# Linear MCMC  

debugonce(lc_linear)
# debugonce(R2jags::jags)
result <- lc_linear(data_test_prep, plot_input = TRUE, plot_norm = TRUE)  

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

