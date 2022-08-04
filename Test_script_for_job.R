
library(dplyr)
library(dplyr)
# devtools::load_all()

# Simulate and estimate  
sim_est <- function(..., mean_y = NULL, sd_y = NULL, seed = 22){
  set.seed(seed)
  sim <- lc_simulate(...)  
  # mean(sim$data$x, na.rm = TRUE) %>% cat("mean x:", ., "\n")
  # mean(sim$data$y_real, na.rm = TRUE) %>% cat("mean y:", ., "\n")
  # sd(sim$data$y_real, na.rm = TRUE) %>% cat("sd y:", ., "\n")
  # - analysis   
  result <- lc_linear_measerror(sim$data, measurement_error = 0.1, mean_y = mean_y, sd_y = sd_y)  
  a <- result$intercept["50%"]
  b <- result$slope["50%"]
  # abline(a, b, col = "green2")
  # Add confidence interval  
  # lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
  # lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")
  invisible(result)
}

# run for k samples with n = 30
if (TRUE){
  seeds <- round(runif(n = 2, min = 1, max = 10000))
  X_list <- seeds %>%
    purrr::map(~sim_est(x = seq(-4,4,length=30), intercept = 0, slope = -0.25, sigma = 0.75, 
                        threshold_1 = -1000, threshold_2 = -1000, threshold_change = 0, 
                        plot = FALSE,
                        seed = .)
    )
}
