#' Get quantiles for y (given x), mean quantiles for all chains
#' 
#' @param x x value
#' @param mcmc MCMC object
#' @param quantiles Quantiles. Default is c(0.025, 0.25, 0.5, 0.75, 0.975)
#' 
#' @return Vector of (with default quantiles) 5 values, containing quantiles for the given x value 
#' 
#' @examples
#' sim <- leftcensored_simulate()
#' result <- lm_linear(sim$data)
#' predicted_y_quantiles(2000, result$model)
#' 
#' @export
predicted_y_quantiles <- function(x, mcmc, 
                                  quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)){
  # Get 
  M <- 1:length(mcmc) %>%
    purrr::map(~predicted_y_quantiles_onechain(x, result$model, ., quantiles = quantiles)) %>% 
    bind_cols() %>% t()
  apply(M, 2, mean)
}

