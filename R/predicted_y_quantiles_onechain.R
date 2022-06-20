#' Get quantiles for y (given x), mean quantiles for all chains
#' 
#' @param x x value
#' @param mcmc MCMC object
#' @param chain_no Chain number
#' @param quantiles Quantiles. Default is c(0.025, 0.25, 0.5, 0.75, 0.975)
#' 
#' @return Matrix of one column and (with default quantiles) 5 rows, containing quantiles for the given x value 
#' 
#' @examples
#' sim <- lc_simulate()
#' result <- lm_linear(sim$data)
#' predicted_y_quantiles_onechain(2000, result$model, 1)
#' 
#' @export
predicted_y_quantiles_onechain <- function(x, mcmc, chain_no,
                                           quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)){
  a <- mcmc[[chain_no]][,"intercept"]
  b <- mcmc[[chain_no]][,"slope"]
  quantile(a + b*x, quantiles) %>% as.matrix(nrow = 1)
}
