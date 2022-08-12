#' Normalization and removal of linear trend  
#' 
#' This function takes two variables x and y as an input, and removes the linear trend between them as well as
#' normalizes the resulting variables (to mean = 0 and standard deviation = 1). Removal of the linear tremd is done using lm().
#' 
#' @param x A vector  
#' @param y A vector  
#' 
#' @keywords Statistics, regression
#' 
#' @return The function returns a list with three elements: \code{x} and \code{y}, the normalized x and y variables, and 
#' \code{parameters}, which contains the parameters needed to denormalize the data (e.g. restore the original values of x and y).
#'   
#' @examples
#' 
#' # Example data
#' n <- 50
#' x2 <- seq(-30,18, length = n)
#' x <- x + 31
#' y <- 50 + (-x2^3 -12*x2^2 + 800*x2)/300 + rnorm(n, sd = 5)
#' 
#' # Normalization
#' norm <- normalize_lm(x, y)
#' plot(norm$x, norm$y, main = "Normalized data")
#' 
#' # De-normalize the data
#' denorm <- denormalize_lm(norm$x, norm$y, norm$parameters)
#' plot(denorm$x, denorm$y, main = "Normalized and denormalized data")
#' 
#' # Compare with the original data  
#' plot(x, y, main = "Original data")
#' 
#' @export
normalize_lm <- function(x, y){
  mod <- lm(y ~ x, data = data.frame(x=x, y=y))
  sd_x <- sd(x)
  mean_x <- mean(x)
  sd_res <- sd(mod$residuals)
  list(
    x = (x-mean_x)/sd_x,
    y = mod$residuals/sd_res,
    parameters = list(
      sd_x=sd_x, mean_x=mean_x, 
      sd_res=sd_res, 
      interc = coef(mod)[1], slope = coef(mod)[2])
  )
}

#' @export
denormalize_lm <- function(x, y, parameters){
  x_denorm = x*parameters$sd_x + parameters$mean_x
  y_denorm = y*parameters$sd_res + parameters$interc + parameters$slope*x_denorm
  list(
    x = x_denorm,
    y = y_denorm
  )
}