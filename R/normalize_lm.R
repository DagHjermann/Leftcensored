#' Normalization and removal of correlation
#' 
#' This function takes two variables x and y as an input, and removes the linear part of the correlation between them, as well as
#' normalizes the resulting variables (to mean = 0 and standard deviation = 1). Removal of the linear trend is done using lm(y ~ x).
#'
#' @param x 
#' @param y 
#'
#' @return List of x, y and parameters
#' @export
#'
#' @examples
#' 
#' # Example data
#' n <- 50
#' x2 <- seq(-30,18, length = n)
#' x <- x2 + 31
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
#' # De-normalize new data
#' x_new <- c(-1.5, 1.5, 1.5, -1.5)
#' y_new <- c(-2, -2, 2, 2)
#' 
#' denorm <- denormalize_lm(norm$x, norm$y, norm$parameters)
#' plot(denorm$x, denorm$y, main = "Normalized and denormalized data")
#' 

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

#' De-normalization of data
#' 
#' De-normalizes data that has  
#' 
#' @param x Vector of x values (NA values not permitted).   
#' @param y Vector of y values (NA values not permitted). Must have same length as x.    
#' @param parameters List of parameters to use for de-normalizaton. Normally, one will use the 
#' output of \code{normalize_lm}, which contains the list element \code{parameters}.  
#'
#' @export
#' 
#' @examples 
#' 
#' # Example data
#' n <- 50
#' x2 <- seq(-30,18, length = n)
#' x <- x2 + 31
#' y <- 50 + (-x2^3 -12*x2^2 + 800*x2)/300 + rnorm(n, sd = 5)
#' 
#' # Normalization
#' norm <- normalize_lm(x, y)
#' plot(norm$x, norm$y, main = "Normalized data")
#' 
#' # Add rectangle
#' x_new <- c(-1.5, 1.5, 1.5, -1.5)
#' y_new <- c(-2, -1, 2, 2)
#' polygon(x_new, y_new, color = "red")
#' 
#' # De-normalize the data
#' denorm <- denormalize_lm(norm$x, norm$y, norm$parameters)
#' plot(denorm$x, denorm$y, main = "Normalized and denormalized data")
#' 
#' # De-normalize the rectangle data, based on the normalization of the data
#' denorm_new <- denormalize_lm(x_new, y_new, norm$parameters)
#' polygon(denorm_new$x, denorm_new$y, border = "red")

denormalize_lm <- function(x, y, parameters){
  x_denorm = x*parameters$sd_x + parameters$mean_x
  y_denorm = y*parameters$sd_res + parameters$interc + parameters$slope*x_denorm
  list(
    x = x_denorm,
    y = y_denorm
  )
}