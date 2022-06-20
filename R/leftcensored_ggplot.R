#' Plot left-justified data
#' 
#' @param data Data frame with columns named 'x', 'y_cens', 'y_aboveLOQ', and (if you use version = 2) 'y_LOQ'
#' @param type Set to either 1 (don't add y_LOQ points) or 2 (add y_LOQ points). Default is 1.
#' 
#' @import ggplot2
#' 
#' @examples
#' # Plot data
#' data(concentrations)
#' data_prepared <- leftcensored_prepare(concentrations)
#' gg <- leftcensored_ggplot(data_prepared, type = 1)
#' gg
#' 
#' # Add regression line to plot
#' result <- lm_linear(data_prepared)
#' plotdata <- regression_ci(data_prepared, result$model)
#' gg +
#'   geom_ribbon(data = plotdata, aes(x = x, ymin = y_lo, ymax = y_hi),
#'   fill = "red", alpha = 0.25) +
#'   geom_path(data = plotdata, aes(x = x, y = y), color = "red3", size = 1)
#' 
#' @export
leftcensored_ggplot <- function(data, type = 1){   
  # Concentrations below LOD
  sel <- data$y_aboveLOQ == 0
  gg <- ggplot(data) +
    geom_point(data = data[sel,], aes(x, y_cens), pch = 21, size = 3, color = "black", fill = "red") +
    geom_point(data = data[!sel,], aes(x, y_cens), pch = 19, size = 3, color = "black")
    if (type == 2){
      gg <- gg + 
        geom_point(aes(x = x, y = y_LOQ), color = "green3")
    }
  gg + theme_bw()
}

