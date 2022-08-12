#' Plot left-justified data
#' 
#' @param data Data frame with columns named 'x', 'y_cens', 'y_aboveLOQ', and (if you use version = 2) 'y_LOQ'
#' @param size Size of points (default is 3). 
#' @param size Width of jittering of points (to reduce overlap between points with same y value). 
#' The default is 0, i.e. no jittering. 
#' 
#' @import ggplot2
#' 
#' @examples
#' # Plot data
#' data_prepared <- lc_prepare(polybrom, 
#'                         x = "year",
#'                         y = "concentration", 
#'                         censored = "LOQ_flag",
#'                         log = TRUE)
#'                         
#' gg <- lc_ggplot(data_prepared, type = 1)
#' 
#' Faceting, themes etc. can be added to the final object
#' gg + facet_wrap(vars(station))
#' gg + facet_wrap(vars(station)) + theme_dark()
#' 
#' # Add regression line to plot
#' data_23B <- data_prepared %>% filter(station == "23B)
#' result <- lc_linear(data_23B)
#' plotdata <- regression_ci(data_prepared, result$model)
#' gg +
#'   geom_ribbon(data = plotdata, aes(x = x, ymin = y_lo, ymax = y_hi),
#'   fill = "red", alpha = 0.25) +
#'   geom_path(data = plotdata, aes(x = x, y = y), color = "red3", size = 1)
#' 
#' @export
lc_ggplot <- function(data, size = 2, width = 0){   
  # Concentrations below LOD
  sel <- data$uncensored == 0
  gg <- ggplot(data) +
    geom_jitter(data = data[sel,], aes(x, threshold), size = size, 
               color = "black", fill = "red", shape = 25, width = width) +
    geom_jitter(data = data[!sel,], aes(x, y), size = size, 
               color = "black", fill = "white", shape = 21, width = width) + 
    theme_bw()
  gg
}

