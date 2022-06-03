#' Concentrations of Acenaphthene in blue mussels, 1995-2011
#' 
#' Data from contaminantmonitoring of contaminants in blue mussls in Norway. The data are concentrations (Conc)
#' of acenaphthene (a PAH com) in blue mussels in a single site, measured in 1-3 batches of blue mussel per year.
#' The column 'Flag' indicates whther concentrations are below or above the limit of quantification (LOQ). When
#' Flag = '<', the data are below the LOQ, i.e., we know that the concentrations are somewhere between zero and LOQ.
#' 
#' @docType data
#' 
#' @usage data(concentrations)
#' 
#' @format A data frame with 22 rows and 6 columns
#' 
#' @keywords datasets
#' 
#' @references Green et al. 2019. Contaminants in coastal waters of Norway 2018. NIVA report 7412-2019. (\link{http://hdl.handle.net/11250/2635080})
#' 
#' @source NIVA database
#' 
#' @examples
#' # Prepare data
#' data(contaminants)
#' data_prepared <- leftcensored_prepare(concentrations)
#' 
#' # Perform analysis
#' result <- leftcensored_lm(data_prepared)
#' 
#' # Check result
#' result$summary
#' plot(result$model)
#' 
#' # Plots
#' regressiondata_leftcens <- regression_ci(data_prepared, result)
#' 
"concentrations"