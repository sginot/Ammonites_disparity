#-------------------------------------------------------------------------------
# Morphological rarefaction plots at biozone and interval levels

# Required packages:
library(scales)
library(viridis)

#-------------------------------------------------------------------------------
# Load related objects
load("Raref_disparity_biozone.RData")
load("Raref_disparity_intervals.RData")
load("variables_and_data_whorl_sections.RData")

#-------------------------------------------------------------------------------

plot.morph.raref <- function(x, split = F) {
  # x must be a list containing for each time bin a table with pseudovalues
  # split defines if curves should be plotted overlaid or in separate plots
  
  maxN <- max(unlist(lapply(x, dim)))
    # Maximum number of observations across all individual biozones
  
  maxdispa <- max(unlist(x))
    # Maximum value of given disparity index of list x
  
    if (!split) {
      
      plot(1, 1,
           type = "n",
           xlim = c(0, maxN),
           ylim = c(0, maxdispa),
           xlab = "Sample size",
           ylab = "Disparity",
           main = "Morphological rarefaction curves for biozones")
      
      for (i in 1:length(x)) {
        
        mat <- x[[i]]
        
        lines()
        
      }
    }
  
}




