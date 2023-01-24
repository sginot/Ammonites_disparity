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
load("SOV_SOR_disparity.RData")

#-------------------------------------------------------------------------------
# Check if output directory is present, if not create it

fdirs <- c('../Figures/Morpho_rarefaction')

for(s in fdirs) if(!dir.exists(s)) dir.create(s)

#-------------------------------------------------------------------------------
# Define plotting function to be re-used with the different lists containing
# rarefaction pseudo values.

plot.morph.raref <- function(x, 
                             split = F, 
                             cols = cols) {
  # x must be a list containing for each time bin a table with pseudovalues
  # split defines if curves should be plotted overlaid or in separate plots
  
  maxN <- max(unlist(lapply(x, dim)))
    # Maximum number of observations across all individual time bins
  
  maxdispa <- max(unlist(x))
    # Maximum value of given disparity index of list x
  
    if (!split) {
      
      plot(1, 1,
           type = "n",
           xlim = c(0, maxN + 1),
           ylim = c(0, maxdispa),
           xlab = "Sample size",
           ylab = "Disparity",
           main = "Morphological rarefaction curves")
        # Empty plot
      
      for (i in 1:length(x)) {
        
        mat <- x[[i]]
          # isolate matrix [i] for practical reasons
        
        n <- ncol(mat)
          # Find maximum number of subsampled observations for this specific
          # time bin matrix
        
        maxmat <- apply(mat, 
                        2, 
                        max)
        minmat <- apply(mat, 
                        2, 
                        min)
        meanmat <- apply(mat, 
                         2, 
                         mean)
          # These will be the values that make up the curve
        
        lines(2:(n + 1), 
              meanmat, 
              lwd = 2.5,
              col = cols[i])
          # Color palette must have a number of colors >= number time bins
        
        lines(2:(n + 1), 
              maxmat, 
              lty = 2,
              col = cols[i])
          # Add a dashed line showing the maximum value among subsamples 
          # iterations
        
        lines(2:(n + 1), 
              minmat,
              lty = 2,
              col = cols[i])
          # Add a dashed line showing the minimum value among subsamples
          # iterations
        
      }
    }
  
  if (split) {
    # If overlaid curves are not readable, make separate plots
    
    layout(mat = matrix(1:length(x), 
                        nrow = 1))
      # Separate plots along a single line
    
    par(mar = c(5,2,1,1))
      # Reduce margins
    
    for (i in 1:length(x)) {
      
      mat <- x[[i]]

      n <- ncol(mat)
 
      maxmat <- apply(mat, 
                      2, 
                      max)
      minmat <- apply(mat, 
                      2, 
                      min)
      meanmat <- apply(mat, 
                       2, 
                       mean)
      plot(2:(n + 1),
           meanmat,
           lwd = 2.5,
           col = cols[i],
           type = "l",
           xlim = c(0, n + 1),
           ylim = c(0, maxdispa),
           xlab = "")
        # A new plot will be generated at each ith iteration
      
      lines(2:(n + 1), 
            maxmat, 
            lty = 2,
            col = cols[i])
      
      lines(2:(n + 1), 
            minmat,
            lty = 2,
            col = cols[i])
      
    }
    
  }
  
}

#-------------------------------------------------------------------------------
# Apply function to biozones

cols <- viridis(n = 30)
  # Define a color palette with sufficient number of colors

plot.morph.raref(x = rarefSOR_biozones[c(1:6, 8:23, 25:30)], 
                 split = F, 
                 cols = cols)
  # Overlaid curves impossible to read unless only few curves are shown

#-------------------------------------------------------------------------------
# New version of function, directly split plots

plot.morph.raref2 <- function(x, 
                              cols = cols, 
                              names.bins = NA,
                              non.raref.values = NA) {
  
  maxN <- max(unlist(lapply(x, dim)))

  maxdispa <- max(unlist(x))

  layout(mat = matrix(1:length(x), 
                        nrow = 1))
    
  par(mar = c(5, 0, 5, 0))
    # No margins on left and right, so that graph makes a chronological line
    
  for (i in 1:length(x)) {
      
    tryCatch({
      
      mat <- x[[i]]
      
      n <- ncol(mat)
      
      maxmat <- apply(mat, 
                      2, 
                      max)
      minmat <- apply(mat, 
                      2, 
                      min)
      meanmat <- apply(mat, 
                       2, 
                       mean)
      
      plot(2:(n + 1),
           meanmat, 
           lwd = 2.5,
           col = cols[i],
           type = "l",
           xlim = c(0, n + 1),
           ylim = c(0, maxdispa),
           xlab = "",
           bty = "n",
           yaxt = "n",
           main = names.bins[i],
           cex.main = 0.5)
        # A new plot will be generated at each ith iteration
      
      lines(2:(n + 1), 
            maxmat, 
            lty = 2,
            col = cols[i])
    
      lines(2:(n + 1), 
            minmat,
            lty = 2,
            col = cols[i])
    
      lines(x = c(0, n + 1),
            y = rep(non.raref.values[i], 2),
            col = cols[i],
            lwd = 3)
        # Make horizontal line to show the real disparity value, measure across
        # all species present in biozone [i]
      
      abline(v = n + 1,
             lty = 2,
             lwd = 2,
             col = "gray")
        # Make vertical line to separate biozones
      },

    error = function(a) {
      plot(1, 1, 
           type = "n",
           xaxt = "n",
           yaxt = "n",
           bty = "n",
           xlab = "")
      
      text(1, 1, "NA")})
      
    }
  
}

#-------------------------------------------------------------------------------
# Apply new function to biozones level to make chronological plots

output_folder <- "../Figures/Morpho_rarefaction/"

cols <- viridis(n = 30)
  # Define a color palette with sufficient number of colors

raref_SOR_ord <- rarefSOR_biozones[ord]
SOR_ord <- SOR_biozones[ord]
  # Order variables in CHRONOLOGICAL ORDER

pdf(paste(output_folder,
          "Morpho_raref_SOR_biozones.pdf"), 
    width = 40, 
    height = 3)

plot.morph.raref2(x = raref_SOR_ord, 
                  cols = cols,
                  names.bins = Biozones[ord],
                  non.raref.values = SOR_ord)

dev.off()

#-------------------------------------------------------------------------------
# Same for interval level

cols <- viridis(n = 7)
  # Define a color palette with sufficient number of colors

pdf(paste(output_folder,
          "Morpho_raref_SOR_intervals.pdf"), 
    width = 10, 
    height = 3)

plot.morph.raref2(x = rarefSOR_intervals, 
                  cols = cols,
                  names.bins = Interval,
                  non.raref.values = SOR_intervals)

dev.off()

#-------------------------------------------------------------------------------
# Similar plots can be done for the sum of variance (SOV)

# At the biozones level

cols <- viridis(n = 30)

raref_SOV_ord <- rarefSOV_biozones[ord]
SOV_ord <- SOV_biozones[ord]
  # Order variables in CHRONOLOGICAL ORDER

pdf(paste(output_folder,
          "Morpho_raref_SOV_biozones.pdf"), 
    width = 40, 
    height = 3)

plot.morph.raref2(x = raref_SOV_ord, 
                  cols = cols,
                  names.bins = Biozones[ord],
                  non.raref.values = SOV_ord)

dev.off()

# At the interval level

cols <- viridis(n = 7)

pdf(paste(output_folder,
          "Morpho_raref_SOV_intervals.pdf"), 
    width = 10, 
    height = 3)

plot.morph.raref2(x = rarefSOV_intervals, 
                  cols = cols,
                  names.bins = Interval,
                  non.raref.values = SOV_intervals)

dev.off()

#-------------------------------------------------------------------------------
# END OF SCRIPT