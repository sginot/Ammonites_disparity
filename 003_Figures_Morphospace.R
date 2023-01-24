#-------------------------------------------------------------------------------
# Script gathering code to produce various figures of the whorl section morpho-
# space. The morphospace is a PCA based on fourier coefficients of the contours
# of whorl aperture section shapes (see script 001).

# Required packages:
library(MASS)
library(car)
library(scales)
library(viridis)

#-------------------------------------------------------------------------------
# For WHORL SECTION ANALYSES, load objects produced in script 001 and 002 

load("PCA_Whorl_sections.RData")
load("taxa_lists_time_bins.RData")
load("variables_and_data_whorl_sections.RData")

#-------------------------------------------------------------------------------
# Check if output directories are present, otherwise, create them.

fdirs <- c('../Figures', 
           '../Figures/PCA_biozones_names_species', 
           '../Figures/PCA_biozones_contour',
           '../Figures/PCA_stage_contour', 
           '../Figures/PCA_interval_contour', 
           '../Figures/PCA_biozone_density',
           '../Figures/PCA_stage_density', 
           '../Figures/PCA_interval_density', 
           '../Figures/PCA_biozone_colors',
           '../Figures/PCA_stage_colors', 
           '../Figures/PCA_interval_colors')

for(s in fdirs) if(!dir.exists(s)) dir.create(s)

#-------------------------------------------------------------------------------
# Visual inspection: should already have been checked at the end of script 001.
# This plot is just to double check that no outlier is present.

plot(pca$x[,1:2])

# If no obvious outlier is detected, following steps can be run.

#-------------------------------------------------------------------------------
# Produce a series of graphs representing the morphospace occupation for each 
# biozone. This will produce 30 files, so make sure to place them in a specific
# Figures folder. Here we created a Figures folder in the main folder, one level
# higher than the repository folder, plus a subfolder within Figures.

output_folder <- "../Figures/PCA_biozones_names_species/"
  # This can be modified to another folder path.

for (i in 1:length(Biozones)) {
  
  pdf(file = paste(output_folder, 
                   Biozones[i], 
                   "_morphospace.pdf"))
  
  plot(pca$x[, 1:2]) 
    #Initiate plot
  
  points(pca$x[, c(1, 2)], 
         col = "grey", 
         pch = 19)
    # Add global morphospace occupation as grey points
  
  y <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)), 1] 
    # Get coordinates of species points present during biozone [i]
  
  label.species <- list.taxa.biozones[[i]] 
  
  points(x = x, 
         y = y, 
         col = "black", 
         pch = 19)
    # Add points for species present in biozone [i]
  
  text(x = x,
       y = y, 
       labels = label.species, 
       pos = 4, 
       cex = 0.5)
    # Overlay the image name for corresponding points
    # This can be removed if plot gets too messy
  
  dev.off()
  }

#-------------------------------------------------------------------------------
# Make "contour plots" for the same figures at the biozone level

output_folder <- "../Figures/PCA_biozones_contour/"

for (i in c(1:6, 8:23, 25:30)) { 
  # Biozones are specified manually here, to avoid errors in the loop due to
  # some biozones without enough samples to produce a contour
  
  pdf(file = paste(output_folder, 
                   Biozones[i],
                   "_morphospace_contours.pdf"))
  
  y <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)), 1]
  
  f1 <- kde2d(x, y, 
              lims = c(min(pca$x[,1])-0.2, 
                     max(pca$x[,1])+0.2, 
                     min(pca$x[,2])-0.2, 
                     max(pca$x[,2]))+0.2, 
              n = 100)
    # Create the contour, with n grid points and limited over the lims range
  
  plot(pca$x[, c(1, 2)], 
       type = "n")
    # Initiate plot without points
  
  color.palette = function(n) hcl.colors(n, 
                                         "YlOrRd", 
                                         rev = TRUE)
    # Define appropriate color palette
  
  .filled.contour(f1$x, 
                  f1$y, 
                  f1$z, 
                  levels = 1:20, 
                  col = color.palette(15))
    # Add color regions of the contour
  
  contour(f1, 
          xlim = c(min(pca$x[,1])-0.2, 
                   max(pca$x[,1]))+0.2, 
          ylim = c(min(pca$x[,2])-0.2, 
                   max(pca$x[,2])+0.2), 
          col = "grey70", 
          add = T, 
          lwd = 2, 
          drawlabels = F)
    # Add contour lines to the graph
  
  points(pca$x[, c(1, 2)], 
         col = "grey", 
         pch = 19)
    # Overlay with all points
  
  points(x = x,
         y = y,
         col = "black", 
         pch = 19)
    # Overlay with points present in biozone [i]
  
  dev.off()
  }


for (i in c(7, 24)) {
  # Biozones with only n=1
  
  pdf(file = paste(output_folder,
                   Biozones[i],
                   "_morphospace_contours.pdf"))
  
  y <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)),1] 
  
  label.species <- list.taxa.biozones[[i]]
  
  plot(pca$x[, c(1, 2)],
       type="n")
  
  points(pca$x[, c(1, 2)],
         col = "grey",
         pch = 19)
  points(x = x, 
         y = y, 
         col = "black",
         pch = 19)
  
  dev.off()
}

#-------------------------------------------------------------------------------
# Make similar contour plots for stages

output_folder <- "../Figures/PCA_stage_contour/"


for (i in 1:length(Stages)) { #Stages
  
  pdf(file = paste(output_folder,
                   Stages[i],
                   "_stage_morphospace_contours.pdf"))
  
  y <- pca$x[match(list.taxa.stages[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.stages[[i]], 
                   levels(Ap.Image)), 1]
  
  f1 <- kde2d(x, y, 
              lims = c(min(pca$x[,1])-0.2, 
                       max(pca$x[,1])+0.2, 
                       min(pca$x[,2])-0.2, 
                       max(pca$x[,2]))+0.2,
              n=100)
  
  label.species <- list.taxa.biozones[[i]]
  
  plot(pca$x[, c(1, 2)], 
       type = "n")
  
  color.palette = function(n) hcl.colors(n, 
                                         "YlOrRd", 
                                         rev = TRUE)
  
  .filled.contour(f1$x, 
                  f1$y, 
                  f1$z, 
                  levels = 1:10, 
                  col = color.palette(6))
  
  contour(f1, 
          xlim = c(min(pca$x[,1])-0.2, 
                   max(pca$x[,1]))+0.2, 
          ylim = c(min(pca$x[,2])-0.2, 
                   max(pca$x[,2])+0.2), 
          col = "grey70", 
          add = T, 
          lwd = 2, 
          drawlabels = F)
  
  points(pca$x[, c(1, 2)], 
         col = "grey", 
         pch = 19)
  
  points(x = x, 
         y = y, 
         col = "black", 
         pch = 19)
  
  #text(x=x, y=y, labels=label.species, pos=4, cex=0.5)
  
  dev.off()
}

#-------------------------------------------------------------------------------
# Make similar contour plots for stages

output_folder <- "../Figures/PCA_interval_contour/"

for (i in 1:length(Interval)) { #Intervals
  
  pdf(file = paste(output_folder, 
                   Interval[i], 
                   "_interval_morphospace_contours.pdf"))
  
  y <- pca$x[match(list.taxa.intervals[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.intervals[[i]],
                   levels(Ap.Image)), 1]
  
  f1 <- kde2d(x, y, 
              lims = c(min(pca$x[,1])-0.2, 
                       max(pca$x[,1])+0.2, 
                       min(pca$x[,2])-0.2, 
                       max(pca$x[,2]))+0.2, 
              n = 100)
  
  label.species <- list.taxa.intervals[[i]]
  
  plot(pca$x[, c(1, 2)], 
       type = "n")
  
  color.palette = function(n) hcl.colors(n, 
                                         "YlOrRd", 
                                         rev = TRUE)
  
  .filled.contour(f1$x, 
                  f1$y, 
                  f1$z, 
                  levels = 1:10, 
                  col = color.palette(6))
  
  contour(f1, 
          xlim = c(min(pca$x[,1])-0.2, 
                   max(pca$x[,1]))+0.2, 
          ylim = c(min(pca$x[,2])-0.2, 
                   max(pca$x[,2])+0.2), 
          col = "grey70", 
          add = T, 
          lwd = 2, 
          drawlabels = F)
  
  points(pca$x[, c(1, 2)], 
         col = "grey", 
         pch = 19)
  
  points(x = x,
         y = y,
         col = "black", 
         pch = 19)
  
  #text(x=x, y=y, labels=label.species, pos=4, cex=0.5)
  
  dev.off()
}

#-------------------------------------------------------------------------------
# Plot the same morphospaces, with density curves on the edges
# First for biozones

output_folder <- "../Figures/PCA_biozone_density/"

for (i in c(1:6, 8:23, 25:30)) { 
  # Biozones 7 and 24 have not enough points (N = 1) for density curves
  
  pdf(file = paste(output_folder,
                   Biozones[i], 
                   "_morphospace_density.pdf"))
  
  y <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)), 1]
    # Get coordinates for points present in Biozone [i]
  
  dsx <- density(x, 
                 from = min(pca$x[,1]), 
                 to = max(pca$x[,1]))
    # Kernel density estimation for x axis (PC1)
  
  dsy <- density(y, 
                 from = min(pca$x[,2]), 
                 to = max(pca$x[,2]))
    # Kernel density estimation for y axis (PC2)
  
  layout(matrix(c(1, 1, 4, 2, 2, 3, 2, 2, 3), 
                ncol = 3, 
                byrow = T))
    # Prepare graphical device layout to leave space around the graph for curves
  
  par(mar = c(0, 4, 0, 0))
    # Adjust margins around the first plot
  
  plot(dsx, 
       xaxt = "n",
       yaxt = "n", 
       xlab = "", 
       ylab = "", 
       main = "", 
       bty = "n",
       zero.line = F)
    # Plot the density curve for x axis
  
  par(mar = c(4, 4, 0, 0))
    # Re-adjust margins for second plot
  
  plot(pca$x[, c(1, 2)], 
       col = "grey", 
       pch = 19)
    # Plot morphospace, ie PCA
  
  points(x = x, 
         y = y, 
         col = "black", 
         pch = 19)
    # Add points which are present in biozone [i]
  
  par(mar = c(4, 0, 0, 0))
    # Re-adjust margins for third plot
  
  plot(-dsy$y, 
       dsy$x, 
       ylim = range(dsy$x), 
       xlim = rev(range(-dsy$y)), 
       type = "l",
       xaxt = "n",
       yaxt = "n", 
       xlab = "", 
       ylab = "", 
       main = "", 
       bty = "n")
    # Plot density curve for y axis (PC2)
  
  par(mar = c(0, 0, 0, 0))
    # Readjust margins for final plot
  
  plot(1, 1, 
       "n", 
       xlab = "", 
       ylab = "",
       xaxt = "n", 
       yaxt = "n", 
       bty = "n")
    # Empty plotting area for title
  
  text(1, 1, 
       labels = Biozones[i], 
       cex=1.2)
    # Add title, ie name of biozone
  
  dev.off()
}

# Do the plots for the two biozones with N = 1


for (i in c(7, 24)) { 
  # biozone 7 and 24 have not enough points for density
  
  pdf(file = paste(output_folder,
                   Biozones[i], 
                   "_morphospace_density.pdf"))
  
  y <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)), 1] 
  
  layout(matrix(c(1, 1, 4,2, 2, 3, 2, 2, 3), 
                ncol=3, 
                byrow=T))
  
  par(mar = c(0, 4, 0, 0))
  plot(1, 1, 
       xaxt = "n",
       yaxt = "n", 
       xlab = "", 
       ylab = "", 
       main = "", 
       bty = "n", 
       type = "n")
  
  par(mar = c(4, 4, 0, 0))
  plot(pca$x[, c(1, 2)], 
       col = "grey", 
       pch = 19)
  points(x = x, 
         y = y, 
         col = "black", 
         pch = 19)
  
  par(mar = c(4, 0, 0, 0))
  plot(1, 1, 
       xaxt = "n",
       yaxt = "n", 
       xlab = "", 
       ylab = "", 
       main = "", 
       bty = "n", 
       type = "n")
  
  par(mar = c(0, 0, 0, 0))
  plot(1, 1, 
       "n", 
       xlab = "", 
       ylab = "",
       xaxt = "n", 
       yaxt = "n", 
       bty = "n")
  text(1, 1, 
       labels = Biozones[i], 
       cex = 1.2)
  
  #text(x=x, y=y, labels=label.species, pos=4, cex=0.5)
  
  dev.off()
}

#-------------------------------------------------------------------------------
# Make same plots at the stage level

output_folder <- "../Figures/PCA_stage_density/"

for (i in c(1:length(Stages))) {
  
  pdf(file = paste(output_folder,
                   Stages[i], 
                   "_stage_morphospace_density.pdf"))
  
  y <- pca$x[match(list.taxa.stages[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.stages[[i]], 
                   levels(Ap.Image)), 1] 
  
  dsx <- density(x, 
                 from = min(pca$x[, 1]), 
                 to = max(pca$x[, 1]))
  
  dsy <- density(y, 
                 from = min(pca$x[, 2]), 
                 to = max(pca$x[, 2]))
  
  layout(matrix(c(1, 1, 4, 2, 2, 3, 2, 2, 3), 
                ncol = 3, 
                byrow = T))
  
  par(mar = c(0, 4, 0, 0))
  
  plot(dsx, 
       xaxt = "n",
       yaxt = "n", 
       xlab = "", 
       ylab = "", 
       main = "", 
       bty = "n",
       zero.line = F)
  
  par(mar = c(4, 4, 0, 0))
  
  plot(pca$x[, c(1, 2)], 
       col = "grey", 
       pch = 19)
  
  points(x = x, 
         y = y, 
         col = "black", 
         pch = 19)
  
  par(mar = c(4, 0, 0, 0))
  
  plot(-dsy$y, 
       dsy$x, 
       ylim = range(dsy$x), 
       xlim = rev(range(-dsy$y)), 
       type = "l",
       xaxt = "n",
       yaxt = "n", 
       xlab = "", 
       ylab = "", 
       main = "", 
       bty = "n")
  
  par(mar = c(0, 0, 0, 0))
  
  plot(1, 1, 
       "n", 
       xlab = "", 
       ylab = "",
       xaxt = "n",
       yaxt = "n", 
       bty = "n")
  
  text(1, 1, 
       labels = Stages[i], 
       cex = 1.2)
  
  dev.off()
}

#-------------------------------------------------------------------------------
# Make same plots at the interval level

output_folder <- "../Figures/PCA_interval_density/"

for (i in c(1:7)) {
  
  pdf(file = paste(output_folder,
                   Interval[i], 
                   "_interval_morphospace_density.pdf"))
  
  y <- pca$x[match(list.taxa.intervals[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.intervals[[i]], 
                   levels(Ap.Image)), 1] 
  
  dsx <- density(x, 
                 from = min(pca$x[, 1]), 
                 to = max(pca$x[, 1]))
  
  dsy <- density(y, 
                 from = min(pca$x[, 2]), 
                 to = max(pca$x[, 2]))
  
  layout(matrix(c(1, 1, 4, 2, 2, 3, 2, 2, 3),
                ncol = 3, 
                byrow = T))
  
  par(mar = c(0, 4, 0, 0))
  
  plot(dsx, 
       xaxt = "n", 
       yaxt = "n", 
       xlab = "", 
       ylab = "", 
       main = "", 
       bty = "n",
       zero.line = F)
  
  par(mar = c(4, 4, 0, 0))
  
  plot(pca$x[, c(1, 2)], 
       col = "grey", 
       pch = 19)
  
  points(x = x, 
         y = y, 
         col = "black", 
         pch = 19)
  
  par(mar = c(4, 0, 0, 0))
  
  plot(-dsy$y, 
       dsy$x, 
       ylim = range(dsy$x), 
       xlim = rev(range(-dsy$y)), 
       type = "l", 
       xaxt = "n",
       yaxt = "n", 
       xlab = "", 
       ylab = "", 
       main = "", 
       bty = "n")
  
  par(mar = c(0, 0, 0, 0))
  
  plot(1, 1, 
       "n", 
       xlab = "", 
       ylab = "",
       xaxt = "n", 
       yaxt = "n", 
       bty = "n")
  
  text(1, 1, 
       labels = Interval[i], 
       cex = 1.2)
  
  dev.off()
}

#-------------------------------------------------------------------------------
# Make the same morphospaces with different colors showing superfamilies
# Add ellipses to the interval-level graphs

palette(value = "ggplot2")
  # Renew color palette for more marked color differences

#-------------------------------------------------------------------------------
# Start with biozone level

output_folder <- "../Figures/PCA_biozone_colors/"

for (i in 1:length(Biozones)) { # number of biozones
  
  pdf(file = paste(output_folder,
                   Biozones[i], 
                   "_morphospace_colours.pdf"))
  
  plot(pca$x[, 1:2],
       col = "grey",
       pch = 19)
    # Plot background morphospace
  
  y <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)), 1] 
  
  label.species <- list.taxa.biozones[[i]] 
  
  points(x = x, 
         y = y, 
         col = c(1:7)[list.superfamily.biozones[[i]]], 
         pch = 19)
  
  text(x = x, 
       y = y, 
       labels = label.species, 
       pos = 4, 
       cex = 0.5)
  
  dev.off()
}

#-------------------------------------------------------------------------------
# Same at the stage level, without species names

output_folder <- "../Figures/PCA_stage_colors/"

for (i in 1:length(Stages)) { # Number of stages
  
  pdf(file = paste(output_folder, 
                   Stages[i], 
                   "_morphospace_stage_colors.pdf"))
  
  plot(pca$x[, 1:2],
       col = "grey",
       pch = 19)
  
  y <- pca$x[match(list.taxa.stages[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.stages[[i]],
                   levels(Ap.Image)), 1] 
  
  label.species <- list.taxa.stages[[i]] 
  
  points(x = x,
         y = y,
         col = c(1:7)[list.superfamily.stages[[i]]], 
         pch = 19)
  
  #text(x=x, y=y, labels=label.species, pos=4, cex=0.5)
  
  dev.off()
}

#-------------------------------------------------------------------------------
# Same at the interval level, with confidence ellipses overlaid

output_folder <- "../Figures/PCA_interval_colors/"

for (i in 1:length(Interval)) { # Number of intervals
  
  pdf(file = paste(output_folder, 
                   Interval[i], 
                   "_morphospace_interval_ell50.pdf"))
    # First for 50% ellipses
  
  plot(pca$x[, 1:2], 
       col = "grey", 
       pch = 19)
  
  y <- pca$x[match(list.taxa.intervals[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.intervals[[i]], 
                   levels(Ap.Image)), 1] 
  
  label.species <- list.taxa.intervals[[i]] 
  
  for (j in 1:length(Superfam)) { # Number of superfamilies
    
    tryCatch({
      # tryCatch avoids the loop from stopping when there are not enough
      # specimens to create an ellipse
      
      dataEllipse(x[which(list.superfamily.interval[[i]] == 
                            levels(list.superfamily.interval[[i]])[j])], 
                  y[which(list.superfamily.interval[[i]] == 
                            levels(list.superfamily.interval[[i]])[j])], 
                  add = T, 
                  levels = 0.5, 
                  col = palette()[j], 
                  fill = T,
                  center.pch = F, 
                  plot.points = F)}, 
        # Add ellipses for each family [j] present in interval [i]
      
      error = function(a) {})
        # If family is not present in this interval, return only message
    
  }
  
  points(x = x,
         y = y, 
         col = c(1:7)[list.superfamily.interval[[i]]], 
         pch = 19)
  
  dev.off()
  
}
  # The loop may issue warning messages, with relate to the drawing of some
  # "line shaped ellipses" when only few specimens are present for a given
  # superfamily


# Make the same graphs, with 75% ellipses

for (i in 1:length(Interval)) {
  
  pdf(file = paste(output_folder,
                   Interval[i], 
                   "_morphospace_interval_ell75.pdf"))
  
  plot(pca$x[, c(1, 2)],
       col = "grey", 
       pch = 19)
  
  y <- pca$x[match(list.taxa.intervals[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.intervals[[i]], 
                   levels(Ap.Image)), 1] 
  
  label.species <- list.taxa.intervals[[i]] 
  
  for (j in 1:length(Superfam)) {
    
    tryCatch({
      
      dataEllipse(x[which(list.superfamily.interval[[i]] == 
                            levels(list.superfamily.interval[[i]])[j])], 
                  y[which(list.superfamily.interval[[i]] == 
                            levels(list.superfamily.interval[[i]])[j])], 
                  add = T, 
                  levels = 0.75, 
                  col = palette()[j], 
                  fill = T, 
                  center.pch = F, 
                  plot.points = F)}, 
      
      error = function(a) {})
    
  }
  
  points(x = x, 
         y = y, 
         col = c(1:7)[list.superfamily.interval[[i]]], 
         pch = 19)
  
  dev.off()
  
}

# Make the same graphs with 99% ellipses

for (i in 1:length(Interval)) { 
  
  pdf(file = paste(output_folder,
                   Interval[i], 
                   "_morphospace_interval_ell100.pdf"))
  
  plot(pca$x[, c(1, 2)], 
       col = "grey", 
       pch = 19)
  
  y <- pca$x[match(list.taxa.intervals[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.intervals[[i]], 
                   levels(Ap.Image)), 1] 
  
  label.species <- list.taxa.intervals[[i]] 
  
  for (j in 1:length(Superfam)) {
    
    tryCatch({
      
      dataEllipse(x[which(list.superfamily.interval[[i]] == 
                            levels(list.superfamily.interval[[i]])[j])], 
                  y[which(list.superfamily.interval[[i]] == 
                            levels(list.superfamily.interval[[i]])[j])], 
                  add = T, 
                  levels = 0.99, 
                  col = palette()[j], 
                  fill = T, 
                  center.pch = F, 
                  plot.points = F)}, 
      
      error = function(a) {})
    
  }
  
  points(x = x, 
         y = y, 
         col = c(1:7)[list.superfamily.interval[[i]]], 
         pch = 19)
  
  dev.off()
  
}

#-------------------------------------------------------------------------------
# END OF SCRIPT