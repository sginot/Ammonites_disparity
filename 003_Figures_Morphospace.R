#-------------------------------------------------------------------------------
# Script gathering code to produce various figures of the whorl section morpho-
# space. The morphospace is a PCA based on fourier coefficients of the contours
# of whorl aperture section shapes (see script 001).

# Required packages:
library(MASS)
library(scales)

#-------------------------------------------------------------------------------
# Load objects produced in script 001 and 002

load("PCA_Whorl_sections.RData")
load("taxa_lists_time_bins.RData")
load("variables_and_data_whorl_sections.RData")

#-------------------------------------------------------------------------------
# Original plot: visual inspection is necessary if script 001 has been run for
# the first time. 
# In some cases an outlier shape may appear, due to slight noise in the 
# interpolation of shapes causing the first landmark not to align
# with those of the other shapes (therefore starting from an non-homologous 
# location). 
# If this happens, simply re-run script 001. Once outlier shape is removed, 
# save the morphospace and use it in further analysis, rather than 
# re-running script 001.

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

for (i in c(1:6,8:23,25:30)) { 
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


for (i in c(7,24)) {
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
