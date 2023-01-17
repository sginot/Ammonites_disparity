#-------------------------------------------------------------------------------
# Script gathering code to produce various figures of the conch ratio morpho-
# space. The morphospace is a PCA based on fourier coefficients of the contours
# of whorl aperture section shapes (see script 001).

# Required packages:
library(MASS)
library(car)
library(scales)
library(viridis)

#-------------------------------------------------------------------------------
# For CONCH RATIO ANALYSES, load objects produced in script 011

load("PCA_conch_ratios.RData.RData")
load("taxa_lists_time_bins_conchs.RData")
load("variables_and_data_conch.RData")

#-------------------------------------------------------------------------------
# Biozone level plots

output_folder <- "../Figures/PCA_biozones_contour/"

for (i in c(1:6, 8:23, 25:30)) {
  
  pdf(file = paste(output_folder,
                   Biozones[i], 
                   "_morphospace_contours.pdf"))
  
  y <- pca$x[which(biozone_factor == Biozones[i]), 2] 
  x <- pca$x[which(biozone_factor == Biozones[i]), 1] 
  
  f1 <- kde2d(x, y, 
              lims = c(min(pca$x[,1])-0.2, 
                       max(pca$x[,1])+0.2, 
                       min(pca$x[,2])-0.2, 
                       max(pca$x[,2]))+0.2, 
              n = 100)
  
  label.species <- list.taxa.biozones[[i]]
  
  plot(pca$x[, c(1, 2)], 
       type = "n")
  
  color.palette = function(n) hcl.colors(n, 
                                         "YlOrRd", 
                                         rev = TRUE)
  
  .filled.contour(f1$x, 
                  f1$y, 
                  f1$z, 
                  levels = 1:20, 
                  col = color.palette(15))
  
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


for (i in c(7, 24)){ # biozones with n=1
  
  pdf(file = paste(output_folder,
                   Biozones[i], 
                   "_morphospace_contours.pdf"))
  
  y <- pca$x[which(biozone_factor == Biozones[i]), 2] 
  x <- pca$x[which(biozone_factor == Biozones[i]), 1] 

  label.species <- list.taxa.biozones[[i]]
  
  plot(pca$x[, c(1,2)], 
       type = "n")

  points(pca$x[, c(1,2)], 
         col = "grey", 
         pch = 19)
  
  points(x = x, 
         y = y, 
         col = "black", 
         pch = 19)

    dev.off()
}

#-------------------------------------------------------------------------------
# Stage level plots

output_folder <- "../Figures/PCA_stage_contour/"

for (i in 1:3){ #Stages
  pdf(file = paste(output_folder, 
                   Stages[i], 
                   "_stage_conch_morphospace_contours.pdf"))
  
  y <- pca$x[which(stage_factor == Stages[i]), 2] 
  x <- pca$x[which(stage_factor == Stages[i]), 1] 
  
  f1 <- kde2d(x, 
              y, 
              lims=c(min(pca$x[,1])-0.2, 
                     max(pca$x[,1])+0.2, 
                     min(pca$x[,2])-0.2, 
                     max(pca$x[,2]))+0.2, 
              n=100)
  
  label.species <- list.taxa.biozones[[i]]
  
  plot(pca$x[,c(1, 2)], type="n")
  
  color.palette = function(n) hcl.colors(n, 
                                         "YlOrRd", 
                                         rev = TRUE)
  
  .filled.contour(f1$x, 
                  f1$y, 
                  f1$z, 
                  levels = 1:10, 
                  col = color.palette(10))
  
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
