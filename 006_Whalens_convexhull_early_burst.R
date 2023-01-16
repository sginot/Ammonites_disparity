#-------------------------------------------------------------------------------
# Compute convex hull area of morphospace per biozones and per interval.
# This is inspired from Whalen's et al. (2020) paper in Science Advances.
# Note that this restricts disparity only to the first TWO PCs of the PCA, due
# to the implementation of function 'chull'.

# In addtion to the convex hull area, we re-implement Whalen's et al approach
# to test for early bust, by producing null distributions of the convex hull
# area across time. To do so, morphologies are randomly assigned to the 
# species, and the sample stays the same.

# Required packages:
library(sp)
library(scales)

#-------------------------------------------------------------------------------
# Load necessary objects from previous scripts
load("PCA_Whorl_sections.RData")
load("taxa_lists_time_bins.RData")

#-------------------------------------------------------------------------------
# Compute convex hull area for all biozones

CH_Area_biozones <- rep(NA, 
                        length(Biozones))
  # Make empty vector to record results

for (i in 1:length(Biozones)) {
  
  y <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.biozones[[i]], 
                   levels(Ap.Image)), 1]
    # Coordinates of all obeservations present in Biozone [i]
  
  CH <- chull(x, y)
    # Indices of points forming the outer edge of distribution, ie convex hull
  
  CH_xcoo <- x[CH]
  CH_ycoo <- y[CH]
    # Coordinates of points forming the outer edge of distribution, ie the 
    # coordinates of the convex hull
  
  sp.chull <- Polygon(coords = cbind(c(CH_xcoo, CH_xcoo[1]), 
                                     c(CH_ycoo, CH_ycoo[1])), 
                      hole = F)
    # Note that to close the convex hull, the coordinates of the first point
    # must be repeated as the last point. Cf. help of function Polygon
  
  CH_Area_biozones[i] <- sp.chull@area
    # Record the convex hull area for biozone [i]
  
}
  # Loop may issue warning messages for biozones in which there are less then
  # two observations.

#-------------------------------------------------------------------------------
# Do the same for intervals

CH_Area_intervals <- rep(NA, 
                        length(Interval))

for (i in 1:length(Interval)) {
  
  y <- pca$x[match(list.taxa.intervals[[i]], 
                   levels(Ap.Image)), 2] 
  x <- pca$x[match(list.taxa.intervals[[i]], 
                   levels(Ap.Image)), 1]

  CH <- chull(x, y)
  
  CH_xcoo <- x[CH]
  CH_ycoo <- y[CH]

  sp.chull <- Polygon(coords = cbind(c(CH_xcoo, CH_xcoo[1]), 
                                     c(CH_ycoo, CH_ycoo[1])), 
                      hole = F)
  
  CH_Area_intervals[i] <- sp.chull@area

}


#-------------------------------------------------------------------------------
# Reassign existing morphologies randomly with replacement to the observations
# which are contained in the data table. Then compute the new values for convex
# hull areas, to produce a null distribution for each time bin.

iter <- 10000
  # Define number of iterqtions of resampling and computing of pseudo-values

null_CH_Area_biozones <- matrix(nrow = iter,
                                ncol = length(Biozones))
colnames(null_CH_Area_biozones) <- Biozones
  # Make empty matrix to record results, with iterations as rows and biozones
  # as columns

for (it in 1:iter) {
  
  Ap.rand <- as.factor(sample(table.aperture$Aperture.Image,
                                    replace = T)) 
    # Make shuffled vector of which images correspond to all observations 
  
  CH_Areas <- rep(NA, length(Biozones))
    # Make empty vector to record results across all biozones for iteration [it] 
  
  for (i in 1:length(Biozones)) {
    
    oo <- match(list.taxa.biozones[[i]], Ap.Image)
      # Indices of taxa present in biozone in within the imqge name vector
    
    taxa.rand <- Ap.rand[oo]
      # Indexing the shuffled image name vector with the same indices keeps the
      # correct position of "real" specimens, but assign them a different image
    
    y <- pca$x[match(taxa.rand, 
                     levels(Ap.Image)), 2] 
    x <- pca$x[match(taxa.rand, 
                     levels(Ap.Image)), 1]
      # The basic PCA remains the same with all 582 morphologies. The shuffled
      # image name vector does not contain all these image names, but all image
      # names that are resampled are present in original vector. The original
      # vector levels correspond to the lines of the PCA. Therefore, the taxa.
      # rand vector will have the right number of taxa for biozone [i], but with
      # positions corresponding to other taxa

    CH <- chull(x, y)
    
    CH_xcoo <- x[CH]
    CH_ycoo <- y[CH]
    
    sp.chull <- Polygon(coords = cbind(c(CH_xcoo, CH_xcoo[1]), 
                                       c(CH_ycoo, CH_ycoo[1])), 
                        hole = F)
    
    CH_Areas[i] <- sp.chull@area
    
  }
  
  null_CH_Area_biozones[it,] <- CH_Areas
  
}

#-------------------------------------------------------------------------------
# Same for intervals

iter <- 10000
# Define number of iterqtions of resampling and computing of pseudo-values

null_CH_Area_intervals <- matrix(nrow = iter,
                                ncol = length(Interval))

for (it in 1:iter) {
  
  Ap.rand <- as.factor(sample(table.aperture$Aperture.Image,
                              replace = T)) 

  CH_Areas <- rep(NA, length(Interval))

  for (i in 1:length(Interval)) {
    
    oo <- match(list.taxa.intervals[[i]], Ap.Image)

    taxa.rand <- Ap.rand[oo]

    y <- pca$x[match(taxa.rand, 
                     levels(Ap.Image)), 2] 
    x <- pca$x[match(taxa.rand, 
                     levels(Ap.Image)), 1]

    CH <- chull(x, y)
    
    CH_xcoo <- x[CH]
    CH_ycoo <- y[CH]
    
    sp.chull <- Polygon(coords = cbind(c(CH_xcoo, CH_xcoo[1]), 
                                       c(CH_ycoo, CH_ycoo[1])), 
                        hole = F)
    
    CH_Areas[i] <- sp.chull@area
    
  }
  
  null_CH_Area_intervals[it,] <- CH_Areas
  
}

#-------------------------------------------------------------------------------
# Bootstrapping data for confidence intervals of convex hull area

bootCHA_biozones <- list()

for (i in 1:length(Biozones)) {
  
  tryCatch({
    
    submorphospace <- pca$x[match(list.taxa.biozones[[i]], 
                                  levels(Ap.Image)),]
    
    CH_Areas <- rep(NA, 1000)
    
    for (j in 1:1000) {
      
      bootspace <- submorphospace[sample(nrow(submorphospace), 
                                         replace = T),]
      
      x <- bootspace[, 1]
      
      y <- bootspace[, 2]
      
      CH <- chull(x, y)
      
      CH_xcoo <- x[CH]
      CH_ycoo <- y[CH]
      
      sp.chull <- Polygon(coords = cbind(c(CH_xcoo, CH_xcoo[1]), 
                                         c(CH_ycoo, CH_ycoo[1])), 
                          hole = F)
      
      CH_Areas[j] <- sp.chull@area
      
    }
    
    bootCHA_biozones[[i]] <- CH_Areas}, 
    
    error = function(a) {return(NA)})
}

sbootCHA_biozones <- lapply(bootCHA_biozones, sort)

upCI_CHA_biozones <- 
  loCI_CHA_biozones <- 
  rep(NA, length(Biozones))

for (i in c(1:6,8:23,25:30)) {
  
  upCI_CHA_biozones[i] <- sbootCHA_biozones[[i]][975]
  loCI_CHA_biozones[i] <- sbootCHA_biozones[[i]][25]
  
}

#-------------------------------------------------------------------------------
# Save objects for further analysis and plots

save(list = c("CH_Area_biozones", 
              "upCI_CHA_biozones",
              "loCI_CHA_biozones",
              "CH_Area_intervals", 
              "null_CH_Area_biozones", 
              "null_CH_Area_intervals",
              "upCI_CHA_biozones",
              "loCI_CHA_biozones"),
     file = "convex_hull_areas.RData")

#-------------------------------------------------------------------------------
# END OF SCRIPT