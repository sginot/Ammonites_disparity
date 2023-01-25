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
load("PCA_conch_ratios.RData")
load("taxa_lists_time_bins_conchs.RData")
load("variables_and_data_conch.RData")
     
#-------------------------------------------------------------------------------
# Compute convex hull area for all biozones

CH_Area_biozones <- rep(NA, 
                        length(Biozones))
# Make empty vector to record results

for (i in 1:length(Biozones)) {
  
  y <- pca$x[which(biozone_factor == Biozones[i]), 2]
  x <- pca$x[which(biozone_factor == Biozones[i]), 1]
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
# Compute null distribution of CHA for BIOZONES

iter <- 10000
# Define number of iterations of resampling and computing of pseudo-values

null_CH_Area_biozones <- matrix(nrow = iter,
                                 ncol = length(Biozones))

for (it in 1:iter) {
  
  rand <- sample(nrow(table.conch), 
                 replace = T)
  # Resampled indices of individuals
  
  conch <- (ncol(table.conch) - 4) : ncol(table.conch)
  # Colums containing conch ratios
  
  conch.rand <- table.conch[rand, conch]
  # Randomized conch ratios
  
  pca.rand <- prcomp(conch.rand, 
                     scale. = T)
  # Randomized global morphospace
  
  CH_Areas <- rep(NA, length(Biozones))
  
  for (i in 1:length(Biozones)) {
    
    y <- pca.rand$x[which(biozone_factor == Biozones[i]), 2]
    x <- pca.rand$x[which(biozone_factor == Biozones[i]), 1]
    # Select non-randomized observations in randomized pca
    
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
# Do the same for intervals

CH_Area_intervals <- rep(NA, 
                         length(Interval))

for (i in 1:length(Interval)) {
  
  y <- pca$x[which(interval_factor == Interval[i]), 2]
  x <- pca$x[which(interval_factor == Interval[i]), 1]
  
  CH <- chull(x, y)
  
  CH_xcoo <- x[CH]
  CH_ycoo <- y[CH]
  
  sp.chull <- Polygon(coords = cbind(c(CH_xcoo, CH_xcoo[1]), 
                                     c(CH_ycoo, CH_ycoo[1])), 
                      hole = F)
  
  CH_Area_intervals[i] <- sp.chull@area
  
}

#-------------------------------------------------------------------------------
# Compute null distribution of CHA for intervals

iter <- 10000
# Define number of iterations of resampling and computing of pseudo-values

null_CH_Area_intervals <- matrix(nrow = iter,
                                 ncol = length(Interval))

for (it in 1:iter) {
  
  rand <- sample(nrow(table.conch), 
                 replace = T)
    # Resampled indices of individuals
  
  conch <- (ncol(table.conch) - 4) : ncol(table.conch)
    # Colums containing conch ratios
  
  conch.rand <- table.conch[rand, conch]
    # Randomized conch ratios
  
  pca.rand <- prcomp(conch.rand, 
                     scale. = T)
    # Randomized global morphospace
  
  CH_Areas <- rep(NA, length(Interval))
  
  for (i in 1:length(Interval)) {
    
    y <- pca.rand$x[which(interval_factor == Interval[i]), 2]
    x <- pca.rand$x[which(interval_factor == Interval[i]), 1]
      # Select non-randomized observations in randomized pca
    
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
# Bootstrapping data at interval level

bootCHA_intervals <- list()

for (i in 1:length(Interval)) {
  
  tryCatch({
    
    submorphospace <- pca$x[which(interval_factor == Interval[i]), ]
    
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
    
    bootCHA_intervals[[i]] <- CH_Areas}, 
    
    error = function(a) {return(NA)})
}

sbootCHA_intervals <- lapply(bootCHA_intervals, sort)

upCI_CHA_intervals <- 
  loCI_CHA_intervals <- 
  rep(NA, length(Interval))

for (i in 1:length(Interval)) {
  
  upCI_CHA_intervals[i] <- sbootCHA_intervals[[i]][975]
  loCI_CHA_intervals[i] <- sbootCHA_intervals[[i]][25]
  
}

#-------------------------------------------------------------------------------
# Bootstrapping data at biozones level

bootCHA_biozones <- list()

for (i in 1:length(Biozones)) {
  
  tryCatch({
    
    submorphospace <- pca$x[which(biozone_factor == Biozones[i]), ]
    
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

for (i in 1:length(Biozones)) {
  
  upCI_CHA_biozones[i] <- sbootCHA_biozones[[i]][975]
  loCI_CHA_biozones[i] <- sbootCHA_biozones[[i]][25]
  
}

#-------------------------------------------------------------------------------
# Save objects

save(list = c("CH_Area_intervals",
              "upCI_CHA_intervals",
              "loCI_CHA_intervals",
              "CH_Area_biozones",
              "upCI_CHA_biozones",
              "loCI_CHA_biozones",
              "null_CH_Area_intervals",
              "null_CH_Area_biozones"),
     file = "convex_hull_areas_conch.RData")
