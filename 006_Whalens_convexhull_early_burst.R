#-------------------------------------------------------------------------------
# Compute convex hull area of morphospace per biozones and per interval.
# This is inspired from Whalen's et al. (2020) paper in Science Advances.
# Note that this restricts disparity only to the first TWO PCs of the PCA, due
# to the implementation of function 'chull'.

# In addtion to the convex hull area, we re-implement Whalen's et al approach
# to test for eqrly bust, by producing null distributions of the convex hull
# area across time. To do so, morphologies are randomly assigned to the 
# species, and the sample stays the same.

# Required packages:
library(sp)

#-------------------------------------------------------------------------------
# Reassign existing morphologies randomly with replacement to the observations
# which are contained in the data table.

table.aperture <- read.csv2("Dataset_whorl_profiles.csv", 
                            h=T, 
                            sep=";", 
                            dec=",",
                            na.strings = "XX")

Ap_Image_rand <- sample(table.aperture$Aperture.Image,
                        replace = T)
  # Assign random image name to each observation

ltb <- list()

for (i in 1:length(Biozones)) { # Fill the list
  
  l <- which(table.aperture$Biozone == Biozones[i])
  ltb[[i]] <- Ap_Image_rand[l]
  
}
