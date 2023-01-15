#-------------------------------------------------------------------------------
# Script dedicated to the computation of average displacement index from 
# Guillerme et al (2020). In addition, the data is bootstrapped to create 
# confidence intervals.

# Required packages:
NA

#-------------------------------------------------------------------------------
# Load objects produced in script 001 and 002
load("PCA_Whorl_sections.RData")
load("taxa_lists_time_bins.RData")
load("variables_and_data_whorl_sections.RData")

#-------------------------------------------------------------------------------
# Compute average displacement for each biozone, and store them in vector

AD_biozones <- rep(NA,
                    length(Biozones))
  # Make empty vector to store results

for (i in 1:length(Biozones)) {
  
  tryCatch({
    # tryCatch avoids the loop from stopping when there are not enough points
    # to compute variance.
    
    submorphospace <- pca$x[match(list.taxa.biozones[[i]], 
                                  levels(Ap.Image)),]
      # Select data to compute index: all species present in biozone [i] x
      # all PCs
    
    centroid <- apply(submorphospace, 
                      2,
                      mean)
      # Compute centroid for all observations of biozone [i]
    
    pos2cent <- t(t(submorphospace) - centroid)
      # Compute positions to the centroid of the submorphospace. The original
      # matrix is transposed so that the correct values of the centroid vector
      # are subtracted. The resulting matrix is transposed again, to get back
      # the PCs as columns.
    
    denominator <- sqrt(sum(pos2cent)^2)
      # See Table 3 of Guillerme et al. 2020
    
    numerator <- sqrt(sum(submorphospace)^2)
      # Sqrt of sum of Position to the center of morphospace (ie 0,0,0,0...)

    AD_biozones[i] <- numerator / denominator},
      # Ratio of positions to centroid over position to center of space
    
    error = function(a) {return(NA)})
  # If variances could not be computed, functions returns NA
  
}
