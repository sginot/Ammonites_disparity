#-------------------------------------------------------------------------------
# Script dedicated to the computation of average displacement index from 
# Guillerme et al (2020). In addition, the data is bootstrapped to create 
# confidence intervals.

# Required packages:
NA

#-------------------------------------------------------------------------------
# Load objects produced in script 001 and 002
load("PCA_conch_ratios.RData")
load("taxa_lists_time_bins_conchs.RData")
load("variables_and_data_conch.RData")

#-------------------------------------------------------------------------------
# Compute distance of centroid to center for intervals.

pos_intervals <- rep(NA,
                     length(Interval))
# Make empty vector to store results

for (i in 1:length(Interval)) {
  
  tryCatch({
    
    submorphospace <- pca$x[which(interval_factor == Interval[i]),]
    
    centroid <- apply(submorphospace, 
                      2,
                      mean)
    # Compute centroid for all observations of biozone [i]
    
    center <- rep(0, 
                  length(centroid))
    # Center of global morphospace is 0,0,0,0,...
    
    pos_intervals[i] <- dist(rbind(centroid, 
                                   center))},
    # Ratio of positions to centroid over position to center of space
    
    error = function(a) {return(NA)})
  # If indices could not be computed, functions returns NA
  
}

#-------------------------------------------------------------------------------
# Bootstrap distance of centroid to center for interval

bootpos_int <- list()

for (i in 1:length(Interval)) {
  
  tryCatch({
    
    submorphospace <- pca$x[which(interval_factor == Interval[i]),]
    
    bootpos <- rep(NA, 1000)
    
    for (j in 1:1000) {
      
      bootspace <- submorphospace[sample(nrow(submorphospace), 
                                         replace = T),]
      
      centroid <- apply(bootspace, 
                        2,
                        mean)
      
      center <- rep(0, 
                    length(centroid))
      
      bootpos[j] <- dist(rbind(centroid, 
                               center))
      
    }
    
    bootpos_int[[i]] <- bootpos}, 
    
    error = function(a) {return(NA)})
}


sbootpos_int <- lapply(bootpos_int, 
                       sort)

upCI_pos_int <- 
  loCI_pos_int <- 
  rep(NA, length(Interval))

for (i in 1:length(Interval)) {
  
  upCI_pos_int[i] <- sbootpos_int[[i]][975]
  loCI_pos_int[i] <- sbootpos_int[[i]][25]
  # To obtain bilateral 95% CI, the 50 most extreme values are excluded,
  # therefore 25 values on each side of distribution
}

#-------------------------------------------------------------------------------
# Save objects for later analyses and plotting

position_intervals <- data.frame(Position = pos_intervals, 
                                 upCI = upCI_pos_int, 
                                 loCI = loCI_pos_int)


save(list = c("position_intervals"),
     file = "position_centroid_conch.RData")

