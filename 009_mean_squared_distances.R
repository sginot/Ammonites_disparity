#-------------------------------------------------------------------------------
# Computation of mean squared distance from centroid (and mean pairwise 
# distances) for biozones and intervals. In addition, bootstrapping is used for 
# computation of confidence intervals.

# Required packages:
NA

#-------------------------------------------------------------------------------
# Load objects produced in script 001 and 002
load("PCA_Whorl_sections.RData")
load("taxa_lists_time_bins.RData")
load("variables_and_data_whorl_sections.RData")

#-------------------------------------------------------------------------------
# Compute mean squared distance from centroid and mean pairwise distance for 
# each biozone


cdist_biozones <- 
  mID_biozones <- 
  rep(NA, length(Biozones))
# Make empty vector to store results

for (i in 1:length(Biozones)) {
  
  tryCatch({
    # tryCatch avoids the loop from stopping when there are not enough points
    
    submorphospace <- pca$x[match(list.taxa.biozones[[i]], 
                                  levels(Ap.Image)),]
      # Select data to compute index: all species present in biozone [i] x
      # all PCs
    
    mID <- mean(dist(submorphospace))
      # Mean pairwise distance
    
    centroid <- apply(submorphospace, 
                      2,
                      mean)
      # Compute centroid for all observations of biozone [i]
    
    mat <- rbind(centroid, submorphospace)
      # Make matrix including centroid to be used in function dist
    
    cdist <- as.matrix(dist(mat))[-1,1]
      # Compute distances across mat. Change dist object into matrix, seclect
      # the first column, which contains distances to centroid, but remove the 
      # first row of first column (distance from centroid to centroid = 0)
    
    mID_biozones[i] <- mID

    cdist_biozones[i] <- mean(cdist^2)},
      # Mean squared distance to centroid
    
    error = function(a) {return(NA)})
 
  
}

#-------------------------------------------------------------------------------
# Same for intervals


cdist_intervals <- 
  mID_intervals <- 
  rep(NA, length(Interval))
# Make empty vector to store results

for (i in 1:length(Interval)) {
  
  tryCatch({
    # tryCatch avoids the loop from stopping when there are not enough points
    
    submorphospace <- pca$x[match(list.taxa.intervals[[i]], 
                                  levels(Ap.Image)),]
    # Select data to compute index: all species present in biozone [i] x
    # all PCs
    
    mID <- mean(dist(submorphospace))
    # Mean pairwise distance
    
    centroid <- apply(submorphospace, 
                      2,
                      mean)
    # Compute centroid for all observations of biozone [i]
    
    mat <- rbind(centroid, submorphospace)
    # Make matrix including centroid to be used in function dist
    
    cdist <- as.matrix(dist(mat))[-1,1]
    # Compute distances across mat. Change dist object into matrix, seclect
    # the first column, which contains distances to centroid, but remove the 
    # first row of first column (distance from centroid to centroid = 0)
    
    mID_intervals[i] <- mID
    
    cdist_intervals[i] <- mean(cdist^2)},
    
    
    error = function(a) {return(NA)})
  
  
}

#-------------------------------------------------------------------------------
# Bootstrap data to obtain confidence intervals for the indices

bootcdist_biozones <- 
  bootmID_biozones <- 
  list()
# Make empty list to be filled with pseudo-values for each time bin

for (i in 1:length(Biozones)) {
  
  tryCatch({
    # tryCatch avoids the loop from crashing when there are not enough points
    # to compute ranges or variances
    
    submorphospace <- pca$x[match(list.taxa.biozones[[i]], 
                                  levels(Ap.Image)),]
    # Data frame containing data for all species present in biozone [i]
    
    bootmID <- rep(NA, 1000)
    bootcdist <- rep(NA, 1000)
    # Empty vectors to record pseudo-values (1000 iterations)
    
    for (j in 1:1000) {
      
      bootspace <- submorphospace[sample(nrow(submorphospace), 
                                         replace = T),]
        # Resample with replacement species among species present within
        # biozone [i]
      
      mID <- mean(dist(bootspace))
        # Mean pairwise distance
      
      centroid <- apply(bootspace, 
                        2,
                        mean)
        # Compute centroid for all observations of biozone [i]
      
      mat <- rbind(centroid, bootspace)
        # Make matrix including centroid to be used in function dist
      
      cdist <- as.matrix(dist(mat))[-1, 1]
        # Compute distances across mat. Change dist object into matrix, seclect
        # the first column, which contains distances to centroid, but remove the 
        # first row of first column (distance from centroid to centroid = 0)
      
      bootmID[j] <- mID
      
      bootcdist[j] <- mean(cdist^2)
        # Mean squared distance to centroid
      
      }
    
    bootcdist_biozones[[i]] <- bootcdist
    bootmID_biozones[[i]] <- bootmID}, 
    
    error = function(a) {return(NA)})
}


sbootcdist_biozones <- lapply(bootcdist_biozones, 
                            sort)
sbootmID_biozones <- lapply(bootmID_biozones, 
                              sort)
  # Sort pseudo-values in increasing order within all biozone

upCI_cdist_biozones <- 
  loCI_cdist_biozones <- 
  upCI_mID_biozones <-
  loCI_mID_biozones <-
  rep(NA, length(Biozones))
  # Make empty vectors to store confidence interval (CI) values


for (i in c(1:6,8:23,25:30)) {

  upCI_cdist_biozones[i] <- sbootcdist_biozones[[i]][975]
  loCI_cdist_biozones[i] <- sbootcdist_biozones[[i]][25]
  
  upCI_mID_biozones[i] <- sbootmID_biozones[[i]][975]
  loCI_mID_biozones[i] <- sbootmID_biozones[[i]][25]
  # To obtain bilateral 95% CI, the 50 most extreme values are excluded,
  # therefore 25 values on each side of distribution
}


#-------------------------------------------------------------------------------
# Same for intervals

bootcdist_intervals <- 
  bootmID_intervals <- 
  list()

for (i in 1:length(Interval)) {
  
  tryCatch({

    submorphospace <- pca$x[match(list.taxa.intervals[[i]], 
                                  levels(Ap.Image)),]

    bootmID <- rep(NA, 1000)
    bootcdist <- rep(NA, 1000)

    for (j in 1:1000) {
      
      bootspace <- submorphospace[sample(nrow(submorphospace), 
                                         replace = T),]
 
      mID <- mean(dist(bootspace))
  
      centroid <- apply(bootspace, 
                        2,
                        mean)
   
      mat <- rbind(centroid, bootspace)
   
      cdist <- as.matrix(dist(mat))[-1, 1]
 
      bootmID[j] <- mID
      
      bootcdist[j] <- mean(cdist^2)
  
    }
    
    bootcdist_intervals[[i]] <- bootcdist
    bootmID_intervals[[i]] <- bootmID}, 
    
    error = function(a) {return(NA)})
}


sbootcdist_int <- lapply(bootcdist_intervals, 
                              sort)
sbootmID_int <- lapply(bootmID_intervals, 
                            sort)

upCI_cdist_int <- 
  loCI_cdist_int <- 
  upCI_mID_int <-
  loCI_mID_int <-
  rep(NA, length(Interval))

for (i in 1:length(Interval)) {
  
  upCI_cdist_int[i] <- sbootcdist_int[[i]][975]
  loCI_cdist_int[i] <- sbootcdist_int[[i]][25]
  
  upCI_mID_int[i] <- sbootmID_int[[i]][975]
  loCI_mID_int[i] <- sbootmID_int[[i]][25]

}


#-------------------------------------------------------------------------------
# Save objects for later analyses and plotting

distances_intervals <- data.frame(Centroid_dist = cdist_intervals, 
                                 upCI_cdist = upCI_cdist_int, 
                                 loCI_cdist = loCI_cdist_int,
                                 Pairwise_dist = mID_intervals,
                                 upCI_mID = upCI_mID_int,
                                 loCI_mID = loCI_mID_int)

distances_biozones <- data.frame(Centroid_dist = cdist_biozones, 
                                 upCI_cdist = upCI_cdist_biozones, 
                                 loCI_cdist = loCI_cdist_biozones,
                                 Pairwise_dist = mID_biozones,
                                 upCI_mID = upCI_mID_biozones,
                                 loCI_mID = loCI_mID_biozones,
                                 row.names = Biozones)

save(list = c("distances_biozones", 
              "distances_intervals"),
     file = "disparity_distances.RData")

#-------------------------------------------------------------------------------
# END OF SCRIPT
