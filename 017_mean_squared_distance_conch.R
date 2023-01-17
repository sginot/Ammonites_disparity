#-------------------------------------------------------------------------------
# Computation of mean squared distance from centroid (and mean pairwise 
# distances) for intervals. In addition, bootstrapping is used for 
# computation of confidence intervals.

# Required packages:
NA

#-------------------------------------------------------------------------------
# Load objects produced in script 011
load("PCA_conch_ratios.RData")
load("taxa_lists_time_bins_conchs.RData")
load("variables_and_data_conch.RData")

#-------------------------------------------------------------------------------
# Compute distance for intervals

cdist_intervals <- 
  mID_intervals <- 
  rep(NA, length(Interval))
# Make empty vector to store results

for (i in 1:length(Interval)) {
  
  tryCatch({

    submorphospace <- pca$x[which(interval_factor == Interval[i]),]

    mID <- mean(dist(submorphospace))

    centroid <- apply(submorphospace, 
                      2,
                      mean)

    mat <- rbind(centroid, submorphospace)

    cdist <- as.matrix(dist(mat))[-1,1]

    mID_intervals[i] <- mID
    
    cdist_intervals[i] <- mean(cdist^2)},
    
    
    error = function(a) {return(NA)})
  
  
}

plot(cdist_intervals, type = "b")

#-------------------------------------------------------------------------------
# Bootstrap data to get confidence intervals

bootcdist_intervals <- 
  bootmID_intervals <- 
  list()

for (i in 1:length(Interval)) {
  
  tryCatch({
    
    submorphospace <-  pca$x[which(interval_factor == Interval[i]),]
    
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

save(list = c("distances_intervals"),
     file = "disparity_distances_conch.RData")

