#-------------------------------------------------------------------------------
# Script dedicated to the computation of  basic disparity indices:
# Sum of ranges (SOR) and Sum of variances (SOV).
# Both are computed at the various temporal resolutions.
# In addition, the data is bootstrapped to create SOR/SOV confidence intervals.

# Required packages:
NA

#-------------------------------------------------------------------------------
# Load objects produced in script 001 and 002
load("PCA_conch_ratios.RData")
load("taxa_lists_time_bins_conchs.RData")
load("variables_and_data_conch.RData")

#-------------------------------------------------------------------------------
# Compute sum of variances (SOV) for each biozone, and store them in vector

SOV_biozones <- rep(NA,
                    length(Biozones))
# Make empty vector to store results

for (i in 1:length(Biozones)) {
  
  tryCatch({
    # tryCatch avoids the loop from stopping when there are not enough points
    # to compute variance.
    
    submorphospace <- pca$x[which(biozone_factor == Biozones[i]),]
    # Select data to compute variance: all species present in biozone [i] x
    # all PCs
    
    vars <- apply(submorphospace,
                  2,
                  var)
    # Compute variances for each individual PC
    
    SOV_biozones[i] <- sum(vars)},
    # Sum up all those variances
    
    error = function(a) {return(NA)})
  # If variances could not be computed, functions returns NA
  
}

#-------------------------------------------------------------------------------
# Compute sum of variances (SOV) for each interval, and store them in vector

SOV_intervals <- rep(NA,
                     length(Interval))

for (i in 1:length(Interval)) {
  
  tryCatch({
    
    submorphospace <- pca$x[which(interval_factor == Interval[i]),]
    
    vars <- apply(submorphospace,
                  2,
                  var)
    
    SOV_intervals[i] <- sum(vars)}, 
    
    error = function(a) {return(NA)})
  
}

#-------------------------------------------------------------------------------
# Compute sum of variances (SOV) for each stage, and store them in vector


SOV_stages <- rep(NA,
                  length(Stages))

for (i in 1:length(Stages)) {
  
  tryCatch({
    
    submorphospace <- pca$x[which(stage_factor == Stages[i]),]
    
    
    vars <- apply(submorphospace,
                  2,
                  var)
    
    SOV_stages[i] <- sum(vars)}, 
    
    error = function(a) {return(NA)})
  
}

#-------------------------------------------------------------------------------
# Do the same computations for sum of ranges (SOR)
# First with biozones

SOR_biozones <- rep(NA,
                    length(Biozones))
# Make empty vector to store results

for (i in 1:length(Biozones)) {
  
  tryCatch({
    # tryCatch avoids the loop from stopping when there are not enough points
    # to compute ranges (N = 1)
    
    submorphospace <- pca$x[which(biozone_factor == Biozones[i]),]
    # Select data to compute rqnge: all species present in biozone [i] x
    # all PCs
    
    ranges <- apply(apply(submorphospace, 
                          2, 
                          range),
                    2,
                    diff)
    # Find minimum and maximum values for eahc PC, then substract them from 
    # each other to obtain the range value
    
    SOR_biozones[i] <- sum(ranges)},
    # Sum up all those ranges
    
    error = function(a) {return(NA)})
  # If ranges could not be computed, functions returns NA
  
}

#-------------------------------------------------------------------------------
# Same for intervals

SOR_intervals <- rep(NA,
                     length(Interval))

for (i in 1:length(Interval)) {
  
  tryCatch({
    
    submorphospace <- pca$x[which(interval_factor == Interval[i]),]
    
    ranges <- apply(apply(submorphospace, 
                          2, 
                          range),
                    2,
                    diff)
    
    SOR_intervals[i] <- sum(ranges)},
    
    error = function(a) {return(NA)})
  
}

#-------------------------------------------------------------------------------
# Same for stages

SOR_stages <- rep(NA,
                  length(Stages))

for (i in 1:length(Stages)) {
  
  tryCatch({
    
    submorphospace <- pca$x[which(stage_factor == Stages[i]),]
    
    ranges <- apply(apply(submorphospace, 
                          2, 
                          range),
                    2,
                    diff)
    
    SOR_stages[i] <- sum(ranges)},
    
    error = function(a) {return(NA)})
  
}

#-------------------------------------------------------------------------------
# Resampling approaches to produce confidence intervals
# First approach is a classical bootstrapping of points, within each time bin
# Another is a rarefaction approach, basically jacknifing within time bins

#-------------------------------------------------------------------------------
# Bootstrapping to obtain 95% confidence intervals, at the biozone level

bootSOR_biozones <- 
  bootSOV_biozones <- 
  list()
# Make empty lists to be filled with pseudo-values for each biozones

for (i in 1:length(Biozones)) {
  
  tryCatch({
    # tryCatch avoids the loop from crashing when there are not enough points
    # to compute ranges or variances
    
    submorphospace <- pca$x[which(biozone_factor == Biozones[i]),]
    # Data frame containing data for all species present in biozone [i]
    
    bootSOR <- 
      bootSOV<- 
      rep(NA, 1000)
    # Empty vectors to record pseudo-values (1000 iterations)
    
    for (j in 1:1000) {
      
      bootspace <- submorphospace[sample(nrow(submorphospace), 
                                         replace = T),]
      # Resample with replacement species among species present within
      # biozone [i]
      
      ranges <- apply(apply(bootspace, 
                            2, 
                            range),
                      2,
                      diff)
      # Compute ranges in bootstrapped biozone PCA
      
      bootSOR[j] <- sum(ranges)
      # Record jth pseudo-value for sum of ranges
      
      vars <- apply(bootspace,
                    2,
                    var)
      # Compute variances in bootstrapped biozone PCA
      
      bootSOV[j] <- sum(vars)
      # Record jth pseudo-value for sum of variances
      
    }
    
    bootSOR_biozones[[i]] <- bootSOR
    
    bootSOV_biozones[[i]] <- bootSOV}, 
    
    error = function(a) {return(NA)})
}

sbootSOR <- lapply(bootSOR_biozones, sort)
sbootSOV <- lapply(bootSOV_biozones, sort)
# Sort pseudo-values in increasing order within all biozone

upCI_SOR <- 
  loCI_SOR <- 
  upCI_SOV <-
  loCI_SOV <-
  rep(NA, length(Biozones))
# Make empty vectors to store confidence interval (CI) values


for (i in c(1:6,8:23,25:30)) {
  # Biozones 7 and 24 do not have any values due to small sample size
  
  upCI_SOR[i] <- sbootSOR[[i]][975]
  loCI_SOR[i] <- sbootSOR[[i]][25]
  # To obtain bilateral 95% CI, the 50 most extreme values are excluded,
  # therefore 25 values on each side of distribution
  
  upCI_SOV[i] <- sbootSOV[[i]][975]
  loCI_SOV[i] <- sbootSOV[[i]][25]
}

#-------------------------------------------------------------------------------
# Bootstrapping to obtain 95% confidence intervals, at the interval level

bootSOR_intervals <- 
  bootSOV_intervals <- 
  list()
# Make empty lists to be filled with pseudo-values for each interval

for (i in 1:length(Interval)) {
  
  tryCatch({
    
    submorphospace <- pca$x[which(interval_factor == Interval[i]),]
    
    bootSOR <- 
      bootSOV<- 
      rep(NA, 1000)
    
    for (j in 1:1000) {
      
      bootspace <- submorphospace[sample(nrow(submorphospace), 
                                         replace = T),]
      
      ranges <- apply(apply(bootspace, 
                            2, 
                            range),
                      2,
                      diff)
      
      bootSOR[j] <- sum(ranges)
      
      vars <- apply(bootspace,
                    2,
                    var)
      
      bootSOV[j] <- sum(vars)
      
    }
    
    bootSOR_intervals[[i]] <- bootSOR
    
    bootSOV_intervals[[i]] <- bootSOV}, 
    
    error = function(a) {return(NA)})
}

sbootSOR_int <- lapply(bootSOR_intervals, sort)
sbootSOV_int <- lapply(bootSOV_intervals, sort)

upCI_SOR_int <- 
  loCI_SOR_int <- 
  upCI_SOV_int <-
  loCI_SOV_int <-
  rep(NA, length(Interval))

for (i in 1:length(Interval)) {
  
  upCI_SOR_int[i] <- sbootSOR_int[[i]][975]
  loCI_SOR_int[i] <- sbootSOR_int[[i]][25]
  
  upCI_SOV_int[i] <- sbootSOV_int[[i]][975]
  loCI_SOV_int[i] <- sbootSOV_int[[i]][25]
  
}

#-------------------------------------------------------------------------------
# Make data frame with SOR SOV values and their 95% CIs at biozone level

disparity_conch_biozones <- data.frame(SOV = SOV_biozones,
                                 SOV_upCI = upCI_SOV,
                                 SOV_loCI = loCI_SOV,
                                 SOR = SOR_biozones,
                                 SOR_upCI = upCI_SOR,
                                 SOR_loCI = loCI_SOR,
                                 row.names = Biozones)

# Same at interval level

disparity_conch_intervals <- data.frame(SOV = SOV_intervals,
                                  SOV_upCI = upCI_SOV_int,
                                  SOV_loCI = loCI_SOV_int,
                                  SOR = SOR_intervals,
                                  SOR_upCI = upCI_SOR_int,
                                  SOR_loCI = loCI_SOR_int,
                                  row.names = Interval)

#-------------------------------------------------------------------------------
# Save objects for further analyses and plots

save(list = c("SOV_intervals", 
              "SOV_biozones", 
              "SOR_intervals", 
              "SOR_biozones",
              "disparity_biozones", 
              "disparity_intervals"),
     file = "SOV_SOR_disparity_conch.RData")

#-------------------------------------------------------------------------------
# END OF SCRIPT
