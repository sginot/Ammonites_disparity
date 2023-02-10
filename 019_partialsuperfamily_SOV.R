#-------------------------------------------------------------------------------
# Calculate SOV per time bin and per superfamily and plot as barplots. 
# Also calculate partial disparity in the fashion used by Foote 1993.
# Note this will apply the script for conch ratio data only.

# Required packages:
library(scales)

#-------------------------------------------------------------------------------
# Load objects from previous scripts

load("PCA_conch_ratios.RData")
load("taxa_lists_time_bins_conchs.RData")
load("convex_hull_areas_conch.RData")
load("SOV_SOR_disparity_conch.RData")
load("disparity_distances_conch.RData")
load("position_centroid_conch.RData")
load("variables_and_data_conch.RData")
#-------------------------------------------------------------------------------
# Define output folder

fdirs <- c('../Figures')

for(s in fdirs) if(!dir.exists(s)) dir.create(s)

output_folder <- '../Figures/'
#-------------------------------------------------------------------------------
# Define color palette

cols <- c("#F7941D", 
          "#26AAE1", 
          "#283991", 
          "#662D91", 
          "#03A753", 
          "#8CC63F", 
          "#ED247E")
#-------------------------------------------------------------------------------
# Load files with ages for time bins

intervals_ages <- read.csv("intervals_dates.csv",
                           h = T,
                           sep= ",",
                           dec = ".")

biozones_ages <- read.csv("biozones_dates.csv",
                          h = T,
                          sep= ",",
                          dec = ".")

#-------------------------------------------------------------------------------
# Compute and plot "Partial SOV" at interval level 

PSOV <- list()

for (i in 1:7){
  # 7 is number of intervals
  
  tryCatch({
    
    submorphospace <- pca$x[which(interval_factor == Interval[i]), ]
    
    parts <- rep(NA, 7)
    
    for (j in 1:7) {
      # 7 is number of superfamilies
      
      super <- list.superfamily.interval[[i]]
      
      m <- submorphospace[which(super == levels(super)[j]), ]
      
      parts[j] <- sum(apply(m, 2, var))
      
    }
    
    PSOV[[i]] <- parts},
    
    error = function(a){return(0)})
}

mat <- matrix(NA, ncol = 7, 
              nrow = 7)

for(i in 1:7) {mat[,i] <- PSOV[[i]]}

mat[which(is.na(mat))] <- 0

pdf(paste(output_folder,
          "Partial_disparity_interval_conch.pdf",
          sep = ""))

barplot(mat, 
        ylab = "SOV",
        width = rev(intervals_ages[,3]),
        space = 0,
        names.arg = rev(intervals_ages[,1]),
        col = cols,
        las = 2)

legend("topleft", 
       lty = 1, 
       lwd = 10, 
       cex = 0.7, 
       col = cols, 
       legend = c("Agoniatitoidea", 
                  "Anarcestoidea", 
                  "Gephuroceratoidea",
                  "Mimagoniatitoidea", 
                  "Mimosphinctoidea",
                  "Pharciceratoidea",
                  "Tornoceratoidea"))
 
dev.off()

#-------------------------------------------------------------------------------
# Compute and plot "Partial SOV" at biozone level 

PSOV <- list()

for (i in 1:30){
  # 30 is number of biozones
  
    submorphospace <- pca$x[which(biozone_factor == Biozones[i]), ]
    
    parts <- rep(NA, 7)
    
    for (j in 1:7) {
      # 7 is number of superfamilies
      
      super <- list.superfamily.biozones[[i]]
      
      tryCatch({
      
      m <- submorphospace[which(super == levels(super)[j]), ]
      
      parts[j] <- sum(apply(m, 2, var))},
      
      error = function(a){return(0)})
      
    }
    
    PSOV[[i]] <- parts
}
    
mat <- matrix(NA, ncol = 30, 
              nrow = 7)

for(i in 1:30) {
  
  tryCatch({
  
  mat[,i] <- PSOV[[i]]},
  
  error = function(a){return(0)})
  
}

mat[which(is.na(mat))] <- 0

mat <- mat[,ord]

pdf(paste(output_folder,
          "Partial_disparity_biozones_conch.pdf",
          sep = ""),
    width = 11,
    height = 14)

par(mar = c(22, 4, 1, 1))
barplot(mat, 
        ylab = "SOV",
        width = rev(biozones_ages[,3]),
        space = 0,
        names.arg = rev(biozones_ages[,1]),
        col = cols,
        las = 2)

legend("topleft", 
       lty = 1, 
       lwd = 10, 
       cex = 1.4, 
       col = cols, 
       legend = c("Agoniatitoidea", 
                  "Anarcestoidea", 
                  "Gephuroceratoidea",
                  "Mimagoniatitoidea", 
                  "Mimosphinctoidea",
                  "Pharciceratoidea",
                  "Tornoceratoidea"))

dev.off()

#-------------------------------------------------------------------------------
# Foote (1993) similar figure of partial disparities
# First at interval level

MD <- rep(NA, 7)

PD <- list()

for (i in 1:7){
  
  tryCatch({
    
    submorphospace <- pca$x[which(interval_factor == Interval[i]),]
    
    N <- nrow(submorphospace)
    
    centroid <- apply(submorphospace, 
                      2, 
                      mean)
    
    mat <- rbind(centroid, 
                 submorphospace)
    
    distcent <- dist(mat)[1:N]
    
    parts <- rep(NA, 7)
    
    for (j in 1:7) {
      
      super <- list.superfamily.interval[[i]]
      
      parts[j] <- sum(distcent[which(super == levels(super)[j])]^2)/(N-1)
      
    }
    
    PD[[i]] <- parts
    
    MD[i] <- sum(distcent^2)/(N-1)}, 
    
    error = function(a){return(NA)})
  
}

mat <- matrix(NA, 
              ncol = 7,
              nrow = 7)

for(i in 1:7) {mat[, i] <- PD[[i]]}

pdf(paste(output_folder,
          "Partial_disparity_Foote_style_interval_conch.pdf",
          sep = ""))

barplot(mat, 
        ylab = "Sum of squared Euclidean distance from centroid",
        width = rev(intervals_ages[,3]),
        space = 0,
        names.arg = rev(intervals_ages[,1]),
        col = cols,
        las = 2)

legend("topleft", 
       lty = 1, 
       lwd = 10, 
       cex = 0.7, 
       col = cols, 
       legend = c("Agoniatitoidea", 
                  "Anarcestoidea", 
                  "Gephuroceratoidea",
                  "Mimagoniatitoidea", 
                  "Mimosphinctoidea",
                  "Pharciceratoidea",
                  "Tornoceratoidea"))

dev.off()

#-------------------------------------------------------------------------------
# Foote (1993) similar figure of partial disparities
# Then at biozone level

MD <- rep(NA, 30)

PD <- list()

for (i in 1:30) {
  
    submorphospace <- pca$x[which(biozone_factor == Biozones[i]),]
    
    N <- nrow(submorphospace)
    
    centroid <- apply(submorphospace, 
                      2, 
                      mean)
    
    mat <- rbind(centroid, 
                 submorphospace)
    
    distcent <- dist(mat)[1:N]
    
    parts <- rep(NA, 7)
    
    for (j in 1:7) {
      
      super <- list.superfamily.biozones[[i]]
      
      tryCatch({
        
      parts[j] <- sum(distcent[which(super == levels(super)[j])]^2)/(N-1)}, 
      
      error = function(a) {return(NA)})
      
    }
    
    PD[[i]] <- parts
    
    MD[i] <- sum(distcent^2)/(N-1)
  
}

mat <- matrix(NA, 
              ncol = 30, 
              nrow = 7)

for (i in 1:30) {mat[,i] <- PD[[i]]}

mat <- mat[,ord]

pdf(paste(output_folder,
          "Partial_disparity_Foote_style_biozones_conch.pdf",
          sep = ""),
    width = 11,
    height = 14)

par(mar = c(22, 4, 1, 1))
barplot(mat, 
        ylab = "Sum of squared Euclidean distance from centroid",
        width = rev(biozones_ages[,3]),
        space = 0,
        names.arg = rev(biozones_ages[,1]),
        col = cols,
        las = 2)

legend("topleft", 
       lty = 1, 
       lwd = 10, 
       cex = 1.4, 
       col = cols, 
       legend = c("Agoniatitoidea", 
                  "Anarcestoidea", 
                  "Gephuroceratoidea",
                  "Mimagoniatitoidea", 
                  "Mimosphinctoidea",
                  "Pharciceratoidea",
                  "Tornoceratoidea"))

dev.off()
