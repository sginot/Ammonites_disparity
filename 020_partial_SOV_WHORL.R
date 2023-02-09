#-------------------------------------------------------------------------------
# Calculate SOV per time bin and per superfamily (ie "partial disparity"), and
# plot as barplots. Note thqt "partial" Sum of variances is perhaps not the best
# way to describe it. This script applies to whorl sections data.

# Required packages:
library(scales)

#-------------------------------------------------------------------------------
# Load objects from previous scripts

load("PCA_Whorl_sections.RData")
load("taxa_lists_time_bins.RData")
load("convex_hull_areas.RData")
load("SOV_SOR_disparity.RData")
load("disparity_distances.RData")
load("position_centroid.RData")
load("variables_and_data_whorl_sections.RData")
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

intervals_ages <- read.csv("intervals_dates.csv",
                           h = T,
                           sep= ",",
                           dec = ".")

PSOV <- list()

for (i in 1:7){
  # 7 is number of intervals
  
  tryCatch({
    
    submorphospace <- pca$x[match(list.taxa.intervals[[i]], 
                                  levels(Ap.Image)), ]
    
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

for (i in 1:7) {mat[,i] <- PSOV[[i]]}

mat[which(is.na(mat))] <- 0

pdf(paste(output_folder,
          "Partial_disparity_interval_whorl.pdf",
          ""))

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
       cex = 1.1, 
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
# Do the same thing at the biozones level 


bizones_ages <- read.csv("biozones_dates.csv",
                         h = T,
                         sep= ",",
                         dec = ".")

ls.sup.ord <- list.superfamily.biozones[ord]
ls.taxa.ord <- list.taxa.biozones[ord]

PSOV <- list()

for (i in 1:30){
  # 30 is number of biozones
  
  submorphospace <- pca$x[match(ls.taxa.ord[[i]], 
                                levels(Ap.Image)), ]
  
  parts <- rep(NA, 7)
  
  for (j in 1:7) {
    # 7 is number of superfamilies
    
    super <- ls.sup.ord[[i]]
    
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

pdf(paste(output_folder,
          "Partial_disparity_biozones_whorl.pdf",
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
