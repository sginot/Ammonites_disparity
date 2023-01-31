#-------------------------------------------------------------------------------
# Calculate SOV per time bin and per superfamily (ie "partial disparity"), and
# plot as barplots.

# Required packages:
library(scales)
#-------------------------------------------------------------------------------

intervals_ages <- read.csv("intervals_dates.csv",
                           h = T,
                           sep= ",",
                           dec = ".")

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

pdf("Partial_disparity_interval_new.pdf")

barplot(mat, 
        ylab = "SOV",
        width = rev(intervals_ages[,3]),
        space = 0,
        names.arg = rev(intervals_ages[,1]),
        col = c("darkorange", "deepskyblue", "darkblue",
                           "purple4",
                           "darkgreen", 
                           "violetred", 
                           "yellowgreen"),
        las = 2)

legend("topleft", 
       lty = 1, 
       lwd = 2, 
       cex = 0.7, 
       col = c("darkgreen",
                     "purple4",
                     "deepskyblue", 
                     "darkorange", 
                     "violetred", 
                     "darkblue",
                     "yellowgreen"), 
       legend = c("Mimosphinctaceae", 
                "Mimagoniatitaceae", 
                "Anarcestaceae", 
                "Agoniatitaceae", 
                "Tornocerataceae", 
                "Gephurocerataceae", 
                "Pharcicerataceae"))
 
dev.off()
