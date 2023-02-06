#-------------------------------------------------------------------------------
# Plotting curves for the various indices of disparity computed in the previous
# scripts: Sum of ranges, sum of variances, convex hull area, distance to
# centroid, pairwise distance, distance to centroid (not disparity but position)

# Required packages:
library(viridis)
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
# Make global data frames containing data for intervals

global_biozones <- data.frame(disparity_conch_biozones,
                              distances_biozones,
                              CH_Area_biozones,
                              upCI_CHA_biozones,
                              loCI_CHA_biozones,
                              position_biozones)

global_intervals <- data.frame(disparity_conch_intervals,
                               distances_intervals,
                               CH_Area_intervals,
                               upCI_CHA_intervals,
                               loCI_CHA_intervals,
                               position_intervals)

#-------------------------------------------------------------------------------
# Order rows chronologically

global_biozones <- global_biozones[ord,]

Biozones <- Biozones[ord]

#-------------------------------------------------------------------------------
# Write as csv file

write.csv(global_biozones, 
          file = "disparities_biozones_conch.csv")

write.csv(global_intervals, 
          file = "disparities_intervals_conch.csv")

#-------------------------------------------------------------------------------
# Load time bin data and produce vectors to use them as x coordinates in the 
# plots

biozones_ages <- read.csv("biozones_dates.csv",
                          h = T,
                          sep= ",",
                          dec = ".")

intervals_ages <- read.csv("intervals_dates.csv",
                           h = T,
                           sep= ",",
                           dec = ".")

middle_biozones <- biozones_ages$age_base - biozones_ages$duration / 2

middle_intervals <- intervals_ages$age_base - intervals_ages$duration / 2

#-------------------------------------------------------------------------------
# Combined curves plots for biozones

g <- global_biozones
# For simplification purposes

x <- -rev(middle_biozones)
x.base <- -rev(biozones_ages$age_base)
x.coo <- c(x, rev(x))
# For simplification purposes

output_folder <- "../Figures/"

pdf(file = paste(output_folder,
                 "all_curves_biozones_conch.pdf"),
    width = 7,
    height = 14)

layout(matrix(c(1:7, 7), 
              nrow = 8))

par(mar = c(0, 5, 0, 5))

# 1 - SOV (all PCs)

plot(x = x,
     y = g$SOV,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "SOV",
     ylim = c(min(g$SOV_loCI, na.rm = T),
              max(g$SOV_upCI, na.rm = T)))

abline(v = x.base, 
       lty = 2,
       col = "gray")

points(x = x,
       y = g$SOV_upCI,
       pch = "-",
       col = "blue")
points(x = x,
       g$SOV_loCI,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$SOV_upCI,
          rev(g$SOV_loCI)),
        border = NA,
        col = alpha("blue", alpha = 0.2))

# 2 - Mean squared distance to centroid (all PCs)

plot(x = x,
     y = g$Centroid_dist,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "Squared distance to centroid",
     ylim = c(min(g$loCI_cdist, na.rm = T),
              max(g$upCI_cdist, na.rm = T)))

abline(v = x.base, 
       lty = 2,
       col = "gray")

points(x = x,
       y = g$upCI_cdist,
       pch = "-",
       col = "blue")
points(g$loCI_cdist,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$upCI_cdist,
          rev(g$loCI_cdist)),
        border = NA,
        col = alpha("blue", alpha = 0.2)) 


# 3 - Mean pairwise distance (all PCs)

plot(x = x,
     y = g$Pairwise_dist,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "Pairwise distance",
     ylim = c(min(g$loCI_mID, na.rm = T),
              max(g$upCI_mID, na.rm = T)))

abline(v = x.base, 
       lty = 2,
       col = "gray")

points(x = x,
       y = g$upCI_mID,
       pch = "-",
       col = "blue")
points(x = x,
       g$loCI_mID,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$upCI_mID,
          rev(g$loCI_mID)),
        border = NA,
        col = alpha("blue", alpha = 0.2)) 


# 4 - SOR (all PCs)

plot(x = x,
     y = g$SOR,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "SOR",
     ylim = c(min(g$SOR_loCI, na.rm = T),
              max(g$SOR_upCI, na.rm = T)))

abline(v = x.base, 
       lty = 2,
       col = "gray")

points(x = x,
       y = g$SOR_upCI,
       pch = "-",
       col = "blue")
points(x = x,
       y = g$SOR_loCI,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$SOR_upCI,
          rev(g$SOR_loCI)),
        border = NA,
        col = alpha("blue", alpha = 0.2))   

# 5 - Convex hull area (PC1-2 only)

plot(x = x,
     y = g$CH_Area_biozones,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "Convex hull area (PCs 1-2)",
     ylim = c(min(g$loCI_CHA_biozones, na.rm = T),
              max(g$upCI_CHA_biozones, na.rm = T)))

abline(v = x.base, 
       lty = 2,
       col = "gray")

points(x = x,
       y = g$upCI_CHA_biozones,
       pch = "-",
       col = "blue")
points(x = x,
       y = g$loCI_CHA_biozones,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$upCI_CHA_biozones,
          rev(g$loCI_CHA_biozones)),
        border = NA,
        col = alpha("blue", alpha = 0.2)) 

# 6 - Distance of centroid to center of plot

plot(x = x,
     y = g$Position,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "Distance of centroid to center",
     ylim = c(min(g$loCI, na.rm = T),
              max(g$upCI, na.rm = T)))

abline(v = x.base, 
       lty = 2,
       col = "gray")

points(x = x,
       y = g$upCI,
       pch = "-",
       col = "blue")
points(x = x,
       y = g$loCI,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$upCI,
          rev(g$loCI)),
        border = NA,
        col = alpha("blue", alpha = 0.2)) 

# 7 - biozones names

par(mar = c(24, 5, 0, 5))

plot(x, rep(0, 30), 
     type = "n",
     bty = "n",
     yaxt = "n",
     xaxt = "n",
     ylab = "",
     xlab = "")

axis(side = 1, 
     at = x.base,
     labels = F,
     lwd = 2)

axis(side = 1, 
     at = x, 
     labels = Biozones,
     tick = F,
     las = 2)

dev.off()

#-------------------------------------------------------------------------------
# Combined plot for intervals

g <- global_intervals

Intervals <- rev(intervals_ages$names)

x <- -rev(middle_intervals)
x.base <- -rev(intervals_ages$age_base)
x.coo <- c(x, rev(x))
# For simplification purposes

output_folder <- "../Figures/"

pdf(file = paste(output_folder,
                 "all_curves_intervals_conch.pdf"),
    width = 7,
    height = 14)

layout(matrix(c(1:7), 
              nrow = 7))

par(mar = c(0, 5, 0, 5))

# 1 - SOV (all PCs)

plot(x = x,
     y = g$SOV,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "SOV",
     ylim = c(min(g$SOV_loCI, na.rm = T),
              max(g$SOV_upCI, na.rm = T)))

abline(v = x.base, 
       lty = 2,
       col = "gray")

points(x = x,
       y = g$SOV_upCI,
       pch = "-",
       col = "blue")
points(x = x,
       y = g$SOV_loCI,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$SOV_upCI,
          rev(g$SOV_loCI)),
        border = NA,
        col = alpha("blue", alpha = 0.2))

# 2 - Mean squared distance to centroid (all PCs)

plot(x = x,
     y = g$Centroid_dist,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "Squared distance to centroid",
     ylim = c(min(g$loCI_cdist, na.rm = T),
              max(g$upCI_cdist, na.rm = T)))

abline(v = x.base, 
       lty = 2,
       col = "gray")

points(x = x,
       y = g$upCI_cdist,
       pch = "-",
       col = "blue")
points(x = x,
       y = g$loCI_cdist,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$upCI_cdist,
          rev(g$loCI_cdist)),
        border = NA,
        col = alpha("blue", alpha = 0.2)) 


# 3 - Mean pairwise distance (all PCs)

plot(x = x,
     y = g$Pairwise_dist,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "Pairwise distance",
     ylim = c(min(g$loCI_mID, na.rm = T),
              max(g$upCI_mID, na.rm = T)))

abline(v = x.base, 
       lty = 2,
       col = "gray")

points(x = x,
       y = g$upCI_mID,
       pch = "-",
       col = "blue")
points(x = x,
       y = g$loCI_mID,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$upCI_mID,
          rev(g$loCI_mID)),
        border = NA,
        col = alpha("blue", alpha = 0.2)) 


# 4 - SOR (all PCs)

plot(x = x,
     y = g$SOR,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "SOR",
     ylim = c(min(g$SOR_loCI, na.rm = T),
              max(g$SOR_upCI, na.rm = T)))

abline(v = x.base, 
       lty = 2,
       col = "gray")

points(x = x,
       y = g$SOR_upCI,
       pch = "-",
       col = "blue")
points(x = x,
       y = g$SOR_loCI,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$SOR_upCI,
          rev(g$SOR_loCI)),
        border = NA,
        col = alpha("blue", alpha = 0.2))   

# 5 - Convex hull area (PC1-2 only)

plot(x = x,
     y = g$CH_Area_intervals,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "Convex hull area (PCs 1-2)",
     ylim = c(min(g$loCI_CHA_intervals, na.rm = T),
              max(g$upCI_CHA_intervals, na.rm = T)))

abline(v = x.base, 
       lty = 2,
       col = "gray")

points(x = x,
       y = g$upCI_CHA_intervals,
       pch = "-",
       col = "blue")
points(x = x,
       y = g$loCI_CHA_intervals,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$upCI_CHA_intervals,
          rev(g$loCI_CHA_intervals)),
        border = NA,
        col = alpha("blue", alpha = 0.2)) 

# 6 - Distance of centroid to center of plot

plot(x = x,
     y = g$Position,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "Distance of centroid to center",
     ylim = c(min(g$loCI, na.rm = T),
              max(g$upCI, na.rm = T)))

abline(v = x.base, 
       lty = 2,
       col = "gray")

points(x = x,
       y = g$upCI,
       pch = "-",
       col = "blue")
points(x = x,
       y = g$loCI,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$upCI,
          rev(g$loCI)),
        border = NA,
        col = alpha("blue", alpha = 0.2)) 

# 7 - Intervals names

par(mar = c(12, 5, 0, 5))

plot(x, rep(0, 7), 
     type = "n",
     bty = "n",
     yaxt = "n",
     xaxt = "n",
     ylab = "",
     xlab = "")

axis(side = 1, 
     at = x.base,
     labels = F,
     lwd = 2)

axis(side = 1, 
     at = x, 
     labels = Intervals,
     tick = F,
     las = 2)


dev.off()

#-------------------------------------------------------------------------------
# END OF SCRIPT