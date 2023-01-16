#-------------------------------------------------------------------------------
# Plotting curves for the various indices of disparity computed in the previous
# scripts: Sum of ranges, sum of variances, convex hull area, distance to
# centroid, pairwise distance, distance to centroid (not disparity but position)

# Required packages:
library(viridis)
library(scales)

#-------------------------------------------------------------------------------
# Load objects from previous scripts

load("PCA_Whorl_sections.RData")
load("taxa_lists_time_bins.RData")
load("convex_hull_areas.RData")
load("SOV_SOR_disparity.RData")
load("disparity_distances.RData")
load("position_centroid.RData")

#-------------------------------------------------------------------------------
# Make global data frames containing data for biozones and data for intervals

global_biozones <- data.frame(disparity_biozones,
                              distances_biozones,
                              CH_Area_biozones,
                              upCI_CHA_biozones,
                              loCI_CHA_biozones,
                              position_biozones)

global_intervals <- data.frame(disparity_intervals,
                              distances_intervals,
                              CH_Area_intervals,
                              position_intervals)

#-------------------------------------------------------------------------------
# Order rows chronologically

global_biozones <- global_biozones[ord,]

#-------------------------------------------------------------------------------
# Plots

g <- global_biozones
  # For simplification purposes

x.coo <- c(1:30, 30:1)

layout(matrix(1:6, nrow = 6))

par(mar = c(0, 5, 0, 5))

# 1 - SOV (all PCs)

plot(g$SOV,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "SOV",
     ylim = c(min(g$SOV_loCI, na.rm = T),
              max(g$SOV_upCI, na.rm = T)))

points(g$SOV_upCI,
       pch = "-",
       col = "blue")
points(g$SOV_loCI,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$SOV_upCI,
          rev(g$SOV_loCI)),
        border = NA,
        col = alpha("blue", alpha = 0.2))
     
# 2 - Mean squared distance to centroid (all PCs)

plot(g$Centroid_dist,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "Squared distance to centroid",
     ylim = c(min(g$loCI_cdist, na.rm = T),
              max(g$upCI_cdist, na.rm = T)))

points(g$upCI_cdist,
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

plot(g$Pairwise_dist,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "Pairwise distance",
     ylim = c(min(g$loCI_mID, na.rm = T),
              max(g$upCI_mID, na.rm = T)))

points(g$upCI_mID,
       pch = "-",
       col = "blue")
points(g$loCI_mID,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$upCI_mID,
          rev(g$loCI_mID)),
        border = NA,
        col = alpha("blue", alpha = 0.2)) 


# 4 - SOR (all PCs)

plot(g$SOR,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "SOR",
     ylim = c(min(g$SOR_loCI, na.rm = T),
              max(g$SOR_upCI, na.rm = T)))

points(g$SOR_upCI,
       pch = "-",
       col = "blue")
points(g$SOR_loCI,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$SOR_upCI,
          rev(g$SOR_loCI)),
        border = NA,
        col = alpha("blue", alpha = 0.2))   

# 5 - Convex hull area (PC1-2 only)

plot(g$CH_Area_biozones,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "Convex hull area (PCs 1-2)",
     ylim = c(min(g$loCI_CHA_biozones, na.rm = T),
              max(g$upCI_CHA_biozones, na.rm = T)))

points(g$upCI_CHA_biozones,
       pch = "-",
       col = "blue", cex = 3)
points(g$loCI_CHA_biozones,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$upCI_CHA_biozones,
          rev(g$loCI_CHA_biozones)),
        border = NA,
        col = alpha("blue", alpha = 0.2)) 

# 6 - Distance of centroid to center of plot

plot(g$Position,
     type = "b",
     lwd = 3,
     xaxt = "n",
     ylab = "Distance of centroid to center",
     ylim = c(min(g$loCI, na.rm = T),
              max(g$upCI, na.rm = T)))

points(g$upCI,
       pch = "-",
       col = "blue")
points(g$loCI,
       pch = "-",
       col = "blue")

polygon(x.coo,
        c(g$upCI,
          rev(g$loCI)),
        border = NA,
        col = alpha("blue", alpha = 0.2))   
