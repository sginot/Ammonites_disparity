#-------------------------------------------------------------------------------
# Script dedicated to the plotting of results from Whalen's et al (2020) 
# approach. The aim is to show the chronological changes in convex hull area
# and the null model values for the same index

# Required packages:
library(scales)
library(viridis)

#-------------------------------------------------------------------------------
# Load objects produced in scripts 002 and 006
load("convex_hull_areas.RData")
load("variables_and_data_whorl_sections.RData")

#-------------------------------------------------------------------------------
# Compute the residual convex hull area, by subtracting the median of null model
# from actual values and from null model distribution

null_median <- apply(null_CH_Area_biozones[,ord], 
                     2, 
                     median)

q95 <- apply(null_CH_Area_biozones[,ord], 
             2, 
             quantile, 
             probs = 0.95)

q5 <- apply(null_CH_Area_biozones[,ord], 
            2, 
            quantile,
            probs = 0.05)

q75 <- apply(null_CH_Area_biozones[,ord], 
             2, 
             quantile, 
             probs = 0.75)

q25 <- apply(null_CH_Area_biozones[,ord], 
             2, 
             quantile,
             probs = 0.25)

CHA <- CH_Area_biozones[ord]

res_CHA <- CHA - null_median

res_q5 <- q5 - null_median
res_q25 <- q25 - null_median
res_q75 <- q75 - null_median
res_q95 <- q95 - null_median

#-------------------------------------------------------------------------------
# Basic plot

output_folder <- "../Figures/Convex_hull_area_disparity_curve/"
  # Define output folder containing figures.

pdf(file = paste(output_folder, 
                 "Convex_hull_area_biozones.pdf"),
    width = 7,
    height = 7)
  # Create image pdf file amd set the size of plotting area

plot(1:30, 
     CHA, 
     type = "b", 
     lwd = 3, 
     pch = 19, 
     ylim = c(-0.1, 1))
  # Plot measured values for convex hull area

lines(1:30,
      null_median,
      col = cols[1],
      lty = 2)
  # Add median of the null model inspired from Whalen's et al paper

polygon(c(1:30, 30:1), 
        c(q5, rev(q95)), 
        col = alpha(cols[1], 
                    alpha = 0.2))
  # Add colored area upper and lower bound are defined as 5 and 95 percentiles 

polygon(c(1:30, 30:1), 
        c(q25, rev(q75)), 
        col = alpha(cols[1], 
                    alpha = 0.2),
        border = NA)
  # Add colored area upper and lower bound are defined as 25 and 75 percentiles 

legend("topleft",
       lty = c(1, 2, 1, 1),
       col = c("black", 
               "black", 
               alpha(cols[1], alpha = 0.2),
               alpha(cols[1], alpha = 0.4)),
       lwd = c(3, 1, 3, 3),
       legend = c("Measured convex hull area",
                  "Median convex hull area of null model",
                  "90% envelope of null model",
                  "50% envelope of null model"))

dev.off()

#-------------------------------------------------------------------------------
# Residual convex hull area plot

pdf(file = paste(output_folder, 
                 "Residual_convex_hull_area_biozones.pdf"),
    width = 7,
    height = 7)

plot(1:30, 
     res_CHA, 
     type = "b", 
     lwd = 3, 
     pch = 19, 
     ylim = c(-0.3, 0.3))

polygon(c(1:30, 30:1), 
        c(res_q5, rev(res_q95)), 
        col = alpha(cols[1], 
                    alpha = 0.2))

polygon(c(1:30, 30:1), 
        c(res_q25, rev(res_q75)), 
        col = alpha(cols[1], 
                    alpha = 0.2),
        border = NA)

legend("topleft",
       lty = c(1, 1, 1),
       col = c("black",
               alpha(cols[1], alpha = 0.2),
               alpha(cols[1], alpha = 0.4)),
       lwd = c(3, 3, 3),
       legend = c("Residual convex hull area",
                  "Residual 90% envelope of null model",
                  "Residual 50% envelope of null model"))

dev.off()

#-------------------------------------------------------------------------------
# Do the same for intervals
# Compute the residual convex hull area, by subtracting the median of null model
# from actual values and from null model distribution

null_median_int <- apply(null_CH_Area_intervals, 
                     2, 
                     median)

q95_int <- apply(null_CH_Area_intervals, 
             2, 
             quantile, 
             probs = 0.95)

q5_int <- apply(null_CH_Area_intervals, 
            2, 
            quantile,
            probs = 0.05)

q75_int <- apply(null_CH_Area_intervals, 
             2, 
             quantile, 
             probs = 0.75)

q25_int <- apply(null_CH_Area_intervals, 
             2, 
             quantile,
             probs = 0.25)

CHA_int <- CH_Area_intervals

res_CHA_int <- CHA_int - null_median_int

res_q5_int <- q5_int - null_median_int
res_q25_int <- q25_int - null_median_int
res_q75_int <- q75_int - null_median_int
res_q95_int <- q95_int - null_median_int

#-------------------------------------------------------------------------------
# Basic plot

pdf(file = paste(output_folder, 
                 "Convex_hull_area_intervals.pdf"),
    width = 7,
    height = 7)

plot(1:7, 
     CHA_int, 
     type = "b", 
     lwd = 3, 
     pch = 19, 
     ylim = c(-0.1, 1))

lines(1:7,
      null_median_int,
      col = cols[1],
      lty = 2)

polygon(c(1:7, 7:1), 
        c(q5_int, rev(q95_int)), 
        col = alpha(cols[1], 
                    alpha = 0.2))

polygon(c(1:7, 7:1), 
        c(q25_int, rev(q75_int)), 
        col = alpha(cols[1], 
                    alpha = 0.2),
        border = NA)

legend("topleft",
       lty = c(1, 2, 1, 1),
       col = c("black", 
               "black", 
               alpha(cols[1], alpha = 0.2),
               alpha(cols[1], alpha = 0.4)),
       lwd = c(3, 1, 3, 3),
       legend = c("Measured convex hull area",
                  "Median convex hull area of null model",
                  "90% envelope of null model",
                  "50% envelope of null model"))

dev.off()

#-------------------------------------------------------------------------------
# Residual plot

pdf(file = paste(output_folder, 
                 "Residual_Convex_hull_area_intervals.pdf"),
    width = 7,
    height = 7)

plot(1:7, 
     res_CHA_int, 
     type = "b", 
     lwd = 3, 
     pch = 19, 
     ylim = c(-0.2, 0.2))

polygon(c(1:7, 7:1), 
        c(res_q5_int, rev(res_q95_int)), 
        col = alpha(cols[1], 
                    alpha = 0.2))

polygon(c(1:7, 7:1), 
        c(res_q25_int, rev(res_q75_int)), 
        col = alpha(cols[1], 
                    alpha = 0.2),
        border = NA)

legend("topleft",
       lty = c(1, 1, 1),
       col = c("black",
               alpha(cols[1], alpha = 0.2),
               alpha(cols[1], alpha = 0.4)),
       lwd = c(3, 3, 3),
       legend = c("Residual convex hull area",
                  "Residual 90% envelope of null model",
                  "Residual 50% envelope of null model"))

dev.off()

#-------------------------------------------------------------------------------
# END OF SCRIPT