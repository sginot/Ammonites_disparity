#-------------------------------------------------------------------------------
# Updated and cleaned up script for importing ammonoid whorl section data from
# JPG images. These are black shapes on white backgrounds, extracted from the
# Primary litterature

# This specific script uses package Momocs to import those shapes as matrices of
# X/Y coordinates.

# Required packages

library(Momocs)
library(rgeos)
library(jpeg)
#-------------------------------------------------------------------------------
# Define raw data folder and files
# Raw data folder is in the upper level folder as the repo

input_folder <- "../images_whorl_sections"

# Files are JPG
fil <- list.files(input_folder, 
                  pattern = "jpg",
                  full.names = T)

#-------------------------------------------------------------------------------
# Import images
WS <- import_jpg(jpg.paths = fil, 
                 auto.notcentered=T)

#-------------------------------------------------------------------------------
# Subsample and superimpose shapes

# Make WS into an "Out" object (Momocs syntax)
OutWS <- Out(WS)

# Subsample all to the same number of points (here 200)
cooWS <- coo_interpolate(OutWS, 
                         n=200)

# Center (to 0,0) and scale (to unit centroid size) all shapes
# Rotation is not necessary as all shapes are schematic and drawn 
# in conventional orientation
centerWS <- coo_center(cooWS)
scaleWS <- coo_scale(centerWS) 

# Find the point of the shape nearest to vertical line starting from center 
ids <- suppressWarnings(coo_intersect_angle(scaleWS, 
                           angle=pi/2)) 

# Slide number of the coordinates so that the first point of each shape is
# homologous (i.e. the upper "tip" of the shape, defined in vector ids)
slideWS <- coo_slide(scaleWS, 
                     id=unlist(ids))

panel(slideWS, cols = "lightgray")

# Run elliptical fourier analysis with default automatic number of harmonics
fou <- efourier(slideWS, 
                norm=F) # 6 harmonics for 99% of variance

#-------------------------------------------------------------------------------
# Use those Fourier coefficients to produce a morphospace (by a pca)

pca <- prcomp(fou$coe)

plot(pca$x[,1:2], 
     asp = 1, 
     pch = 21, 
     bg = "lightgray")
# Visual inspection is necessary if script 001 has been run for
# the first time.
# In some cases an outlier shape may appear, due to slight noise in the 
# interpolation of shapes causing the first landmark not to align
# with those of the other shapes (therefore starting from an non-homologous 
# location). 
# If this happens, simply re-run script 001. Once outlier shape is removed, 
# save the morphospace and load it in further analysis, rather than 
# re-running script 001.

#-------------------------------------------------------------------------------
# Save that PCA object for loading in further scripts

save(pca, file = "PCA_Whorl_sections.RData")

#-------------------------------------------------------------------------------
# END OF SCRIPT