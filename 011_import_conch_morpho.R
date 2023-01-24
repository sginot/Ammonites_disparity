#-------------------------------------------------------------------------------
# Script for importing conch measurements, and producing a global morphospace.
# By using the same object name "pca" as in script 001, and the same vector name 
# for unique occurences "Ap.Image", one should be able to simply make the 
# morphospace for conchs (this script), and then run scripts 003 to 010, to 
# obtain the desired results. 

# HOWEVER, note that if one does not wish to overwrite all figure / tables /
# R objects produced for the whorl section dataset (i.e. when starting from 
# sript 001), one must change the names of those previous files, or move them to
# a separate folder, that will not be affected by running the scripts again.
# For example, one can create a new "whorl" subfolder in the global folder, 
# and place all result file for the whorl section sbefore running the current 
# script and the following scripts 002 to 010.

# Required packages
NA

#-------------------------------------------------------------------------------
# Read table
table.conch <- read.csv2("Dataset_conch_ratios.csv", 
                            h=T, 
                            sep=";", 
                            dec=",",
                            na.strings = "XX")

Species <- table.conch$Species

Ap.Image <- as.factor(table.conch$Specimen)
  # Use the same name for the vector as in script 002, to allow running the same
  # following scripts

Superfam <- as.factor(table.conch$Superfamily)

#-------------------------------------------------------------------------------
# Define which taxa are present in each time bin (BIOZONE level)

# Create vector of unique biozones names
biozone_vector <- table.conch$Biozone
biozone_factor <- as.factor(biozone_vector)
Biozones <- levels(biozone_factor)
# NOTA BENE: in "Biozones", biozones are ordered ALPHABETICALLY, 
# not chronologically.

# Create empty list to record taxa (more precisely images names)
list.taxa.biozones <- list()

for (i in 1:length(Biozones)) { # Fill the list
  
  l <- which(table.conch$Biozone == Biozones[i])
  list.taxa.biozones[[i]] <- Ap.Image[l]
  
}

# Create empty list to record superfamilies
list.superfamily.biozones <- list()

for (i in 1:length(Biozones)) {  
  
  l <- which(table.conch$Biozone == Biozones[i])
  list.superfamily.biozones[[i]] <- Superfam[l]
  
}

#-------------------------------------------------------------------------------
# Define proper chronological order of biozones. Not used in this script but
# Useful for future scripts, especially when plotting!

ord <- c(24, 11, 7, 15, 18, 13, 16, 26, 5, 6, 12, 23, 28, 9, 10, 3, 20, 4, 14,
         8, 27, 25, 2, 1, 22, 19, 17, 29, 30, 21)
# This vector gathers the indices of biozones to reorder them from
# alphabetical order to chronological order

#-------------------------------------------------------------------------------
# Define which taxa are present in each time bin (INTERVAL level)

# Create vector of unique interval names
interval_factor <- as.factor(table.conch$Interval)
Interval <- levels(interval_factor)

# Make empty list 
list.taxa.intervals <- list()

# Fill the list
for (i in 1:length(Interval)) {  
  
  l <- which(table.conch$Interval == Interval[i])
  list.taxa.intervals[[i]] <- Ap.Image[l]
  
}

# Same for superfamilies
list.superfamily.interval <- list()

for (i in 1:length(Interval)) {
  
  l <- which(table.conch$Interval == Interval[i])
  list.superfamily.interval[[i]] <- Superfam[l]
  
}

#-------------------------------------------------------------------------------
# Define which taxa are present in each time bin (STAGE level)

stage_factor <- as.factor(table.conch$Stage)
Stages <- levels(stage_factor)

# Create empty list
list.taxa.stages <- list()

# Fill with vectors of taxa names for each time bin
for (i in 1:length(Stages)) {
  
  l <- which(table.conch$Stage == Stages[i])
  list.taxa.stages[[i]] <- Ap.Image[l]
  
}

# Same for superfamilies
list.superfamily.stages <- list()

for (i in 1:3){
  l <- which(table.conch$Stage == Stages[i])
  list.superfamily.stages[[i]] <- Superfam[l]
}

#-------------------------------------------------------------------------------
# Create morphospace (PCA of conch ratios)

conch.ratios <- table.conch[, (ncol(table.conch) - 4) : ncol(table.conch)]
  # Select only colums containing the conch ratios values (final five columns)

pca <- prcomp(conch.ratios, 
              scale. = T)

#-------------------------------------------------------------------------------
# Save objects

save(list = ls(pattern = "list."), 
     file = "taxa_lists_time_bins_conchs.RData")

save(list = "pca",
     file = "PCA_conch_ratios.RData")

save(list = c("table.conch", 
              "Species", 
              "Superfam", 
              "Ap.Image", 
              "Interval",
              "Biozones",
              "Stages",
              "interval_factor",
              "biozone_factor",
              "stage_factor",
              "ord"),
     file = "variables_and_data_conch.RData")

#-------------------------------------------------------------------------------
# END OF SCRIPT