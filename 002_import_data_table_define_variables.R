#-------------------------------------------------------------------------------
# The data table used in the study records data associated with all specimens of
# Devonian Ammonoid from Morocco. Data variables include the systematics of 
# specimens, their time bins, and corresponding illustration.

#-------------------------------------------------------------------------------
# Read table
table.aperture <- read.csv2("Dataset_whorl_profiles.csv", 
                            h=T, 
                            sep=";", 
                            dec=",",
                            na.strings = "XX")

# Change some character variables to factors
Ap.Image <- as.factor(table.aperture$Aperture.Image)
  # Number of factor levels should match number of images (N = 582)

Species <- table.aperture$Species
Superfam <- as.factor(table.aperture$Superfamily)

#-------------------------------------------------------------------------------
# Define which taxa are present in each time bin (BIOZONE level)

# Create vector of unique biozones names
biozone_vector <- table.aperture$Biozone
biozone_factor <- as.factor(biozone_vector)
Biozones <- levels(biozone_factor)
  # NOTA BENE: in "Biozones", biozones are ordered ALPHABETICALLY, 
  # not chronologically.

# Create empty list to record taxa (more precisely images names)
list.taxa.biozones <- list()

for (i in 1:length(Biozones)) { # Fill the list
  
  l <- which(table.aperture$Biozone == Biozones[i])
  list.taxa.biozones[[i]] <- Ap.Image[l]
  
  }

# Create empty list to record superfamilies
list.superfamily.biozones <- list()

for (i in 1:length(Biozones)) {  
  
  l <- which(table.aperture$Biozone == Biozones[i])
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
interval_factor <- as.factor(table.aperture$Interval)
Interval <- levels(interval_factor)

# Make empty list 
list.taxa.intervals <- list()

# Fill the list
for (i in 1:length(Interval)) {  
  
  l <- which(table.aperture$Interval == Interval[i])
  list.taxa.intervals[[i]] <- Ap.Image[l]
  
  }

# Same for superfamilies
list.superfamily.interval <- list()

for (i in 1:length(Interval)) {
  
  l <- which(table.aperture$Interval == Interval[i])
  list.superfamily.interval[[i]] <- Superfam[l]
  
  }

#-------------------------------------------------------------------------------
# Define which taxa are present in each time bin (STAGE level)

stage_factor <- as.factor(table.aperture$Stage)
Stages <- levels(stage_factor)

# Create empty list
list.taxa.stages <- list()

# Fill with vectors of taxa names for each time bin
for (i in 1:length(Stages)) {
  
  l <- which(table.aperture$Stage == Stages[i])
  list.taxa.stages[[i]] <- Ap.Image[l]
  
  }

# Same for superfamilies
list.superfamily.stages <- list()

for (i in 1:3){
  l <- which(table.aperture$Stage == Stages[i])
  list.superfamily.stages[[i]] <- Superfam[l]
}

#-------------------------------------------------------------------------------
# Save objects

save(list = ls(pattern = "list."), 
     file = "taxa_lists_time_bins.RData")

save(list = c("table.aperture", 
              "Species", 
              "Superfam", 
              "Ap.Image", 
              "Interval",
              "Biozones",
              "Stages",
              "ord"),
     file = "variables_and_data_whorl_sections.RData")

#-------------------------------------------------------------------------------
# END OF SCRIPT