#-------------------------------------------------------------------------------
# Script dedicated to the computation of disparity index mean squared euclidean
# distance. In addition, the data is bootstrapped to create confidence intervals

# Required packages:
NA

#-------------------------------------------------------------------------------
# Load objects produced in script 001 and 002
load("PCA_Whorl_sections.RData")
load("taxa_lists_time_bins.RData")
load("variables_and_data_whorl_sections.RData")

#-------------------------------------------------------------------------------
# Compute for each biozone, and store them in vector
