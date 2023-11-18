# sim-main.R
#
# Main file for conducting simulations for the me-dependent-sampling project
#-------------------------------------------------------------------------------
rm(list=ls())
source('data-functions.R')
source('estimation-functions.R')
#-------------------------------------------------------------------------------
# Ensure all required packages are loaded

packages <- c("tidyverse","mvtnorm",'SuperLearner')
package_list <- rownames(installed.packages())
for (pkg in packages) {  # This should output the actual package names
  if (!(pkg %in% package_list)) {
    stop("Package '", pkg, "' not installed")
  }
}
#-------------------------------------------------------------------------------