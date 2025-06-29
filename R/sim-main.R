# sim-main.R -- cluster version
#
# Main file for conducting simulations for the me-dependent-sampling project
# Execute this from the command line with the syntax
# sbatch --array=1-<# of arrays> sim-main-sbatch.sh <# of sims per array>
#-------------------------------------------------------------------------------
rm(list=ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('data-functions.R')
source('estimation-functions.R')
source('sim-functions.R')
source('sampling_funcs.R')
#-------------------------------------------------------------------------------
# Check if on cluster to determine library path

if (Sys.info()['sysname'] == 'Linux') { # if on cluster, specify library path
  .libPaths("~/R_packages")
}
#-------------------------------------------------------------------------------
# Ensure all required packages are installed

packages <- c("tidyverse","mvtnorm",'SuperLearner','ggridges',
              'doParallel','foreach','AIPW')
package_list <- rownames(installed.packages())
for (pkg in packages) {  # This should output the actual package names
  if (!(pkg %in% package_list)) {
    stop("Package '", pkg, "' not installed")
  }
}
#-------------------------------------------------------------------------------
# Handle command line arguments

command_line_args <- commandArgs(trailingOnly = TRUE)

array_index <- as.numeric(command_line_args[1]) # ids which job this is in array
nsim_per_job <- as.numeric(command_line_args[2]) # number of sims per job
full_job_id <- as.character(command_line_args[3]) # id for the whole array
#-------------------------------------------------------------------------------
# Need a few of the packages loaded
library(tidyverse)
library(SuperLearner)
library(ggridges)
library(doParallel)
library(foreach)
library(AIPW)
# ------------------------------------------------------------------------------
# Set up relevant paths for output
flnm <- paste0('sim_results_job',array_index,'_array',full_job_id,'.csv')

simdir <- '../output/' # directory
fullpath <- paste(simdir,flnm,sep='') # construct path to final file name
# ------------------------------------------------------------------------------
# Setup parameters that remain constant across simulations

p <- 3 # number of covariates
delta <-  c(-0.5,-0.5,-0.5) # c(-0.1,-0.6,-0.9) # rep(-0.5,p) #
nu <- c(0.1,-0.1,0.1)
beta <- c(1,2,-2) # outcome model main effects
tau <- 1 # constant treatment effect
gamma <- rep(1,p) # interaction effects with tmt c(0.5,2.1,-1.2)
eta <- c(0.5,-0.25,-0.25,0.75,0.25)  # c(0.6,-0.2,0.8,0.1,-0.3) c(0.5,-0.25,-0.25,1.5,0.2) # c(0.5,-0.25,-0.25,1,0.2) 
#-------------------------------------------------------------------------------
# Setup parameters that vary across simulations

rho <- c(0.1, 0.2, 0.3, 0.4, 0.5) # relative size of the validation data
n <- c(1000, 2500, 5000) # sample size
omega <- 0.15
#-------------------------------------------------------------------------------
Y <- 'Y' ; A <- 'A' ; Ystar <- 'Ystar' ; Astar <- 'Astar' ; X <- c('X1','X2','X3')
sl.lib <- c('SL.glm','SL.glm.interaction','SL.gam')
est_r <- c(TRUE) # true if want to estimate P(R=1|Z) (known by researcher)
#-------------------------------------------------------------------------------
# Main simulation
nsim <- nsim_per_job # number of simulation iterations

# Vectorize every combination of rho, n, and omega and loop over these combos
param_combos <- expand.grid(rho=rho, n=n, omega=omega,est_r=est_r)
options(error = traceback)
res_list <- list()
for (s in 1:nrow(param_combos)) {
  print(s)

  # Set up list of params for this iteration
  params <- list(n=param_combos[s,'n'], omega=param_combos[s,'omega'],
                 rho=param_combos[s,'rho'],
                 rho_s=param_combos[s,'rho']*0.5,
                 p=p,
                 delta=delta,
                 nu=nu,
                 beta=beta,
                 gamma=gamma,
                 eta=eta,
                 tau=tau,
                 est_r=param_combos[s,'est_r'],
                 sl.lib=sl.lib)

  # For current combination of parameters, run simulation
  res <- main_sim(nsim,
                  params,parallel=T)

  # Store results
  res_list[[s]] <- res

}
final_res <- do.call(rbind, res_list)
#------------------------------------------------------------------------------
# Save output to csv
write.csv(final_res,fullpath,row.names = FALSE)
