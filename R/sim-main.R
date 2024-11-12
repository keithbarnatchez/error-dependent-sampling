# sim-main.R
#
# Main file for conducting simulations for the me-dependent-sampling project
#-------------------------------------------------------------------------------
rm(list=ls())
source('data-functions.R') 
source('estimation-functions.R')
source('sim-functions.R')
source('sampling_funcs.R') 
#-------------------------------------------------------------------------------
# Check if on cluster to determine library path

par_opt <- FALSE # true if want to parallelize
#-------------------------------------------------------------------------------
# Ensure all required packages are installed

packages <- c("tidyverse","mvtnorm",'SuperLearner','ggridges',
              'doParallel','foreach')
package_list <- rownames(installed.packages())
for (pkg in packages) {  # This should output the actual package names
  if (!(pkg %in% package_list)) {
    stop("Package '", pkg, "' not installed")
  }
}
#-------------------------------------------------------------------------------
# Need a few of the packages loaded
library(tidyverse)
library(SuperLearner)
library(ggridges)
library(doParallel)
library(foreach)
# ------------------------------------------------------------------------------
# Set up relevant paths for output
flnm <- gsub(':','-',Sys.time()) # filename suffixed with date/time
flnm <- paste('sim_results_',gsub(' ','_',flnm),'.csv',sep = '')

simdir <- 'output/' # directory
fullpath <- paste(simdir,flnm,sep='') # construct path to final file name 
# ------------------------------------------------------------------------------
Y <- 'Y' ; A <- 'A' ; Ystar <- 'Ystar' ; Astar <- 'Astar' ; X <- c('X1','X2','X3')
sl.lib <- c('SL.glm.interaction','SL.earth')
# ------------------------------------------------------------------------------
# Setup parameters that remain constant across simulations

n <- 5000 # sample size
p <- 3 # covaraites
delta <- c(-0.1,-0.6,-0.9)
nu <- c(0.1,-0.1,0.1)
beta <- c(1,2,-2) # outcome model main effects
tau <- 1 # constant portion of treatment effect

gamma <- rep(1,3) # interaction effects with tmt c(0.5,2.1,-1.2)
eta <- c(0.6,-0.2,0.8,0.1,-0.3)
omega <- 0.1 
rho <- 0.2 # relative size of the validation data
rho_s <- 0.1
#-------------------------------------------------------------------------------
# Setup parameters that vary across simulations

rho <- 0.2 # c(0.2,0.3) # c(0.1, 0.2, 0.3, 0.4, 0.5) # relative size of the validation data
n <- 5000 # c(500,1000) # c(500, 1000, 2500, 5000) # sample size
omega <- 0.1 # c(0.1,0.11) # c(0.05,0.1,0.15,0.2) # severity of the measurement error
#-------------------------------------------------------------------------------
Y <- 'Y' ; A <- 'A' ; Ystar <- 'Ystar' ; Astar <- 'Astar' ; X <- c('X1','X2','X3')
S <- 'S' ; R <- 'R'
#-------------------------------------------------------------------------------
# Main simulation 

nsim <- 25 # number of simulation iterations

# Vectorize every combination of rho, n, and omega and loop over these combos
param_combos <- expand.grid(rho=rho, n=n, omega=omega)

res_list <- list()
for (s in 1:nrow(param_combos)) {

  # Set up list of params for this iteration
  params <- list(n=param_combos[s,'n'], omega=param_combos[s,'omega'],
                 rho=param_combos[s,'rho'],
                 p=p,
                 delta=delta,
                 nu=nu,
                 beta=beta,
                 gamma=gamma,
                 eta=eta,
                 tau=tau, 
                 sl.lib=sl.lib)
  
  # For current combination of parameters, run simulation
  res <- main_sim(nsim,params,parallel=FALSE)
  
  # Store results 
  res_list[[s]] <- res 

}
final_res <- do.call(rbind, res_list)
#------------------------------------------------------------------------------
# Save output to csv

write.csv(res,fullpath,row.names = FALSE)

# Move off of cluster onto desktop

#--------------------------------------------------------------------------------
# Make plots

# ridge plot of the one-step and plug-in
# 
# # keep only the plug-in and one-step estimates and reshape to long format
# psi_hat_df_long <- psi_hat_df %>% select(onestep,plugin) %>%
#   pivot_longer(cols = everything(), 
#                names_to = "estimator", values_to = "psi_hat") %>%
#   mutate(estimator = case_when(estimator == 'onestep' ~ 'One step',
#                                estimator == 'plugin' ~ 'Plug in'))
# 
# # Form the plot (distributions of estimates across sim iterations, by estimator)
# tmt_eff_plot <- psi_hat_df_long  %>%
#   ggplot(aes(x = psi_hat,
#              y = as.factor(estimator),
#              fill=as.factor(estimator))) +
#   geom_density_ridges(alpha=0.8) + theme_ridges() +
#   labs(title = "Distribution of treatment effect estimates",
#        subtitle = "n = 2500") +
#   labs(x='Estimate',
#        y='Method') +
#   theme(legend.position = 'none') +
#   geom_vline(xintercept=2.5,color='black',linetype='dashed') ; tmt_eff_plot
# 
# # save the plot 
# ggsave('../figures/ridge_plot.pdf',tmt_eff_plot,width=6,height=4)
# 
