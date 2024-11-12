# sim-main.R
#
# Main file for conducting simulations for the me-dependent-sampling project
#-------------------------------------------------------------------------------
rm(list=ls())
source('data-functions.R')
source('estimation-functions.R')
#-------------------------------------------------------------------------------
# Check if on cluster to determine library path

if (Sys.info()['sysname'] == 'Linux') { # if on cluster, specify library path
  .libPaths("~/R_packages")
}
#-------------------------------------------------------------------------------
# Ensure all required packages are installed

packages <- c("tidyverse","mvtnorm",'SuperLearner','ggridges')
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
#-------------------------------------------------------------------------------
# Setup parameters

n <- 2500 ; p <- 3 # sample size, number of covariates
delta <- c(-0.1,-0.6,-0.9)
nu <- c(0.1,-0.1,0.1)
beta <- c(1,2,-2) # outcome model main effects
tau <- 1 # constant treatment effect

gamma <- rep(1,3) # interaction effects with tmt c(0.5,2.1,-1.2)
eta <- c(0.6,-0.2,0.8,0.1,-0.3)
omega <- 0.1 
rho <- 0.1 # relative size of the validation data
#-------------------------------------------------------------------------------
# Setup grids for parameters we want to vary

n_grid <- c(1000, 2500, 5000, 1e4)
rho_grid <- c(0.1, 0.2, 0.3, 0.4, 0.5)



#-------------------------------------------------------------------------------
# Main simulation 

nsim <- 250 # number of iterations

psi_hat_df <- matrix(ncol = 8, nrow = nsim)
for (s in 1:nsim) {
  
  # print every 100 iterations
  if (s %% 100 == 0) {
    print(paste('On iteration',s,'of',nsim,sep=' '))
  }
  
  # Simulate data
  df <- gen_data(n,p,delta,nu,beta,gamma,eta,omega,tau)
  
  # Get tmt effect estimates / SEs
  psi_hat <- est_psi_a(df)
  
  psi_hat_df[s,] <- unlist(psi_hat)
}

colnames(psi_hat_df) <- colnames(data.frame(t(unlist(psi_hat))))
psi_hat_df <- as.data.frame(psi_hat_df)
#------------------------------------------------------------------------------
# Save output to csv

write.csv(psi_hat_df,'psi_hat_df.csv',row.names = FALSE)

# Move off of cluster onto desktop

#--------------------------------------------------------------------------------
# Make plots

# ridge plot of the one-step and plug-in

# keep only the plug-in and one-step estimates and reshape to long format
psi_hat_df_long <- psi_hat_df %>% select(onestep,plugin) %>%
  pivot_longer(cols = everything(), 
               names_to = "estimator", values_to = "psi_hat") %>%
  mutate(estimator = case_when(estimator == 'onestep' ~ 'One step',
                               estimator == 'plugin' ~ 'Plug in'))

# Form the plot (distributions of estimates across sim iterations, by estimator)
psi_hat_df_long  %>%
  ggplot(aes(x = psi_hat,
             y = as.factor(estimator),
             fill=as.factor(estimator))) +
  geom_density_ridges(alpha=0.8) + theme_ridges() +
  labs(title = "Distribution of treatment effect estimates",
       subtitle = "n = 2500") +
  labs(x='Estimate',
       y='Method') +
  theme(legend.position = 'none') +
  geom_vline(xintercept=2.5,color='black',linetype='dashed') 

ggsave
