# test-funcs.R
#
# Testing code for data generating/estimation/simulation functions
#-------------------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(SuperLearner)
source('estimation-functions.R')
source('data-functions.R')
#-------------------------------------------------------------------------------
# Turn off warning messages

options(warn=-1)
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
# Setup parameters

n <- 500 ; p <- 3
delta <- c(-0.1,-0.6,-0.9)
nu <- c(0.1,-0.1,0.1)
beta <- c(1,2,-2)
tau <- 1
gamma <- rep(1,3) # c(0.5,2.1,-1.2)
eta <- c(0.6,-0.2,0.8,0.1,-0.3)
rho <- 0.1 # relative size of the validation data
#-------------------------------------------------------------------------------
# Test data gen funcs

X <- gen_covariates(n,p)
A <- gen_treatment(X,delta)
Astar <- treatment_meas_model(A, X, omega)
Y <- outcome_model(A,X,beta,tau,gamma)
Ystar <- outcome_meas_model(A,X,nu)
R <- gen_selection(X,Astar,Ystar,eta)

df <- gen_data(n,p,delta,nu,beta,gamma,eta,omega,tau)
#-------------------------------------------------------------------------------
# Test nuisance estimation funcs

mu_hat <- est_mu_a(df)
lambda_hat <- est_lambda_a(df)
kappa_hat <- est_kappa(df)
  
df_aug <- cbind(df,mu_hat,lambda_hat,kappa_hat)
eta_hat <- est_eta_a(df_aug)
pi_hat <- est_pi_a(df_aug)

df_aug_aug <- cbind(df_aug,eta_hat,pi_hat)
#-------------------------------------------------------------------------------
# Test treatment effect estimation funcs
df <- gen_data(n,p,delta,nu,beta,gamma,eta,omega,tau)
output <- est_psi_a(df) ; output
#-------------------------------------------------------------------------------
# Simulate 25 datasets and estimate treatment effect for each

nsim <- 100
# create an empty dataframe with 4 columns and nsim rows
psi_hat_df <- matrix(ncol = 8, nrow = nsim)
for (s in 1:nsim) {
  print(s)
  # Simulate data
  df <- gen_data(n,p,delta,nu,beta,gamma,eta,omega,tau)
  
  # Get tmt effect estimate
  psi_hat <- est_psi_a(df)
 
  # convert psi_hat into a dataframe
  psi_hat_df[s,] <- unlist(psi_hat)
}


colnames(psi_hat_df) <- colnames(data.frame(t(unlist(psi_hat))))
psi_hat_df <- as.data.frame(psi_hat_df)

# keep only the plugin and onstep estimates and reshape to long format
psi_hat_df_long <- psi_hat_df %>% select(onestep,plugin) %>%
  pivot_longer(cols = everything(), names_to = "estimator", values_to = "psi_hat")

# make a new version of dataframe that filters observations if they exceed 10 in abs
# value
psi_hat_df_long_windsor <- psi_hat_df_long %>% filter(abs(psi_hat) < 5) %>%
  mutate(estimator = case_when(estimator == 'onestep' ~ 'One step',
                               estimator == 'plugin' ~ 'Plug in'))

# make a ridgeplot of the estimates

psi_hat_df_long_windsor  %>% filter(estimator=='One step') %>%
  ggplot(aes(x = psi_hat,
             y = as.factor(estimator),
             fill=as.factor(estimator))) +
  geom_density_ridges(alpha=0.8) + theme_ridges() +
  labs(title = "Distribution of treatment effect estimates",
       subtitle = "n = 500") +
  labs(x='Estimate',
       y='Method') +
  theme(legend.position = 'none') +
  geom_vline(xintercept=2.5,color='black',linetype='dashed') 

