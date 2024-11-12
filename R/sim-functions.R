# sim-functions.R
#
# Code for implementing simulations in the error-dependent sampling project
#-------------------------------------------------------------------------------

sim_iter <- function(n,p,delta,nu,beta,gamma,eta,omega,tau,rho,sl.lib,
                     control_rate = TRUE) {
  #' Function for simulating data
  #' 
  #' Arguments:
  #' - n: Sample size
  #' - p: Number of covariates
  #' - delta: Treatment effect
  #' - nu: Outcome measurement error variance
  #' - beta, tau, gamma: Outcome model coefficients
  #' 
  #' Outputs:
  #' - A dataframe containing the sim data
  #' 
  #' 

  # Simulate data
 df <- gen_data(n,p,delta,nu,beta,gamma,eta,omega,tau,rho)
 
  # Estimate psi with the ad-hoc rule
  psi_hat <- est_psi_a(df,sl.lib) 
  
  # Estimate psi with the adaptive rule
  df_adapt <- two_stage_main(df,
                             Y,A,Astar,Ystar,X,
                             rho_s,rho,
                             sl.lib) 
  psi_hat_adapt <- est_psi_a(df_adapt, sl.lib) 
  
  # Turn into dfs
  psi_hat <- data.frame(t(unlist(psi_hat)),stringsAsFactors = FALSE)
  
  psi_hat_adapt <- data.frame(t(unlist(psi_hat_adapt)),stringsAsFactors = FALSE)
  names(psi_hat_adapt) <- paste0(names(psi_hat_adapt),"_adapt")
 
 # Unpack output to store in df
 psi_hat <- cbind(psi_hat,psi_hat_adapt)
 
  # Return data
  return(psi_hat)
}

main_sim <- function(nsim, params,
                     parallel = TRUE) {
  #' Main function for carrying out simulations in the error-dependent sampling
  #' project
  #' 
  #' Arguments:
  #' - nsim: Number of simulations to run
  #' - params: A list comntaining all simulation parameters
  #' 
  #' Outputs:
  #' - A dataframe containing the simulation results
  #' 
  
  # Extract names from params
  list2env(params, envir = environment())
  
  if (parallel==TRUE) {
    # Register parallel backend
    registerDoParallel(cores = detectCores())
    
    # Carry out the simulation 
    sim_results <- foreach(i = 1:nsim, .combine = rbind) %dopar% {
      sim_iter(n,p,delta,nu,beta,gamma,eta,omega,tau,rho,sl.lib)
    }
    
    stopImplicitCluster()
  } else {
    # Carry out the simulation 
    sim_results <- foreach(i = 1:nsim, .combine = rbind) %do% {
      sim_iter(n,p,delta,nu,beta,gamma,eta,omega,tau,rho,sl.lib)
    }
  }
  
  # add parameters to results dataframe
  params_df <- as.data.frame(t(unlist(params)))
  
  # replicate params_df so that it has n total identical rows
  params_df <- params_df[rep(seq_len(nrow(params_df)), each = nsim), ]
  
  # combine params_df and sim_results
  sim_results <- cbind(params_df, sim_results)
  
  # Return results
  return(sim_results)
   
}
