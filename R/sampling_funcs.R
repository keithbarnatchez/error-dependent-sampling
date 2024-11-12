#' Function for performing two-phase sampling scheme, following the variance
#' minimizing method presented in Wang et al. (2024). Given a dataframe and
#' sampling probs for the pilot wave + overall, function returns a new dataframe 
#' that adds on sampling indicators for 
#' 
#' INPUTS:
#' - data: dataframe created by gen_data()
#' - Y: name of outcome variable
#' - A: name of treatment variable
#' - Astar: name of error-prone tmt variable 
#' - Ystar: name of error-prone outcome variable 
#' - X: vector of covariate names
#' - rho_s: sampling probability for pilot wave
#' - rho: overall sampling probability
#' - sl.lib: superlearner library to use for estimating the conditional variance
#' 
#'
#'
#'
#' OUTPUTS:
#' - new_data: dataframe with new sampling indicators. S=wave 1, R=wave 2
#'
two_stage_main <- function(data,
                           Y,A,Astar,Ystar,X,
                           rho_s,rho,
                           sl.lib) {
  
  # Stage 1: Get initial pilot sample, of relative size rho_s
  S <- rbinom(nrow(data),1,rho_s)
  pilot_data <- data[S == 1,] %>% mutate(R=0,S=1) 
  
  # Stage 2: Get full data influence function estimate with the pilot data
  full_data_eif_mod <- AIPW::AIPW$new(Y = pilot_data[,Y],
                              A = pilot_data[,A],
                              W = subset(pilot_data,select=X),
                              Q.SL.library = sl.lib,
                              g.SL.library = sl.lib,
                              k_split = 1,
                              verbose = FALSE)$fit()$summary()
  full_data_eif <- full_data_eif_mod$obs_est$aipw_eif1 - full_data_eif_mod$obs_est$aipw_eif0
  
  # Get conditional variance estimate of full_data_eif, conditional on X, Astar and Ystar
  # First, estimate E[fulldataeif|X,Astar,Ystar] with superlearner
  rhs <- pilot_data %>% select(all_of(c(X,Astar,Ystar)))
  full_data_eif_reg <- SuperLearner::SuperLearner(Y = full_data_eif,
                                                 X = rhs,
                                                 SL.library = sl.lib)
 
  # Get the residuals
  full_data_eif_res <- full_data_eif - full_data_eif_reg$SL.predict
  
  # Estimate the conditional variance by regressing squared resids
  full_data_eif_var <- SuperLearner::SuperLearner(Y = full_data_eif_res^2,
                                                  X = rhs,
                                                  SL.library = 'SL.ranger')
  
  # Next, predict conditional variance over the full data
  sigma2_pred <- predict(full_data_eif_var,newdata=data)$pred
  
  # Define sampling func for generic t
  est_eq <- function(t) {
           mean(
             (1-S)*(as.numeric(sigma2_pred > t) + (sigma2_pred/t)*as.numeric(sigma2_pred <= t)) - (rho - rho_s)
           )
  }
  
  # Solve for 0 and plug in to get the sampling func
  thresh <- uniroot(est_eq, c(1e-6,1e5*max(sigma2_pred)))$root
  sampling_func <- (1-S)*((sigma2_pred > thresh) + (sigma2_pred/thresh)*(sigma2_pred <= thresh))
  
  # Update data with sampling indicators
  data$S <- S
  data$V <- ifelse(S==1,0,rbinom(nrow(data),1,sampling_func))
  data$R <- data$S + data$V
  
  # Also add on the final sampling probs
  data$final_probs <- rho_s + (1-rho_s)*sampling_func
  
  return(data)
  
}


two_stage_adapt <- function(data,
                            Y, A, Astar, Ystar, X, 
                            Sprobs,
                            rho_s, rho,
                            sl.lib) {
  
  ##### Stage 1: Extract the initial sampling probabilities, and make them mean rho_s
  samp_probs <- data[,Rprobs]
  est_eq_s <- function(t) {
    mean( (samp_probs > t) + samp_probs/t*(samp_probs <= t) - rho_s ) 
  }
  tstar <- uniroot(est_eq_s, c(1e-3,100*max(samp_probs)))$root
  updated_samp_probs <- est_eq_s(tstar) + rho_s
  
  # Use these to get the first stage sample
  S <- rbinom(nrow(data),1,updated_samp_probs) ; data$S <- S
  
  ##### Stage 2: Get full-data EIF for this subset and get cond var est
  full_data_eif <- with(data,
                        (eta_hat1/pi_hat1 - eta_hat0/pi_hat0) + 
                          A/pi_hat1*(Y-eta_hat1/pi_hat1) - 
                          (1-A)/(1-pi_hat0)*(Y-eta_hat0/pi_hat0)
  )
  data$full_data_eif <- full_data_eif
  
  # Among those with R=1, regress the full data EIF on X, Astar and Ystar
  df_reg <- data %>% filter(R==1)
  full_data_eif_mod <- SuperLearner::SuperLearner(Y=df_reg$full_data_eif,
                                                  X=df_reg %>% select(X,Astar,Ystar),
                                                  SL.library='SL.ranger',
                                                  verbose=FALSE)
  
  # Predict the full data EIF for all observations
  full_data_eif_proj <- predict(full_data_eif_mod,newdata=data)$pred
  
  # Get the conditional variance estimate
  full_data_eif_res <- full_data_eif - full_data_eif_proj
  full_data_eif_var <- SuperLearner::SuperLearner(Y = full_data_eif_res^2,
                                                  X = df_reg %>% select(X,Astar,Ystar),
                                                  SL.library = 'SL.ranger')
  sigma2_pred <- predict(full_data_eif_var,newdata=data)$pred
  
  ##### Step 3: With the conditional variance, get the new sampling func
  est_eq <- function(t) {
    mean(
      (1-S)*(as.numeric(sigma2_pred > t) + (sigma2_pred/t)*as.numeric(sigma2_pred <= t)) - (rho - rho_s)
    )
  }
  thresh <- uniroot(est_eq, c(1e-6,100*max(sigma2_pred)))$root
  sampling_func <- (1-S)*((sigma2_pred > thresh) + (sigma2_pred/thresh)*(sigma2_pred <= thresh))
  
  # Update data with sampling indicators
  data$V <- ifelse(S==1,0,rbinom(nrow(data),1,sampling_func))
  data$R <- max(data$S,data$V)
  
  # Also add on the final sampling probs
  data$final_probs <- rho_s + (1-rho_s)*sampling_func
  
  return(data)
  
}