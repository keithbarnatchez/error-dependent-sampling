# estimation-functions.R
#
# Functions for performing estimation (plug-in and one-step) in the measurement
# error-dependent sampling project
#-------------------------------------------------------------------------------

est_mu_a <- function(data, sl.lib) {
  #' Function for estimating E(Y|A=a,X,Astar,Ystar,R=1), i.e. the imputation
  #' model for Y  
  #'
  #' Inputs:
  #' - Y: Outcome
  #' - A: Treatment assignment
  #' - X: Covariate matrix
  #' - Astar: Treatment assignment in validation study
  #' - Ystar: Outcome in validation study
  #' 
  #' Outputs:
  #' - A vector of estimated conditional means

    # Make a df with all the RHS variables for superlearner
  rhsvars <- c(names(data)[grep("^X",names(data))],'A','Astar','Ystar')
  
  # Unpack data
  est_data <- data %>% filter(R==1)
  Y <- est_data$Y
  covs <- est_data %>% select(all_of(rhsvars))
  
  # Estimate mu_a via superlearner
  mu_a_mod <- SuperLearner::SuperLearner(Y=Y, X=covs, 
              SL.library = sl.lib, 
              family = gaussian(), method = "method.NNLS", verbose = FALSE)
  
  # Get predicted values under treatment
  XA1 <- data %>% mutate(A=1) 
  mu_hat1 <- predict(mu_a_mod,
                     XA1,
                     onlySL=TRUE)$pred
  
  # Get predicted values under no treatment
  XA0 <- data %>% mutate(A=0)
  mu_hat0 <- predict(mu_a_mod,
                     XA0,
                     onlySL=TRUE)$pred
  
  return(data.frame(mu_hat1,mu_hat0))
  
}

est_lambda_a <- function(data, sl.lib) {
  #' Function for estimating P(A=a|X,Astar,Ystar,R=1), i.e. the imputation model
  #' for A
  #'
  #' Inputs:
  #' - A: Treatment assignment
  #' - X: Covariate matrix
  #' - Astar: Treatment assignment in validation study
  #' - Ystar: Outcome in validation study
  #'  
  #' Outputs:
  #' - A vector of estimated conditional probabilities
  #' 
  #'   
  
  # Make a df with all the RHS variables for superlearner
  rhsvars <- c(names(data)[grep("^X",names(data))],'Astar','Ystar')
  
  # Unpack data
  est_data <- data %>% filter(R==1)
  A <- est_data$A
  covs <- est_data %>% select(all_of(rhsvars))
  
  # Estimate lambda_a via superlearner
  lambda_a_mod <- SuperLearner::SuperLearner(Y=A, X=covs, 
              SL.library = c("SL.mean","SL.glmnet","SL.glm","SL.step"), 
              family = binomial(), method = "method.NNLS", verbose = FALSE)
  
  # Get predicted values for the whole dataset
  lambda_hat <- predict(lambda_a_mod,
                        data,
                        onlySL=TRUE)$pred
  
  lamdf <- data.frame(lambda_hat1 = lambda_hat,
                      lambda_hat0 = 1-lambda_hat)
  
  return(lamdf)
  
}

est_kappa <- function(data, sl.lib) {
  #' Function for estimating probability of selection into validation study, ie
  #' P(R=1|X,Astar,Ystar)
  #'
  #' Inputs:
  #' - R: Selection indicator
  #' - X: Covariate matrix
  #' - Astar: Error prone treatment assignment
  #' - Ystar: Error prone outcome
  #' 
  #' Outputs:
  #' - A vector of estimated conditional probabilities
  #' 
  
  # Make a df with all the RHS variables for superlearner
  rhsvars <- c(names(data)[grep("^X",names(data))],'Astar','Ystar')
  covs <- data %>% select(all_of(rhsvars))
  
  # Estimate kappa via superlearner
  kappa_mod <- SuperLearner::SuperLearner(Y=data$R, X=covs, 
              SL.library = sl.lib, 
              family = binomial(), verbose = FALSE)
  
  # Extract the predictions from superlearner model
  kappa <- predict(kappa_mod,
                   data,
                   onlySL=TRUE)$pred
  
  return(kappa)
}

est_eta_a <- function(data,sl.lib) {
  #' Function for predicting
  #'        E(E(Y|A=a,X,Ystar,Astar,R=1)*E(A|X,Ystar,Astar,R=1) | X)
  #' Recall this function effectively marginalizes out Ystar
  #' 
  #' Inputs
  #' - data: A dataframe from gen_data() that has been augmented with the 
  #'.        predicted mu_hat and lambda_hat values
  #' 
  #' Outputs:
  #' - A vector of predictions
  
  # Get beta_a, the product of mu_a and lambda_a
  beta_hat1 <- data$mu_hat1*data$lambda_hat1
  beta_hat0 <- data$mu_hat0*data$lambda_hat0
  
  rhsvars <- c(names(data)[grep("^X",names(data))])
  covs <- data %>% select(all_of(rhsvars))
  
  # Regress beta_hat on X with superlearner
  # Need separate models for A=0 and A=1
  eta_a_mod1 <- SuperLearner::SuperLearner(Y=beta_hat1, X=covs, 
              SL.library = sl.lib, 
              family = gaussian(), verbose = FALSE)
  eta_a_mod0 <- SuperLearner::SuperLearner(Y=beta_hat0, X=covs, 
                SL.library = sl.lib, 
                family = gaussian(), verbose = FALSE)
  
  # Get the predicted values for the whole dataset
  eta_hat1 <- predict(eta_a_mod1,
                      data,
                      onlySL=TRUE)$pred
  eta_hat0 <- predict(eta_a_mod0,
                      data,
                      onlySL=TRUE)$pred
  
  return(data.frame(eta_hat1,eta_hat0))
  
}

est_pi_a <- function(data,sl.lib) {
  #' Function for estimating
  #'       E( P(A=a | X, Astar, Ystar, R=1) | X, Astar)
  #' 
  #' Inputs
  #' - data: A dataframe from gen_data() that has been augmented with the
  #'        predicted lambda_hat values
  #'        
  #' Outputs:
  #' - A vector of predictions
  
  # Regress lambda_hat on X and Astar with superlearner
  rhsvars <- c(names(data)[grep("^X",names(data))])
  covs <- data %>% select(all_of(rhsvars))
  
  # Need separate models for A=0 and A=1
  # pi_1 = E(lambda_1 | X, Astar)
  # Note: setting family to binomial since pi_a should be in (0,1)
  pi_a_mod1 <- SuperLearner::SuperLearner(Y=data$lambda_hat1, X=covs, 
              SL.library = sl.lib, 
              family = binomial(), verbose = FALSE)
  
  # pi_1 = E(lambda_0 | X, Astar) = E(1-lambda_1 | X, Astar)
  pi_a_mod0 <- SuperLearner::SuperLearner(Y=data$lambda_hat0, X=covs, 
                SL.library = sl.lib, 
                family = binomial(), verbose = FALSE)
  
  # Get the predicted values for the whole dataset
  pi_hat1 <- predict(pi_a_mod1,
                      data,
                      onlySL=TRUE)$pred
  pi_hat0 <- predict(pi_a_mod0,
                      data,
                      onlySL=TRUE)$pred
  
  return(data.frame(pi_hat1,pi_hat0))
}

est_psi_a <- function(data,sl.lib) {
  #' Function for estimating 
  #'    E(Y(a)) = E(eta_a(X)/pi_a(X))
  #' via efficient one-step estimation (and plug-in estimation)
  #' 
  #' Inputs:
  #' - data: Dataframe from the gen_data function
  #' - cross_fit: Logical indicating whether to use cross-fitting
  #' 
  #' Outputs:
  #' - An output object with an effect estimate, 95% confidence interval and SE
  #'   for both the one-step and plug-in estimators
  
  # Estimate mu_a and lambda_a with the validation data
  mu_a <- est_mu_a(data,sl.lib)
  lambda_a <- est_lambda_a(data,sl.lib)
  
  # Append new predictions onto the dataframe
  data <- cbind(data,mu_a,lambda_a)
  
  # Regress mu_a and lambda_a predictions on X and Astar across the whole dataset
  # to get eta_a predictions
  eta_a <- est_eta_a(data,sl.lib)
  data <- cbind(data,eta_a)
  
  # Regress lambda_a predictions on X and Astar across the whole dataset
  pi_a <- est_pi_a(data, sl.lib)
  data <- cbind(data,pi_a)
  
  # Get P(R=1|X,Ystar,Astar)
  if('final_probs' %in% names(data)) {
    kappa_hat <- data$final_probs
  } else {
    kappa_hat <- data$Rprobs # est_kappa(data,sl.lib)
  }
  
  # Get the plug-in estimate and its SE
  plugin <- mean(data$eta_hat1/data$pi_hat1) - mean(data$eta_hat0/data$pi_hat0)
  pluginse <- sd(data$eta_hat1/data$pi_hat1 - data$eta_hat0/data$pi_hat0)/sqrt(nrow(data))
  
  # Get the one-step bias corrected estimate and its SE
  # First for A=1
  psi1hatold <- with(data,
          (eta_hat1/pi_hat1 + 
              A*R*(Y-mu_hat1)/(kappa_hat*pi_hat1) +
              (R*(A-lambda_hat1)/kappa_hat + lambda_hat1-pi_hat1)*
              ((mu_hat1-eta_hat1/pi_hat1)/pi_hat1)
          )     
  )
  psi1hat <- with(data,
                  eta_hat1/pi_hat1 +
                  (lambda_hat1/pi_hat1)*(mu_hat1 - eta_hat1/pi_hat1) +
                  (A*R/(kappa_hat*pi_hat1))*(Y - eta_hat1/pi_hat1) -
                  (R*lambda_hat1)/(kappa_hat*pi_hat1)*(mu_hat1 - eta_hat1/pi_hat1)   
                    )
                    
  # Now with A=0
  psi0hatold <- with(data,
               (eta_hat0/pi_hat0 + 
                      (1-A)*R*(Y-mu_hat0)/(kappa_hat*pi_hat0) +
                      (R*((1-A)-(1-lambda_hat0))/kappa_hat + (1-lambda_hat0)-pi_hat0)*
                      ((mu_hat0-eta_hat0/pi_hat0)/pi_hat0)
               )     
  )
  psi0hat <- with(data,
                  eta_hat0/pi_hat0 +
                    (lambda_hat0/pi_hat0)*(mu_hat0 - eta_hat0/pi_hat0) +
                    ((1-A)*R/(kappa_hat*pi_hat0))*(Y - eta_hat0/pi_hat0) -
                    (R*lambda_hat0)/(kappa_hat*pi_hat0)*(mu_hat0 - eta_hat0/pi_hat0)   
  )
  
  # Form the final one-step ATE estimate
  onestep <- mean(psi1hat - psi0hat)
  onestepse <- sd(psi1hat - psi0hat)/sqrt(nrow(data))
  
  # Also get the alternate one-step estimate that projects the full-data EIF 
  # onto the always-observed data
  onestepalt_res <- get_alt_onestep(data,kappa_hat,sl.lib)
  onestepalt <- onestepalt_res$psihat ; onestepaltse <- onestepalt_res$psihatse
  
  # browser()
  # And get re-weighted TMLE
  tmle_res <- tmle_est(data,sl.lib)
  tml <- tmle_res$tml ; tmlse <- tmle_res$tmlse
  
  # And get the re-weighted one-step
  
  
  
  return(list(plugin=plugin,pluginse=pluginse,
              onestep=onestep,onestepse=onestepse,
              tml=tml,tmlse=tmlse))
              # maxl1=max(data$lambda_hat1),maxpi1=max(data$pi_hat1),
              # mink=min(kappa_hat),maxk=max(kappa_hat)))
}

#' Function for implementing the two-phase sampling strategy outlined in 
#' Wang et al. (2023). Main idea: obtain a completely random val. subset
#' of relative size rho_S, use this to obtain a full-data EIF estimate,
#' and then use this full-data EIF to obtain the final sampling probs
#'
two_phase_est <- function(data,
                          rho_S,rho,sl.lib) {
  
  # Part 1: get sampling indicators
  data <- two_stage_main(data,Y='Y',A='A',Astar='Astar',
                         Ystar='Ystar',X='X',rho_S=rho_S,
                         rho=rho,sl.lib=sl.lib)
  
  pilot_data <- data %>% filter(S==1)
  non_pilot_data <- data %>% filter(S==0)
  
  # Part 2: get full data EIF estimate. Estimate E[Y|A,X] with the pilot data
  mu_hat_mod <- SuperLearner::SuperLearner(Y=pilot_data[,Y],
                                           A=subset(pilot_data, select=c(X,A)),
                                           SL.library = sl.lib)
  pi_hat_mod <- SuperLearner::SuperLearner(Y=pilot_data[,A],
                                           A=subset(pilot_data, select=X),
                                           SL.library = sl.lib)
  mu1_hat_p <- predict(mu_hat_mod,
                           newdata=pilot_data %>% mutate(A=1))$pred
  mu0_hat_p <- predict(mu_hat_mod,
                           newdata=pilot_data %>% mutate(A=0))$pred
  pi_hat_p <- predict(pi_hat_mod)$pred
  full_data_eif <- with(pilot_data,
                      (mu1_hat_p - mu0_hat_p) + Y*A/pi_hat_p - Y*(1-A)/(1-pi_hat_p))
  
  # Part 3: Regress the full_data_eif on X, Astar and Ystar with superlearner
  full_data_eif_mod <- SuperLearner::SuperLearner(Y=full_data_eif,
                                        X=pilot_data %>% select(X,Astar,Ystar),
                                        SL.library=sl.lib,
                                        family=binomial(),
                                        verbose=FALSE)
  
  # Get the final ingredients for estimators
  pilot_data$full_data_eif <- full_data_eif
  non_pilot_data$full_data_eif <- 0 
  pilot_data$full_data_eif_proj <- predict(full_data_eif_mod,newdata=pilot_data)
  non_pilot_data$full_data_eif_proj <- predict(full_data_eif_mod,newdata=non_pilot_data)
  
  # Form the final estimate
  final_est <- with(pilot_data,
                    mean())

}

get_alt_onestep <- function(data,
                            Y,A,W,X,
                            kappa_hat,
                            sl.lib) {
  #'
  #' 
  #' INPUTS:
  #'  - A df which has been augmented with eta and pi estimates
  #'  - Sampling probabilities
  
  # if ('Rprobs' %in% colnames(data)) {
  #   aipw_res <- AIPW::AIPW$new(Y=data$Y,
  #                A=data$A,
  #                W=subset(data,select=X),
  #                Q.SL.library = sl.lib,
  #                g.SL.library = sl.lib,
  #                k_split = 1,
  #                verbose=FALSE)$fit()$summary()
  # }
  
  # Get the full data EIF estimate
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
                                        SL.library=sl.lib,
                                        verbose=FALSE)
  
  # Predict the full data EIF for all observations
  full_data_eif_proj <- predict(full_data_eif_mod,newdata=data)$pred
  
  # Get the observed data eif
  R <- data$R
  obs_data_eif <- full_data_eif*R/kappa_hat - 
                       (R/kappa_hat - 1)*full_data_eif_proj
  
  # Get the one-step estimate
  psihat <- mean(obs_data_eif)
  psihatse <- sd(obs_data_eif)/sqrt(nrow(data))
  
  return(list(psihat=psihat,psihatse=psihatse))
  
}

shift_n_scale <- function(x) {
  #' Shifts and scales a vector to be between 0 and 1
  #' 
  #' INPUTS:
  #' - x: a numeric vector
  #' 
  #' OUTPUTS:
  #' - A numeric vector that has been shifted and scaled
  #'
  
  if (min(x)>=0 & max(x)<=1) {
    return(x)
  }
  
  return((x - min(x))/(max(x) - min(x)))
}

logit_transform <- function(x) {
  
  # First, check if x is bounded between 0 and 1. If not, transform it
  if (min(x) < 0 | max(x) > 1) {
    x <- shift_n_scale(x)
  }
  
  return(log(x/(1-x)))
  
}

tmle_evm <- function(data,Z,sl.lib) {
  #' Implements a TML estimator using the ideas from Hejazi (2020) and Rose & 
  #' van der Laan (2011), with empirical variance minimization
  #' 
  #' INPUTS:
  #' - data:
  #' - sl.lib: string vector of superlearner libraries
  #' 
  #' Outputs:
  #' - Point estimate and estimated SE
  #' 
  
  data=gen_data(n,p,delta,nu,beta,gamma,eta,omega,tau,rho)
  
  Z <- c(X,Astar,Ystar)

  # First, get the initial reweighted TMLE 
  rvdl_tmle <- with(data,
                    tmle::tmle(Y=data$Y,A=data$A,W=data %>% select(X),
                               obsWeights = data$R/data$Rprobs,
                               g.SL.library=sl.lib,
                               Q.SL.library=sl.lib))
  
  # Extract the IC
  IC <- rvdl_tmle$estimates$IC$IC.ATE + rvdl_tmle$estimates$ATE$psi
  
  # Re-scale this IC by (R/Rprobs - 1)^{-1}
  IC_rescaled <- with(data,(R/Rprobs - 1)^(-1)*IC)

  
  # Now regress IC on the rescaled Z with splines, using (weights R/Rprobs-1)^2
  varphi_mod <- SuperLearner::SuperLearner(Y=IC_rescaled,
                                         X=data %>% select(Z),
                                         obsWeights = (data$R/data$Rprobs - 1)^2,
                                         SL.library='SL.earth',
                                         verbose=FALSE)
  
  # Get the predicted IC
  varphi_Z <- predict(varphi_mod,newdata=data)$pred
  
  # Get final estimate
  tml_evm <- with(data, mean(IC - (R/Rprobs-1)*varphi_Z))
  tml_evmse <- with(data, sd(IC - (R/Rprobs-1)*varphi_Z)/sqrt(nrow(data)))
  
  return(list(tml_evm=tml_evm,tml_evmse=tml_evmse))
  
  
}

tmle_cv <- function(data, sl.lib) {
  
  # First, get initial unbiased estimate
  rvdl_tmle <- with(data,
                    tmle::tmle(Y=data$Y,A=data$A,W=data %>% select(X),
                               obsWeights = data$R/data$Rprobs,
                               g.SL.library=sl.lib,
                               Q.SL.library=sl.lib))
  tml_val <- rvdl_tmle$estimates$ATE$psi
  tmlse_val <- sqrt(rvdl_tmle$estimates$ATE$var.psi)
  tml_eic <- rvdl_tmle$estimates$IC$IC.ATE + tml_val
  
  
  
  # Next, get error-prone versions
  rvdl_tmle_epval <-  with(data,
                           tmle::tmle(Y=data$Ystar,A=data$Astar,W=data %>% select(X),
                                      obsWeights = data$R/data$Rprobs,
                                      g.SL.library=sl.lib,
                                      Q.SL.library=sl.lib))
  
  tml_epval <- rvdl_tmle_epval$estimates$ATE$psi
  tmlse_epval <- sqrt(rvdl_tmle_epval$estimates$ATE$var.psi)
  tml_eicval <- rvdl_tmle_epval$estimates$IC$IC.ATE + tml_epval
  
  # Finally, with full data
  rvdl_tmle_epfull <-  with(data,
                           tmle::tmle(Y=data$Ystar,A=data$Astar,W=data %>% select(X),
                                      g.SL.library=sl.lib,
                                      Q.SL.library=sl.lib))
  tml_epfull <- rvdl_tmle_epfull$estimates$ATE$psi
  tmlse_epfull <- sqrt(rvdl_tmle_epfull$estimates$ATE$var.psi)
  tml_eicfull <- rvdl_tmle_epfull$estimates$IC$IC.ATE + tml_epfull
  
  # Now, form the CV estimator
  gamma_hat <- cov(tml_eic, (R/Rprobs -1)*tml_eicfull)
  V_hat <- var( (R/Rprobs -1)*tml_eicfull )
  
  tau_cv <- tml_val - gamma_hat/V_hat * (tml_epval - tml_epfull)
  
  return(list(tau_cv=tau_cv))
  
}

tmle_est <- function(data, sl.lib) {
  #' Implements a TML estimator using the ideas from Hejazi (2020) and Rose & 
  #' van der Laan (2011). 
  #'
  #' INPUTS:
  #' - data:
  #' - sl.lib: string vector of superlearner libraries
  #' 
  #' Outputs:
  #' - Point estimate and estimated SE
  #'

  # Y <- data$Y ; A <- data$A ; X <- data$X
  # R <- data$R ; Ystar <- data$Ystar ; Astar <- data$Astar
  
  # Do Rose & van der Laan (2011) method if we have known probs
  if ('Rprobs' %in% colnames(data)) { 
    
    # Idea: just do a weighted TMLE 
    rvdl_tmle <- with(data,
                      tmle::tmle(Y=data$Y,A=data$A,W=data %>% select(X),
                            obsWeights = data$R/data$Rprobs,
                            g.SL.library=sl.lib,
                            Q.SL.library=sl.lib))
    
    # Extract results
    tml <- rvdl_tmle$estimates$ATE$psi
    tmlse <- sqrt(rvdl_tmle$estimates$ATE$var.psi)
    
    # Return in a list
    return(
      list(tml=tml,tmlse=tmlse) 
    )
    
  }
  
  # Shift and scale Y (and keep track of scaling param)
  data <- data %>% mutate(Ysc = shift_n_scale(Y))
  ylow <- min(data$Y) ; yhigh <- max(data$Y) ; sc <- yhigh-ylow
  
  # First, estimate the sampling probs
  kappa_hat <- est_kappa(data,sl.lib)
  data$kappa_hat <- kappa_hat
  
  # Otherwise, do the Hejazi (2021) method
  data_reg <- data %>% filter(R==1)

  # get estimate of E[Y|A=a,X] and P(A=a|X) via weighted regressions
  outcome_mod <- SuperLearner::SuperLearner(Y=data_reg$Ysc,
                             X=data_reg %>% select(X, A),
                             SL.library=c('SL.earth'),
                             obsWeights=1/data_reg$kappa_hat,
                             verbose=FALSE)
  
  m_hat1 <- predict(outcome_mod,
                    newdata=data %>% mutate(A=1))$pred
  m_hat0 <- predict(outcome_mod,
                    newdata=data %>% mutate(A=0))$pred
  
  # Similarly, get estimates of P(A=a|X) with weighted regressions
  pi_mod <- SuperLearner::SuperLearner(Y=data$A,
                             X=data %>% select(X),
                             SL.library='SL.earth',
                             obsWeights=data$R/kappa_hat,
                             verbose=FALSE)
  pi_hat1 <- predict(pi_mod,new_data=data)$pred ; pi_hat0 <- 1-pi_hat1
  
  # Form the plug-in estimates
  plugin <- m_hat1 - m_hat0
  
  
  # Use these to form full-data-EIF pseudo-outcomes
  psuedo_outcome <- with(data,
                         (m_hat1 - m_hat0) + Ysc*A/pi_hat1 - Ysc*(1-A)/(1-pi_hat0))
  
  
  # Regress the pseudo-outcomes on Z to get projection estimate of EIF
  pseudo_mod <- SuperLearner::SuperLearner(Y=psuedo_outcome,
                             X=data %>% select(X,Astar,Ystar),
                             SL.library=sl.lib,
                             verbose=FALSE)
  eif_proj <- predict(pseudo_mod,newdata=data)$pred
  
  # Get the full data EIF
  
  # The TML is carried out in two updates. First, update kappa estimate
  kappa_new <- glm(R ~  -1 + I(eif_proj/kappa_hat),
        data=data,
        offset=logit_transform(kappa_hat),
        family=binomial())
  kappa_new <- predict(kappa_new,newdata=data,type='response')
  
  # Next, update the outcome regression
  
  # Construct clever covariates
  H_A <- with(data, A/pi_hat1 - (1-A)/(1-pi_hat0))
  H_1 <- with(data, 1/pi_hat1)
  H_0 <- with(data, -1/(1-pi_hat0))
  
  # Get eps estimate
  m_mod_new <- glm(Ysc ~ -1 + H_A,
                   data=data,
                   family=binomial(),
                   weights=data$R/kappa_new,
                   offset=logit_transform(plugin))
  eps <- coef(m_mod_new)
  
  # With eps, we can update m1 and m0
  m1_new <- plogis(qlogis(m_hat1) + eps*H_1)
  m0_new <- plogis(qlogis(m_hat0) + eps*H_0)
  
  mean(m1_new - m0_new)
  
  # Then update the outcome regression estimate
  return(1)
  
  # Finally, do a weighted marginalization of the outcome regression
  
}
