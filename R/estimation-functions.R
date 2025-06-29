# estimation-functions.R (CLUSTER VERSION)
#
# Functions for performing estimation (plug-in and one-step) in the measurement
# error-dependent sampling project
#-------------------------------------------------------------------------------

#' @title Get oracle estimate
#'
#' @description Function for estimating the oracle estimate of ATE
#'
#' @param data A dataframe containing data simulated from gen_data()
#' @param sl.lib A list of superlearner algorithms to use
#'
#' @return The oracle estimate of the ATE
#' @export
get_oracle_est <- function(data,sl.lib) {

  # estimate via AIPW
  oracle_est <-  AIPW$new(Y=data$Y,
                          A=data$A,
                          W=subset(data,select=c("X1","X2","X3")),
                          Q.SL.library = sl.lib,
                          g.SL.library = sl.lib,
                          k_split = 1,
                          verbose=FALSE)$fit()$summary()
  
  # Extract estimate and SE
  oracleest <- oracle_est$estimates$RD['Estimate']
  oraclese <- oracle_est$estimates$RD['SE']

  return(list(oracle=unname(oracleest),
              oraclese=unname(oraclese)))
}


#' @title Get naive estimate
#'
#' @description Function for estimating the naive (ignoring measurement error,
#' treating Ystar and Astar as the truth) estimate of ATE
#'
#' @param data A dataframe containing data simulated from gen_data()
#' @param sl.lib A list of superlearner algorithms to use
#'
#' @return The naive estimate of the ATE
#' @export
get_naive_est <- function(data,sl.lib) {

  # estimate via AIPW
  naive_est <-  AIPW$new(Y=data$Ystar,
                         A=data$Astar,
                         W=subset(data,select=c("X1","X2","X3")),
                         Q.SL.library = sl.lib,
                         g.SL.library = sl.lib,
                         k_split = 1,
                         verbose=FALSE)$fit()$summary()
  
  # Extract estimate and SE
  naiveest <- naive_est$estimates$RD['Estimate']
  naivese <- naive_est$estimates$RD['SE']

  return(list(naive=unname(naiveest),
              naivese=unname(naivese)))
}

#' @title Estimate outcome imputation model
#'
#' @description Function for estimating the outcome imputation model
#'
#' @param data A dataframe containing data simulated from gen_data()
#' @param sl.lib A list of superlearner algorithms to use
#'
#' @return A vector of estimated conditional means
#' @export
est_mu_a <- function(data, sl.lib) {

    # Make a df with all the RHS variables for superlearner
  rhsvars <- c(names(data)[grep("^X",names(data))],'A','Astar','Ystar')

  # Unpack data
  est_data <- data %>% filter(R==1)
  Y <- est_data$Y
  covs <- est_data %>% select(all_of(rhsvars))

  # Estimate mu_a via superlearner
  mu_a_mod <- SuperLearner::SuperLearner(Y=Y, X=covs,
              SL.library = sl.lib,
              family = gaussian(), verbose = FALSE)

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

  # Truncate predictions to not exceed max / min of Y
  mu_hat1 <- pmin(pmax(mu_hat1,min(data$Y)),max(data$Y))
  mu_hat0 <- pmin(pmax(mu_hat0,min(data$Y)),max(data$Y))

  return(data.frame(mu_hat1,mu_hat0))

}

#' @title Estimate treatment imputation model
#'
#' @description Function for estimating the treatment imputation model
#'
#' @param data A dataframe containing data simulated from gen_data()
#' @param sl.lib A list of superlearner algorithms to use
#'
#' @return A vector of estimated conditional means
#' @export
est_lambda_a <- function(data, sl.lib) {

  # Make a df with all the RHS variables for superlearner
  rhsvars <- c(names(data)[grep("^X",names(data))],'Astar','Ystar')

  # Unpack data
  est_data <- data %>% filter(R==1)
  A <- est_data$A
  covs <- est_data %>% select(all_of(rhsvars))

  # Estimate lambda_a via superlearner
  # lambda_a_mod <- SuperLearner::SuperLearner(Y=A, X=covs,
  #             SL.library = sl.lib, # c("SL.mean","SL.glmnet","SL.glm","SL.step"),
  #             family = binomial(), method = "method.NNLS", verbose = FALSE)

  lambda_a_mod <- SuperLearner::SuperLearner(
    Y=A, X=covs,
   SL.library = sl.lib, # method = "method.NNloglik",
   family = binomial(), verbose = FALSE)
  # Get predicted values for the whole dataset
  lambda_hat <- predict(lambda_a_mod,
                        data,
                        onlySL=TRUE)$pred
  
  # trim to (.01, .99)
  lambda_hat <- pmin(pmax(lambda_hat,.025),.975)

  lamdf <- data.frame(lambda_hat1 = lambda_hat,
                      lambda_hat0 = 1-lambda_hat)
  
  

  return(lamdf)

}

#' @title Estimate validation probabilities
#'
#' @description Function for estimating validation probabilities
#'
#' @param data A dataframe containing data simulated from gen_data()
#' @param sl.lib A list of superlearner algorithms to use
#'
#' @return A vector of estimated validation probabilities
#' @export
est_kappa <- function(data, sl.lib) {

  # Make a df with all the RHS variables for superlearner
  rhsvars <- c(names(data)[grep("^X",names(data))],'Astar','Ystar')
  covs <- data %>% select(all_of(rhsvars))

  # Estimate kappa via superlearner
  kappa_mod <- SuperLearner::SuperLearner(Y=data$R, X=covs,
              SL.library = sl.lib,
              family = binomial(),
              verbose = FALSE)

  # Extract the predictions from superlearner model
  kappa <- predict(kappa_mod,
                   data,
                   onlySL=TRUE)$pred

  return(kappa)
}

#' @title Estimate E[lambda_a(Z)*mu_a(Z)|X]
#'
#' @description Function for estimating E[mu_a(Z)|X]
#'
#' @param data A dataframe containing data simulated from gen_data()
#' @param sl.lib A list of superlearner algorithms to use
#'
#' @return A vector of estimated conditional means
#' @export
est_eta_a <- function(data,sl.lib) {

  # Get beta_a, the product of mu_a and lambda_a
  beta_hat1 <- data$mu_hat1*data$lambda_hat1
  beta_hat0 <- data$mu_hat0*data$lambda_hat0

  rhsvars <- c(names(data)[grep("^X",names(data))],'Astar')
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


#' @title Estimate E[lambda_a(Z)|X]
#'
#' @description Function for estimating E[mu_a(Z)|X]
#'
#' @param data A dataframe containing data simulated from gen_data()
#' @param sl.lib A list of superlearner algorithms to use
#'
#' @return A vector of estimated conditional means
#' @export
est_pi_a <- function(data,sl.lib) {

  # Regress lambda_hat on X and Astar with superlearner
  rhsvars <- c(names(data)[grep("^X",names(data))],'Astar')
  covs <- data %>% select(all_of(rhsvars))

  # Need separate models for A=0 and A=1
  # pi_1 = E(lambda_1 | X, Astar)
  # Note: setting family to binomial since pi_a should be in (0,1)
  pi_a_mod1 <- SuperLearner::SuperLearner(Y=data$lambda_hat1, X=covs,
              SL.library = sl.lib, verbose = FALSE)

  # pi_1 = E(lambda_0 | X, Astar) = E(1-lambda_1 | X, Astar)
  pi_a_mod0 <- SuperLearner::SuperLearner(Y=data$lambda_hat0, X=covs,
                SL.library = sl.lib,
                verbose = FALSE)

  # Get the predicted values for the whole dataset
  pi_hat1 <- predict(pi_a_mod1,
                      data,
                      onlySL=TRUE)$pred
  pi_hat0 <- predict(pi_a_mod0,
                      data,
                      onlySL=TRUE)$pred

  # Trim pi_hat 1 with pmax pmin
  pi_hat1 <- pmin(0.975, pmax(0.025, pi_hat1))
  pi_hat0 <- 1-pi_hat1

  return(data.frame(pi_hat1,pi_hat0))
}

#' @title Function for implementing Approach 1 one-step estimator
#'
#' @description Function for implementing the one-step estimator for the
#'
#' @param data A dataframe containing data simulated from gen_data()
#' @param sl.lib A list of superlearner algorithms to use
#'
#' @return An output object with an effect estimate, 95% confidence interval and SE
#' for both the one-step and plug-in estimators
#' @export
est_psi_a <- function(data,
                      est_r,sl.lib) {

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
  kappa_hat <- data$Rprobs
  if (est_r) {
    kappa_hat <- est_kappa(data,sl.lib)
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

  # Oracle and naive
  oracle <- get_oracle_est(data,sl.lib)
  naive <- get_naive_est(data,sl.lib)

  return(
    list(plugin=plugin,pluginse=pluginse,
              onestep=onestep,onestepse=onestepse,
              oracle=oracle$oracle,naive=naive$naive,
              oraclese=oracle$oraclese,naivese=naive$naivese,
              eif=psi1hat-psi0hat)
    )
              # maxl1=max(data$lambda_hat1),maxpi1=max(data$pi_hat1),
              # mink=min(kappa_hat),maxk=max(kappa_hat)))
}


# # Do Rose & van der Laan (2011) method if we have known probs
# if (FALSE) { # ('Rprobs' %in% colnames(data)) {
#
#   # Idea: just do a weighted TMLE
#   rvdl_tmle <- with(data,
#                     tmle::tmle(Y=data$Y,A=data$A,W=data %>% select(X),
#                                obsWeights = data$R/data$Rprobs,
#                                g.SL.library=sl.lib,
#                                Q.SL.library=sl.lib))
#
#   # Extract results
#   tml <- rvdl_tmle$estimates$ATE$psi
#   tmlse <- sqrt(rvdl_tmle$estimates$ATE$var.psi)
#
#   # Return in a list
#   return(
#     list(tml=tml,tmlse=tmlse)
#   )
#
# }  # Otherwise, do the Hejazi (2021) method

tmle_est <- function(data,
                     est_r,sl.lib,
                     eem=FALSE) {
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
  #'

  # First, estimate the sampling probs
  kappa_hat <- data$Rprobs
  if (est_r) {
    kappa_hat <- est_kappa(data,sl.lib)
  }
  data$kappa_hat <- as.double(kappa_hat)

  # Scale outcome to [0,1] interval if its continuous
  data$Ysc = shift_n_scale(data$Y)
  data_reg <- data # %>% filter(R==1

  # Get subset of variables: anything starting with X, or A
  Xvars <- colnames(data_reg)[grepl("^X",colnames(data_reg))]

  # get estimate of E[Y|A=a,X] and P(A=a|X) via weighted regressions
  outcome_mod <- SuperLearner::SuperLearner(Y=data_reg$Ysc,
                                            X=data_reg %>% select(c(Xvars,'A')),
                                            SL.library= sl.lib,
                                            method="method.CC_nloglik",
                                            obsWeights=as.double(data$R/data$kappa_hat),
                                            verbose=FALSE)
  m_hat_init <- predict(outcome_mod,newdata=data)$pred
  m_hat1 <- predict(outcome_mod,
                    newdata=data %>% mutate(A=1))$pred
  m_hat0 <- predict(outcome_mod,
                    newdata=data %>% mutate(A=0))$pred

  # Similarly, get estimates of P(A=a|X) with weighted regressions
  pi_mod <- SuperLearner::SuperLearner(Y=data$A,
                                       X=data %>% select(X),
                                       SL.library= sl.lib, #'SL.earth',
                                       obsWeights=as.double(data$R/kappa_hat),
                                       family=binomial(),
                                       verbose=FALSE)
  pi_hat1 <- predict(pi_mod,new_data=data)$pred ; pi_hat0 <- 1-pi_hat1

  # Form the plug-in estimates
  plugin <- m_hat1 - m_hat0

  # Use these to form full-data-EIF pseudo-outcomes
  psuedo_outcome <- with(data,
                         (m_hat1 - m_hat0) + (Ysc-m_hat1)*A/pi_hat1 -
                           (Ysc-m_hat0)*(1-A)/(1-pi_hat1) - mean(plugin))
  # Regress the pseudo-outcomes on Z to get projection estimate of EIF
  pseudo_mod <- SuperLearner::SuperLearner(Y=psuedo_outcome,
                                           X=data %>% select(c(Xvars,'Astar','Ystar')),
                                           SL.library=sl.lib,
                                           obsWeights=data$R,
                                           verbose=FALSE)
  eif_proj <- predict(pseudo_mod,newdata=data)$pred

  # Rescale the EIF
  psuedo_rsc <- with(data,
                     (R/kappa_hat - 1)^(-1)*(R/kappa_hat)*psuedo_outcome)
  pseudo_mod_eem <- SuperLearner::SuperLearner(Y=psuedo_rsc,
                                           X=data %>% select(c(Xvars,'Astar','Ystar')),
                                           SL.library=sl.lib,
                                           obsWeights=(data$R/data$kappa_hat - 1)^2,
                                           verbose=FALSE)
  eif_proj_eem <- predict(pseudo_mod_eem,newdata=data)$pred

  # The TML is carried out in two updates. First, update kappa estimate
  kappa_new <- glm(R ~  -1 + I(eif_proj/kappa_hat), #+ I(eif_proj_eem/kappa_hat),
                   data=data,
                   offset=logit_transform(kappa_hat),
                   family=binomial())
  kappa_new <- predict(kappa_new,type='response')
  kappa_new <- trim(kappa_new,0.01)

  kappa_eem <- glm(R ~  -1 + I(eif_proj/kappa_hat) + I(eif_proj_eem/kappa_hat),
                   data=data,
                   offset=logit_transform(kappa_hat),
                   family=binomial())
  kappa_eem <- predict(kappa_eem,type='response')
  kappa_eem<- trim(kappa_eem,0.01)

  # Next, update the outcome regression

  # Construct clever covariates
  H_A <- with(data, A/pi_hat1 - (1-A)/(1-pi_hat1))
  H_1 <- with(data, 1/pi_hat1)
  H_0 <- with(data, -1/(1-pi_hat1))

  # Get eps estimates
  m_mod_new <- glm(Ysc ~ -1 + H_A,
                   data=data,
                   family=binomial(),
                   weights=data$R/kappa_new,
                   offset=logit_transform(m_hat_init))
  eps <- coef(m_mod_new)

  m_mod_eem <- glm(Ysc ~ -1 + H_A,
                   data=data,
                   family=binomial(),
                   weights=data$R/kappa_eem,
                   offset=logit_transform(m_hat_init))
  eps_eem <- coef(m_mod_eem)

  # With eps, we can update m1 and m0
  m1_new <- plogis(qlogis(shift_n_scale(m_hat1)) + eps*H_1)
  m0_new <- plogis(qlogis(shift_n_scale(m_hat0)) + eps*H_0)
  tml_hat_uns <- mean(m1_new - m0_new)
  tml_hat <- (max(data$Y) - min(data$Y))*tml_hat_uns

  # Do same with eem
  m1_eem <- plogis(qlogis(shift_n_scale(m_hat1)) + eps_eem*H_1)
  m0_eem <- plogis(qlogis(shift_n_scale(m_hat0)) + eps_eem*H_0)
  tml_hat_uns_eem <- mean(m1_eem - m0_eem)
  tml_hat_eem <- (max(data$Y) - min(data$Y))*tml_hat_uns_eem

  nuisances <- data.frame(
    Ysc=data$Ysc,
    kappa_hat=kappa_hat,
    m_hat1=m_hat1,
    m_hat0=m_hat0,
    pi_hat1=pi_hat1,
    pi_hat0=pi_hat0,
    varphi_hat_diff=eif_proj,
    kappa_new=kappa_new,
    m1_new=m1_new,
    m0_new=m0_new,
    phi_hat_diff=psuedo_outcome
  )

  # Then update the outcome regression estimate
  return(list(estimates=list(tml_hat=tml_hat,tml_hat_eem=tml_hat_eem),
              nuisances=nuisances))

  # Finally, do a weighted marginalization of the outcome regression

}

logit_transform <- function(x) {

  # First, check if x is bounded between 0 and 1. If not, transform it
  if (min(x) < 0 | max(x) > 1) {
    x <- shift_n_scale(x)
  }

  return(log(x/(1-x)))

}

shift_n_scale <- function(x,
                          buffer=FALSE) {
  #' Shifts and scales a vector to be between 0 and 1
  #'
  #' INPUTS:
  #' - x: a numeric vector
  #'
  #' OUTPUTS:
  #' - A numeric vector that has been shifted and scaled
  #'

  # No need to do anything if already on [0,1] interval
  if (min(x)>=0 & max(x)<=1) {
    return(x)
  }

  x <- (x - min(x))/(max(x) - min(x))
  if (buffer) {
    x[x==0] <- .Machine$double.neg.eps
    x[x==1] <- 1-.Machine$double.neg.eps
  }
  return(x)
}


tmle_evm <- function(data,X,Astar,Ystar,sl.lib) {
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

  Z <- cbind(X,Astar,Ystar)

  # First, get the initial re-weighted TMLE
  rvdl_tmle <- with(data,
                    tmle::tmle(Y=data$Y,A=data$A,W=X,
                               obsWeights = data$R/data$Rprobs,
                               g.SL.library=sl.lib,
                               Q.SL.library=sl.lib))

  # Extract the IC
  # Note: by construction IC is already weighted by R/Rprobs (tmle weights it)
  IC <- rvdl_tmle$estimates$IC$IC.ATE + rvdl_tmle$estimates$ATE$psi

  # Re-scale this IC by (R/Rprobs - 1)^{-1} to form the pseudo-outcome
  IC_rescaled <- with(data,(R/Rprobs - 1)^(-1)*IC)

  # Now regress IC on the rescaled Z with splines
  varphi_mod <- with(data,SuperLearner::SuperLearner(Y=IC_rescaled,
                                           X=Z, #data %>% select(Z),
                                           SL.library=sl.lib, # sl.lib,
                                           obsWeights = (data$R/data$Rprobs - 1)^2,
                                           verbose=FALSE))

  # Get the predicted IC
  varphi_Z <- predict(varphi_mod,newdata=data)$pred

  # Get final estimate
  tml_evm <- with(data, mean(IC - (R/Rprobs-1)*varphi_Z))
  tml_evmse <- with(data, sd(IC - (R/Rprobs-1)*varphi_Z)/sqrt(nrow(data)))

  # and get the version that fits the varphi conventionally
  # note: multiplying by Rprobs bc we want to reg the unweighted EIC on Z
  varphi_mid <- with(data,SuperLearner::SuperLearner(Y=IC*data$Rprobs,
                                                     X=Z,
                                                     SL.library=sl.lib,# sl.lib,
                                                     obsWeights = data$R,
                                                     verbose=FALSE))
  varphi_Z_mid <- predict(varphi_mid,newdata=data)$pred
  tml_mid <- with(data, mean(IC - (R/Rprobs-1)*varphi_Z_mid))
  tml_midse <- with(data, sd(IC - (R/Rprobs-1)*varphi_Z_mid)/sqrt(nrow(data)))

  return(list(tml_evm=tml_evm,tml_evmse=tml_evmse,
              tml_mid=tml_mid,tml_midse=tml_midse))


}



SL.glmnet_int <- function(Y, X, newX, family, obsWeights, id, alpha = 1,
                          nfolds = 10, nlambda = 100, useMin = TRUE, loss = "deviance", ...) {
  
  # Generate interaction terms
  X_interact <- model.matrix(~ .^2 - 1, data = X)  # Include all pairwise interactions, drop intercept
  newX_interact <- model.matrix(~ .^2 - 1, data = newX)  # Ensure same transformations
  
  # Fit cross-validated glmnet
  fitCV <- glmnet::cv.glmnet(x = X_interact, y = Y, weights = obsWeights,
                             lambda = NULL, type.measure = loss, nfolds = nfolds,
                             family = family$family, alpha = alpha, nlambda = nlambda,
                             ...)
  
  # Store column names to enforce consistency during prediction
  fit <- list(object = fitCV, useMin = useMin, Xnames = colnames(X_interact))
  class(fit) <- "SL.glmnet_int"
  
  # Predict using best lambda
  pred <- predict(fitCV, newx = newX_interact, type = "response",
                  s = ifelse(useMin, "lambda.min", "lambda.1se"))
  
  return(list(pred = pred, fit = fit))
}

# Prediction function for the custom learner
predict.SL.glmnet_int <- function(object, newdata, ...) {
  # Recreate interaction terms with same column names
  newX_interact <- model.matrix(~ .^2 - 1, data = newdata)[, object$Xnames, drop = FALSE]
  
  # Predict using stored glmnet model
  predict(object$object, newx = newX_interact, type = "response",
          s = ifelse(object$useMin, "lambda.min", "lambda.1se"))
}

SL.lightgbm <- function(Y, X, newX, family, obsWeights, id = NULL, ...) {
  require(lightgbm)
  require(Matrix)
  
  # LightGBM expects numeric matrices
  X_mat <- as.matrix(X)
  newX_mat <- as.matrix(newX)
  
  # Convert to LightGBM Dataset
  dtrain <- lgb.Dataset(data = X_mat, label = Y, weight = obsWeights)
  
  # Choose objective
  if (family$family == "gaussian") {
    objective <- "regression"
    metric <- "l2"
  } else if (family$family == "binomial") {
    objective <- "binary"
    metric <- "binary_logloss"
  } else {
    stop("SL.lightgbm only supports gaussian and binomial families.")
  }
  
  # Train LightGBM
  params <- list(objective = objective,
                 metric = metric,
                 verbosity = -1)
  
  model <- lgb.train(params = params,
                     data = dtrain,
                     nrounds = 100,
                     verbose = -1)
  
  # Predict
  pred <- predict(model, newX_mat)
  
  # If binomial, predictions are probabilities
  if (family$family == "binomial") {
    pred <- pmin(pmax(pred, 1e-6), 1 - 1e-6)  # clamp for stability
  }
  
  fit <- list(object = model)
  class(fit) <- "SL.lightgbm"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.lightgbm <- function(object, newdata, ...) {
  newX <- as.matrix(newdata)
  predict(object$object, newX)
}
