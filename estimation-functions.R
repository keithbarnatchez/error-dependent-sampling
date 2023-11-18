# estimation-functions.R
#
# Functions for performing estimation (plug-in and one-step) in the measurement
# error-dependent sampling project
#-------------------------------------------------------------------------------

est_mu_a <- function(data) {
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
              SL.library = c("SL.mean","SL.glmnet","SL.glm","SL.step"), 
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

est_lambda_a <- function(data) {
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

est_kappa <- function(data) {
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
              SL.library = c("SL.mean","SL.glmnet","SL.glm"), 
              family = binomial(), method = "method.NNLS", verbose = FALSE)
  
  # Extract the predictions from superlearner model
  kappa <- predict(kappa_mod,
                   data,
                   onlySL=TRUE)$pred
  
  return(kappa)
}

est_eta_a <- function(data,sl.lib='') {
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
  
  rhsvars <- c(names(data)[grep("^X",names(data))],'Astar')
  covs <- data %>% select(all_of(rhsvars))
  
  # Regress beta_hat on X with superlearner
  # Need separate models for A=0 and A=1
  eta_a_mod1 <- SuperLearner::SuperLearner(Y=beta_hat1, X=covs, 
              SL.library = c("SL.mean","SL.glmnet","SL.glm","SL.step"), 
              family = gaussian(), method = "method.NNLS", verbose = FALSE)
  eta_a_mod0 <- SuperLearner::SuperLearner(Y=beta_hat0, X=covs, 
                SL.library = c("SL.mean","SL.glmnet","SL.glm","SL.step"), 
                family = gaussian(), method = "method.NNLS", verbose = FALSE)
  
  # Get the predicted values for the whole dataset
  eta_hat1 <- predict(eta_a_mod1,
                      data,
                      onlySL=TRUE)$pred
  eta_hat0 <- predict(eta_a_mod0,
                      data,
                      onlySL=TRUE)$pred
  
  return(data.frame(eta_hat1,eta_hat0))
  
}

est_pi_a <- function(data) {
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
  rhsvars <- c(names(data)[grep("^X",names(data))],'Astar')
  covs <- data %>% select(all_of(rhsvars))
  
  # Need separate models for A=0 and A=1
  # pi_1 = E(lambda_1 | X, Astar)
  # Note: setting family to binomial since pi_a should be in (0,1)
  pi_a_mod1 <- SuperLearner::SuperLearner(Y=data$lambda_hat1, X=covs, 
              SL.library = c("SL.mean","SL.glm","SL.step"), 
              family = binomial(), method = "method.NNLS", verbose = FALSE)
  # pi_1 = E(lambda_0 | X, Astar) = E(1-lambda_1 | X, Astar)
  pi_a_mod0 <- SuperLearner::SuperLearner(Y=data$lambda_hat0, X=covs, 
                SL.library = c("SL.mean","SL.glm","SL.step"), 
                family = binomial(), method = "method.NNLS", verbose = FALSE)
  
  # Get the predicted values for the whole dataset
  pi_hat1 <- predict(pi_a_mod1,
                      data,
                      onlySL=TRUE)$pred
  pi_hat0 <- predict(pi_a_mod0,
                      data,
                      onlySL=TRUE)$pred
  
  return(data.frame(pi_hat1,pi_hat0))
}

est_psi_a <- function(data,cross_fit=FALSE) {
  #' Function for estimating 
  #'.    E(Y(a)) = E(eta_a(X)/pi_a(X))
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
  mu_a <- est_mu_a(data)
  lambda_a <- est_lambda_a(data)
  
  # Append new predictions onto the dataframe
  data <- cbind(data,mu_a,lambda_a)
  
  # Regress mu_a and lambda_a predictions on X and Astar across the whole dataset
  # to get eta_a predictions
  eta_a <- est_eta_a(data)
  data <- cbind(data,eta_a)
  
  # Regress lambda_a predictions on X and Astar across the whole dataset
  pi_a <- est_pi_a(data)
  data <- cbind(data,pi_a)
  
  # Get P(R=1|X,Ystar,Astar)
  kappa_hat <- est_kappa(data)
  
  # Get the plug-in estimate and its SE
  plugin <- mean(data$eta_hat1/data$lambda_hat1) - mean(data$eta_hat0/data$lambda_hat0)
  pluginse <- sd(data$eta_hat1/data$pi_hat1 - data$eta_hat0/data$pi_hat0)/sqrt(nrow(data))
  
  # Get the one-step bias corrected estimate and its SE
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
  onestep <- mean(psi1hat - psi0hat)
  onestepse <- sd(psi1hat - psi0hat)/sqrt(nrow(data))
  
  return(list(plugin=plugin,pluginse=pluginse,
              onestep=onestep,onestepse=onestepse,
              maxl1=max(data$lambda_hat1),maxpi1=max(data$pi_hat1),
              mink=min(kappa_hat),maxk=max(kappa_hat)))
}

