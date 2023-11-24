# data-functions.R
#
# Code for simulating data in the measurement error-dependent sampling project
#-------------------------------------------------------------------------------
expit <- function(o) {
  return(exp(o)/(1+exp(o)))
}

trim <- function(p,cutoff=0.05) {
  #' Function for trimming probabilities
  #' 
  #' Arguments:
  #' - p: A vector of probabilities
  #' - cutoff: The cutoff value
  #' 
  #' Outputs:
  #' - A vector of probabilities, with values below/above the cutoff set to the 
  #'  cutoff
  
  # trim lower piece
  p <- ifelse(p<cutoff,cutoff,p)
  
  # trim upper piece
  p <- ifelse(p>(1-cutoff),1-cutoff,p)
  
  return(p)
}

gen_covariates <- function(n,p) {
  #' Function for simulating covariates
  #'
  #' Simulates p covariates from a multivariate uniform(0,1) distribution
  #' 
  #' Arguments:
  #' - n: Sample size
  #' - p: Number of covariates
  #' 
  #' Outputs:
  #' - A dataframe containing the covariates
  
  X <- matrix(runif(n*p),nrow=n,ncol=p)
  return(X)
  
}

gen_treatment <- function(X, delta) {
  #' Function for simulating treatment assignment
  #'
  #' Arguments:
  #' - X: Covariate matrix
  #' - delta: Treatment effect
  #' 
  #' Outputs:
  #' - A vector of treatment assignments
  
  probs <- trim(expit(as.matrix(X)%*%delta))
  A <- rbinom(nrow(X),size=1,prob = probs)
  return(A)
}

gen_selection <- function(X,Astar,Ystar,eta) {
  #' Function for simulating selection into the validation study 
  #' 
  #' Arguments:
  #' - X: Covariate matrix
  #' - eta: Selection coefficients
  #' 
  #' Outputs:
  #' - A vector of selection indicators
  #' 
  
  probs <- trim(expit(as.matrix(cbind(X,Astar,Ystar))%*%eta))
  R <- rbinom(nrow(X),size=1,prob = probs)
  
  return(R)
  
}

treatment_meas_model <- function(A, X, omega) {
  #' Function for simulating treatment measurement
  #' 
  #' Arguments:
  #' - A: Treatment assignment
  #' - X: Covariate matrix
  #' - omega: Measurement error variance
  #' 
  #' Outputs:
  #' - A vector of treatment measurements
  #' 
  
  # Simulate treatment measurements using a binary misclassification model,
  # where the misclassification probabilities depend on A and X
  Astar <- ifelse(runif(length(A)) < 0.85, A, 1 - A)
  
  return(Astar)
  
}

outcome_meas_model <- function(Y, A, X, nu) {
  #' Function for simulating outcome measurement
  #' 
  #' Arguments:
  #' - A: Treatment assignment
  #' - X: Covariate matrix
  #' - nu: Measurement error variance
  #' 
  #' Outputs:
  #' - A vector of outcome measurements
  #' 
  
  # Simulate outcome measurements 
  Ystar <- Y + X%*%nu + rnorm(nrow(X),mean=0,sd=1)
  
  return(Ystar)
  
}

outcome_model <- function(A, X, beta, tau, gamma) {
  #' Function for simulating outcome measurements
  #' 
  #' Arguments:
  #' - A: Treatment assignment
  #' - X: Covariate matrix
  #' - beta, tau, gamma: Outcome model coefficients
  #' - scenario: Simulation scenario (1=sharp causal null, 2=const tmt effect, 3=het tmt effect)
  #' 
  #' Outputs:
  #' - A vector of outcome measurements
  #' 
  
  # Simulate outcomes 
  Y <- X%*%beta + tau*A + A*(X%*%gamma) + rnorm(nrow(X),mean=0,sd=1)
  
  return(Y)
  
}

gen_data <- function(n,p,delta,nu,beta,gamma,eta,omega,tau) {
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
  
  # Simulate covariates
  X <- gen_covariates(n,p)
  
  # Simulate treatment
  A <- gen_treatment(X,delta)
  
  # Treatment measurement
  Astar <- treatment_meas_model(A, X, omega)
  
  # Outcome
  Y <- outcome_model(A,X,beta,tau,gamma)
  
  # Outcome measurement
  Ystar <- outcome_meas_model(Y,A,X,nu)
  
  # Selection into validation data
  R <- gen_selection(X,Astar,Ystar,eta)
  
  return(data.frame(X,A,Astar,Y,Ystar,R))
}



