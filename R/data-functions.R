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

#' @title Simulate covariates
#'
#' @description Simulates p covariates from a multivariate uniform(0,1) distribution
#'
#' @param n Sample size
#' @param p Number of covariates
#'
#' @return A dataframe containing the covariates
#' @export
gen_covariates <- function(n,p) {

  X <- matrix(runif(n*p),nrow=n,ncol=p)
  return(X)

}

#' @title Generate treatment
#'
#' @description Function for simulating treatment assignment from logistic
#' propensity score model
#'
#' @param X Covariate matrix
#' @param delta Coefficients of propensity score model
#'
#' @return A vector of treatment assignments
#' @export
gen_treatment <- function(X, delta) {

  probs <- trim(expit(as.matrix(X)%*%delta))
  A <- rbinom(nrow(X),size=1,prob = probs)
  return(A)
}

#' @title Generate validation indicators
#'
#' @description Function for simulating selection into validation data
#'
#' @param X Covariate matrix
#' @param Astar Error-prone treatment measurement
#' @param Ystar Error-prone outcome measurement
#' @param eta Selection coefficients
#' @param rho Selection probability
#' @return A list containing the selection indicators and selection probabilities
#' @export
gen_selection <- function(X,Astar,Ystar,eta,
                          rho) {

  probs <- rho*trim(expit(as.matrix(cbind(X,Astar,Ystar))%*%eta))/mean(trim(expit(as.matrix(cbind(X,Astar,Ystar))%*%eta)))
  R <- rbinom(nrow(X),size=1,prob = probs)

  return(list(R=R,Rprobs=probs))

}

#' @title Treatment measurement model
#'
#' @description Function for simulating treatment measurement
#'
#' @param A True treatment values
#' @param X Covariate matrix
#' @param omega Misclassification rate
#'
#' @return A vector of error-prone treatment measurements
#' @export
treatment_meas_model <- function(A, X, omega=0.95) {

  # Simulate treatment measurements using a binary misclassification model,
  # where the misclassification probabilities depend on A and X
  Astar <- ifelse(runif(length(A)) < omega, A, 1 - A)

  return(Astar)

}


#' @title Outcome measurement model
#'
#' @description Function for simulating outcome measurement
#'
#' @param Y True outcome values
#' @param A True treatment assignment
#' @param X Covariate matrix
#' @param nu Coefficientson X (allows for systematic measurement error)
#'
#' @return A vector of error-prone outcome measurements
#' @export
outcome_meas_model <- function(Y, A, X, nu) {

  # Simulate outcome measurements
  Ystar <- Y + X%*%nu + rnorm(nrow(X),mean=0,sd=1)

  return(Ystar)

}

#' @title Outcome model
#'
#' @description Function for simulating outcome measurements
#'
#' @param A Treatment assignment vector
#' @param X Covariate matrix
#' @param beta Coefficients on X
#' @param tau Treatment effect
#' @param gamma Coefficients on X for treatment effect
#'
#' @return A vector of outcome measurements
#' @export
outcome_model <- function(A, X, beta, tau, gamma) {

  # Simulate outcomes
  Y <- X%*%beta + tau*A + A*(X%*%gamma) + rnorm(nrow(X),mean=0,sd=1)

  return(Y)

}

#' @title Generate data
#'
#' @description Main function for simulating data
#'
#' @param n Sample size
#' @param p Number of covariates
#' @param delta Treatment effect
#' @param nu Outcome measurement error variance
#' @param beta Outcome model covariate coefficients
#' @param gamma Outcome model interaction coefficients
#' @param eta Selection model coefficients
#' @param tau Outcome model treatment effect coefficient
#' @param omega Treatment measurement error rate
#' @param rho Selection probability
#'
#' @return A dataframe containing the simulated data
#' @export
gen_data <- function(n,p,delta,nu,beta,gamma,eta,omega,tau,
                     rho) {

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
  Rres <- gen_selection(X,Astar,Ystar,eta,rho)
  R <- Rres$R ; Rprobs <- Rres$Rprobs

  return(data.frame(X,A,Astar,Y,Ystar,R,Rprobs))
}
