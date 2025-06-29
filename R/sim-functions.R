# sim-functions.R (CLUSTER VERSION)
#
# Code for implementing simulations in the error-dependent sampling project
#-------------------------------------------------------------------------------

#' @title Run single simulation iterations
#'
#' @description Function for running a single simulation iteration called by main_sim
#'
#' @param n Sample size
#' @param p Number of covariates (all iid Uniform(0,1))
#' @param delta
#' @param nu Outcome measurement error variance
#' @param beta Outcome model main covariate effect coefficients
#' @param gamma Outcome model interaction coefficients
#' @param eta
#' @param omega
#' @param tau
#' @param rho
#' @param rho_s
#' @param sl.lib
#'
#' @return A dataframe containing the simulation results
#' @export
sim_iter <- function(n,p,delta,nu,beta,gamma,eta,omega,tau,rho,rho_s,
                     est_r,
                     sl.lib) {

 # df <- gen_data(n,p,delta,nu,beta,gamma,eta,omega,tau,rho)
  X <- matrix(runif(n*p),nrow=n,ncol=p)

  # treatment
  probs <- trim(expit(as.matrix(X)%*%delta))
  A <- rbinom(nrow(X),size=1,prob = probs)

  # outcome
  Y <- X%*%beta + tau*A + A*(X%*%gamma) + rnorm(nrow(X),mean=0,sd= 1+rowSums(X)) # exp(rowSums(X)))
  # Y <- rbinom(n,1,
  #             plogis(-rowSums(X) - A))

  # Y <- rbinom(n,1,0.5)
  # measurements
  Ystar <- Y + X%*%nu + rnorm(nrow(X),mean=0,sd=1)
  Astar <- ifelse(runif(length(A)) < 1-omega, A, 1 - A)
  # Ystar <- ifelse(runif(length(Y)) < 1-omega, Y, 1 - Y)

  Rprobs <- rho*trim(expit(as.matrix(cbind(X,Astar,Ystar))%*%eta))/mean(trim(expit(as.matrix(cbind(X,Astar,Ystar))%*%eta)))
  # Rprobs <- rep(rho,n) # rho*expit(as.matrix(cbind(X,Astar,Ystar))%*%eta)/mean(expit(as.matrix(cbind(X,Astar,Ystar))%*%eta))
  R <- rbinom(nrow(X),size=1,prob = Rprobs)
  X <- as.data.frame(X)
  colnames(X) <- c('X1','X2','X3')

  df <- cbind(X,
              data.frame(A=A,Y=Y,Astar=Astar,Ystar=Ystar,R=R,Rprobs=Rprobs))

  # Estimate psi with the ad-hoc rule
  psi_hat <- est_psi_a(df,
                       est_r,sl.lib)
  a1_eif <- as.double(psi_hat$eif)
  psi_hat$eif <- NULL

  # Turn into df
  psi_hat <- data.frame(t(unlist(psi_hat)),stringsAsFactors = FALSE)

  # Finally, implement approach 2 methods with the drcmd package
  # drcmd requires missing values be entered as NA, so set those up
  df <- df %>% mutate(Y = ifelse(R==0,NA,Y),
                      A = ifelse(R==0,NA,A))

  if (est_r) {
    drcmd_evm_res <- drcmd::drcmd(Y=df$Y,A=df$A,X=X,
                                 W=df %>% select(Ystar, Astar),
                                 default_learners=sl.lib,
                                 r_learners = sl.lib,
                                 eem_ind=TRUE)
    drcmd_mid_res <- drcmd::drcmd(Y=df$Y,A=df$A,X=X,
                                  W=df %>% select(Ystar, Astar),
                                  default_learners=sl.lib,
                                  r_learners = sl.lib,
                                  eem_ind=FALSE)

    psi_hat_tml_drcmd <- drcmd::drcmd(Y=df$Y,A=df$A,X=X,
                                      W=df %>% select(Ystar, Astar),
                                      default_learners=sl.lib,
                                      r_learners = sl.lib,
                                      eem_ind=FALSE,tml=T)
  } else {
      drcmd_evm_res <- drcmd::drcmd(Y=df$Y,A=df$A,X=X,
                                    W=df %>% select(Ystar, Astar),
                                    default_learners=sl.lib,
                                    Rprobs=df$Rprobs,
                                    eem_ind=TRUE)
      drcmd_mid_res <- drcmd::drcmd(Y=df$Y,A=df$A,X=X,
                                    W=df %>% select(Ystar, Astar),
                                    default_learners=sl.lib,
                                    Rprobs=df$Rprobs,
                                    eem_ind=FALSE)

      psi_hat_tml_drcmd <-  drcmd::drcmd(Y=df$Y,A=df$A,X=X,
                                        W=df %>% select(Ystar, Astar),
                                        default_learners=sl.lib,
                                        Rprobs=df$Rprobs,
                                        eem_ind=FALSE,tml=T)
  }

  drcmd_tml <- data.frame(drcmd_tml = psi_hat_tml_drcmd$results$estimates$psi_hat_ate)
  
  # Form CV estimator
  # get eif from onestep
  IF1 <- mean(drcmd_evm_res$results$nuis$m_1_hat) + drcmd_evm_res$results$nuis$varphi_1_hat + df$R/df$Rprobs*(drcmd_evm_res$results$nuis$phi_1_hat - drcmd_evm_res$results$nuis$varphi_1_hat)
  IF0 <- mean(drcmd_evm_res$results$nuis$m_0_hat) + drcmd_evm_res$results$nuis$varphi_0_hat + df$R/df$Rprobs*(drcmd_evm_res$results$nuis$phi_0_hat - drcmd_evm_res$results$nuis$varphi_0_hat)
  a2_eif <- IF1-IF0
  est_cv_mod <- lm(a2_eif ~ I(a2_eif-a1_eif))
  est_cv <- est_cv_mod$coefficients[1]
  est_cvse <- vcov(est_cv_mod)[1,1]
  
  # Get the optimal combo
  opt_weight <- (var(a2_eif) - cov(a2_eif,a1_eif))/(var(a2_eif) + var(a1_eif) - 2*cov(a2_eif,a1_eif))
  wgtd_avg <- opt_weight*mean(a1_eif) + (1-opt_weight)*mean(a2_eif)
  wgtd_avg_eif <-  opt_weight*a1_eif + (1-opt_weight)*a2_eif
  
  # Compile results
  drcmd_evm <- data.frame(drcmd_evm=drcmd_evm_res$results$estimates$psi_hat_ate,
                          drcmd_evmse=drcmd_evm_res$results$ses$psi_hat_ate,
                          drcmd_evm_direct=drcmd_evm_res$results$estimates$psi_hat_ate_direct,
                          drmd_evm_directse=drcmd_evm_res$results$ses$psi_hat_ate_direct,
                          drcmd_cv=est_cv,
                          wgtd_avg=wgtd_avg,
                          wgtd_avgse=sd(wgtd_avg_eif)/sqrt(length(wgtd_avg_eif)))

  drcmd_mid <- data.frame(drcmd_mie=drcmd_mid_res$results$estimates$psi_hat_ate,
                          drcmd_midse=drcmd_mid_res$results$ses$psi_hat_ate)
  
  # Unpack output to store in df
  psi_hat <- cbind(psi_hat,
                   drcmd_evm, drcmd_mid, drcmd_tml)

  # Return data
  return(psi_hat)
}

#' @title Main simulation function
#'
#' @description Function for running simulation iterations, for a fixed set of
#' parameters, for <paper title here>
#'
#' @param nsim Number of simulations to run
#' @param params A list containing all simulation parameters, collected in the
#' file sim_main.R
#' @param parallel Logical indicating whether to run simulations in parallel
#'
#' @return A dataframe containing the simulation results
#' @export
main_sim <- function(nsim, params,
                     parallel = TRUE) {

  # Extract names from params
  list2env(params, envir = environment())

  if (parallel==TRUE) {
    # Register parallel backend
    registerDoParallel(cores = 8)#  detectCores())

    # Carry out the simulation
    # SL fails very, very rarely due to numerical issues which would crash whole
    # simulation, so catch those instances and re-run that iteration
    sim_results <- foreach(i = 1:nsim, .combine = rbind) %dopar% {
      attempt <- 1
      max_attempts <- 3
      success <- FALSE
      result <- NULL
      while (attempt <= max_attempts & !success) { #
        tryCatch({
          result <- sim_iter(n, p, delta, nu, beta, gamma, eta, omega, tau, rho, rho_s,
                             est_r, sl.lib)
          success <- TRUE  # Mark as successful if no error occurs
        }, error = function(e) {
          attempt <- attempt + 1
          if (attempt > max_attempts) {
            stop("Simulation failed after 3 attempts")
          }
        })
      }
      result
    }
    stopImplicitCluster()
  } else {
    # Carry out the simulation
    sim_results <- foreach(i = 1:nsim, .combine = rbind) %do% {
      sim_iter(n,p,delta,nu,beta,gamma,eta,omega,tau,rho,rho_s,
               est_r,sl.lib)
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
