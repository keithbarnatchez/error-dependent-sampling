# Plots for the error-dependent sampling project
#-------------------------------------------------------------------------------
# Packages
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(reshape2)
library(ggridges)
library(ggtext)
rm(list=ls())
#-------------------------------------------------------------------------------
# Specify path to data

path <- '../output/sim_results-2025-04-23-07-16-47.487575.csv' #                                                                                                                                 '../output/sim_results-2025-02-03-11-23-37.csv'  # '../output/sim_results-2025-01-20-08-07-44.csv'   #  '../output/sim_results-2024-09-05-07-34-53.csv'

res <- read.csv(path)
#-------------------------------------------------------------------------------
# Main code

# Rename the plugin and onestep variables in res to pluginest and onestepest
res <- res %>% rename(pluginest = plugin, onestepest = onestep,
                      oracleest = oracle, naiveest = naive,
                      drcmdevmest=drcmd_evm,
                      onestepse=onestepse,
                      drcmdmidest=drcmd_mie,
                      drcmdevmdirectest=drcmd_evm_direct,
                      drcmdevmse = drcmd_evmse,
                      drcmdmidse = drcmd_midse,
                      drcmdtmlest = drcmd_tml,
                      wgtdest = wgtd_avg,
                      wgtdse=wgtd_avgse,
                      cvest= drcmd_cv) 

# Reshape to long format on the est vars
res_long <- res %>% 
  pivot_longer(cols = ends_with('est'),
               names_to = 'method',
               values_to = 'tau_est') %>%
  mutate(method=case_when(method == 'onestepest' ~ 'Approach 1 One-step',
                          method == 'pluginest' ~ 'Approach 1 Plug-in',
                          method == 'oracleest' ~ 'Oracle',
                          method == 'naiveest' ~ 'Naive',
                          method == 'onestep_adaptest' ~ 'Adaptive one-step',
                          method == 'plugin_adaptest' ~ 'Adaptive plug-in',
                          method == 'cvest' ~ 'TML-CV',
                          method == 'drcmdevmest' ~ 'Approach 2 One-step (EEM)',
                          method == 'drcmdmidest' ~ 'Approach 2 One-step',
                          method == 'drcmdevmdirectest' ~ 'drcmd EVM direct',
                          method == 'drcmdtmlest' ~ 'Approach 2 TML',
                          method == 'wgtdest' ~ 'Weighted',
                          method == 'cvest' ~ 'CV'
                          ))
#-------------------------------------------------------------------------------
# Set up dfs for plotting

refval = mean(res$oracleest) # actual treatment effect

# Summary stats for mean, SE, RMSE, things of that nature
grp_vec <- c('method','omega','n','rho','est_r')
res_summ <- res_long %>% group_by(across(all_of(grp_vec))) %>%
  summarize(per_bias = mean(100*((tau_est-refval)/abs(refval)),na.rm=T),
            avg_est=mean(tau_est,na.rm=T),
            sd_est=sd(tau_est,na.rm=T),
            rmse = sqrt(mean( (tau_est-refval)^2 , na.rm=T) ) ) %>%
  mutate(rho_desc = paste('rho =',rho))

res_ci<- res %>%
  mutate(
    # tmlcov = ((tmlest -  qnorm(0.975)*tmlse)  <= 2.5 ) & ((tmlest +  qnorm(0.975)*tmlse)  >= 2.5 ),
    onestepcov = ((onestepest -  qnorm(0.975)*onestepse)  <= refval ) & ((onestepest +  qnorm(0.975)*onestepse)  >= refval ),
    plugincov = ((pluginest -  qnorm(0.975)*pluginse)  <= refval ) & ((pluginest +  qnorm(0.975)*pluginse)  >= refval ),
    drcmdevmcov = ((drcmdevmest -  qnorm(0.975)*drcmdevmse)  <= refval ) & ((drcmdevmest +  qnorm(0.975)*drcmdevmse)  >= refval ),
    drcmdmidcov = ((drcmdmidest -  qnorm(0.975)*drcmdmidse)  <= refval ) & ((drcmdmidest +  qnorm(0.975)*drcmdmidse)  >= refval ),
    wgtdcov = ((wgtdest -  qnorm(0.975)*wgtdse)  <= refval ) & ((wgtdest +  qnorm(0.975)*wgtdse)  >= refval ),
    oraclecov = ((oracleest -  qnorm(0.975)*oraclese)  <= refval ) & ((oracleest +  qnorm(0.975)*oraclese)  >= refval ),
    naivecov = ((naiveest -  qnorm(0.975)*naivese)  <= refval ) & ((naiveest +  qnorm(0.975)*naivese)  >= refval )
  ) %>%
  pivot_longer(cols = ends_with('cov'),
               names_to = 'method',
               values_to = 'cov') %>%
  mutate(method=case_when(method == 'onestepcov' ~ 'Approach 1 One-step',
                          method == 'plugincov' ~ 'Approach 1 Plug-in',
                          method == 'oraclecov' ~ 'Oracle',
                          method == 'naivecov' ~ 'Naive',
                          method == 'drcmdevmcov' ~ 'Approach 2 One-step (EEM)',
                          method == 'drcmdmidcov' ~ 'Approach 2 One-step',
                          method=='wgtdcov' ~ 'Weighted' ))

res_summ_ci <- res_ci %>% group_by(across(all_of(grp_vec))) %>%
  summarize(cov = mean(cov,na.rm=T)) %>%
  mutate(rho_desc = paste('rho =',rho))
#-------------------------------------------------------------------------------
# Operating characteristics plot

opchar_df <- res_summ_ci %>% 
  right_join(res_summ, by=c(grp_vec,'rho_desc')) %>%
  pivot_longer(cols=c('per_bias','rmse','cov','avg_est','sd_est'),
               names_to = 'outcome',
               values_to = 'value') %>%
  filter(outcome %in% c('per_bias','rmse','cov')) %>%
  filter(!(method=='Approach 1 Plug-in' & outcome=='cov') ) %>%
  filter(method %in% c('Approach 2 One-step',
                       'Approach 2 One-step (EEM)',
                       'Oracle',
                       'Approach 1 One-step',
                       'Approach 1 Plug-in',
                       'Naive',
                       'Weighted',
                       'Approach 2 TML'
                       )) %>%
  filter(!(method=='Naivhe' & outcome=='rmse')) %>%
  mutate(outcome = case_when(outcome == 'per_bias' ~ 'Bias',
                             outcome == 'rmse' ~ 'RMSE',
                             outcome == 'cov' ~ 'Coverage'))

opchar_plot <- opchar_df %>% filter(method != 'Approach 2 TML') %>%
  filter(est_r==FALSE) %>%
  mutate(method = ifelse(method=='Weighted','Ensemble',method)) %>%
  ggplot(aes(x=rho,y=value,group=method,color=method)) +
  geom_line(size = 0.8, linetype = "solid") +  # Slightly thicker, solid lines
  geom_point(size = 2.5, shape = 6, stroke=1,
             position=position_jitter(w=0.0075,h=0)) +
  facet_grid(outcome ~ paste0('n=',n), scales='free') +
  geom_hline(data=subset(opchar_df,outcome=='Coverage'),aes(yintercept = 0.95), linetype = "dashed", color = "black") +
  geom_hline(data=subset(opchar_df,outcome=='Coverage'),aes(yintercept = 1), linetype = "dashed",alpha=0, color = "white") +
  theme_bw() +
  labs(title='**Operating characteristics of proposed estimators**',
       subtitle = 'Varying (1) overall sample size (2) relative size of validation data',
       x='Relative size of validation data, P(R=1)',
       y='Operating characteristics',
       color='Method') +
  theme(legend.position = 'bottom',
        plot.subtitle = element_markdown(),
        plot.title = element_markdown(),
        panel.grid.major = element_line(size = 0.40),  # Lighter gridlines
        panel.grid.minor = element_line(size=0),  # Remove minor gridlines for a cleaner look
        strip.background = element_blank(),  # Remove facet background for a sleeker look
        strip.text = element_text(size = 13, face = "bold"),  # Bold facet titles for emphasis
  ) ; opchar_plot
ggsave('../figures/opchars.pdf',opchar_plot,
       width=8,height=6,units = 'in')
#-------------------------------------------------------------------------------
# Supplement: estimating P(R=1|Z)

opchar_plot <- opchar_df %>% filter(method != 'Approach 2 TML') %>%
  mutate(method = ifelse(method=='Weighted','Ensemble',method)) %>%
  filter(est_r==TRUE) %>%
  ggplot(aes(x=rho,y=value,group=method,color=method)) +
  geom_line(size = 0.8, linetype = "solid") +  # Slightly thicker, solid lines
  geom_point(size = 2.5, shape = 6, stroke=1,
             position=position_jitter(w=0.0075,h=0)) +
  facet_grid(outcome ~ paste0('n=',n), scales='free') +
  geom_hline(data=subset(opchar_df,outcome=='Coverage'),aes(yintercept = 0.95), linetype = "dashed", color = "black") +
  geom_hline(data=subset(opchar_df,outcome=='Coverage'),aes(yintercept = 1), linetype = "dashed",alpha=0, color = "white") +
  theme_bw() +
  labs(title='**Operating characteristics of proposed estimators**',
       subtitle = 'Varying (1) overall sample size (2) relative size of validation data',
       x='Relative size of validation data, P(R=1)',
       y='Operating characteristics',
       color='Method') +
  theme(legend.position = 'bottom',
        plot.subtitle = element_markdown(),
        plot.title = element_markdown(),
        panel.grid.major = element_line(size = 0.40),  # Lighter gridlines
        panel.grid.minor = element_line(size=0),  # Remove minor gridlines for a cleaner look
        strip.background = element_blank(),  # Remove facet background for a sleeker look
        strip.text = element_text(size = 13, face = "bold"),  # Bold facet titles for emphasis
  ) ; opchar_plot
ggsave('../figures/opchars_estr.pdf',opchar_plot,
       width=8,height=6,units = 'in')

#-------------------------------------------------------------------------------
# Supplement: version with TML

opchar_plot_tml <- opchar_df %>% filter(!method %in% c('Approach 2 One-step',
                                                       'Approach 1 Plug-in')) %>%
  mutate(method = ifelse(method=='Weighted','Ensemble',method)) %>%
  filter(est_r==FALSE) %>%
  ggplot(aes(x=rho,y=value,group=method,color=method)) +
  geom_line(size = 0.8, linetype = "solid") +  # Slightly thicker, solid lines
  geom_point(size = 2.5, shape = 6, stroke=1, 
             position=position_jitter(w=0.0075,h=0)) +
  facet_grid(outcome ~ paste0('n=',n), scales='free') +
  geom_hline(data=subset(opchar_df,outcome=='Coverage'),aes(yintercept = 0.95), linetype = "dashed", color = "black") +
  geom_hline(data=subset(opchar_df,outcome=='Coverage'),aes(yintercept = 1), linetype = "dashed",alpha=0, color = "white") +
  theme_bw() +
  labs(title='**Operating characteristics of proposed estimators**',
       subtitle = 'Varying (1) overall sample size (2) relative size of validation data',
       x='Relative size of validation data, P(R=1)',
       y='Operating characteristics',
       color='Method') +
  theme(legend.position = 'bottom',
        plot.subtitle = element_markdown(),
        plot.title = element_markdown(),
        panel.grid.major = element_line(size = 0.40),  # Lighter gridlines
        panel.grid.minor = element_line(size=0),  # Remove minor gridlines for a cleaner look
        strip.background = element_blank(),  # Remove facet background for a sleeker look
        strip.text = element_text(size = 13, face = "bold"),  # Bold facet titles for emphasis
  ) ; opchar_plot_tml
ggsave('../figures/opchars_tml.pdf',opchar_plot_tml,
       width=8,height=6,units = 'in')
