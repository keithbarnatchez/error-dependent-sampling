# Plots for the error-dependent sampling project
#-------------------------------------------------------------------------------
# Packages
library(tidyverse)
library(reshape2)
library(ggridges)
library(ggtext)
#-------------------------------------------------------------------------------
# Specify path to data

path <- '../output/sim_results-2024-09-05-07-34-53.csv'

res <- read.csv(path) 
#-------------------------------------------------------------------------------
# Main code

# Rename the plugin and onestep variables in res to pluginest and onestepest
res <- res %>% rename(pluginest = plugin, onestepest = onestep,
                      oracleest = oracle, naiveest = naive,
                      onestep_adaptest=onestep_adapt,plugin_adaptest=plugin_adapt,
                      tmlest=tml,cvest=psi_hat_cv,tmlevmest=tml_evm)

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
                          method == 'tmlest' ~ 'TML',
                          method == 'tmlevmest' ~ 'Approach 2 One-step (EEM)',
                          method == 'cvest' ~ 'TML-CV')) 

# First, for tmt effect estimates
grp_vec <- c('method','omega','n','rho')
res_summ <- res_long %>% group_by(across(all_of(grp_vec))) %>%
  summarize(per_bias = mean(100*(tau_est-2.5),na.rm=T),
            rmse = sqrt(mean( (tau_est-2.5)^2 , na.rm=T) ) ) %>% 
  mutate(rho_desc = paste('rho =',rho))

ofyesplot <- res_long %>% filter(n %in% c(1000, 2500, 5000)) %>%
  filter(!(method %in% c('Adaptive one-step','Adaptive plug-in'))) %>%
  mutate(rho_desc = paste('P(R=1) =',rho)) %>%
  mutate(n_desc = paste('n =',n)) %>%
  filter(omega==0.2) %>%  ggplot(aes(x=tau_est,y=as.factor(method),
                         fill=as.factor(method))) + 
  xlim(c(-3,8)) +
  geom_density_ridges(alpha = 0.5,
                      scale = 1,
                      rel_min_height = 0.02,
                      quantile_lines = T, 
                      quantiles = 0.5, 
                      panel_scaling = F) + 
  facet_grid(fct_reorder(n_desc,n) ~ as.factor(rho_desc),scales='free') + theme_bw() +
  theme(legend.position = 'bottom',
        axis.text.y = element_blank(),  # Hide y-axis labels
        axis.ticks.y = element_blank(),
        plot.subtitle = element_markdown(),
        plot.title = element_markdown()) + 
  geom_vline(xintercept = 2.5, color='black',linetype='dashed') +
  labs(title='**Performance of one-step estimators**',
       subtitle = 'Varying (1) overall sample size (2) relative size of validation data<br>Distribution of treatment effect estimates across 2500 iterations',
       y='Method',
       x='ATE estimate',
       fill='Method') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) ; ofyesplot
ggsave('../figures/onestep-plugin-grid.pdf',ofyesplot,
       width=8,height=6,units = 'in')

ofyesplot <- res_long %>% filter(n %in% c(1000, 2500, 5000)) %>%
  filter(!(method %in% c('Adaptive one-step','Adaptive plug-in',
                         'TML','TML-CV'))) %>%
  mutate(rho_desc = paste('P(R=1) =',rho)) %>%
  mutate(n_desc = paste('n =',n)) %>%
  filter(omega==0.2) %>%  ggplot(aes(x=tau_est,y=as.factor(method),
                                     fill=as.factor(method))) + 
  xlim(c(-3,8)) +
  geom_density_ridges(alpha = 0.5,
                      scale = 1,
                      rel_min_height = 0.02,
                      quantile_lines = T, 
                      quantiles = 0.5, 
                      panel_scaling = F) + 
  facet_grid(fct_reorder(n_desc,n) ~ as.factor(rho_desc),scales='free') + theme_bw() +
  theme(legend.position = 'bottom',
        axis.text.y = element_blank(),  # Hide y-axis labels
        axis.ticks.y = element_blank(),
        plot.subtitle = element_markdown(),
        plot.title = element_markdown()) + 
  geom_vline(xintercept = 2.5, color='black',linetype='dashed') +
  labs(title='**Performance of one-step estimators**',
       subtitle = 'Varying (1) overall sample size (2) relative size of validation data<br>Distribution of treatment effect estimates across 2500 iterations',
       y='Method',
       x='ATE estimate',
       fill='Method') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) ; ofyesplot
ggsave('../figures/onestep-plugin-grid.pdf',ofyesplot,
       width=8,height=6,units = 'in')


# Using res_summ, plot rmse as a function of rho, with facetsf for sample size
# Omit the naive method to prevent scale from being distorted
rmseplot <- res_summ %>% filter(n %in% c(1000, 2500, 5000)) %>%
  filter(!(method %in% c('Adaptive one-step','Adaptive plug-in','Naive',
                         'TML','TML-CV'))) %>%
  mutate(n_desc = paste('n =',n)) %>%
  ggplot(aes(x=rho,y=rmse,group=method,color=method)) + 
  geom_line(size=0.7) + geom_point() +
  facet_wrap(~n_desc) + 
  labs(title='**RMSE of treatment effect estimates**',
       subtitle = 'Varying (1) overall sample size (2) relative size of validation data',
       x='P(R=1)',
       y='RMSE',
       color='Method') + 
  theme_bw() + 
  theme(legend.position = 'bottom',
        plot.subtitle = element_markdown(),
        plot.title = element_markdown()) ; rmseplot
ggsave('../figures/rmse-plot-n-grid.pdf',width=8,height=6,units = 'in') 

# Plot RMSEs with naive + oracle methods included
rmseplot <- res_summ %>% filter(n %in% c(1000, 2500, 5000)) %>%
  filter(!(method %in% c('Adaptive one-step','Adaptive plug-in','Naive',
                         'TML','TML-CV'))) %>%
  ggplot(aes(x=rho,y=rmse,group=method,color=method)) + 
  geom_line(size = 0.8, linetype = "solid") +  # Slightly thicker, solid lines
    geom_point(size = 2.5, shape = 21, fill = "white", stroke = 1.2) +
  facet_wrap(~paste0('n=',n)) + 
  labs(title='**RMSE of treatment effect estimates**',
       subtitle = 'Varying (1) overall sample size (2) relative size of validation data',
       x='P(R=1)',
       y='RMSE',
       color='Method') + 
  theme_bw() + 
  theme(legend.position = 'bottom',
        plot.subtitle = element_markdown(),
        plot.title = element_markdown(),
              panel.grid.major = element_line(size = 0.40),  # Lighter gridlines
              panel.grid.minor = element_line(size=0),  # Remove minor gridlines for a cleaner look
              strip.background = element_blank(),  # Remove facet background for a sleeker look
              strip.text = element_text(size = 13, face = "bold"),  # Bold facet titles for emphasis
            )  ; rmseplot
ggsave('../figures/rmse-plot-n-grid-with-naive.pdf',width=8,height=6,units = 'in') 


# res_df_plot %>%
#   filter(method %in% c('oracle','tmlevm'), outcome=='Average standard error') %>%
#   mutate(value=value/oracle_se) %>%
#   ggplot(aes(x = rho, y = value, group = method, color = method)) +
#   geom_line(size = 0.8, linetype = "solid") +  # Slightly thicker, solid lines
#   geom_point(size = 2.5, shape = 21, fill = "white", stroke = 1.2) +  
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = 'bottom',
#     panel.grid.major = element_line(size = 0.40),  # Lighter gridlines
#     panel.grid.minor = element_line(size=0),  # Remove minor gridlines for a cleaner look
#     strip.background = element_blank(),  # Remove facet background for a sleeker look
#     strip.text = element_text(size = 13, face = "bold"),  # Bold facet titles for emphasis
#     plot.title = element_text(face = "bold", size = 16),  # Slightly larger and bold title
#     plot.subtitle = element_text(size = 12, margin = margin(b = 10))  # Add spacing for subtitle
#   ) + 
#   labs(
#     x = 'Relative size of the validation data, P(R=1)', 
#     y = 'Average relative efficiency', 
#     title = 'Outcome: 3-year ADE risk', 
#     subtitle = 'Exposure: early ART initiation', 
#     color = 'Method'
#   ) 





         