library(tidyverse)
library(patchwork)
library(brms)
library(orthopolynom)
library(SobolSequence)
library(kableExtra)
library(latex2exp)
source("PCE_helpers.R")
source("ishigami_helpers.R")
source("plot_helpers.R")
theme_set(theme_bw())

# plots over training points
ishigami_summary <- readRDS("ishigami_summary_noise.rds") %>%
  bind_rows(
    readRDS("ishigami_summary.rds") %>%
      mutate(noise = 0)
  ) %>%
  mutate(
    bias_mean = PCE_mean_estimate - ishigami_mean,
    bias_sd = PCE_sd_estimate - ishigami_sd,
    L2AE_mean = abs(PCE_mean_estimate - ishigami_mean),
    L2AE_sd = abs(PCE_sd_estimate - ishigami_sd),
    type = factor(type, levels = c("full", "b_sel", "varsel"),
                  labels = c("Full", "Sel-Sobol", "Sel-Projpred")),
    noise_sd = paste0("sigma[noise]~'='~", noise),
    time_per_chain_minutes = time_total / nchains / 60
  ) %>%
  group_by(mc, type, noise_sd) %>%
  summarize(
    bias_mean_avg = mean(bias_mean),
    bias_mean_min = min(bias_mean),
    bias_mean_max = max(bias_mean),
    bias_sd_avg = mean(bias_sd),
    bias_sd_min = min(bias_sd),
    bias_sd_max = max(bias_sd),
    L2AE_mean_avg = mean(L2AE_mean),
    L2AE_mean_min = min(L2AE_mean),
    L2AE_mean_max = max(L2AE_mean),
    L2AE_sd_avg = mean(L2AE_sd),
    L2AE_sd_min = min(L2AE_sd),
    L2AE_sd_max = max(L2AE_sd),
    RMSE_oos_mean_avg = mean(RMSE_oos_mean),
    RMSE_oos_mean_min = min(RMSE_oos_mean),
    RMSE_oos_mean_max = max(RMSE_oos_mean),
  ) %>%
  ungroup() %>%
  rename(Model = type) %>%
  # D2 selection is almost equivalent to Sobol selection
  # no need to introduce and discuss this method
  filter(Model != "Sel-D2")


xbreaks <- sort(unique(ishigami_summary$mc))
xbreaks <- setdiff(xbreaks, 286)
xlims <- c(8, 850)


# bias of Mean estimate
gg_ishigami_bias_mean_summary <- ishigami_summary %>%
  ggplot(aes(mc, bias_mean_avg, ymin = bias_mean_min, ymax = bias_mean_max, 
             color = Model, fill = Model)) +
  geom_smooth(size = 0.8, stat = "identity") +
  facet_grid(cols = vars(noise_sd), labeller=label_parsed) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  ylab(expression(mu(hat(y))-mu(y))) +
  xlab("T") +
  # log scale makes for a quite ugly plot
  # scale_y_continuous(
  #   trans=make_lal_trans(
  #     'trexp',
  #     threshold=0.00001,
  #     exponent=10,
  #     force_thresholds_in_breaks = TRUE
  #   ),
  #   breaks = c(-0.1, -0.01, -0.001, -0.0001, 0.0001, -0.00001,
  #             0, 0.00001, 0.0001, 0.001, 0.01, 0.1),
  #   limits = c(-2.5, 2.5)
  # ) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                      limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_ishigami_bias_mean_summary


# absolute bias of Mean estimate
gg_ishigami_L2AE_mean_summary <- ishigami_summary %>%
  ggplot(aes(mc, L2AE_mean_avg, ymin = L2AE_mean_min, ymax = L2AE_mean_max, 
             color = Model, fill = Model)) +
  geom_smooth(size = 0.8, stat = "identity") +
  facet_grid(cols = vars(noise_sd), labeller=label_parsed) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  ylab(expression("|"~mu(hat(y))-mu(y)~"|")) +
  xlab("T") +
  scale_y_continuous(trans='log10', breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1),
                     limits = c(0.00001, 10)) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                   limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_ishigami_L2AE_mean_summary


# bias of SD estimate
gg_ishigami_bias_sd_summary <- ishigami_summary %>%
  # truncate data for better plotting (affects Sel-Sobol only)
  # mutate(bias_sd_min = ifelse(bias_sd_min < -3, -3, bias_sd_min),
  #        bias_sd_max = ifelse(bias_sd_max > 3, 3, bias_sd_max)) %>%
  ggplot(aes(mc, bias_sd_avg, ymin = bias_sd_min, ymax = bias_sd_max, 
             color = Model, fill = Model)) +
  geom_smooth(size = 0.8, stat = "identity") +
  facet_grid(cols = vars(noise_sd), labeller=label_parsed) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  ylab(expression(sigma(hat(y))-sigma(y))) +
  xlab("T") +
  # log scale makes for a quite ugly plot
  # scale_y_continuous(
  #   trans=make_lal_trans(
  #     'trexp',
  #     threshold=0.00001,
  #     exponent=10,
  #     force_thresholds_in_breaks = TRUE
  #   ),
  #   breaks = c(-0.1, -0.01, -0.001, -0.0001, 0.0001, -0.00001,
  #             0, 0.00001, 0.0001, 0.001, 0.01, 0.1),
  #   limits = c(-2.5, 2.5)
  # ) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_ishigami_bias_sd_summary


# absolute bias of SD estimate
gg_ishigami_L2AE_sd_summary <- ishigami_summary %>%
  # truncate data for better plotting (affects Sel-Sobol only)
  # mutate(L2AE_sd_max = ifelse(L2AE_sd_max > 10, 10, L2AE_sd_max)) %>%
  ggplot(aes(mc, L2AE_sd_avg, ymin = L2AE_sd_min, ymax = L2AE_sd_max, 
             color = Model, fill = Model)) +
  geom_smooth(size = 0.8, stat = "identity") +
  facet_grid(cols = vars(noise_sd), labeller=label_parsed) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  ylab(expression("|"~sigma(hat(y))-sigma(y)~"|")) +
  xlab("T") +
  scale_y_continuous(trans='log10', breaks = c(0.001, 0.01, 0.1, 1),
                     limits = c(0.001, 20)) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_ishigami_L2AE_sd_summary


# out of sample RMSE
gg_ishigami_RMSE_oos_summary <- ishigami_summary %>%
  ggplot(aes(mc, RMSE_oos_mean_avg, ymin = RMSE_oos_mean_min, 
             ymax = RMSE_oos_mean_max, color = Model, fill = Model)) +
  geom_smooth(size = 0.8, stat = "identity") +
  facet_grid(cols = vars(noise_sd), labeller=label_parsed) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  ylab(expression("Out-of-sample RMSE")) +
  xlab("T") +
  scale_y_continuous(trans='log10', breaks = c(0.01, 0.1, 1)) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_ishigami_RMSE_oos_summary

# combine the four plots
gg_ishigami_L2AE_mean_summary /
  gg_ishigami_L2AE_sd_summary /
  gg_ishigami_RMSE_oos_summary /
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom") 

ggsave("plots/ishigami_summaries_noise.jpg", width = 10, height = 8)

