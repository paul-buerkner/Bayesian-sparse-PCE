library(tidyverse)
library(patchwork)
library(brms)
library(orthopolynom)
library(SobolSequence)
library(kableExtra)
source("PCE_helpers.R")
source("sobol_helpers.R")
source("plot_helpers.R")
theme_set(theme_bw())

a <- c(1, 2, 5, 10, 20, 50, 100, 500)

fit_ref <- brm(file = "models/fit_sobol_mc300_M8")
var_ref <- variables(fit_ref)[c(1, 2, 2508, 1386)+1]
gg_pairs <- plot(pairs(fit_ref, variable = var_ref))
plot(gg_pairs)
ggsave("plots/sobol_pairs_mc300_M8.jpg", gg_pairs, width = 10, height = 7)


mcmc_plot(fit_ref, variable = var_ref, type = "trace",
          facet_args = list(nrow = 1), size = 2)
ggsave("plots/sobol_trace_mc300_M8.jpg", width = 10, height = 3)


# plots for some sparse projpred models
fit_M8_varsel <- brm(file = "models/fit_sobol_mc300_M8_varsel")

plot(conditional_effects_sobol(fit_M8_varsel, "x1s", a = a)) +
plot(conditional_effects_sobol(fit_M8_varsel, "x2s", a = a)) +
plot(conditional_effects_sobol(fit_M8_varsel, "x3s", a = a)) +
plot(conditional_effects_sobol(fit_M8_varsel, "x4s", a = a)) +
plot(conditional_effects_sobol(fit_M8_varsel, "x5s", a = a)) +
plot(conditional_effects_sobol(fit_M8_varsel, "x6s", a = a)) +
plot(conditional_effects_sobol(fit_M8_varsel, "x7s", a = a)) +
plot(conditional_effects_sobol(fit_M8_varsel, "x8s", a = a)) +
  plot_layout(nrow = 2)

ggsave("plots/sobol_ce_mc300_M8_varsel.jpg", width = 8, height = 5)


plot(plot_poly_degree(fit_M8_varsel, p = 6, M = 8) + 
       scale_x_continuous(breaks = 1:6, limits = c(1, 6.3))) +
  plot(plot_sobol_coef_degree(fit_M8_varsel, yintercept = sobol_sd_num(a)))

ggsave("plots/sobol_poly_mc300_M8_varsel.jpg", width = 8, height = 3)

# table for the most important polynomials and the degrees
df_M8_varsel <- df_poly_degree(fit_M8_varsel, p = 6, M = 8, digits = 3)
df_M8_varsel %>%
  filter(row_number() %in% 1:14) %>%
  kbl(format= "latex", align=c("l", "l" ,"r", "r", "r", "r")) %>%
  kable_classic(full_width = F, html_font = "helvetica")



fit_M4_varsel <- brm(file = "models/fit_sobol_mc300_M4_varsel")

# conditional effects plots are not very informative here
plot(conditional_effects_sobol(fit_M8_varsel, "x1s", a = a)) +
  plot(conditional_effects_sobol(fit_M8_varsel, "x2s", a = a)) +
  plot(conditional_effects_sobol(fit_M8_varsel, "x3s", a = a)) +
  plot(conditional_effects_sobol(fit_M8_varsel, "x4s", a = a)) +
  plot_layout(nrow = 2)
# ggsave("plots/sobol_ce_mc300_M4_varsel.jpg", width = 8, height = 5)

plot(plot_poly_degree(fit_M4_varsel, p = 10, M = 4) +
       scale_x_continuous(breaks = 1:10, limits = c(1, 10.3))) +
  plot(plot_sobol_coef_degree(fit_M4_varsel, yintercept = sobol_sd_num(a)))
ggsave("plots/sobol_poly_mc300_M4_varsel.jpg", width = 8, height = 3)

# table for the most important polynomials and the degrees
df_M4_varsel <- df_poly_degree(fit_M4_varsel, p = 10, M = 4, digits = 3)
df_M4_varsel %>%
  filter(row_number() %in% 1:14) %>%
  kbl(format= "latex", align=c("l", "l" ,"r", "r", "r", "r")) %>%
  kable_classic(full_width = F, html_font = "helvetica")


# plots over training points
sobol_summary <- readRDS("sobol_summary.rds") %>%
  mutate(
    bias_mean = PCE_mean_estimate - sobol_mean,
    bias_sd = PCE_sd_estimate - sobol_sd,
    L2AE_mean = abs(PCE_mean_estimate - sobol_mean),
    L2AE_sd = abs(PCE_sd_estimate - sobol_sd),
    type = factor(type, levels = c("full", "b_sel", "phi_sel", "varsel"),
                  labels = c("Full", "Sel-Sobol", "Sel-D2", "Sel-Projpred")),
    time_per_chain_minutes = time_total / nchains / 60,
    poly = factor(ifelse(M == 4, "N = 4, d = 10", "N = 8, d = 6")) %>%
      relevel("N = 8, d = 6")
  ) %>%
  rename(Model = type) %>%
  # D2 selection is almost equivalent to Sobol selection
  # no need to introduce and discuss this method
  filter(Model != "Sel-D2")

xminor_breaks <- sort(unique(sobol_summary$mc))
xbreaks <- setdiff(xminor_breaks, c(1001, 3003))
xlims <- c(80, 9000)

# Bias of mean
gg_sobol_bias_mean_summary <- sobol_summary %>%
  ggplot(aes(mc, bias_mean, color = Model, linetype = poly)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression(mu(hat(y))-mu(y))) +
  xlab("T") +
  scale_y_continuous(
    trans=make_lal_trans(
      'trexp',
      threshold=0.00001,
      exponent=10,
      force_thresholds_in_breaks = TRUE
    ),
    breaks = c(-0.1, -0.01, -0.001, -0.0001, 0, 0.0001),
    limits = c(-0.12, 0.0005)
  ) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     minor_breaks = xminor_breaks,
                     limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dotted")

gg_sobol_bias_mean_summary

# absolute bias of mean
gg_sobol_L2AE_mean_summary <- sobol_summary %>%
  ggplot(aes(mc, L2AE_mean, color = Model, linetype = poly)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression("|"~mu(hat(y))-mu(y)~"|")) +
  xlab("T") +
  scale_y_continuous(trans="log10", breaks = c(0.0001, 0.001, 0.01, 0.1),
                     limits = c(0.0001, 0.1)) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     minor_breaks = xminor_breaks,
                     limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dotted")

gg_sobol_L2AE_mean_summary

# bias of SD
gg_sobol_bias_sd_summary <- sobol_summary %>%
  ggplot(aes(mc, bias_sd, color = Model, linetype = poly)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression(sigma(hat(y))-sigma(y))) +
  xlab("T") +
  scale_y_continuous(
    trans=make_lal_trans(
      'trexp',
      threshold=0.0001,
      exponent=10,
      force_thresholds_in_breaks = TRUE
    ),
    breaks = c(-0.1, -0.01, -0.001, 0),
    limits = c(-0.1, 0.00)
  ) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     minor_breaks = xminor_breaks,
                     limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dotted")

gg_sobol_bias_sd_summary

# absolute bias of SD
gg_sobol_L2AE_sd_summary <- sobol_summary %>%
  ggplot(aes(mc, L2AE_sd, color = Model, linetype = poly)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression("|"~sigma(hat(y))-sigma(y)~"|")) +
  xlab("T") +
  scale_y_continuous(trans='log10', breaks = c(0.001, 0.01, 0.1),
                     limits = c(0.001, 0.1)) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     minor_breaks = xminor_breaks,
                     limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dotted")

gg_sobol_L2AE_sd_summary

# out of sample RMSE
gg_sobol_RMSE_oos_summary <- sobol_summary %>%
  ggplot(aes(mc, RMSE_oos_mean, color = Model, linetype = poly)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression("Out-of-sample RMSE")) +
  xlab("T") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     minor_breaks = xminor_breaks,
                     limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dotted")

gg_sobol_RMSE_oos_summary


# estimation time
gg_sobol_time_summary <- sobol_summary %>%
  ggplot(aes(mc, time_per_chain_minutes, color = Model, linetype = poly)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab("Estimation time per chain (minutes)") +
  xlab("T") +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     minor_breaks = xminor_breaks,
                     limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dotted")

gg_sobol_time_summary

# combine the four plots
(gg_sobol_bias_mean_summary + gg_sobol_bias_sd_summary) /
  (gg_sobol_RMSE_oos_summary + gg_sobol_time_summary) /
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom") 

ggsave("plots/sobol_summaries.jpg", width = 10, height = 6)

