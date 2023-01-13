library(tidyverse)
library(patchwork)
library(brms)
library(orthopolynom)
library(SobolSequence)
library(kableExtra)
source("PCE_helpers.R")
source("ishigami_helpers.R")
source("plot_helpers.R")
theme_set(theme_bw())

p <- 10
M <- 3
a <- 7
b <- 0.1

fit_ref <- brm(file = "models/fit_ishigami_mc100")
var_ref <- variables(fit_ref)[c(1, 2, 38, 66)+1]
gg_pairs <- plot(pairs(fit_ref, variable = var_ref))
plot(gg_pairs)
ggsave("plots/ishigami_pairs_small.jpg", gg_pairs, width = 10, height = 7)


mcmc_plot(fit_ref, variable = var_ref, type = "trace",
          facet_args = list(nrow = 1)) + 
  xlab("Post-warmup iteration")
ggsave("plots/ishigami_trace_small.jpg", width = 10, height = 3)


# plots for the projpred (winning) model for 100 training points
fit_varsel <- brm(file = "models/fit_ishigami_mc100_varsel")

ce1 <- conditional_effects_ishigami(fit_varsel, "x1s", a = a, b = b) +
  labs(x = expression(omega[1]~"(scaled)"))
ce2 <- conditional_effects_ishigami(fit_varsel, "x2s", a = a, b = b) +
  labs(x = expression(omega[2]~"(scaled)"))
ce3 <- conditional_effects_ishigami(fit_varsel, "x3s", a = a, b = b) +
  labs(x = expression(omega[3]~"(scaled)"))
ce4 <- conditional_effects_ishigami(fit_varsel, "x1s:x3s", a = a, b = b) +
  scale_color_viridis_d() + 
  scale_fill_viridis_d() + 
  labs(x = expression(omega[1]~"(scaled)"),
       color = expression(omega[3]~"(scaled)"), 
       fill = expression(omega[3]~"(scaled)"))
ce5 <- conditional_effects_ishigami(fit_varsel, "x3s:x1s", a = a, b = b) +
  scale_color_viridis_d() + 
  scale_fill_viridis_d() + 
  labs(x = expression(omega[3]~"(scaled)"),
       color = expression(omega[1]~"(scaled)"), 
       fill = expression(omega[1]~"(scaled)"))

(ce1 + ce2 + ce3) / (ce4 + ce5)
ggsave("plots/ishigami_ce_small_varsel.jpg", width = 8, height = 5)


plot(plot_poly_degree(fit_varsel, p = p, M = M)) +
  plot(plot_sobol_coef_degree(fit_varsel, yintercept = ishigami_sd(a, b)))

ggsave("plots/ishigami_poly_small_varsel.jpg", width = 8, height = 3)

# table for the most important polynomials and the orders
df_varsel <- df_poly_degree(fit_varsel, p = p, M = M)
df_varsel %>%
  filter(row_number() %in% 1:10) %>%
  kbl(format= "latex", align=c("l", "l" ,"r", "r", "r", "r")) %>%
  kable_classic(full_width = F, html_font = "helvetica")


# plots for the maximal-mean model
# fit_b_sel <- brm(file = "models/fit_ishigami_mc100_b_sel")
# 
# plot(plot_poly_degree(fit_b_sel, p = p, M = M)) +
#   plot(plot_sobol_coef_degree(fit_b_sel, yintercept = ishigami_sd(a, b)))
# 
# ggsave("plots/ishigami_poly_small_b_sel.jpg", width = 8, height = 3)


# plots over training points
ishigami_summary <- readRDS("ishigami_summary.rds") %>%
  mutate(
    bias_mean = PCE_mean_estimate - ishigami_mean,
    bias_sd = PCE_sd_estimate - ishigami_sd,
    L2AE_mean = abs(PCE_mean_estimate - ishigami_mean),
    L2AE_sd = abs(PCE_sd_estimate - ishigami_sd),
    type = factor(type, levels = c("full", "b_sel", "phi_sel", "varsel"),
                  labels = c("Full", "Sel-Sobol", "Sel-D2", "Sel-Projpred")),
    time_per_chain_minutes = time_total / nchains / 60
  ) %>%
  rename(Model = type) %>%
  # D2 selection is almost equivalent to Sobol selection
  # no need to introduce and discuss this method
  filter(Model != "Sel-D2")


xbreaks <- sort(unique(ishigami_summary$mc))
xlims <- c(8, 850)


# bias of Mean estimate
gg_ishigami_bias_mean_summary <- ishigami_summary %>%
  ggplot(aes(mc, bias_mean, color = Model)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression(mu(hat(y))-mu(y))) +
  xlab("T") +
  scale_y_continuous(
    trans=make_lal_trans(
      'trexp',
      threshold=0.000001,
      exponent=10,
      force_thresholds_in_breaks = TRUE
    ),
    breaks = c(-0.1, -0.01, -0.001, -0.0001, 0.0001, -0.00001, 
               0, 0.00001, 0.0001, 0.001, 0.01, 0.1),
    limits = c(-0.01, 1)
  ) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_ishigami_bias_mean_summary


# absolute bias of Mean estimate
gg_ishigami_L2AE_mean_summary <- ishigami_summary %>%
  ggplot(aes(mc, L2AE_mean, color = Model)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression("|"~mu(hat(y))-mu(y)~"|")) +
  xlab("T") +
  scale_y_continuous(trans='log10', breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1),
                      limits = c(0.00001, 1)) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                      limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_ishigami_L2AE_mean_summary

# bias of SD estimate
gg_ishigami_bias_sd_summary <- ishigami_summary %>%
  ggplot(aes(mc, bias_sd, color = Model)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression(sigma(hat(y))-sigma(y))) +
  xlab("T") +
  scale_y_continuous(
    trans=make_lal_trans(
      'trexp',
      threshold=0.001,
      exponent=10,
      force_thresholds_in_breaks = TRUE
    ),
    breaks = c(-1, -0.1, -0.01, 0, 0.01, 0.1, 1),
    limits = c(-1.5, 1.5)
  ) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_ishigami_bias_sd_summary

# absolute bias of SD estimate
gg_ishigami_L2AE_sd_summary <- ishigami_summary %>%
  ggplot(aes(mc, L2AE_sd, color = Model)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression("|"~sigma(hat(y))-sigma(y)~"|")) +
  xlab("T") +
  scale_y_continuous(trans='log10', breaks = c(0.001, 0.01, 0.1),
                    limits = c(0.01, 1.2)) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_ishigami_L2AE_sd_summary

# out of sample RMSE
gg_ishigami_RMSE_oos_summary <- ishigami_summary %>%
  ggplot(aes(mc, RMSE_oos_mean, color = Model)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression("Out-of-sample RMSE")) +
  xlab("T") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_ishigami_RMSE_oos_summary
ggsave("plots/ishigami_RMSE_oos_summary.jpg", width = 7, height = 3)


# estimation time
gg_ishigami_time_summary <- ishigami_summary %>%
  ggplot(aes(mc, time_per_chain_minutes, color = Model)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab("Estimation time per chain (minutes)") +
  xlab("T") +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                      limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_ishigami_time_summary

# combine the four plots
(gg_ishigami_bias_mean_summary + gg_ishigami_bias_sd_summary) /
  (gg_ishigami_RMSE_oos_summary + gg_ishigami_time_summary) /
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom") 

ggsave("plots/ishigami_summaries.jpg", width = 10, height = 6)

