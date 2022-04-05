library(tidyverse)
library(patchwork)
library(brms)
library(orthopolynom)
library(SobolSequence)
library(kableExtra)
source("PCE_helpers.R")
source("signum_helpers.R")
source("plot_helpers.R")
theme_set(theme_bw())

# plot conditional_effects
fit_SoSeq_11_R2D2 <- brm(file = "models/fit_signum_SoSeq_11_R2D2")
fit_SoSeq_11_flat <- brm(file = "models/fit_signum_SoSeq_11_flat")
fit_GaInt_11_R2D2 <- brm(file = "models/fit_signum_GaInt_11_R2D2")
fit_GaInt_11_flat <- brm(file = "models/fit_signum_GaInt_11_flat")

newdata <- data.frame(x = seq(-1, 1, length.out = 100))
thres <- 5

data_ce_plot <- bind_rows(
  data.frame(
    Estimate = predict_standard_PCE(newdata$x, p = 10, strategy = "SoSeq"),
    strategy = "SoSeq", model = "Standard"
  ),
  data.frame(
    Estimate = predict_standard_PCE(newdata$x, p = 10, strategy = "GaInt"),
    strategy = "GaInt", model = "Standard"
  ),
  cbind(
    as.data.frame(fitted(fit_SoSeq_11_R2D2, newdata = newdata)),
    strategy = "SoSeq", model = "Bayesian-R2D2"
  ),
  cbind(
    as.data.frame(fitted(fit_GaInt_11_R2D2, newdata = newdata)),
    strategy = "GaInt", model = "Bayesian-R2D2"
  ),
  cbind(
    as.data.frame(fitted(fit_SoSeq_11_flat, newdata = newdata)),
    strategy = "SoSeq", model = "Bayesian-Flat"
  ),
  cbind(
    as.data.frame(fitted(fit_GaInt_11_flat, newdata = newdata)),
    strategy = "GaInt", model = "Bayesian-Flat"
  )
) %>%
  mutate(
    x = rep(newdata$x, 6),
    truth = sign(x),
    strategy = factor(strategy, levels = c("GaInt", "SoSeq"),
                      labels = c("Gaussian-Integration", "Sobol-Sequence")),
    model = factor(
      model, 
      levels = c("Standard", "Bayesian-Flat", "Bayesian-R2D2")
    )
  ) %>%
  group_by(model, strategy) %>%
  mutate(
    index_first_large = which(
      abs(Estimate) > thres | 
        !is.na(Q2.5) & (abs(Q2.5) > thres | abs(Q97.5) > thres)
    )[1],
    index_first_large = ifelse(is.na(index_first_large), Inf,
                               index_first_large),
    Estimate = ifelse(row_number() >= index_first_large, 
                        NA, Estimate)
  ) %>%
  ungroup() %>%
  filter(!is.na(Estimate)) %>%
  rename(
    `Training-Strategy` = strategy,
    `PCE-Model` = model
  )

data_training <- bind_rows(
  tibble(
    x = fit_SoSeq_11_R2D2$data$x,
    `Training-Strategy` = "Sobol-Sequence"
  ),
  tibble(
    x = fit_GaInt_11_R2D2$data$x,
    `Training-Strategy` = "Gaussian-Integration"
  )
)

ggplot(data_ce_plot, aes(x, Estimate, ymin = Q2.5, ymax = Q97.5,
                         color = `PCE-Model`, fill = `PCE-Model`)) +
  geom_smooth(stat = "identity") +
  geom_line(aes(y = truth), size = 1.5, linetype = "dashed", color = "black") + 
  geom_rug(aes(x), inherit.aes = FALSE, data = data_training, size = 1) +
  facet_wrap("`Training-Strategy`") +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  labs(y = "y", x = expression(omega)) +
  theme(legend.position = "bottom") 

ggsave("plots/signum_ce.jpg", width = 10, height = 4)

# (conditional_effects_signum(fit_SoSeq_11_R2D2) + 
#     ggtitle("Sobol-Sequence: Bayesian-RD2D2") +
#   conditional_effects_signum(fit_SoSeq_11_flat) +
#     ggtitle("Sobol-Sequence: Bayesian-Flat")) /
#   (conditional_effects_signum(fit_GaInt_11_R2D2) + 
#      ggtitle("Gaussian Integration: Bayesian-RD2D2") +
#      conditional_effects_signum(fit_GaInt_11_flat) +
#      ggtitle("Gaussian Integration: Bayesian-Flat"))


# reconstruct standard PCE results
set.seed(5235426)
N_test <- 500
test <- data.frame(x = runif(N_test, -1, 1)) %>%
  mutate(y = sign(x))

mc <- 1:11
strategy <- c("GaInt", "SoSeq")
design <- expand.grid(mc = mc, strategy = strategy) %>%
  mutate(
    signum_mean = signum_mean(),
    signum_sd = signum_sd()
  ) %>% 
  rowwise() %>%
  mutate(
    RMSE_oos_mean = 
      rmse_mean2(test$y, predict_standard_PCE(test$x, p = mc - 1, strategy = strategy))
  )

signum_summary_std_GaInt <- design %>% 
  filter(strategy == "GaInt") %>%
  bind_cols(read.table("data/signum_function/GaussIntegration/OutputMeanStd_GaussIntegration.txt")) %>%
  rename(PCE_mean_estimate = V1, PCE_sd_estimate = V2) %>%
  mutate(model = "Standard", time_total = 0, nchains = 1)

signum_summary_std_SoSeq <- design %>% 
  filter(strategy == "SoSeq") %>%
  bind_cols(read.table("data/signum_function/SobolSequence/OutputMeanStd_SobolSequence.txt")) %>%
  rename(PCE_mean_estimate = V1, PCE_sd_estimate = V2)  %>%
  mutate(model = "Standard", time_total = 0, nchains = 1)

# truncate error to make plots easier to read
trunc <- 1e-05

signum_summary <- readRDS("signum_summary.rds") %>%
  mutate(model = paste0("Bayesian-", prior)) %>%
  bind_rows(signum_summary_std_GaInt, signum_summary_std_SoSeq) %>%
  mutate(
    bias_mean = PCE_mean_estimate - signum_mean,
    bias_sd = PCE_sd_estimate - signum_sd,
    L2AE_mean = abs(PCE_mean_estimate - signum_mean),
    L2AE_sd = abs(PCE_sd_estimate - signum_sd),
    L2AE_mean_trunc = ifelse(L2AE_mean < trunc, trunc, L2AE_mean),
    L2AE_sd_trunc = ifelse(L2AE_sd < trunc, trunc, L2AE_sd),
    strategy = factor(strategy, levels = c("GaInt", "SoSeq"),
                  labels = c("Gaussian-Integration", "Sobol-Sequence")),
    model = factor(
      model, 
      levels = c("Standard", "Bayesian-flat", "Bayesian-R2D2"),
      labels = c("Standard", "Bayesian-Flat", "Bayesian-R2D2")
    ),
    time_per_chain_minutes = time_total / nchains / 60
  ) %>%
  rename(
    `Training-Strategy` = strategy,
    `PCE-Model` = model
  ) %>%
  # N=1 training point is a bit pointless
  filter(mc > 1)

# plots over training points
xbreaks <- 2:11
xlims <- c(2, 11)

# bias of Mean estimate
# gg_signum_bias_mean_summary <- signum_summary %>%
#   ggplot(aes(mc, bias_mean, color = `PCE-Model`, linetype = `Training-Strategy`)) +
#   geom_line(size = 0.8) +
#   scale_color_viridis_d() +
#   ylab(expression(mu(hat(y))-mu(y))) +
#   xlab("T") +
#   scale_y_continuous(
#     trans=make_lal_trans(
#       'trexp',
#       threshold=0.000001,
#       exponent=10,
#       force_thresholds_in_breaks = TRUE
#     ),
#     breaks = c(-0.1, -0.01, -0.001, -0.0001, 0.0001, -0.00001,
#                0, 0.00001, 0.0001, 0.001, 0.01, 0.1),
#     limits = c(-0.01, 1)
#   ) +
#   scale_x_continuous(breaks = xbreaks, limits = xlims) +
#   geom_hline(yintercept = 0, linetype = "dashed")
# 
# gg_signum_bias_mean_summary


# absolute bias of Mean estimate
gg_signum_L2AE_mean_summary <- signum_summary %>%
  ggplot(aes(mc, L2AE_mean_trunc, color = `PCE-Model`, linetype = `Training-Strategy`)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression("|"~mu(hat(y))-mu(y)~"|")) +
  xlab("T") +
  scale_y_continuous(
    trans='log10', 
    breaks = c(1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
    limits = c(1e-5, 1e2)
  ) +
  scale_x_continuous(breaks = xbreaks, limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_signum_L2AE_mean_summary

# bias of SD estimate
# gg_signum_bias_sd_summary <- signum_summary %>%
#   ggplot(aes(mc, bias_sd, color = `PCE-Model`, linetype = `Training-Strategy`)) +
#   geom_line(size = 1) +
#   scale_color_viridis_d() +
#   ylab(expression(sigma(hat(y))-sigma(y))) +
#   xlab("T") +
#   scale_y_continuous(
#     trans=make_lal_trans(
#       'trexp',
#       threshold=0.001,
#       exponent=10,
#       force_thresholds_in_breaks = TRUE
#     ),
#     breaks = c(-1, -0.1, -0.01, 0, 0.01, 0.1, 1),
#     limits = c(-1.5, 1.5)
#   ) +
#   scale_x_continuous(breaks = xbreaks, limits = xlims) +
#   geom_hline(yintercept = 0, linetype = "dashed")
# 
# gg_signum_bias_sd_summary

# absolute bias of SD estimate
gg_signum_L2AE_sd_summary <- signum_summary %>%
  ggplot(aes(mc, L2AE_sd_trunc, color = `PCE-Model`, linetype = `Training-Strategy`)) +
  geom_line(size = 1) +
  scale_color_viridis_d() +
  ylab(expression("|"~sigma(hat(y))-sigma(y)~"|")) +
  xlab("T") +
  scale_y_continuous(
    trans='log10', 
    breaks = c(1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
    limits = c(1e-5, 200)
  ) +
  scale_x_continuous(breaks = xbreaks, limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  guides(color = "none", linetype = "none")

gg_signum_L2AE_sd_summary

# out of sample RMSE
gg_signum_RMSE_oos_summary <- signum_summary %>%
  ggplot(aes(mc, RMSE_oos_mean, color = `PCE-Model`, linetype = `Training-Strategy`)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression("Out-of-sample RMSE")) +
  xlab("T") +
  scale_y_continuous(trans='log10',
                     breaks = c(1e2, 1e1, 1e0, 1e-1),
                     limits = c(1e-1, 1e2+200)) +
  scale_x_continuous(breaks = xbreaks, limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_signum_RMSE_oos_summary
ggsave("plots/signum_RMSE_oos_summary.jpg", width = 7, height = 3)

# estimation time
gg_signum_time_summary <- signum_summary %>%
  ggplot(aes(mc, time_per_chain_minutes, color = `PCE-Model`, linetype = `Training-Strategy`)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab("Estimation time per chain (minutes)") +
  xlab("T") +
  scale_x_continuous(breaks = xbreaks, limits = xlims) +
  geom_hline(yintercept = 0, linetype = "dashed")

gg_signum_time_summary

# combine the four plots
(gg_signum_L2AE_mean_summary + gg_signum_L2AE_sd_summary) /
  (gg_signum_RMSE_oos_summary + gg_signum_time_summary) /
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom") 

ggsave("plots/signum_summaries.jpg", width = 10, height = 6)


# plot illustrative example
data_illus <- data.frame(x = seq(-1, 1, length.out = 1000))
data_illus$truth <- sign(data_illus$x)
for (p in c(10, 15, 25)) {
  path_train <- paste0("data/signum_function/gibbs_illustration/pce_No", p, ".dat")
  path_poly <- paste0("data/signum_function/gibbs_illustration/pol_cfs_No", p, "_new.dat")
  path_roots <- paste0("data/signum_function/gibbs_illustration/roots_No", p, ".dat")
  data_train <- read.table(path_roots)
  data_illus[[paste0("pred_", p)]] <- 
    predict_standard_PCE_lm(data_train$V1, p = p, new_x = data_illus$x, path_poly = path_poly)
}
data_illus <- data_illus %>%
  gather(key = "p", value = "pred", pred_10:pred_25) %>%
  mutate(p = sub("pred_", "", p))

ggplot(data_illus, aes(x, pred, color = p)) +
  geom_line(stat = "identity") +
  geom_line(aes(y = truth), size = 1.2, linetype = "dashed", color = "black") + 
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  labs(y = "y", x = expression(omega), color = "d") +
  theme(legend.position = "bottom") 

ggsave("plots/signum_illus.jpg", width = 7, height = 2.8)
