library(ggplot2)
library(patchwork)
library(dplyr)
source("plot_helpers.R")
theme_set(theme_bw())

CO2_summary <- readRDS("CO2_summary.rds") %>%
  mutate(Prior = "R2D2(0.5, 2)")

CO2_summary_high_R2 <- readRDS("CO2_summary_high_R2.rds") %>%
  mutate(Prior = "R2D2(0.9, 10)")

CO2_summary <- CO2_summary %>%
  bind_rows(CO2_summary_high_R2) %>%
  mutate(
    type = factor(type, levels = c("full", "b_sel", "b_sel_50", "varsel", "varsel_50"),
                  labels = c("Full", "Sel-Sobol-25", "Sel-Sobol-50", 
                             "Sel-Projpred-25", "Sel-Projpred-50")),
    diff_mean = PCE_mean_estimate - CO2_mean,
    diff_mean_lower = PCE_mean_lower - CO2_mean,
    diff_mean_upper = PCE_mean_upper - CO2_mean,
    diff_sd = PCE_sd_estimate - CO2_sd,
    diff_sd_lower = PCE_sd_lower - CO2_sd,
    diff_sd_upper = PCE_sd_upper - CO2_sd,
    time_total_minutes = time_total / 60,
    mc_points = (d + 1)^3
  ) %>%
  rename(Model = type)

CO2_summary_mean = CO2_summary %>% 
  group_by(Model, mc_points, Prior) %>%
  summarise(
    L2AE_mean = sqrt(sum(diff_mean^2)) / n(),
    L2AE_sd = sqrt(sum(diff_sd^2)) / n(),
    MAE_mean = mean(abs(diff_mean)),
    MAE_sd = mean(abs(diff_sd)),
    RMSE_ins_mean = mean(RMSE_ins_mean),
    RMSE_oos_mean = mean(RMSE_oos_mean),
    time_mean = mean(time_total_minutes),
    time_lower = quantile(time_total_minutes, 0.05),
    time_upper = quantile(time_total_minutes, 0.95)
  )


xbreaks <- sort(unique(CO2_summary$mc))
xlims <- c(25, 1500)


# mean over locations
CO2_summary %>%
  ggplot(aes(coord, PCE_mean_estimate, color = Model, fill = Model,
             ymin = PCE_mean_lower, ymax = PCE_mean_upper,
             linetype = Prior)) +
  facet_wrap("mc_points", ncol = 1, scales = "free_y") +
  geom_smooth(size = 0.8, stat = "identity") +
  geom_line(aes(coord, CO2_mean), inherit.aes = FALSE, size = 0.5) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  ylab(expression(mu(hat(y)))) +
  xlab("Location")
ggsave("plots/CO2_mean.jpg", width = 8, height = 6)


# sd over locations
CO2_summary %>%
  ggplot(aes(coord, PCE_sd_estimate, color = Model, fill = Model,
             ymin = PCE_sd_lower, ymax = PCE_sd_upper,
             linetype = Prior)) +
  facet_wrap("mc_points", ncol = 1, scales = "free_y") +
  geom_smooth(size = 0.8, stat = "identity") +
  geom_line(aes(coord, CO2_sd), inherit.aes = FALSE, size = 0.5) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  ylab(expression(sigma(hat(y)))) +
  xlab("Location")
ggsave("plots/CO2_sd.jpg", width = 8, height = 6)


# mean bias over locations
CO2_summary %>%
  ggplot(aes(coord, diff_mean, color = Model, fill = Model,
             ymin = diff_mean_lower, ymax = diff_mean_upper,
             linetype = Prior)) +
  facet_wrap("mc_points", ncol = 1, scales = "free_y") +
  geom_smooth(size = 0.8, stat = "identity") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  ylab(expression(mu(hat(y))-mu(y))) +
  xlab("Location")
ggsave("plots/CO2_diff_mean.jpg", width = 8, height = 6)

# sd bias over locations
CO2_summary %>%
  ggplot(aes(coord, diff_sd, color = Model, fill = Model,
             ymin = diff_sd_lower, ymax = diff_sd_upper,
             linetype = Prior)) +
  facet_wrap("mc_points", ncol = 1, scales = "free_y") +
  geom_smooth(size = 0.8, stat = "identity") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  ylab(expression(sigma(hat(y))-sigma(y))) +
  xlab("Location")
ggsave("plots/CO2_diff_sd.jpg", width = 8, height = 6)

# OOS RMSE over locations
CO2_summary %>%
  ggplot(aes(coord, RMSE_oos_mean, color = Model, fill = Model,
             ymin = RMSE_oos_lower, ymax = RMSE_oos_upper,
             linetype = Prior)) +
  facet_wrap("mc_points", ncol = 1, scales = "free_y") +
  geom_smooth(size = 0.8, stat = "identity") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  ylab("Out-of-sample RMSE") +
  xlab("Location")
ggsave("plots/CO2_RMSE_oos.jpg", width = 8, height = 6)


# plots over number of training points
gg_CO2_L2AE_mean_summary <- CO2_summary_mean %>%
  ggplot(aes(mc_points, L2AE_mean, color = Model, linetype = Prior)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression("||"~mu(hat(y))-mu(y)~"|| / L")) +
  xlab("T") +
  scale_y_continuous(trans='log10', breaks = c(0.0001, 0.001, 0.01, 0.1),
                     limits = c(0.0001, 0.1)) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     limits = xlims)

gg_CO2_L2AE_mean_summary
#ggsave("plots/CO2_L2AE_mean_summary.jpg", width = 10, height = 4)

gg_CO2_L2AE_sd_summary <- CO2_summary_mean %>%
  ggplot(aes(mc_points, L2AE_sd, color = Model, linetype = Prior)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression("||"~sigma(hat(y))-sigma(y)~"|| / L")) +
  xlab("T") +
  scale_y_continuous(trans='log10', breaks = c(0.0001, 0.001, 0.01, 0.1),
                     limits = c(0.0001, 0.1)) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     limits = xlims)

gg_CO2_L2AE_sd_summary
#ggsave("plots/CO2_L2AE_sd_summary.jpg", width = 10, height = 4)

gg_CO2_RMSE_oos_summary <- CO2_summary_mean %>%
  ggplot(aes(mc_points, RMSE_oos_mean, color = Model, linetype = Prior)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab(expression("Out-of-sample"~bar(RMSE))) +
  xlab("T") +
  scale_y_continuous(trans='log10', limits = c(0.03, 0.5)) +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     limits = xlims)

gg_CO2_RMSE_oos_summary
#ggsave("plots/CO2_RMSE_oos_summary.jpg", width = 10, height = 4)


# estimation time
gg_CO2_time_summary <- CO2_summary_mean %>%
  ggplot(aes(mc_points, time_mean, color = Model, linetype = Prior)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab("Estimation time (minutes)") +
  xlab("T") +
  scale_x_continuous(trans='log10', breaks = xbreaks,
                     limits = xlims)

gg_CO2_time_summary
#ggsave("plots/CO2_time_summary.jpg", width = 10, height = 4)

# combine the four plots
(gg_CO2_L2AE_mean_summary + gg_CO2_L2AE_sd_summary) /
  (gg_CO2_RMSE_oos_summary + gg_CO2_time_summary) /
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom") 

ggsave("plots/CO2_summaries.jpg", width = 10, height = 6)

