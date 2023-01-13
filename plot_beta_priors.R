library(tidyverse)
library(patchwork)
library(latex2exp)
theme_set(theme_bw())

dbeta2 <- function(x, mu, nu, ...) {
  shape1 <- mu * nu
  shape2 <- (1-mu) * nu
  dbeta(x, shape1, shape2, ...)
}

dbeta2_prime <- function(x, mu, nu) {
  shape1 <- mu * nu
  shape2 <- (1-mu) * nu
  x^(shape1-1) * (1 + x)^(-shape1 - shape2) / beta(shape1, shape2)
}

# beta prior
df_beta <- data.frame(x = seq(0.001, 0.999, 0.001)) %>%
  mutate(
    `y(0.5, 2)` = dbeta2(x, 0.5, 2),
    `y(0.4, 5)` = dbeta2(x, 0.4, 5),
    `y(0.66, 15)` = dbeta2(x, 0.66, 15),
  ) %>%
  gather("Shapes", "dens", starts_with("y")) %>%
  mutate(
    Shapes = str_remove(Shapes, "^y") %>%
      factor(levels = c("(0.66, 15)", "(0.4, 5)", "(0.5, 2)"))
  )

df_beta_shade <- bind_rows(
  data.frame(x = 0, dens = 0, Shapes = unique(df_beta$Shapes)),
  df_beta,
  data.frame(x = 1, dens = 0, Shapes = unique(df_beta$Shapes))
)

gg_beta <- 
  ggplot(df_beta, aes(x, dens, fill = Shapes)) +
  geom_line(size = 1) +
  geom_polygon(data = df_beta_shade, alpha = 0.7) +
  scale_fill_viridis_d() +
  xlab(expression(R^2)) +
  ylab("Density") +
  theme(
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.line.y = element_blank(),
    legend.position = "none"
  )


# beta prime prior
df_beta_prime <- data.frame(x = seq(0.001, 5, 0.001)) %>%
  mutate(
    `y(0.5, 2)` = dbeta2_prime(x, 0.5, 2),
    `y(0.4, 5)` = dbeta2_prime(x, 0.4, 5),
    `y(0.66, 15)` = dbeta2_prime(x, 0.66, 15),
  ) %>%
  gather("Shapes", "dens", starts_with("y")) %>%
  mutate(
    Shapes = str_remove(Shapes, "^y") %>%
      factor(levels = c("(0.66, 15)", "(0.4, 5)", "(0.5, 2)"))
  )

df_beta_prime_shade <- bind_rows(
  data.frame(x = 0, dens = 0, Shapes = unique(df_beta_prime$Shapes)),
  df_beta_prime,
  data.frame(x = 5, dens = 0, Shapes = unique(df_beta_prime$Shapes))
)

gg_beta_prime <- 
  ggplot(df_beta_prime, aes(x, dens, fill = Shapes)) +
  geom_line(size = 1) +
  geom_polygon(data = df_beta_prime_shade, alpha = 0.7) +
  scale_fill_viridis_d() +
  xlab(expression(tau^2)) +
  ylab("Density") +
  labs(fill = TeX("($zeta$, $nu$):")) +
  theme(
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.line.y = element_blank(),
    legend.position = "bottom",
    legend.margin=margin(l = -7, t = -0.3, unit='cm')
  )

gg_beta + gg_beta_prime 
ggsave("plots/beta-priors.jpg", height = 2.5, width = 6)
