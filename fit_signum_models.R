library(dplyr)
library(tidyr)
# library(ggplot2)
library(foreach)
library(doParallel)
library(orthopolynom)
# library(SobolSequence)
# options(width = 100)
# theme_set(theme_default())
# cmdstanr::set_cmdstan_path("/scratch/work/burknep1/cmdstanr/cmdstan-2.26.1")
source("PCE_helpers.R")
source("signum_helpers.R")

set.seed(5235426)
N_test <- 500
test <- data.frame(x = runif(N_test, -1, 1)) %>%
  mutate(y = sign(x))

mc <- 1:11
strategy <- c("GaInt", "SoSeq")
prior <- c("flat", "R2D2")

design <- expand.grid(mc = mc, strategy = strategy, prior = prior)

data_list <- vector("list", nrow(design))
for (i in seq_along(data_list)) {
  mc_i <- design$mc[i]
  strategy_i <- design$strategy[i] 
  prior <- design$prior[i]
  if (strategy_i == "SoSeq") {
    path_train <- "data/signum_function/SobolSequence/TrainingPoints_SobolSequence.txt"
    training <- read.table(path_train) %>%
      rename(x = V1) %>%
      filter(row_number() %in% 1:mc_i)
  } else if (strategy_i == "GaInt") {
    path_train <- 
      paste0("data/signum_function/GaussIntegration/TrainingPoints_GaussIntegration_", mc_i, ".txt")
    training <- read.table(path_train) %>%
      rename(x = V1)
  }
  training$y <- sign(training$x)
  
  # poly <- with(training, PCE_CO2(x, p = mc_i))
  # colnames(poly) <- paste0("P", seq_len(ncol(poly)))
  # training <- bind_cols(training, as.data.frame(poly))
  
  data_list[[i]] <- training
}

# fit the models in parallel
# cl <- makeForkCluster(detectCores())
# registerDoParallel(cl)

exports <- NULL

I <- seq_len(nrow(design))

results <- foreach(i = I, .export = exports, .combine = bind_rows) %dopar% {
  library(dplyr)
  library(brms)
  # library(projpred)
  library(SobolSequence)
  source("PCE_helpers.R")
  source("signum_helpers.R")
  env <- new.env()
  
  mc_i <- design$mc[i]
  strategy_i <- design$strategy[i]
  prior_i <- design$prior[i]
  training <- data_list[[i]]
  
  suffix <- ifelse(prior_i == "flat", "_flat", "_R2D2")
  bprior <- prior(normal(0, 0.05), class = "sigma")
  if (mc_i > 1) {
    if (prior_i == "R2D2") {
      bprior <- bprior + prior(R2D2(mean_R2 = 0.9, prec_R2 = 10)) 
    }
    bform <- y ~ PCE_signum(x, p = mc_i - 1)
  } else {
    bform <- y ~ 1
  }
  
  fit <- brm(
    formula = bform,
    data = training,
    prior = bprior,
    data2 = list(mc_i = mc_i),
    chains = 2, cores = 2,
    iter = 11000, warmup = 1000,
    control = list(adapt_delta = 0.999),
    backend = "cmdstanr",
    file = paste0("models/fit_signum_", strategy_i, "_", mc_i, suffix)
  )
  
  out <- design[i, ]
  # MCMC info
  out$nchains <- fit$fit@sim$chains
  out$niterations <- fit$fit@sim$iter
  out$nwarmup <- fit$fit@sim$warmup
  
  # sampling times
  time_matrix <- rstan::get_elapsed_time(fit$fit)
  out$time_warmup <- sum(time_matrix[, "warmup"])
  out$time_sample <- sum(time_matrix[, "sample"])
  out$time_total <- sum(time_matrix)
  
  pce_mean <- PCE_mean(fit)
  out$signum_mean <- signum_mean()
  out$PCE_mean_estimate <- pce_mean[1,1]
  out$PCE_mean_se <- pce_mean[1,2]
  out$PCE_mean_lower <- pce_mean[1,3]
  out$PCE_mean_upper <- pce_mean[1,4]
  
  pce_sd <- PCE_sd(fit)
  out$signum_sd <- signum_sd()
  out$PCE_sd_estimate <- pce_sd[1,1]
  out$PCE_sd_se <- pce_sd[1,2]
  out$PCE_sd_lower <- pce_sd[1,3]
  out$PCE_sd_upper <- pce_sd[1,4]
  
  # mean after centering PCE predictors (not sensible)
  # pce_mean2 <- PCE_mean2(fit)
  
  # yname <- fit$formula$resp
  # yrep_ins <- posterior_epred(fit)
  # rmse_ins <- posterior_summary(rmse(test[[yname]], yrep_ins))
  # out$RMSE_ins_mean = rmse_ins[1, 1]
  # out$RMSE_ins_se = rmse_ins[1, 2]
  # out$RMSE_ins_lower = rmse_ins[1, 3]
  # out$RMSE_ins_upper = rmse_ins[1, 4]
  # out$RMSE_ins_mean = rmse_mean(fit$data[[yname]], yrep_ins)
  
  yname <- fit$formula$resp
  yrep_oos <- posterior_epred(fit, test)
  # rmse_oos <- posterior_summary(rmse(test[[yname]], yrep_oos))
  # out$RMSE_oos_mean = rmse_oos[1, 1]
  # out$RMSE_oos_se = rmse_oos[1, 2]
  # out$RMSE_oos_lower = rmse_oos[1, 3]
  # out$RMSE_oos_upper = rmse_oos[1, 4]
  out$RMSE_oos_mean = rmse_mean(test[[yname]], yrep_oos)
  
  out
}

# stopCluster(cl)

saveRDS(results, "signum_summary.rds")

