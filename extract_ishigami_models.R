library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(orthopolynom)
library(SobolSequence)
source("PCE_helpers.R")
source("ishigami_helpers.R")


# get training data sets
mc <- c(10, 25, 50, 100, 200, 286, 400, 800)
nterms <- c(15, 25, 25, 25, 25, 25, 25, 25)

a <- 7
b <- 0.1

types <- c("full", "b_sel", "phi_sel", "varsel")
design <- expand.grid(mc = mc, type = types)
design$nterms <- nterms[match(design$mc, mc)]

set.seed(5235426)
N_test <- 500
test <- data.frame(
  x1 = runif(N_test, -pi, pi),
  x2 = runif(N_test, -pi, pi),
  x3 = runif(N_test, -pi, pi)
) %>%
  mutate(
    y = ishigami(x1, x2, x3, a, b),
    x1s = scale_to_1(x1),
    x2s = scale_to_1(x2),
    x3s = scale_to_1(x3)
  )

pol <- with(test, PCE(x1s, x2s, x3s, p = 10)) %>% as.data.frame() 
names(pol) <- paste0("P", 1:ncol(pol))
test <- bind_cols(test, pol)

cl <- makePSOCKcluster(6)
registerDoParallel(cl)

packages <- c("dplyr", "brms", "projpred", "orthopolynom")
results <- foreach(i = seq_len(nrow(design)),
                   .combine = bind_rows,
                   .packages = packages) %dopar% {
  print(i)
  # library(dplyr)
  # library(brms)
  # library(projpred)
  # library(orthopolynom)
  source("PCE_helpers.R")
  source("ishigami_helpers.R")
  
  out <- design[i, ]
  suffix <- paste0("mc", out$mc)
  if (out$type != "full") {
   suffix <- paste0(suffix, "_", out$type)
  }
  fit <- readRDS(paste0("models/fit_ishigami_", suffix, ".rds"))
  
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
  out$ishigami_mean <- ishigami_mean(a, b)
  out$PCE_mean_estimate <- pce_mean[1,1]
  out$PCE_mean_se <- pce_mean[1,2]
  out$PCE_mean_lower <- pce_mean[1,3]
  out$PCE_mean_upper <- pce_mean[1,4]
  
  pce_sd <- PCE_sd(fit)
  out$ishigami_sd <- ishigami_sd(a, b)
  out$PCE_sd_estimate <- pce_sd[1,1]
  out$PCE_sd_se <- pce_sd[1,2]
  out$PCE_sd_lower <- pce_sd[1,3]
  out$PCE_sd_upper <- pce_sd[1,4]
  
  # mean after centering PCE predictors (not sensible)
  # pce_mean2 <- PCE_mean2(fit)
  
  yname <- fit$formula$resp
  yrep_ins <- posterior_epred(fit)
  rmse_ins <- posterior_summary(rmse(test[[yname]], yrep_ins))
  out$RMSE_ins_mean = rmse_ins[1, 1]
  out$RMSE_ins_se = rmse_ins[1, 2]
  out$RMSE_ins_lower = rmse_ins[1, 3]
  out$RMSE_ins_upper = rmse_ins[1, 4]
  # out$RMSE_ins_mean = rmse_mean(fit$data[[yname]], yrep)
  
  yname <- fit$formula$resp
  yrep_oos <- posterior_epred(fit, test)
  rmse_oos <- posterior_summary(rmse(test[[yname]], yrep_oos))
  out$RMSE_oos_mean = rmse_oos[1, 1]
  out$RMSE_oos_se = rmse_oos[1, 2]
  out$RMSE_oos_lower = rmse_oos[1, 3]
  out$RMSE_oos_upper = rmse_oos[1, 4]
  # out$RMSE_oos_mean = rmse_mean(test[[yname]], yrep)
  
  out
}

stopCluster(cl)

saveRDS(results, "ishigami_summary.rds")
