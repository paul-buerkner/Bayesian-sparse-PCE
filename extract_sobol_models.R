library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(orthopolynom)
library(SobolSequence)
source("PCE_helpers.R")
source("sobol_helpers.R")

# get training data sets
mc <- c(100, 300, 900, 1001, 2700, 3003, 8100)
nterms <- c(25, 25, 25, 25, 25, 25, 25)

a <- c(1, 2, 5, 10, 20, 50, 100, 500)

types <- c("full", "b_sel", "phi_sel", "varsel")
design <- expand.grid(mc = mc, M = c(4, 8), type = types)
design$nterms <- nterms[match(design$mc, mc)]

set.seed(5235426)
N_test <- 500
test <- data.frame(
  x1 = runif(N_test, 0, 1),
  x2 = runif(N_test, 0, 1),
  x3 = runif(N_test, 0, 1),
  x4 = runif(N_test, 0, 1),
  x5 = runif(N_test, 0, 1),
  x6 = runif(N_test, 0, 1),
  x7 = runif(N_test, 0, 1),
  x8 = runif(N_test, 0, 1)
) %>%
  mutate(
    y = sobol(x1, x2, x3, x4, x5, x6, x7, x8, a = a),
    x1s = scale_to_1(x1, 0, 1),
    x2s = scale_to_1(x2, 0, 1),
    x3s = scale_to_1(x3, 0, 1),
    x4s = scale_to_1(x4, 0, 1),
    x5s = scale_to_1(x5, 0, 1),
    x6s = scale_to_1(x6, 0, 1),
    x7s = scale_to_1(x7, 0, 1),
    x8s = scale_to_1(x8, 0, 1)
  )

pol8 <- with(test, PCE(x1s, x2s, x3s, x4s, x5s, x6s, x7s, x8s, p = 6)) %>% 
  as.data.frame() 
names(pol8) <- paste0("PE", 1:ncol(pol8))

# some older M = 8 models still use "P" rather than "PE"
pol8_2 <- pol8
names(pol8_2) <- paste0("P", 1:ncol(pol8_2))


pol4 <- with(test, PCE(x1s, x2s, x3s, x4s, p = 10)) %>% 
  as.data.frame() 
names(pol4) <- paste0("PF", 1:ncol(pol4))

test <- bind_cols(test, pol4, pol8, pol8_2)

wd <- getwd()

# cl <- makeForkCluster(4)
# registerDoParallel(cl)

packages <- c("dplyr", "brms", "projpred", "orthopolynom")
results <- foreach(i = seq_len(nrow(design)),
                   .combine = bind_rows,
                   .packages = packages) %dopar% {
  print(i)
  # to ensure that all files can be found by the parallel workers
  setwd(wd)
  # library(dplyr)
  # library(brms)
  # library(projpred)
  # library(orthopolynom)
  source("PCE_helpers.R")
  source("sobol_helpers.R")
  
  out <- design[i, ]
  suffix <- paste0("mc", out$mc, "_M", out$M)
  if (out$type != "full") {
    suffix <- paste0(suffix, "_", out$type)
  }
  file <- paste0("models/fit_sobol_", suffix, ".rds")
  if (file.exists(file)) {
    fit <- readRDS(file)
  } else {
    return(out)
  }
  
  
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
  out$sobol_mean <- sobol_mean(a)
  out$PCE_mean_estimate <- pce_mean[1,1]
  out$PCE_mean_se <- pce_mean[1,2]
  out$PCE_mean_lower <- pce_mean[1,3]
  out$PCE_mean_upper <- pce_mean[1,4]
  
  pce_sd <- PCE_sd(fit)
  out$sobol_sd <- sobol_sd_num(a)
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

#stopCluster(cl)

saveRDS(results, "sobol_summary.rds")
