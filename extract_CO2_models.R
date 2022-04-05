library(brms)
library(dplyr)
library(readr)
library(orthopolynom)
library(foreach)
library(doParallel)
source("PCE_helpers.R")
source("CO2_helpers.R")

# create test data
response_full <- read.table("data/CO2_Response/samples_10k__S_0010.dat")
names(response_full) <- paste0("y", seq_len(ncol(response_full)))

full <- read.table("data/CO2_Response/InputParameters.txt") %>%
  select(-V2) %>%
  rename(x1 = V1, x2 = V3, x3 = V4) %>%
  bind_cols(response_full)

set.seed(12345)
test <- full %>% sample_n(500)

poly_test <- with(test, PCE_CO2(x1, x2, x3, p = 10))
colnames(poly_test) <- paste0("P", seq_len(ncol(poly_test)))
test <- bind_cols(test, as.data.frame(poly_test))


d <- c(2, 3, 4, 5, 10)
coords <- 1:250
types <- c("full", "b_sel", "b_sel_50", "varsel", "varsel_50")
design <- expand.grid(d = d, type = types, coord = coords, 
                      stringsAsFactors = FALSE)

wd <- getwd()

cl <- makePSOCKcluster(10)
registerDoParallel(cl)

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
  source("CO2_helpers.R")
  
  out <- design[i, ]
  suffix <- paste0("d", out$d, "_co", out$coord)
  if (out$type != "full") {
    suffix <- paste0(suffix, "_", out$type)
  }
  fit <- readRDS(paste0("models/fit_CO2_", suffix, ".rds"))
  
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
  out$CO2_mean <- CO2_mean(out$coord)
  out$PCE_mean_estimate <- pce_mean[1,1]
  out$PCE_mean_se <- pce_mean[1,2]
  out$PCE_mean_lower <- pce_mean[1,3]
  out$PCE_mean_upper <- pce_mean[1,4]
  
  pce_sd <- PCE_sd(fit)
  out$CO2_sd <- CO2_sd(out$coord)
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

write_rds(results, "CO2_summary.rds")
