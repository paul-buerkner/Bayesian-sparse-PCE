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
nterms <- c(25, 25, 25, 25, 25, 25, 25, 25)

a <- 7
b <- 0.1

set.seed(1234)
data_list <- vector("list", length(mc))
names(data_list) <- mc
for (i in seq_along(data_list)) {
  tmp <- sobolSequence.points(3, count = mc[i]) %>%
    as.data.frame() %>%
    rename(x1 = V1, x2 = V2, x3 = V3) %>%
    mutate_all(~scale_from_1(scale_to_1(., lb = 0, ub = 1))) %>%
    mutate(
      y = ishigami(x1, x2, x3, a, b),
      x1s = scale_to_1(x1),
      x2s = scale_to_1(x2),
      x3s = scale_to_1(x3)
    )
  
  pol <- with(tmp, PCE(x1s, x2s, x3s, p = 10)) %>% as.data.frame() 
  names(pol) <- paste0("P", 1:ncol(pol))
  tmp <- bind_cols(tmp, pol)
  
  data_list[[i]] <- tmp
}

# model fitting preparations
noise <- c(0.1, 0.3, 0.5)
ntrials <- 10
design <- expand.grid(mc = mc, noise = noise, trial = 1:ntrials)
design$nterms <- nterms[match(design$mc, mc)]

exports <- NULL

# fit the models in parallel
# cl <- makeForkCluster(detectCores())
# registerDoParallel(cl)

J <- seq_len(nrow(design))

foreach(j = J, .export = exports) %dopar% {
  library(dplyr)
  library(brms)
  library(projpred)
  library(SobolSequence)
  source("PCE_helpers.R")
  source("ishigami_helpers.R")
  env <- new.env()
  
  mc <- as.character(design$mc[j]) 
  nterms <- design$nterms[j]
  noise <- design$noise[j]
  trial <- design$trial[j]
  suffix <- paste0("mc", mc, "_noise", noise, "_trial", trial)
  
  training <- data_list[[mc]]
  # add gaussian noise
  training$y <- training$y + rnorm(nrow(training), 0, noise)
  training_sel <- training %>% select("y", starts_with("P"))
  
  fit_ishigami <- brm(
    y ~ ., 
    data = training_sel,
    prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2)) + 
      prior(normal(0, 0.05), class = "sigma"),
    chains = 2, cores = 2,
    control = list(adapt_delta = 0.95),
    backend = "cmdstanr",
    file = paste0("models/noise/fit_ishigami_", suffix)
  )
  
  b2_mean <- colMeans(PCE_sobol_coef(fit_ishigami))
  idx_b_sel <- order(b2_mean, decreasing = TRUE)[1:nterms]
  idx_phi_sel <- select_poly_R2D2(fit_ishigami, thres = 0, nterms = nterms)
  file_fit_varsel <- paste0("models/noise/fit_ishigami_", suffix, "_varsel")
  if (!file.exists(file_fit_varsel)) {
    idx_varsel <- projpred_L1_search_path(fit_ishigami, nterms_max = nterms)
  }
  rm(fit_ishigami)
  
  fit_ishigami_b_sel <- brm(
    y ~ PCE(x1s, x2s, x3s, p = 10, idx = idx_b_sel), 
    data = training,
    prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2)) + 
      prior(normal(0, 0.05), class = "sigma"),
    data2 = list(idx_b_sel = idx_b_sel),
    chains = 2, cores = 2,
    control = list(adapt_delta = 0.95),
    backend = "cmdstanr",
    file = paste0("models/noise/fit_ishigami_", suffix, "_b_sel")
  )
  rm(fit_ishigami_b_sel)
  
  fit_ishigami_varsel <- brm(
    y ~ PCE(x1s, x2s, x3s, p = 10, idx = idx_varsel), 
    data = training,
    prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2)) + 
      prior(normal(0, 0.05), class = "sigma"),
    data2 = list(idx_varsel = idx_varsel),
    chains = 2, cores = 2,
    control = list(adapt_delta = 0.95),
    backend = "cmdstanr",
    file = file_fit_varsel
  )
  rm(fit_ishigami_varsel)
  NULL
}

# stopCluster(cl)

