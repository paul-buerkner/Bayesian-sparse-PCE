library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(orthopolynom)
library(SobolSequence)
source("PCE_helpers.R")
source("sobol_helpers.R")

sparse_formula <- function(M, idx) {
  stopifnot(M %in% c(4, 8))
  if (M == 4) {
    out <- paste("y ~ PCE(x1s, x2s, x3s, x4s, p = 10, idx = ", idx, ")")
  } else if (M == 8) {
    out <- paste("y ~ PCE(x1s, x2s, x3s, x4s, x5s, x6s, x7s, x8s, p = 6, idx = ", idx, ")")
  }
  as.formula(out, env = new.env())
}

# get training data sets
mc <- c(100, 300, 900, 1001, 2700, 3003, 8100)
nterms <- c(25, 25, 25, 25, 25, 25, 25)

a <- c(1, 2, 5, 10, 20, 50, 100, 500)

set.seed(1234)
data_list <- vector("list", length(mc))
names(data_list) <- mc
for (i in seq_along(data_list)) {
  tmp <- sobolSequence.points(8, count = mc[i]) %>%
    as.data.frame() %>%
    rename(
      x1 = V1, x2 = V2, x3 = V3, x4 = V4, 
      x5 = V5, x6 = V6, x7 = V7, x8 = V8
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
  
  pol8 <- with(tmp, PCE(x1s, x2s, x3s, x4s, x5s, x6s, x7s, x8s, p = 6)) %>% 
    as.data.frame() 
  names(pol8) <- paste0("PE", 1:ncol(pol8))
  
  pol4 <- with(tmp, PCE(x1s, x2s, x3s, x4s, p = 10)) %>% 
    as.data.frame() 
  names(pol4) <- paste0("PF", 1:ncol(pol4))
  
  tmp <- bind_cols(tmp, pol8, pol4)
  
  data_list[[i]] <- tmp
}

# model fitting preparations
design <- expand.grid(mc = mc, M = c(4, 8))
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
  source("sobol_helpers.R")
  
  mc <- as.character(design$mc[j]) 
  M <- design$M[j]
  nterms <- design$nterms[j]
  
  training <- data_list[[mc]]
  P_names <- ifelse(M == 4, "PF", "PE")
  training_sel <- training %>% select("y", starts_with(P_names))
  
  fit_sobol <- brm(
    y ~ ., 
    data = training_sel,
    prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2)) + 
      prior(normal(0, 0.05), class = "sigma"),
    chains = 2, cores = 2,
    control = list(adapt_delta = 0.99),
    backend = "cmdstanr",
    file = paste0("models/fit_sobol_mc", mc, "_M", M)
  )
  
  # extract selected polynomials to remove the main model before
  # fitting the sparse sub-models. This avoids the file sizes of the
  # sub-models to blow up (exact reason for this still unknown)
  b2_mean <- colMeans(PCE_sobol_coef(fit_sobol))
  idx_b_sel <- order(b2_mean, decreasing = TRUE)[1:nterms]
  idx_phi_sel <- select_poly_R2D2(fit_sobol, thres = 0, nterms = nterms)
  file_fit_varsel <- paste0("models/fit_sobol_mc", mc, "_M", M, "_varsel")
  if (!file.exists(file_fit_varsel)) {
    # no need to run full varsel here; getting the search path is enough
    # varsel <- run_varsel(fit_sobol, nterms_max = nterms, method = "L1")
    # idx_varsel <- solution_terms(varsel)
    # idx_varsel <- as.numeric(sub("^[^[:digit:]]+", "", idx_varsel))
    # rm(varsel)
    idx_varsel <- projpred_L1_search_path(fit_sobol, nterms_max = nterms)
  }
  rm(fit_sobol)
  
  fit_sobol_b_sel <- brm(
    sparse_formula(M = M, idx = "idx_b_sel"),
    # y ~ PCE(x1s, x2s, x3s, x4s, x5s, x6s, x7s, x8s, p = 6, idx = idx_b_sel),
    data = training,
    prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2)) + 
      prior(normal(0, 0.05), class = "sigma"),
    data2 = list(idx_b_sel = idx_b_sel),
    chains = 2, cores = 2,
    control = list(adapt_delta = 0.95),
    backend = "cmdstanr",
    file = paste0("models/fit_sobol_mc", mc, "_M", M, "_b_sel")
  )
  rm(fit_sobol_b_sel)
  
  fit_sobol_phi_sel <- brm(
    sparse_formula(M = M, idx = "idx_phi_sel"),
    # y ~ PCE(x1s, x2s, x3s, x4s, x5s, x6s, x7s, x8s, p = 6, idx = idx_phi_sel), 
    data = training,
    prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2)) + 
      prior(normal(0, 0.05), class = "sigma"),
    data2 = list(idx_phi_sel = idx_phi_sel),
    chains = 2, cores = 2,
    control = list(adapt_delta = 0.95),
    backend = "cmdstanr",
    file = paste0("models/fit_sobol_mc", mc, "_M", M, "_phi_sel")
  )
  rm(fit_sobol_phi_sel)
  
  fit_sobol_varsel <- brm(
    sparse_formula(M = M, idx = "idx_varsel"),
    # y ~ PCE(x1s, x2s, x3s, x4s, x5s, x6s, x7s, x8s, p = 6, idx = idx_varsel), 
    data = training,
    prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2)) + 
      prior(normal(0, 0.05), class = "sigma"),
    data2 = list(idx_varsel = idx_varsel),
    chains = 2, cores = 2,
    control = list(adapt_delta = 0.95),
    backend = "cmdstanr",
    file = file_fit_varsel
  )
  rm(fit_sobol_varsel)
  NULL
}

# stopCluster(cl)

