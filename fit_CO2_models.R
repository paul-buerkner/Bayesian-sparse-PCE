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
source("CO2_helpers.R")

# get training data sets
d <- c(2, 3, 4, 5, 10)
mc <- (d + 1)^3

# currently unused variable
nterms <- c(25, 50)


data_list <- vector("list", length(d))
names(data_list) <- d
for (i in seq_along(data_list)) {
  path_resp <- paste0(
    "data/CO2_Response/d", d[i], "/Tad_Fkt5023_MR0_P0_MC", mc[i], 
    "_N250_T8640000.000000_tP-1.000000_tMR-1.000000_alpha0.500000__S_0000.dat"
  )
  response <- read.table(path_resp) 
  names(response) <- paste0("y", seq_len(ncol(response)))
  
  path_train <- paste0("data/CO2_Response/d", d[i], "/IntegrationPoints.txt")
  training <- read.table(path_train) %>%
    select(-V2) %>%
    rename(x1 = V1, x2 = V3, x3 = V4) %>%
    bind_cols(response)
  
  poly <- with(training, PCE_CO2(x1, x2, x3, p = 10))
  colnames(poly) <- paste0("P", seq_len(ncol(poly)))
  training <- bind_cols(training, as.data.frame(poly))
  
  data_list[[i]] <- training
  saveRDS(training, paste0("data/CO2_Response/training_d", d[i], ".rds"))
}

# model fitting preparations
coords <- 1:250
design <- expand.grid(coord = coords, d = d)
# design$nterms <- nterms[match(design$d, d)]

# fit the models in parallel
cl <- makeForkCluster(detectCores())
registerDoParallel(cl)
exports <- NULL

foreach(j = seq_len(nrow(design)), .export = exports) %dopar% {
  library(dplyr)
  library(brms)
  library(projpred)
  library(orthopolynom)
  source("PCE_helpers.R")
  source("CO2_helpers.R")
  env <- new.env()
  
  i <- design$coord[j]
  d <- as.character(design$d[j])
  nterms <- design$nterms[j]

  training <- data_list[[d]]
  training_sel <- training %>% select(paste0("y", i), starts_with("P"))
  
  # sdy <- max(sd(training[[paste0("y", i)]]), 0.05)
  # intercept_prior <- paste0("normal(0, ", sdy, ")")
  fit_CO2 <- brm(
    as.formula(paste0("y", i, "~ ."), env = env), 
    data = training_sel,
    prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2)) + 
      prior(normal(0, 1), class = "Intercept") +
      prior(normal(0, 0.1), class = "sigma"),
    chains = 1, cores = 1, warmup = 1000, iter = 3000,
    control = list(adapt_delta = 0.95),
    backend = "rstan", stan_model_args = list(save_dso = FALSE),
    # backend = "cmdstanr",
    file = paste0("models/fit_CO2_d", d, "_co", i)
  )
  
  b2_mean <- colMeans(PCE_sobol_coef(fit_CO2))
  idx_b_sel_50 <- order(b2_mean, decreasing = TRUE)[1:50]
  idx_b_sel_25 <- idx_b_sel_50[1:25]

  file_fit_varsel_50 <- paste0("models/fit_CO2_d", d, "_co", i, "_varsel_50")  
  # not named after "nterms = 25" to avoid refitting already existing models
  file_fit_varsel_25 <- paste0("models/fit_CO2_d", d, "_co", i, "_varsel")
  if (!file.exists(file_fit_varsel_50) || !file.exists(file_fit_varsel_25)) {
    # don't store varsel but still avoid rerunning it
    # varsel <- run_varsel(fit_CO2, nterms_max = 50, method = "L1")
    # plot(varsel, stats = c("elpd", "rmse"))
    # idx_varsel_50 <- solution_terms(varsel)
    # idx_varsel_50 <- as.numeric(sub("^[^[:digit:]]+", "", idx_varsel_50))
    # idx_varsel_25 <- idx_varsel_50[1:25]
    # rm(varsel)
    idx_varsel_50 <- projpred_L1_search_path(fit_CO2, nterms_max = 50)
    idx_varsel_25 <- idx_varsel_50[1:25]
  }
  
  # we decided to not include this method in the paper because
  # it yielded almost identical results to b_sel
  # idx_phi_sel <- select_poly_R2D2(fit_CO2, thres = 0, nterms = 50)
  # idx_phi_sel_25 <- idx_phi_sel[1:25]
  rm(fit_CO2)
  
  
  formula_b_sel_50 <- paste0("y", i, " ~ PCE_CO2(x1, x2, x3, p = 10, idx = idx_b_sel_50)")
  fit_CO2_b_sel_50 <- brm(
    as.formula(formula_b_sel_50, env = env),
    data = training,
    prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2)) + 
      prior(normal(0, 1), class = "Intercept") +
      prior(normal(0, 0.1), class = "sigma"),
    data2 = list(idx_b_sel_50 = idx_b_sel_50),
    chains = 1, cores = 1, warmup = 1000, iter = 3000,
    control = list(adapt_delta = 0.95),
    backend = "rstan", stan_model_args = list(save_dso = FALSE),
    # backend = "cmdstanr",
    file = paste0("models/fit_CO2_d", d, "_co", i, "_b_sel_50")
  )
  rm(fit_CO2_b_sel_50)
  
  
  formula_b_sel_25 <- paste0("y", i, " ~ PCE_CO2(x1, x2, x3, p = 10, idx = idx_b_sel_25)")
  fit_CO2_b_sel_25 <- brm(
    as.formula(formula_b_sel_25, env = env),
    data = training,
    prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2)) + 
      prior(normal(0, 1), class = "Intercept") +
      prior(normal(0, 0.1), class = "sigma"),
    data2 = list(idx_b_sel_25 = idx_b_sel_25),
    chains = 1, cores = 1, warmup = 1000, iter = 3000,
    control = list(adapt_delta = 0.95),
    backend = "rstan", stan_model_args = list(save_dso = FALSE),
    # backend = "cmdstanr",
    # not named after "nterms = 25" to avoid refitting already existing models
    file = paste0("models/fit_CO2_d", d, "_co", i, "_b_sel")
  )
  rm(fit_CO2_b_sel_25)

  
  formula_varsel_50 <- paste0("y", i, " ~ PCE_CO2(x1, x2, x3, p = 10, idx = idx_varsel_50)")
  fit_CO2_varsel_50 <- brm(
    as.formula(formula_varsel_50, env = env),
    data = training,
    prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2)) + 
      prior(normal(0, 1), class = "Intercept") +
      prior(normal(0, 0.1), class = "sigma"),
    data2 = list(idx_varsel_50 = idx_varsel_50),
    chains = 1, cores = 1, warmup = 1000, iter = 3000,
    control = list(adapt_delta = 0.95),
    backend = "rstan", stan_model_args = list(save_dso = FALSE),
    # backend = "cmdstanr",
    file = file_fit_varsel_50
  )
  rm(fit_CO2_varsel_50)
  
  
  formula_varsel_25 <- paste0("y", i, " ~ PCE_CO2(x1, x2, x3, p = 10, idx = idx_varsel_25)")
  fit_CO2_varsel_25 <- brm(
    as.formula(formula_varsel_25, env = env),
    data = training,
    prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2)) + 
      prior(normal(0, 1), class = "Intercept") +
      prior(normal(0, 0.1), class = "sigma"),
    data2 = list(idx_varsel_25 = idx_varsel_25),
    chains = 1, cores = 1, warmup = 1000, iter = 3000,
    control = list(adapt_delta = 0.95),
    backend = "rstan", stan_model_args = list(save_dso = FALSE),
    # backend = "cmdstanr",
    file = file_fit_varsel_25
  )
  rm(fit_CO2_varsel_25)
  
  
  # formula_phi_sel <- paste0("y", i, " ~ PCE_CO2(x1, x2, x3, p = 10, idx = idx_phi_sel)")
  # fit_CO2_phi_sel <- brm(
  #   as.formula(formula_phi_sel, env = env),
  #   data = training,
  #   prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2)) + 
  #     prior(normal(0, 1), class = "Intercept") +
  #     prior(normal(0, 0.1), class = "sigma"),
  #   data2 = list(idx_phi_sel = idx_phi_sel),
  #   chains = 1, cores = 1, warmup = 1000, iter = 3000,
  #   control = list(adapt_delta = 0.95),
  #   backend = "rstan", stan_model_args = list(save_dso = FALSE),
  #   # backend = "cmdstanr",
  #   file = paste0("models/fit_CO2_d", d, "_co", i, "_phi_sel")
  # )
  # rm(fit_CO2_phi_sel)
  
  NULL
}

stopCluster(cl)

