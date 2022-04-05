# TODO: rename arguments for consistency with notation in the paper

Rcpp::sourceCpp("PCE_helpers.cpp")

scale2 <- function(x, na.rm = TRUE) {
  (x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)
}

center <- function(x, na.rm = TRUE) {
  x - mean(x, na.rm = na.rm)
}

# scale to [-1, 1] from [lb, ub]
scale_to_1 <- function(x, lb = -pi, ub = pi) {
  2 / (ub - lb) * (x - lb) - 1
}

# scale to [lb, ub] from [-1, 1]
scale_from_1 <- function(x, lb = -pi, ub = pi) {
  (ub - lb) / 2 * (x+1) + lb
}

# rmse over data points per posterior draw
# returns a vector of length ndraws
rmse <- function(y, yrep) {
  yrep <- as.matrix(yrep)
  y <- matrix(y, byrow = TRUE, nrow = nrow(yrep), ncol = ncol(yrep))
  sqrt(rowMeans((y - yrep)^2))
}

# rmse over data points ignoring posterior uncertainty
# returns a scalar
rmse_mean <- function(y, yrep) {
  yrep <- as.matrix(yrep)
  yrep_mean <- colMeans(yrep)
  sqrt(mean((y - yrep_mean)^2))
}

# rmse over data points for point predictions
# returns a scalar
rmse_mean2 <- function(y, yrep) {
  sqrt(mean((y - yrep)^2))
}

legendre_polynomials <- function(p) {
  orthopolynom::legendre.polynomials(p, normalized = TRUE)
}

# generates multi-indices of multivariate polynomial base
# uses graded lexicographic ordering (P. 156, Sullivan)
# adapted from Matlab code of Ilja Kröker
# faster version in C++ is poly_idx_cpp
# @param p maximal order of polynomials (named 'd' in the paper)
# @param M number of input variables
poly_idx <- function(p, M) {
  P = choose(p + M, M)
  out = matrix(0, P, M)
  tA = rep(0, M)
  l = 2
  pmax = (p + 1)^M
  for (i in 1:pmax) {
    ri = i
    for (d in 1:M) {
      md = (p + 1)^(M - d)
      val = floor(ri / md)
      tA[d] = val 
      ri = ri - val * md
    }
    if (sum(tA) <= p) {
      out[l, ] = tA
      l = l + 1
    }
  }
  out
}

PCE <- function(..., p = 10, idx = NULL, scale = TRUE,
                poly = legendre_polynomials) {
  dots <- list(...)
  N <- length(dots[[1]])
  M <- length(dots)
  comb <- poly_idx_cpp(p, M)
  # first column is the constant polynomial f0
  comb <- comb[-1, , drop = FALSE]
  rownames(comb) <- seq_len(nrow(comb))
  if (!is.null(idx)) {
    # select only desired polynomials
    comb <- comb[idx, , drop = FALSE]
  }
  if (is.function(poly)) {
    poly <- poly(p)
    poly <- replicate(M, poly, simplify = FALSE)
  } 
  stopifnot(is.list(poly) && length(poly) == M)
  out <- matrix(1, N, nrow(comb))
  for (i in seq_len(nrow(comb))) {
    for (j in seq_len(ncol(comb))) {
      out[, i] <- out[, i] * predict(poly[[j]][[comb[i,j] + 1]], dots[[j]])
    }
  }
  if (scale) {
    # TODO: make more general
    out <- out / sqrt(polynomial_variance_legendre(M)) 
  }
  colnames(out) <- rownames(comb)
  out
}

PCE_mean <- function(fit) {
  # "Intercept" is always the mean of the outcome variable
  as.data.frame(fit, variable = "b_Intercept") %>%
    posterior_summary()
}

PCE_mean2 <- function(fit, variable = "^b_P", regex = TRUE) {
  # "Intercept" is always the mean of the outcome variable
  b_Intercept <- as.matrix(fit, variable = "b_Intercept")[, 1]
  b <- as.matrix(fit, variable = variable, regex = regex)
  means_X <- colMeans(standata(fit)$X[, -1, drop = FALSE])
  means_X <- matrix(means_X, nrow = nrow(b), ncol = ncol(b), byrow = TRUE)
  Intercept <- b_Intercept + rowMeans(means_X * b)
  posterior_summary(Intercept)
}


# variance of evaluated polynomials from a number of variables
polynomial_variance_legendre <- function(M) {
  out <- 0.4900143 / (2 ^ (M - 1))
  # if (!is.null(cols)) {
  #   out <- rep(out, length(cols))
  # }
  out
}

polynomial_variance_empirical <- function(fit, cols = NULL) {
  mm <- standata(fit)$X[, -1, drop = FALSE]
  if (!is.null(cols)) {
    mm <- mm[, cols, drop = FALSE]    
  }
  mm %>% apply(2, var)
}

PCE_sobol_coef <- function(fit, variable = "^b_P", regex = TRUE) {
  if (!any(grepl(variable, variables(fit), fixed = !regex))) {
    return(data.frame())
  }
  ps <- as.data.frame(fit, variable = variable, regex = regex)
  cols <- sub("^b_", "", colnames(ps))
  # use the empirical approach to compute variances?
  # vars <- polynomial_variance(M, cols = cols)
  ps %>%
    mutate_all(~.^2) %>%
    as.matrix() %>%
    # sweep(2, vars, "*") %>%
    identity()
}

PCE_sd <- function(fit, variable = "^b_P", regex = TRUE) {
  out <- PCE_sobol_coef(fit, variable = variable, regex = regex)
  if (!NCOL(out)) {
    return(data.frame(0, 0, 0, 0))
  }
  out %>%
    rowSums() %>%
    sqrt() %>%
    posterior_summary()
}

# select coefficients using parameters of the R2D2 prior
select_poly_R2D2 <- function(fit, thres = 0.005, nterms = -1) {
  phi <- as.data.frame(fit, variable = "^R2D2_phi", regex = TRUE)
  R2 <- as.data.frame(fit, variable = "R2D2_R2")[[1]]
  # use the empirical approach to compute variances?
  # vars <- polynomial_variance(M, cols = colnames(phi))
  # variances are only constant in the infinite sample case
  # phi <- sweep(phi, 2, vars, "*")
  phi <- phi / rowSums(phi)
  ncoef <- ncol(phi)
  out <- NULL
  phi_sum <- 0
  for (i in seq_len(ncoef)) {
    phi_mean_max <- 0
    for (j in setdiff(seq_len(ncoef), out)) {
      phi_sum_tmp <- phi_sum + phi[, j]
      phi_mean_tmp <- mean(phi_sum_tmp)
      if (phi_mean_tmp > phi_mean_max) {
        phi_mean_max <- phi_mean_tmp
        sel <- j
      }
    }
    if (mean(R2 * phi[, sel]) < thres) {
      break
    } else {
      phi_sum <- phi_sum + phi[, sel]
      out <- c(out, sel)
      if (length(out) == nterms) {
        break
      }
    }
  }
  out
}

plot_poly_degree <- function(fit, p, M, variable = "^b_P") {
  comb <- poly_idx_cpp(p, M)[-1, , drop = FALSE]
  idx <- fit$data2$idx
  if (!is.null(idx)) {
    # select only desired polynomials
    comb <- comb[idx, , drop = FALSE]
  }
  order <- rowSums(comb)
  # order <- factor(order, levels = min(order):max(order))
  sobol <- PCE_sobol_coef(fit, variable = variable) 
  sobol_mean <- colMeans(sobol)
  sobol_lower <- apply(sobol, 2, quantile, probs = 0.025)
  sobol_upper <- apply(sobol, 2, quantile, probs = 0.975)
  data.frame(order, sobol_mean, sobol_lower, sobol_upper) %>%
    ggplot(aes(order, sobol_mean, ymin = sobol_lower, ymax = sobol_upper)) +
    geom_pointrange(position = position_jitter(width = 0.1)) +
    ylab("Sobol index") +
    xlab("Joint polynomial degree")
}

plot_sobol_coef_degree <- function(fit, yintercept, variable = "^b_P") {
  sobol <- PCE_sobol_coef(fit, variable = variable) 
  sobol_mean <- colMeans(sobol)
  sobol <- sobol[, order(sobol_mean, decreasing = TRUE), drop = FALSE]
  sobol_cumsum <- t(apply(sobol, 1, cumsum)) %>% 
    unname() %>%
    as.data.frame() %>%
    rename_all(~sub("^V", "", .)) %>%
    mutate_all(as.numeric) %>%
    gather() %>%
    mutate(
      key = as.numeric(as.character(key)),
      value = sqrt(value)
    ) %>%
    group_by(key) %>%
    summarise(
      mean = mean(value), 
      lower = quantile(value, probs = 0.025),
      upper = quantile(value, probs = 0.975),
      prob = mean(value > yintercept)
    ) %>%
    identity()
  
  xintercept <- which(sobol_cumsum$prob > 0.2)[1]
  if (is.na(xintercept)) {
    seql <- seq_len(nrow(sobol_cumsum))
    l <- c(0, sobol_cumsum$mean[seql[-length(seql)]])
    u <- c(sobol_cumsum$mean[seql[-1]], max(sobol_cumsum$mean))
    xintercept <- which(abs(u - l) / yintercept < 0.005)[1]
  }
  if (is.na(xintercept)) {
    xintercept <- ncol(sobol_cumsum)
  }
  
  sobol_cumsum %>%
    ggplot(aes(key, mean, ymin = lower, ymax = upper)) +
    geom_pointrange() +
    geom_hline(yintercept = yintercept) +
    geom_vline(xintercept = xintercept) +
    ylab("Total Sobol index") +
    xlab("Polynomial (sorted)") +
    ylim(c(0, NA))
}

df_poly_degree <- function(fit, p, M, variable = "^b_P", digits = 2) {
  comb <- poly_idx_cpp(p, M)[-1, , drop = FALSE]
  idx <- fit$data2$idx
  if (!is.null(idx)) {
    # select only desired polynomials
    comb <- comb[idx, , drop = FALSE]
  }
  colnames(comb) <- paste0("O(X", 1:M, ")")
  sobol <- PCE_sobol_coef(fit, variable = variable) 
  sobol_mean <- colMeans(sobol)
  sobol_lower <- apply(sobol, 2, quantile, probs = 0.025)
  sobol_upper <- apply(sobol, 2, quantile, probs = 0.975)
  df <- data.frame(sobol_mean, sobol_lower, sobol_upper) %>%
    bind_cols(as.data.frame(comb)) %>%
    arrange(desc(sobol_mean)) %>%
    mutate_all(~round(., digits)) %>%
    rename(Estimate = sobol_mean) %>%
    mutate(
      `95%-CI` = paste0("[", sobol_lower, ", ", sobol_upper, "]")
    ) %>%
    select(-sobol_lower, -sobol_upper) %>%
    select(Estimate, `95%-CI`, everything()) %>%
    identity()
  
  rownames(df) <- NULL
  # if (!is.null(nterms)) {
  #   df <- df[1:nterms, , drop = FALSE]
  # }
  # df %>%
  #   ggplot(aes(order1, order2, color = sobol_mean)) +
  #   geom_point(size = 3, position = position_jitter(width = 0.1)) + 
  #   scale_color_viridis_c() + 
  #   scale_x_continuous(breaks = 0:p) +
  #   scale_y_continuous(breaks = 0:p) +
  #   ylab("2nd Input Dimension") +
  #   xlab("1st Input Dimension")
  df
}

run_varsel <- function(fit, cv = FALSE, file = NULL, ...) {
  if (!is.null(file) && file.exists(file)) {
    return(read_rds(file))
  }
  if (cv) {
    out <- cv_varsel(fit, ...) 
  } else {
    out <- varsel(fit, ...) 
  }
  if (!is.null(file)) {
    out <- strip_varsel_env(out)
    saveRDS(out, file)
  }
  out
}

# remove huge (unused) environments
strip_varsel_env <- function(vsel) {
  is_fun <- unlist(lapply(vsel$refmodel, is.function))
  vsel$refmodel[is_fun] <- NULL
  vsel$refmodel$family$mu_fun <- NULL
  vsel$family$mu_fun <- NULL
  vsel
}
