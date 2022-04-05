sobol <- function(..., a = c(1, 2, 5, 10, 20, 50, 100, 500)) {
  X <- list(...)
  stopifnot(length(a) == length(X))
  out <- 1
  for (i in seq_along(X)) {
    out <- out * (abs(4 * X[[i]] - 2) + a[i]) / (1 + a[i])
  }
  out
} 

sobol_for_int <- function(x, a = c(1, 2, 5, 10, 20, 50, 100, 500)) {
  stopifnot(length(a) == length(x))
  prod((abs(4 * x - 2) + a) / (1 + a))
}

sobol_mean <- function(a) {
  1
}

sobol_mean_num <- function(a) {
  ndim <- length(a)
  cubature::adaptIntegrate(sobol_for_int, rep(0, ndim), rep(1, ndim), a = a)$integral
}

sobol_sd_num <- function(a) {
  sobol2_for_int <- function(x, a) {
    sobol_for_int(x, a)^2
  }
  ndim <- length(a)
  sqrt(cubature::adaptIntegrate(sobol2_for_int, rep(0, ndim), rep(1, ndim), a = a)$integral -
    sobol_mean(a)^2)
}

sobol_mean_sim <- function(a, N = 100000) {
  X <- lapply(a, function(z) runif(N, 0, 1))
  args <- c(X, list(a = a))
  y <- do.call(sobol, args)
  mean(y)
}

sobol_sd_sim <- function(a, N = 100000) {
  X <- lapply(a, function(z) runif(N, 0, 1))
  args <- c(X, list(a = a))
  y <- do.call(sobol, args)
  sd(y)
}

PCE_summary_sobol <- function(fit, summary = TRUE, ceffects = TRUE,
                              corder = TRUE, newdata = NULL,
                              p = NULL, M = NULL, 
                              a = c(1, 2, 5, 10, 20, 50, 100, 500)) {
  require(ggplot2)
  if (summary) {
    summary(fit)
  }
  
  plot(fit, pars = parnames(fit)[1:5], fixed = TRUE)
  plot(fit, pars = "^sigma$")
  plot(pairs(fit, pars = parnames(fit)[1:5], fixed = TRUE))
  plot(mcmc_plot(fit, pars = parnames(fit)[1:24], type = "dens", fixed = TRUE))
  
  if (ceffects) {
    plot(conditional_effects_sobol(fit, "x1s", a = a))
    plot(conditional_effects_sobol(fit, "x2s", a = a))
    plot(conditional_effects_sobol(fit, "x3s", a = a))
    plot(conditional_effects_sobol(fit, "x4s", a = a))
    plot(conditional_effects_sobol(fit, "x5s", a = a))
    plot(conditional_effects_sobol(fit, "x6s", a = a))
    plot(conditional_effects_sobol(fit, "x7s", a = a))
    plot(conditional_effects_sobol(fit, "x8s", a = a))
  }
  
  if (corder) {
    plot(plot_poly_degree(fit, p = p, M = M))
    plot(plot_sobol_coef_degree(fit, yintercept = sobol_sd_num(a)))
  }
  
  cat("In-sample R^2\n")
  print(bayes_R2(fit))
  if (!is.null(newdata)) {
    cat("Out-of-sample R^2\n")
    print(bayes_R2(fit, newdata = newdata))
  }
  cat("\n")
  
  cat("In-sample RMSEs\n")
  yrep <- posterior_epred(fit)
  print(mean(rmse(fit$data$y, yrep)))
  print(rmse_mean(fit$data$y, yrep))
  if (!is.null(newdata)) {
    cat("Out-of-sample RMSEs\n")
    yrep <- posterior_epred(fit, newdata = newdata)
    print(mean(rmse(newdata$y, yrep))) 
    print(rmse_mean(newdata$y, yrep))
  }
  cat("\n")
  
  cat("Approximate PCE mean\n")
  print(PCE_mean(fit))
  cat("\n")
  
  cat("Approximate PCE SD\n")
  print(PCE_sd(fit))
  cat("\n")
  
  invisible(fit)
}

conditional_effects_sobol <- function(fit, effects, 
                                      a = c(1, 2, 5, 10, 20, 50, 100, 500), 
                                      ...) {
  require(ggplot2)
  ce <- conditional_effects(fit, effects, ...)
  stopifnot(length(ce) == 1L)
  ce[[1]] <- ce[[1]] %>%
    mutate(
      x1 = scale_from_1(x1s, 0, 1),
      x2 = scale_from_1(x2s, 0, 1),
      x3 = scale_from_1(x3s, 0, 1),
      x4 = scale_from_1(x4s, 0, 1),
      x5 = scale_from_1(x5s, 0, 1),
      x6 = scale_from_1(x6s, 0, 1),
      x7 = scale_from_1(x7s, 0, 1),
      x8 = scale_from_1(x8s, 0, 1),
      truth = sobol(x1, x2, x3, x4, x5, x6, x7, x8, a = a)
    )
  out <- plot(ce, plot = FALSE)[[1]]
  out + geom_line(aes(y = truth), size = 1.5, linetype = "dashed")
}
