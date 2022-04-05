ishigami <- function(x1, x2, x3, a, b) {
  a * sin(x2)^2 + (1+b*(x3^4)) * sin(x1)
}

ishigami_mean <- function(a, b) {
  1/2 * a
}

ishigami_sd <- function(a, b) {
  sqrt(a^2 / 8 + b * pi^4 / 5 + b^2 * pi^8 / 18 + 1/2)
}

ishigami_mean_sim <- function(a, b, N = 100000) {
  dat <- data.frame(
    x1 = runif(N, -pi, pi),
    x2 = runif(N, -pi, pi),
    x3 = runif(N, -pi, pi)
  ) %>%
    mutate(
      y = ishigami(x1, x2, x3, a, b),
    )
  mean(dat$y)
}

ishigami_sd_sim <- function(a, b, N = 100000) {
  dat <- data.frame(
    x1 = runif(N, -pi, pi),
    x2 = runif(N, -pi, pi),
    x3 = runif(N, -pi, pi)
  ) %>%
    mutate(
      y = ishigami(x1, x2, x3, a, b),
    )
  sd(dat$y)
}

PCE_summary_ishigami <- function(fit, summary = TRUE, ceffects = TRUE,
                                 corder = TRUE, newdata = NULL,
                                 p = 10, M = 3, a = 7, b = 0.1) {
  if (summary) {
    summary(fit)
  }
  
  plot(fit, pars = parnames(fit)[1:5], fixed = TRUE, N = 6)
  plot(pairs(fit, pars = parnames(fit)[1:5], fixed = TRUE))
  plot(mcmc_plot(fit, pars = parnames(fit)[1:24], type = "dens", fixed = TRUE))
  
  if (ceffects) {
    plot(conditional_effects_ishigami(fit, "x1s", a = a, b = b))
    plot(conditional_effects_ishigami(fit, "x2s", a = a, b = b))
    plot(conditional_effects_ishigami(fit, "x3s", a = a, b = b))
    plot(conditional_effects_ishigami(fit, "x1s:x3s", a = a, b = b))
    plot(conditional_effects_ishigami(fit, "x3s:x1s", a = a, b = b))
  }
  
  if (corder) {
    plot(plot_poly_degree(fit, p = p, M = M))
    plot(plot_sobol_coef_degree(fit, yintercept = ishigami_sd(a, b)))
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
    cat("Out-of-sample RMSE\n")
    yrep <- posterior_epred(fit, newdata = newdata)
    print(mean(rmse(newdata$y, yrep))) 
    print(rmse_mean(newdata$y, yrep))
  }
  cat("\n")
  
  cat("Approximate PCE mean\n")
  print(PCE_mean(fit))
  cat("\n")
  
  cat("Approximate PCE SD (assuming independence)\n")
  print(PCE_sd(fit))
  cat("\n")
  
  invisible(fit)
}

conditional_effects_ishigami <- function(fit, effects, a = 7, b = 0.1, ...) {
  ce <- conditional_effects(fit, effects, ...)
  stopifnot(length(ce) == 1L)
  ce[[1]] <- ce[[1]] %>%
    mutate(
      x1 = scale_from_1(x1s),
      x2 = scale_from_1(x2s),
      x3 = scale_from_1(x3s),
      truth = ishigami(x1, x2, x3, a = a, b = b)
    )
  out <- plot(ce, plot = FALSE)[[1]]
  out + geom_line(aes(y = truth), size = 1.5, linetype = "dashed")
}
