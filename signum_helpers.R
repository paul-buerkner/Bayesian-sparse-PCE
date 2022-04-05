source("PCE_helpers.R")
library(dplyr)

signum_mean <- function() {
  0
}

signum_sd <- function()  {
  1
}

poly_signum <- function(path_poly = "data/signum_function/OrthonormalBasis.txt", ...) {
  read.table(path_poly) %>%
    apply(1, as.polynomial)
}

PCE_signum <- function(x, p = 10, idx = NULL, ...) {
  PCE(
    x, p = p, idx = idx, scale = FALSE, 
    poly = list(poly_signum(...))
  )
}

PCE_signum_coefficients <- function(p, strategy, path_train = NULL, ...) {
  if (is.null(path_train)) {
    if (strategy == "SoSeq") {
      path_train <- 
        paste0("data/signum_function/SobolSequence/Coefficients_", p, ".txt")
    } else if (strategy == "GaInt") {
      path_train <- 
        paste0("data/signum_function/GaussIntegration/Coefficients_", p, ".txt")
    }
  } 
  as.vector(as.matrix(read.table(path_train)))
}

conditional_effects_signum <- function(fit, thres = 5, ...) {
  int_cond <- list(x = seq(-1, 1, length.out = 100))
  ce <- conditional_effects(fit, effects = "x", 
                            int_conditions = int_cond, ...)
  stopifnot(length(ce) == 1L)
  ce[[1]] <- ce[[1]] %>%
    mutate(
      truth = sign(x),
      index_first_large = which(abs(estimate__) > thres)[1],
      index_first_large = ifelse(is.na(index_first_large), Inf,
                                 index_first_large),
      estimate__ = ifelse(row_number() >= index_first_large, 
                          NA, estimate__)
    )
  
  plot(ce, plot = FALSE)[[1]] + 
    geom_line(aes(y = truth), size = 1.5, linetype = "dashed") +
    xlim(c(-1, 1)) +
    ylim(c(-thres, thres))
}


predict_standard_PCE <- function(x, p, ...) {
  mc <- length(x)
  poly_eval <- cbind(1, PCE_signum(x, p = p, ...))
  coef <- PCE_signum_coefficients(p, ...)
  coef <- matrix(coef, nrow = mc, ncol = length(coef), byrow = TRUE)
  rowSums(coef * poly_eval)
} 

fit_standard_PCE <- function(x, p, ...) {
  lm(sign(x) ~ 1 + PCE_signum(x, p = p, ...))
}

predict_standard_PCE_lm <- function(x, p, new_x, ...) {
  fit <- fit_standard_PCE(x, p, ...)
  predict(fit, newdata = data.frame(x = new_x))
}
