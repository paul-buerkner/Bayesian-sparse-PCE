CO2_mean <- function(s = NA) {
  read.table("data/CO2_Response/mu.dat")[s, 1]
}

CO2_sd <- function(s = NA)  {
  read.table("data/CO2_Response/sigma.dat")[s, 1]
}


poly1 <- read.table("data/CO2_Response/npc_0_10.dat") %>%
  apply(1, polynom::as.polynomial)
poly2 <- read.table("data/CO2_Response/npc_2_10.dat") %>%
  apply(1, polynom::as.polynomial)
poly3 <- read.table("data/CO2_Response/npc_3_10.dat") %>%
  apply(1, polynom::as.polynomial)


PCE_CO2 <- function(x1, x2, x3, p = 10, idx = NULL, scale = FALSE) {
  out <- PCE(
    x1, x2, x3, p = p, idx = idx, scale = FALSE, 
    poly = list(poly1, poly2, poly3)
  )
  if (scale) {
    # TODO: do we want to scale at all?
    # scale based on the richest data source
    full <- PCE(
      full$x1, full$x2, full$x3,
      p = p, idx = idx, scale = FALSE,
      poly = list(poly1, poly2, poly3)
    )
    means <- colMeans(full)
    sds <- sqrt(diag(cov(full)))
    out <- sweep(out, 2, means, "-")
    out <- sweep(out, 2, sds, "/")
  }
  out
}
