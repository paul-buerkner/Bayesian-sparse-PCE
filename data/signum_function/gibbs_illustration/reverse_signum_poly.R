# puts a reversed polynomial table back in a reasonable structure
reverse_signum_poly <- function(p) {
  path <- paste0("data/signum_function/gibbs_illustration/pol_cfs_No", p, ".dat")
  out <- read.table(path)
  .reverse_signum_poly <- function(x, o) {
    # o = maximum polynomial order
    y <- rev(x[1:(o+1)])
    return(c(y, rep(0, p + 1 - length(y))))
  }
  for (i in seq_len(nrow(out))) {
    out[i, ] <- .reverse_signum_poly(out[i, ], o = i - 1)
  }
  new_path <- paste0("data/signum_function/gibbs_illustration/pol_cfs_No", p, "_new.dat")
  write.table(out, new_path)
  out
}

reverse_signum_poly(10)
reverse_signum_poly(15)
reverse_signum_poly(25)
