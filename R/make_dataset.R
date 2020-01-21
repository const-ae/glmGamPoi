

make_dataset <- function(n_genes = 1000, n_samples = 30){

  sf <- exp(rnorm(n = n_samples, mean = 1.7, sd = 0.7))
  overdispersion <- rgamma(n = n_genes, shape = 3, rate = 6)
  X <- matrix(1, nrow = n_samples)  # Only intercept
  Betas <- matrix(rnorm(n = n_genes, mean = 0, sd = 3), nrow = n_genes)

  Q <- Betas %*% t(X)

  Y <- t(sapply(seq_len(n_genes), function(i){
    distraltparam::raltnbinom(n_samples, mean = sf * exp(Q[i, ]), dispersion = overdispersion[i])
  }))

  row_all_zero <- which(rowSums(Y) == 0)
  for(idx in row_all_zero){
    Y[idx, sample(seq_len(n_samples), 1)] <- 1
  }

  list(
    Y = Y,
    X = X,
    Q = Q,
    overdispersion = overdispersion,
    size_factor = sf,
    Betas = Betas
  )
}
