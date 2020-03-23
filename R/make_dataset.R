

make_dataset <- function(n_genes = 1000, n_samples = 30){

  sf <- exp(rnorm(n = n_samples, mean = 1.7, sd = 0.7))
  dispersions <- rgamma(n = n_genes, shape = 3, rate = 6)
  X <- matrix(1, nrow = n_samples)  # Only intercept
  Betas <- matrix(rnorm(n = n_genes, mean = 0, sd = 3), nrow = n_genes)

  Q <- Betas %*% t(X)

  Y <- t(vapply(seq_len(n_genes), function(i){
    rnbinom(n_samples, mu = sf * exp(Q[i, ]), size = 1/dispersions[i])
  }, FUN.VALUE = rep(0.0, n_samples)))

  row_all_zero <- which(rowSums(Y) == 0)
  for(idx in row_all_zero){
    Y[idx, sample(seq_len(n_samples), 1)] <- 1
  }

  list(
    Y = Y,
    X = X,
    Q = Q,
    overdispersion = dispersions,
    size_factor = sf,
    Betas = Betas
  )
}
