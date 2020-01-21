// #include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;



double compute_gp_deviance (double y, double mu, double theta) {
  if(theta < 1e-6){
    // If theta is so small, calculate Poisson deviance
    if(y == 0){
      return 2.0 * mu;
    }else{
      return 2.0 * (y * std::log(y/mu) - (y - mu));
    }
  }else{
    // Otherwise calculate Gamma-Poisson deviance
    if(y == 0){
      return 2.0/theta * std::log((1 + mu * theta));
    } else {
      double s1 = y * std::log((mu + y * mu * theta) / (y +  y * mu * theta));
      double s2 = 1.0/theta * std::log((1 + mu * theta) / (1 + y * theta));
      return -2.0 * (s1 - s2);
    }
  }
}





//--------------------------------------------------------------------------------------------------//
// The following code was originally copied from https://github.com/mikelove/DESeq2/blob/master/src/DESeq2.cpp
// I adapted it to the needs of this project by:
//  * remove lambda / ridge penality
//  * remove weights
//  * Calculate actual deviance (2 * (log(f_NB(y | mu, theta)) - log(f_NB(y | y, theta))))
//    instead of just 2 * log(f_NB(y | mu, theta)),


// fit the Negative Binomial GLM.
// note: the betas are on the natural log scale
//
// [[Rcpp::export]]
List fitBeta(const arma::mat& y, const arma::mat& x, const arma::mat& nf, SEXP alpha_hatSEXP, SEXP beta_matSEXP,
             SEXP tolSEXP, SEXP maxitSEXP, SEXP minmuSEXP) {

  // arma::mat nf = as<arma::mat>(nfSEXP);
  // arma::mat x = as<arma::mat>(xSEXP);
  int y_n = y.n_rows;
  int y_m = y.n_cols;
  int x_p = x.n_cols;
  arma::vec alpha_hat = as<arma::vec>(alpha_hatSEXP);
  arma::mat beta_mat = as<arma::mat>(beta_matSEXP);
  arma::mat beta_var_mat = arma::zeros(beta_mat.n_rows, beta_mat.n_cols);
  arma::mat contrast_num = arma::zeros(beta_mat.n_rows, 1);
  arma::mat contrast_denom = arma::zeros(beta_mat.n_rows, 1);
  arma::mat hat_diagonals = arma::zeros(y.n_rows, y.n_cols);
  int maxit = as<int>(maxitSEXP);
  arma::colvec yrow, nfrow, beta_hat, mu_hat, z;
  arma::mat sigma;
  arma::vec w_vec, w_sqrt_vec;

  arma::colvec gamma_hat, big_z;
  arma::vec big_w_diag;
  arma::mat weighted_x, q, r, big_w_sqrt;
  // deviance, convergence and tolerance
  double dev, dev_old, conv_test;
  double tol = as<double>(tolSEXP);
  // bound the estimated count, as weights include 1/mu
  double minmu = as<double>(minmuSEXP);
  double large = 30.0;
  NumericVector iter(y_n);
  NumericVector deviance(y_n);
  for (int i = 0; i < y_n; i++) {
    if (i % 100 == 0) checkUserInterrupt();
    nfrow = nf.row(i).t();
    yrow = y.row(i).t();
    beta_hat = beta_mat.row(i).t();
    mu_hat = nfrow % exp(x * beta_hat);
    for (int j = 0; j < y_m; j++) {
      mu_hat(j) = fmax(mu_hat(j), minmu);
    }
    dev = 0.0;
    dev_old = 0.0;

    // make an orthonormal design matrix
    for (int t = 0; t < maxit; t++) {
      iter(i)++;
      w_vec = mu_hat/(1.0 + alpha_hat(i) * mu_hat);
      w_sqrt_vec = sqrt(w_vec);

      // prepare matrices
      weighted_x = x.each_col() % w_sqrt_vec;
      qr_econ(q, r, weighted_x);
      big_w_diag = arma::ones(y_m);
      big_w_diag(arma::span(0, y_m - 1)) = w_vec;
      // big_w_sqrt = diagmat(sqrt(big_w_diag));
      z = arma::log(mu_hat / nfrow) + (yrow - mu_hat) / mu_hat;
      arma::vec w_diag = w_vec;
      arma::mat z_sqrt_w = z.each_col() % sqrt(w_diag);
      arma::colvec big_z_sqrt_w = arma::zeros(y_m);
      big_z_sqrt_w(arma::span(0,y_m - 1)) = z_sqrt_w;
      // IRLS with Q matrix for X
      gamma_hat = q.t() * big_z_sqrt_w;
      solve(beta_hat, r, gamma_hat);
      if (sum(abs(beta_hat) > large) > 0) {
        iter(i) = maxit;
        break;
      }
      mu_hat = nfrow % exp(x * beta_hat);
      for (int j = 0; j < y_m; j++) {
        mu_hat(j) = fmax(mu_hat(j), minmu);
      }
      dev = 0.0;
      for (int j = 0; j < y_m; j++) {
        // note the order for Rf_dnbinom_mu: x, sz, mu, lg
        // dev = dev + -2.0 * Rf_dnbinom_mu(yrow(j), 1.0/alpha_hat(i), mu_hat(j), 1);
        dev = dev + compute_gp_deviance(yrow(j), mu_hat(j), alpha_hat(i));
      }
      conv_test = fabs(dev - dev_old)/(fabs(dev) + 0.1);
      if (std::isnan(conv_test)) {
        iter(i) = maxit;
        break;
      }
      if ((t > 0) & (conv_test < tol)) {
        break;
      }
      dev_old = dev;
    }

    deviance(i) = dev;
    beta_mat.row(i) = beta_hat.t();
    // recalculate w so that this is identical if we start with beta_hat
    w_vec = mu_hat/(1.0 + alpha_hat(i) * mu_hat);
    w_sqrt_vec = sqrt(w_vec);

    arma::vec hat_matrix_diag = arma::zeros(x.n_rows);
    arma::mat xw = x.each_col() % w_sqrt_vec;
    arma::mat xtwxr_inv = (x.t() * (x.each_col() % w_vec)).i();

    // Verbose, but fast way to get diagonal of:
    // hat_matrix = xw * xtwxr_inv * xw.t() ;
    for(int jp = 0; jp < y_m; jp++){
      for(int idx1 = 0; idx1 < x_p; idx1++){
        for(int idx2 = 0; idx2 < x_p; idx2++){
          hat_matrix_diag(jp) += xw(jp, idx1) * (xw(jp, idx2) * xtwxr_inv(idx2, idx1));
        }
      }
    }
    hat_diagonals.row(i) = hat_matrix_diag.t();
    // sigma is the covariance matrix for the betas
    sigma = (x.t() * (x.each_col() % w_vec)).i() * x.t() * (x.each_col() % w_vec) * (x.t() * (x.each_col() % w_vec)).i();
    beta_var_mat.row(i) = diagvec(sigma).t();
  }

  return List::create(Named("beta_mat",beta_mat),
                      Named("beta_var_mat",beta_var_mat),
                      Named("iter",iter),
                      Named("hat_diagonals",hat_diagonals),
                      Named("deviance",deviance));
}




// If there is only one group, there is no need to do the full Fisher-scoring
// Instead a simple Newton-Raphson algorithm will do
// [[Rcpp::export]]
List fitBeta_one_group(NumericMatrix Y, NumericMatrix log_offsets,
                       NumericVector thetas, NumericVector beta_start_values,
                       double tolerance, int maxIter) {
  int n_samples = Y.ncol();
  int n_genes = Y.nrow();
  NumericVector result(n_genes);
  IntegerVector iterations(n_genes);

  for(int gene_idx = 0; gene_idx < n_genes; gene_idx++){
    if (gene_idx % 100 == 0) checkUserInterrupt();
    NumericMatrix::Row counts = Y(gene_idx, _);
    NumericMatrix::Row log_off = log_offsets(gene_idx, _);
    const double& theta = thetas(gene_idx);
    double beta = beta_start_values(gene_idx);
    // Newton-Raphson
    int iter = 0;
    for(; iter < maxIter; iter++){
      double dl = 0.0;
      double ddl = 0.0;
      for(int sample_iter = 0; sample_iter < n_samples; sample_iter++){
        const int count = counts[sample_iter];
        const double mu = std::exp(beta + log_off[sample_iter]);
        const double denom = 1.0 + mu * theta;
        dl += (count - mu) / denom;
        ddl += mu * (1.0 + count * theta) / denom / denom;
        // ddl += mu / denom;           // This is what edgeR is using
      }
      const double step = dl / ddl;
      beta += step;
      if(abs(step) < tolerance){
        break;
      }
    }
    result(gene_idx) = beta;
    iterations(gene_idx) = iter;
  }
  return List::create(
    Named("beta", result),
    Named("iter", iterations)
  );
}



