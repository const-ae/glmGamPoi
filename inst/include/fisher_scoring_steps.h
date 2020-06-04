#ifndef FISHER_SCORING_STEPS_H
#define FISHER_SCORING_STEPS_H

#include <RcppArmadillo.h>

using namespace Rcpp;

template<class NumericType>
arma::vec fisher_scoring_qr_step(const arma::mat& model_matrix, const arma::Col<NumericType>& counts,
                                 const arma::colvec& mu, const arma::colvec& theta_times_mu){
  // The QR decomposition of the model_matrix
  arma::mat q, r;
  arma::vec w_vec = (mu/(1.0 + theta_times_mu));
  arma::vec w_sqrt_vec = sqrt(w_vec);
  // prepare matrices
  arma::mat weighted_model_matrix = model_matrix.each_col() % w_sqrt_vec;
  qr_econ(q, r, weighted_model_matrix);
  // Not actually quite the score vec, but related
  // See Dunn&Smyth GLM Book eq. 6.16
  arma::vec score_vec = (q.each_col() % w_sqrt_vec).t() * ((counts - mu) / mu);
  arma::vec step = solve(arma::trimatu(r), score_vec);
  return step;
}


template<class NumericType>
arma::vec fisher_scoring_qr_ridge_step(const arma::mat& model_matrix, const arma::Col<NumericType>& counts, const arma::colvec& mu,
                                       const arma::colvec& theta_times_mu, const double lambda1){
  // The QR decomposition of the model_matrix
  arma::mat q, r;
  int p = model_matrix.n_cols;
  arma::vec w_vec = (mu/(1.0 + theta_times_mu));
  arma::vec w_sqrt_vec = sqrt(w_vec);
  // Add rows for Ridge Regularization (see https://math.stackexchange.com/a/299508/492945)
  arma::vec extended_w_sqrt_vec = arma::join_cols(w_sqrt_vec, arma::ones(p));
  arma::mat ridge_helper = arma::eye(p, p) * sqrt(lambda1);
  arma::mat extended_model_matrix = arma::join_cols(model_matrix, ridge_helper);
  // prepare matrices
  arma::mat weighted_extended_model_matrix = extended_model_matrix.each_col() % extended_w_sqrt_vec;
  qr_econ(q, r, weighted_extended_model_matrix);
  // The last p rows of q are not needed because y' is zero anyway for them
  q.shed_rows(q.n_rows - p, q.n_rows - 1);
  // Not actually quite the score vec, but related
  // See Dunn&Smyth GLM Book eq. 6.16
  arma::vec score_vec = (q.each_col() % w_sqrt_vec).t() * ((counts - mu) / mu);
  arma::vec step = solve(arma::trimatu(r), score_vec);
  return step;
}



template<class NumericType>
arma::vec fisher_scoring_diagonal_step(const arma::mat& model_matrix, const arma::Col<NumericType>& counts,
                                       const arma::colvec& mu, const arma::colvec& theta_times_mu){
  arma::vec w_vec = (mu/(1.0 + theta_times_mu));
  // prepare matrices
  arma::mat weighted_model_matrix = model_matrix.each_col() % w_vec;
  arma::vec score_vec = weighted_model_matrix.t() * ((counts - mu) / mu);
  // This calculates the diag(Xˆt W X) efficiently. arma::sum(mat, 0) = colSums()
  arma::vec info_vec = arma::sum(arma::mat(arma::pow(model_matrix, 2)).each_col() % w_vec, 0).t();
  arma::vec step = score_vec / info_vec;
  return step;
}

#endif
