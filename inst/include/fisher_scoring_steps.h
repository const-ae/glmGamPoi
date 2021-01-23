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


/**
 * Ridge ression penalizes large values of beta:
 *    1/N \Sum (y - X b)^2 + lambda^2 * \Sum b^2
 *
 * Finding the optimal b that balances those two errors can be expressed using the normal equations
 *    b = (X^t X + diag(lambda^2))^-1 X^t y
 *
 * However, this function does not compute b directly, but the step s to go from b^{(r)}
 * to b^{(r+1)} = b^{(r)} + s.
 * We find s by solving
 *   argmin { y' - (X' b^{r} + X' s) },
 * where y' = [y 0]^t and X' = [X diag(lambda)]^t.
 * We can rearrange the equation above for a given b
 *   argmin { [y-(Xb) 0-lambda*b] - X' s}.
 *
 * For numerically stability we apply X' = Q R and solve
 *   Q^t [y-(Xb) 0-lambda*b] = R s
 *
 * For the actual implementation below, we need to keep the weighting w = mu / (1 + mu * theta)
 * in mind. However, for clarity, I skipped w in the above derivation.
 */
template<class NumericType>
arma::vec fisher_scoring_qr_ridge_step(const arma::mat& model_matrix, const arma::Col<NumericType>& counts, const arma::colvec& mu,
                                       const arma::colvec& theta_times_mu, const arma::mat& ridge_penalty,  const arma::colvec& ridge_target, const arma::colvec& beta){
  // The QR decomposition of the model_matrix
  arma::mat q, r;
  int extra = ridge_penalty.n_rows;
  arma::vec w_vec = (mu/(1.0 + theta_times_mu));
  arma::vec w_sqrt_vec = sqrt(w_vec);
  // Add rows for Ridge Regularization (see https://math.stackexchange.com/a/299508/492945)
  arma::mat ridge_helper = sqrt(model_matrix.n_rows) *  ridge_penalty;
  arma::mat extended_model_matrix = arma::join_cols(model_matrix, ridge_helper);

  arma::vec extended_w_sqrt_vec = arma::join_cols(w_sqrt_vec, arma::ones(extra));
  arma::vec extended_working_resid =  arma::join_cols((counts - mu) / mu, - ridge_helper * (beta - ridge_target));
  // prepare matrices
  arma::mat weighted_extended_model_matrix = extended_model_matrix.each_col() % extended_w_sqrt_vec;
  qr_econ(q, r, weighted_extended_model_matrix);

  // Not actually quite the score vec, but related
  // See Dunn&Smyth GLM Book eq. 6.16
  arma::vec score_vec = (q.each_col() % extended_w_sqrt_vec).t() * extended_working_resid;
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
