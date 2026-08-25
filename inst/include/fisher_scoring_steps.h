#ifndef FISHER_SCORING_STEPS_H
#define FISHER_SCORING_STEPS_H

#include <Eigen/Dense>
template <class D> using EMB = Eigen::MatrixBase<D>;
using Eigen::ArrayXd;
using Eigen::HouseholderQR;
using Eigen::MatrixXd;
using Eigen::TriangularView;
using Eigen::VectorX;
using Eigen::VectorXd;

template <class D1, class D2, class D3, class D4>
inline VectorXd fisher_scoring_qr_step(const EMB<D1> &model_matrix, const EMB<D2> &counts, const EMB<D3> &mu, const EMB<D4> &theta_times_mu) {
  // The QR decomposition of the model_matrix
  const ArrayXd w_sqrt_vec = (mu.array() / (1.0 + theta_times_mu.array())).sqrt();
  const HouseholderQR qr = (model_matrix.array().colwise() * w_sqrt_vec).matrix().householderQr();

  // materialize q and r matrix views
  // HouseholderSequence defines the multiplication operator, so we extract a "thin" Q by multiplying with the identity matrix;
  const MatrixXd q = qr.householderQ() * MatrixXd::Identity(model_matrix.rows(), model_matrix.cols());
  // r is stored in the top left corner of the matrixQR() matrix, of which we only need the triangular view
  const TriangularView r = qr.matrixQR().topLeftCorner(model_matrix.cols(), model_matrix.cols()).template triangularView<Eigen::Upper>();

  // Not actually quite the score vec, but related
  // See Dunn&Smyth GLM Book eq. 6.16
  const VectorXd score_vec = (q.array().colwise() * w_sqrt_vec).matrix().transpose() * (counts - mu).cwiseQuotient(mu);
  return r.solve(score_vec);
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

template <class D1, class D2, class D3, class D4, class D5, class D6, class D7>
inline VectorXd fisher_scoring_qr_ridge_step(const EMB<D1> &model_matrix, const EMB<D2> &counts, const EMB<D3> &mu, const EMB<D4> &theta_times_mu,
                                             const EMB<D5> &ridge_penalty, const EMB<D6> &ridge_target, const EMB<D7> &beta_hat) {
  const size_t extra = ridge_penalty.rows();

  // the sqrt(n) is important to scale the ridge_penalty by the number of samples
  // i.e. pen = sum dev(y, mu) + n * b^t Lambda^t Lambda b^t
  const MatrixXd ridge_helper = std::sqrt(double(model_matrix.rows())) * ridge_penalty;

  // Add rows for Ridge Regularization (see https://math.stackexchange.com/a/299508/492945)
  MatrixXd extended_model_matrix(model_matrix.rows() + extra, model_matrix.cols());
  extended_model_matrix << model_matrix, ridge_helper;
  ArrayXd extended_w_sqrt_vec(mu.size() + extra);
  extended_w_sqrt_vec << (mu.array() / (1.0 + theta_times_mu.array())).sqrt(), VectorXd::Constant(extra, 1.0);
  VectorXd extended_working_resid(counts.size() + extra);
  extended_working_resid << (counts - mu).cwiseQuotient(mu), ridge_helper * (ridge_target - beta_hat);

  // The QR decomposition of the model_matrix
  const HouseholderQR qr = (extended_model_matrix.array().colwise() * extended_w_sqrt_vec).matrix().householderQr();
  // materialize q and r matrix views
  // HouseholderSequence defines the multiplication operator, so we extract a "thin" Q by multiplying with the identity matrix;
  const MatrixXd q = qr.householderQ() * MatrixXd::Identity(extended_model_matrix.rows(), extended_model_matrix.cols());
  // r is stored in the top left corner of the matrixQR() matrix, of which we only need the triangular view
  const TriangularView r =
      qr.matrixQR().topLeftCorner(extended_model_matrix.cols(), extended_model_matrix.cols()).template triangularView<Eigen::Upper>();

  // Not actually quite the score vec, but related
  // See Dunn&Smyth GLM Book eq. 6.16
  const VectorXd score_vec = (q.array().colwise() * extended_w_sqrt_vec).matrix().transpose() * extended_working_resid;
  return r.solve(score_vec);
}

template <class D1, class D2, class D3, class D4>
inline VectorXd fisher_scoring_diagonal_step(const EMB<D1> &model_matrix, const EMB<D2> &counts, const EMB<D3> &mu, const EMB<D4> &theta_times_mu) {
  const ArrayXd w_vec = mu.array() / (1.0 + theta_times_mu.array());

  return ((model_matrix.array().colwise() * w_vec).transpose().matrix() * (counts - mu).cwiseQuotient(mu))
      .cwiseQuotient((model_matrix.array().square().colwise() * w_vec).colwise().sum().matrix().transpose());
}

#endif
