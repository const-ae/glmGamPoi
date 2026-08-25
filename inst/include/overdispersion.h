#ifndef OVERDISP_H
#define OVERDISP_H

#include <calc_helpers.h>

#include <cfloat>  // for DBL_EPSILON
#include <utility> // for std::as_const

#include <Eigen/Dense>
template <class D> using EMB = Eigen::MatrixBase<D>;
using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::Vector;
using Eigen::VectorXd;

// This correction factor is necessary to avoid estimates of
// theta that are basically +Inf. The problem is that for
// some combination of the y, mu, and X the term
// lgamma(1/theta) and the log(det(t(X) %*% W %*% X))
// with W = diag(1/(1/mu + theta)) canceled each other
// exactly out for large theta.
const double cr_correction_factor = 0.99;

//--------------------------------------------------------------------------------------------------//
// The following code was originally copied from https://github.com/mikelove/DESeq2/blob/master/src/DESeq2.cpp
// I adapted it to the needs of this project by:
//  * renaming alpha -> theta for consitency
//  * removing the part for the prior on theta
//  * renaming x -> model_matrix
//  * additional small changes
//  * adding capability to calculate digamma/trigamma only
//    on unique counts

/*
 * DESeq2 C++ functions
 *
 * Author: Michael I. Love, Constantin Ahlmann-Eltze
 * Last modified: May 21, 2020
 * License: LGPL (>= 3)
 *
 * Note: The canonical, up-to-date DESeq2.cpp lives in
 * the DESeq2 library, the development branch of which
 * can be viewed here:
 *
 * https://github.com/mikelove/DESeq2/blob/master/src/DESeq2.cpp
 */

// this function returns the log posterior of dispersion parameter alpha, for negative binomial variables
// given the counts y, the expected means mu, the design matrix x (used for calculating the Cox-Reid adjustment),
// and the parameters for the normal prior on log alpha
template <class D1, class D2, class D3, class D4, class D5>
inline double conventional_loglikelihood_fast_impl(const EMB<D1> &y, const EMB<D2> &mu, const double log_theta, const EMB<D3> &model_matrix,
                                                   const bool do_cr_adj, const EMB<D4> &unique_counts, const EMB<D5> &count_frequencies) {
  const double theta = std::exp(log_theta);
  double cr_term = 0.0;
  if (do_cr_adj) {
    const auto c = model_matrix.cols();
    MatrixXd b(c, c);
    b.noalias() = model_matrix.transpose() * (model_matrix.array().colwise() * (mu.cwiseInverse().array() + theta).cwiseInverse()).matrix();
    // cr_term = -0.5 * log(det(b)) * cr_correction_factor;

    // compute det(B) using LU decomposition

    // eigen computes LU decomposition w/ a _unit_ lower triangular L and upper triangular U,
    // representing the resulting L/U matrices in a single square matrix (returned by the matrixLU method)
    // e.g. leaving the diagonal of L represented implicitly (since it is known to be unit values)
    //
    // as such, when computing the determinant ( det(X) = prod(diag(L)) * prod(diag(U)) )
    // we ignore the diagonal of L given that it has unit values, and use the diagonal of matrixLU (which is in turn the diagonal of U)
    // moreover, we clamp at 1e-50 to avoid generating infinities in the final sum
    cr_term = -0.5 * b.partialPivLu().matrixLU().diagonal().cwiseAbs().cwiseMax(1e-50).array().log().sum() * cr_correction_factor;
  }

  const double theta_neg1 = 1.0 / theta;

  double lgamma_term;
  // If summarized counts are available use those to calculate sum(lgamma(y + theta_neg1))
  if (unique_counts.size() > 0 && unique_counts.size() == count_frequencies.size()) {
    lgamma_term = (count_frequencies.array() * (unique_counts.array() + theta_neg1).unaryExpr(std::ref(lgamma_impl))).sum();
  } else {
    lgamma_term = (y.array() + theta_neg1).unaryExpr(std::ref(lgamma_impl)).sum();
  }
  lgamma_term -= y.size() * lgamma_impl(theta_neg1);

  const double ll_part = (-(y.array() + theta_neg1) * (mu.array() + theta_neg1).log()).sum() - (y.size() * theta_neg1 * log_theta);
  return lgamma_term + ll_part + cr_term;
}

// this function returns the derivative of the log posterior with respect to the log of the
// dispersion parameter alpha, given the same inputs as the previous function
template <class D1, class D2, class D3, class D4, class D5>
inline double conventional_score_function_fast_impl(const EMB<D1> &y, const EMB<D2> &mu, const double log_theta, const EMB<D3> &model_matrix,
                                                    const bool do_cr_adj, const EMB<D4> &unique_counts, const EMB<D5> &count_frequencies) {
  const double theta = std::exp(log_theta);
  double cr_term = 0.0;
  if (do_cr_adj) {
    const ArrayXd w_diag = (mu.cwiseInverse().array() + theta).cwiseInverse();
    const auto c = model_matrix.cols();
    MatrixXd b(c, c), db(c, c);
    b.noalias() = model_matrix.transpose() * (model_matrix.array().colwise() * w_diag).matrix();
    db.noalias() = model_matrix.transpose() * (model_matrix.array().colwise() * (-w_diag.square())).matrix();

    // The diag(1e-6) protects against singular matrices
    MatrixXd b_inv = MatrixXd::Identity(b.rows(), b.cols());
    b.diagonal().array() += 1e-6;
    b.llt().solveInPlace(b_inv);

    cr_term = -0.5 * (b_inv * db).trace() * cr_correction_factor;
  }

  const double theta_neg1 = 1.0 / theta;
  double digamma_term, max_y, sum_y, sum_prod_y;
  // If summarized counts are available use those to calculate sum(digamma(y + theta_neg1))
  if (unique_counts.size() > 0 && unique_counts.size() == count_frequencies.size()) {
    digamma_term = (count_frequencies.array() * (unique_counts.array() + theta_neg1).unaryExpr(std::ref(digamma_impl))).sum();
    max_y = unique_counts.maxCoeff();
    const ArrayXd prod_freqs = (count_frequencies.array() * unique_counts.array());
    sum_y = prod_freqs.sum();
    sum_prod_y = (prod_freqs * (unique_counts.array() - 1)).sum();
  } else {
    digamma_term = (y.array() + theta_neg1).unaryExpr(std::ref(digamma_impl)).sum();
    max_y = y.maxCoeff();
    sum_y = y.sum();
    sum_prod_y = ((y.array() - 1) * y.array()).sum();
  }
  // This approximation is based on the fact that for large x
  // (sum(digamma(y + x))  - length(y) * digamma(x)) * x \approx sum(y)
  // Due to numerical imprecision the digamma_term reaches sum(y) sometimes
  // quicker than the ll_term, thus I subtract the first term of the
  // Laurent series expansion at x -> inf
  const double corr = theta_neg1 > 1e5 ? sum_prod_y / (2 * theta_neg1) : 0.0;
  if (max_y * 1e6 < theta_neg1) {
    digamma_term = sum_y - corr;
  } else {
    digamma_term -= y.size() * digamma_impl(theta_neg1);
    digamma_term *= theta_neg1;
    const double x = sum_y - corr;
    digamma_term = std::min(digamma_term, x);
  }
  const double ll_part = (mu.unaryExpr([&theta = std::as_const(theta), &theta_neg1 = std::as_const(theta_neg1)](double mu) -> double {
                              const double mu_theta = mu * theta;
                              if (mu_theta < 1e-10) {
                                return mu_theta * mu_theta * (1 / (1 + mu_theta) - 0.5);
                              } else if (mu_theta < 1e-4) {
                                // The bounds are based on the Taylor expansion of log(1 + x) for x = 0.
                                const double inv = 1 / (1 + mu_theta);
                                const double upper_bound = mu_theta * mu_theta * inv;
                                const double lower_bound = mu_theta * mu_theta * (inv - 0.5);
                                const double suggest = (std::log1p(mu_theta) - mu / (mu + theta_neg1));
                                return std::clamp(suggest, lower_bound, upper_bound);
                              } else {
                                return std::log1p(mu_theta) - mu / (mu + theta_neg1);
                              }
                            }).array() +
                          (y.array() / (mu.array() + theta_neg1)))
                             .sum() *
                         theta_neg1;

  return ll_part - digamma_term + cr_term * theta;
}

// this function returns the second derivative of the log posterior with respect to the log of the
// dispersion parameter alpha, given the same inputs as the previous function
template <class D1, class D2, class D3, class D4, class D5>
inline double conventional_deriv_score_function_fast_impl(const EMB<D1> &y, const EMB<D2> &mu, const double log_theta, const EMB<D3> &model_matrix,
                                                          const bool do_cr_adj, const EMB<D4> &unique_counts, const EMB<D5> &count_frequencies) {
  const double theta = std::exp(log_theta);
  double cr_term = 0.0;
  double cr_term2 = 0.0;
  if (do_cr_adj) {
    const ArrayXd w_diag = (mu.cwiseInverse().array() + theta).cwiseInverse();

    const auto c = model_matrix.cols();
    MatrixXd b(c, c), db(c, c), d2b(c, c);
    b.noalias() = model_matrix.transpose() * (model_matrix.array().colwise() * w_diag).matrix();
    db.noalias() = model_matrix.transpose() * (model_matrix.array().colwise() * (-w_diag.square())).matrix();
    d2b.noalias() = model_matrix.transpose() * (model_matrix.array().colwise() * (2 * w_diag.cube())).matrix();

    // The diag(1e-6) protects against singular matrices
    MatrixXd b_inv = MatrixXd::Identity(b.rows(), b.cols());
    b.diagonal().array() += 1e-6;
    b.llt().solveInPlace(b_inv);

    const MatrixXd d_i_db = b_inv * db;
    const double ddetb = d_i_db.trace();
    const double ddetb_pow2 = ddetb * ddetb;
    const double d2detb = ((ddetb_pow2 - (d_i_db * d_i_db).trace() + (b_inv * d2b).trace()));
    cr_term = (0.5 * ddetb_pow2 - 0.5 * d2detb) * cr_correction_factor;
    cr_term2 = -0.5 * ddetb * cr_correction_factor;
  }

  const double theta_neg1 = 1.0 / theta;
  const double theta_neg2 = theta_neg1 * theta_neg1;

  double digamma_term, trigamma_term;
  // If summarized counts are available use those to calculate sum(digamma()) and sum(trigamma())
  if (unique_counts.size() > 0 && unique_counts.size() == count_frequencies.size()) {
    const ArrayXd counts_p_theta_neg1 = unique_counts.array() + theta_neg1;
    digamma_term = (count_frequencies.array() * counts_p_theta_neg1.unaryExpr(std::ref(digamma_impl))).sum();
    trigamma_term = (count_frequencies.array() * counts_p_theta_neg1.unaryExpr(std::ref(trigamma_impl))).sum();
    trigamma_term *= theta_neg2;
  } else {
    const ArrayXd counts_p_theta_neg1 = y.array() + theta_neg1;
    digamma_term = counts_p_theta_neg1.unaryExpr(std::ref(digamma_impl)).sum();
    trigamma_term = theta_neg2 * counts_p_theta_neg1.unaryExpr(std::ref(trigamma_impl)).sum();
  }
  digamma_term -= y.size() * digamma_impl(theta_neg1);
  trigamma_term -= theta_neg2 * y.size() * trigamma_impl(theta_neg1);

  const ArrayXd mu_theta = mu.array() * theta;
  const double ll_part_1 = (mu_theta.log1p() + ((y - mu).array() / (mu.array() + theta_neg1))).sum();
  const double ll_part_2 = ((mu.array() >= 1e100).any() && !mu.isZero())
                               // original expression was: (mu^2 * theta + y) / ((1 + mu * theta) * (1 + mu * theta))
                               // however, to avoid computing mu^2 when mu is very large (which can generate infinities/NANs)
                               // we divide both sides of the quotient by mu, yielding: (mu * theta + y / mu) / ((1 / mu + theta) * (1 + mu * theta))
                               ? ((mu_theta + (y.array() / mu.array())) / ((mu_theta + 1) * (theta + mu.cwiseInverse().array()))).sum()
                               // if mu values are outside the bounds where problems can occur, use the original expression
                               : ((mu.array() * mu_theta + y.array()) / (1 + mu_theta).array().square()).sum();

  const double ll_part = -2 * theta_neg1 * (ll_part_1 - digamma_term) + (ll_part_2 + trigamma_term);

  double res = ll_part + cr_term * (theta * theta) + (ll_part_1 - digamma_term) * theta_neg1 + cr_term2 * theta;
  return res;
}

// not using Eigen::MatrixBase as base class here to have access to resize method
template <class V1, class V2, class V3> inline void make_map_if_small(V1 &unique_counts_out, V2 &freqs_out, const V3 &x, const size_t stop_if_larger) {
  std::unordered_map<double, double> count_tab;
  count_tab.reserve(stop_if_larger);
  for (auto v : x) {
    ++count_tab[v];
    if (count_tab.size() > stop_if_larger) {
      return;
    }
  }

  // assuming whatever vector type V1 and V2 are, they have a resize method
  unique_counts_out.resize(count_tab.size());
  freqs_out.resize(count_tab.size());

  size_t i = 0;
  for (auto p : count_tab) {
    unique_counts_out[i] = p.first;
    freqs_out[i] = p.second;
    i++;
  }
}

const double NR_STEP_MIN = -32.;
const double NR_STEP_MAX = 32.;
const double GA_EPS = SQRT_DBL_EPS;
const double GA_STEP_MAX = 2.;
const double GA_STEP_EXP_MAX = 2048.;
template <class D1, class D2, class D3, class D4, class D5>
inline bool estimate_theta(double &log_theta_out, int &iters_out, const EMB<D1> &y, const EMB<D2> &mu_vec, const EMB<D3> &model_matrix,
                           const bool do_cox_reid_adjustment, const double tol, const EMB<D4> &unique_counts, const EMB<D5> &count_frequencies) {
  const auto max_iter = iters_out;

  double step = NAN;
  bool use_analytic_grad = true;

  double cur_val =
      conventional_loglikelihood_fast_impl(y, mu_vec, log_theta_out, model_matrix, do_cox_reid_adjustment, unique_counts, count_frequencies);
  double suff_inc = tol * std::abs(cur_val);

  const auto ls_inc = [&step, &cv = std::as_const(cur_val), &lt = std::as_const(log_theta_out), &y = std::as_const(y), &mu = std::as_const(mu_vec),
                       &mm = std::as_const(model_matrix), &cr = std::as_const(do_cox_reid_adjustment), &uc = std::as_const(unique_counts),
                       &cf = std::as_const(count_frequencies)](double &cand, double &new_val) -> bool {
    if (new_val < cv) {
      bool switch_meth = false;
      do {
        step *= 0.5;

        if (std::abs(step) < DBL_EPSILON) {
          return true;
        }

        cand = lt + step;
        new_val = conventional_loglikelihood_fast_impl(y, mu, cand, mm, cr, uc, cf);
      } while (new_val < cv);
      step *= 0.95;
    }
    return false;
  };

  iters_out = 0;
  for (; iters_out < max_iter; iters_out++) {
    const double grad =
        conventional_score_function_fast_impl(y, mu_vec, log_theta_out, model_matrix, do_cox_reid_adjustment, unique_counts, count_frequencies);
    const double abs_hess = std::abs(conventional_deriv_score_function_fast_impl(y, mu_vec, log_theta_out, model_matrix, do_cox_reid_adjustment,
                                                                                 unique_counts, count_frequencies));
    step = std::clamp(grad / abs_hess, NR_STEP_MIN, NR_STEP_MAX);

    if (std::isnan(step)) {
      return false;
    }

    double cand = log_theta_out + step;
    double new_val = conventional_loglikelihood_fast_impl(y, mu_vec, cand, model_matrix, do_cox_reid_adjustment, unique_counts, count_frequencies);

    const bool ls_failed = ls_inc(cand, new_val);
    if (ls_failed) {
      // initializing the step value for the numeric derivative driven steps
      // we "cheat" by just using the existing hessian value
      // (with a clamp since this should be an adjustment)
      step = std::abs(2. / abs_hess);
      step = std::min(GA_STEP_MAX, step);
      break;
    }

    log_theta_out = cand;
    if ((new_val - cur_val) < suff_inc) {
      return true;
    }
    cur_val = new_val;
    suff_inc = tol * std::abs(cur_val);
  }
  for (; iters_out < max_iter; iters_out++) {
    const double step_bwd = conventional_loglikelihood_fast_impl(y, mu_vec, log_theta_out - GA_EPS, model_matrix, do_cox_reid_adjustment,
                                                                 unique_counts, count_frequencies);
    const double step_fwd = conventional_loglikelihood_fast_impl(y, mu_vec, log_theta_out + GA_EPS, model_matrix, do_cox_reid_adjustment,
                                                                 unique_counts, count_frequencies);

    // addition to insure sufficient increase (otherwise some steps can be ill-behaved)
    const double cur_val_inc = cur_val + suff_inc;

    if (step_bwd > cur_val_inc) {
      step = -std::abs(step);
    } else if (step_fwd > cur_val_inc) {
      step = std::abs(step);
      const double exp_clamp = std::log1p(GA_STEP_EXP_MAX * std::exp(-log_theta_out));
      step = std::min(step, exp_clamp);
    } else {
      // no sufficient increase in either direction, bail
      return true;
    }

    if (std::isnan(step)) {
      return false;
    }

    double cand = log_theta_out + step;
    double new_val = conventional_loglikelihood_fast_impl(y, mu_vec, cand, model_matrix, do_cox_reid_adjustment, unique_counts, count_frequencies);

    const bool ls_failed = ls_inc(cand, new_val);
    if (ls_failed) {
      return false;
    }

    log_theta_out = cand;
    if ((new_val - cur_val) < suff_inc) {
      return true;
    }
    cur_val = new_val;
    suff_inc = tol * std::abs(cur_val);
  }

  return false;
}

template <class D1, class D2, class D3, class CharStarAssignable>
inline void overdispersion_mle_NR_impl(double &est_out, int &iters_out, CharStarAssignable &msg_out, const EMB<D1> &y, const EMB<D2> &mu_vector,
                                       const EMB<D3> &model_matrix, const bool do_cox_reid_adjustment, const int max_iter, const double tol) {
  if (y.isZero()) {
    est_out = 0.;
    iters_out = 0;
    msg_out = "All counts y are 0";
    return;
  }

  // replace zeroes w/ 1e-6 because zeroes can cause issues
  const VectorXd mean_vec_clamp = mu_vector.unaryExpr([](const auto e) { return e < 1e-128 ? 1e-6 : e; });

  VectorXd unique_counts, count_frequencies;
  make_map_if_small(unique_counts, count_frequencies, y, y.size() / 2);

  const double far_left_value =
      conventional_score_function_fast_impl(y, mean_vec_clamp, -18.42, model_matrix, do_cox_reid_adjustment, unique_counts, count_frequencies);

  if (far_left_value < 0) {
    est_out = 0.;
    iters_out = 0;
    msg_out = "Even for very small theta, no maximum identified";
    return;
  }

  const double mu = y.mean();
  double theta_init = ((y.array() - mu).square().mean() - mu) / (mu * mu);
  if (std::isnan(theta_init) || (theta_init < 1e-8) || (theta_init > 1e8)) {
    theta_init = 0.5;
  }
  const double log_theta_init = std::log(theta_init);

  double log_theta_out = log_theta_init;
  iters_out = max_iter;

  double succ =
      estimate_theta(log_theta_out, iters_out, y, mean_vec_clamp, model_matrix, do_cox_reid_adjustment, tol, unique_counts, count_frequencies);

  // if failed and cox-reid adjustment is used: try w/o cox-reid adjustment
  if (do_cox_reid_adjustment && (!succ)) {
    msg_out = "Estimated overdispersion w/o cox-reid adjustment";
    log_theta_out = log_theta_init;
    iters_out = max_iter;

    succ = estimate_theta(log_theta_out, iters_out, y, mean_vec_clamp, model_matrix, false, tol, unique_counts, count_frequencies);
  }

  est_out = std::exp(log_theta_out);
}

#endif