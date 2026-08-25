#ifndef FIT_BETA_H
#define FIT_BETA_H

#include <calc_helpers.h>
#include <deviance.h>
#include <fisher_scoring_steps.h>
#include <opt_max.h>

#include <utility> // for std::as_const

#include <Eigen/Dense>
template <class D> using EMB = Eigen::MatrixBase<D>;
using Eigen::ArrayXd;
using Eigen::VectorXd;

#include <LBFGSpp/include/LBFGS.h>

template <class D1, class D2>
inline double compute_pen_sum(const EMB<D1> &beta_min_ridge_t, const EMB<D2> &ridge_penalty_sq, const double n_samples) {
  return n_samples * (beta_min_ridge_t.transpose() * ridge_penalty_sq * beta_min_ridge_t)(0, 0);
}

/**
 * This method takes in a proposal for a step and checks if it actually
 * decreases the deviance of the model. If does not add it tries again
 * with half the step size, then a quarter and so on.
 *
 * If even after 100 steps the deviance (0.5^100 = 7.9e-31) has not decreased
 * it returns NaN.
 *
 * Note that the first two parameters are changed: beta_hat and mu_hat
 *
 * The function returns the new deviance.
 *
 */
template <class D1, class D2, class D3, class D4, class D5, class D6, class FSConf>
inline double decrease_deviance(
    // in-out params
    EMB<D1> &beta_hat, EMB<D2> &mu_hat,
    // const params
    const EMB<D3> &step, const EMB<D4> &model_matrix,
    // callback for computing gp deviance
    const FSConf &conf,
    // more const params
    const EMB<D5> &exp_off, const EMB<D6> &counts,
    // misc params
    const double theta, const double dev_old, const double tolerance, const double max_rel_mu_change) {
  beta_hat += step;

  double speeding_factor = 1.0;
  double dev = 0;
  const VectorXd mu_old = mu_hat;
  for (auto line_iter = 0;; line_iter++) {
    mu_hat = calculate_mu_mult<D2>(model_matrix, beta_hat, exp_off);
    dev = conf.gpd_sum(counts, mu_hat, theta, beta_hat);
    const double conv_test = std::abs(dev - dev_old) / (std::abs(dev) + 0.1);
    const double mu_rel_change = mu_hat.cwiseQuotient(mu_old).maxCoeff();
    if ((dev < dev_old && mu_rel_change < max_rel_mu_change) || conv_test < tolerance) {
      break; // while loop
    } else if (line_iter >= 100) {
      // speeding factor is very small, something is going wrong here
      dev = std::numeric_limits<double>::quiet_NaN();
      break; // while loop
    } else {
      // Halfing the speed
      speeding_factor = speeding_factor / 2.0;
      beta_hat = beta_hat - (step * speeding_factor);
    }
  }
  return dev;
}

class FisherScoreQR {
public:
  template <class D1, class D2, class D3>
  inline double gpd_sum(const EMB<D1> &y, const EMB<D2> &mu, const double theta, const EMB<D3> &beta_hat) const {
    return compute_gp_deviance_sum(y, mu, theta);
  }
  template <class D1, class D2, class D3, class D4, class D5>
  inline VectorXd fs_step(const EMB<D1> &model_matrix, const EMB<D2> &counts, const EMB<D3> &mu, const EMB<D4> &theta_times_mu,
                              const EMB<D5> &beta_hat) const {
    return fisher_scoring_qr_step(model_matrix, counts, mu, theta_times_mu);
  }
};

template <class D1, class D2> class FisherScoreQRwRidge {
private:
  const EMB<D1> &ridge_target_;
  const EMB<D2> &ridge_penalty_;
  const double &n_samples_;
  const MatrixXd ridge_penalty_sq_;

public:
  FisherScoreQRwRidge(const EMB<D1> &ridge_target, const EMB<D2> &ridge_penalty, const double &n_samples)
      : ridge_target_(ridge_target), ridge_penalty_(ridge_penalty), n_samples_(n_samples),
        ridge_penalty_sq_(ridge_penalty.transpose() * ridge_penalty) {}

  template <class D3, class D4, class D5>
  inline double gpd_sum(const EMB<D3> &y, const EMB<D4> &mu, const double theta, const EMB<D5> &beta_hat) const {
    return compute_gp_deviance_sum(y, mu, theta) + compute_pen_sum(beta_hat - ridge_target_, ridge_penalty_sq_, n_samples_);
  }
  template <class D4, class D5, class D6, class D7, class D8>
  inline VectorXd fs_step(const EMB<D4> &model_matrix, const EMB<D5> &counts, const EMB<D6> &mu, const EMB<D7> &theta_times_mu,
                              const EMB<D8> &beta_hat) const {
    return fisher_scoring_qr_ridge_step(model_matrix, counts, mu, theta_times_mu, ridge_penalty_, ridge_target_, beta_hat);
  }
};

class FisherScoreDiagApprox {
public:
  template <class D1, class D2, class D3>
  inline double gpd_sum(const EMB<D1> &y, const EMB<D2> &mu, const double theta, const EMB<D3> &beta_hat) const {
    return compute_gp_deviance_sum(y, mu, theta);
  }
  template <class D1, class D2, class D3, class D4, class D5>
  inline VectorXd fs_step(const EMB<D1> &model_matrix, const EMB<D2> &counts, const EMB<D3> &mu, const EMB<D4> &theta_times_mu,
                              const EMB<D5> &beta_hat) const {
    return fisher_scoring_diagonal_step(model_matrix, counts, mu, theta_times_mu);
  }
};

template <class D1, class D2, class D3, class FSConf> class BetaOptim {
private:
  const FSConf &conf_;
  const EMB<D1> &model_matrix_;
  const EMB<D2> &counts_;
  const EMB<D3> &exp_off_;
  const double &theta_;
  const double &eps_;
  inline double fn_(const VectorXd &beta_row) const {
    return conf_.gpd_sum(counts_, calculate_mu_mult<VectorXd>(model_matrix_, beta_row, exp_off_), theta_, beta_row);
  }
  inline void grad_(VectorXd &x, VectorXd &grad) const {
    auto x_begin = x.begin();
    const auto x_end = x.end();
    auto grad_begin = grad.begin();
    for (; x_begin != x_end; ++x_begin, ++grad_begin) {
      const auto x_curr = *x_begin;
      *x_begin = x_curr + eps_;
      const auto v1 = fn_(x);
      *x_begin = x_curr - eps_;
      const auto v2 = fn_(x);
      *grad_begin = (v1 - v2) / (2 * eps_);
      *x_begin = x_curr;
    }
  }

public:
  BetaOptim(const FSConf &conf, const EMB<D1> &model_matrix, const EMB<D2> &counts, const EMB<D3> &exp_off, const double &theta, const double &eps)
      : conf_(conf), model_matrix_(model_matrix), counts_(counts), exp_off_(exp_off), theta_(theta), eps_(eps) {};
  inline double operator()(const VectorXd &x, VectorXd &grad) const {
    // we need to make a copy of input since it is provided as a const reference
    VectorXd x_c = x;
    grad_(x_c, grad);
    return fn_(x);
  }
};

template <class D1, class D2, class D3, class FSConf>
inline void fitBeta_FS_optim_step(
    // in-out params
    VectorXd &beta_out, double &dev_out, int &iters_out,
    // const params
    const EMB<D1> &model_matrix, const EMB<D2> &counts, const EMB<D3> &exp_off,
    // callback for optim
    const FSConf &conf,
    // other params
    const double theta, const int max_iter) {
  LBFGSpp::LBFGSParam params;
  params.max_iterations = max_iter;
  
  // using delta-based stopping condition instead of gradient value
  params.epsilon = 0;
  params.epsilon_rel = 0;
  params.past = 1;
  params.delta = 1e-12;

  params.max_step = 32.; 

  // nocedal-wright linesearch often fails to find a sufficiently decreasing value, using a different one instead
  // moreover, we don't want to enforce the wolfe condition here given potential instability
  // of our numerically approximated derivative, and nocedal-wright is only compatible w/ the wolfe condition
  // as such, we enforce only the armijo (sufficient decrease) condition, use a plain bracketing line search, and 
  // set a very generous sufficient decrease line search tolerance
  params.linesearch = LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
  params.max_linesearch = 1000;
  params.ftol = 1e-8;
  LBFGSpp::LBFGSSolver<double, LBFGSpp::LineSearchBracketing> solver(params);

  // borrowing the default used for this by the (unstable) Eigen::NumericalDiff module
  const auto eps = SQRT_DBL_EPS;
  const BetaOptim f(conf, model_matrix, counts, exp_off, theta, eps);

  double fx;
  iters_out = solver.minimize(f, beta_out, fx);

  if (iters_out == max_iter) {
    beta_out.fill(NAN);
    dev_out = NAN;
    return;
  }
  dev_out = fx;
}

const double FS_CLAMP_HI = 1e8;
template <class D1, class D2, class D3, class D4, class FSConf>
inline void fitBeta_FS_internal_step(
    // in-out params
    EMB<D1> &beta_out, double &dev_out, int &iters_out,
    // const params
    const EMB<D2> &model_matrix, const EMB<D3> &counts, const EMB<D4> &exp_off,
    // callbacks for steps
    const FSConf &conf,
    // other params
    const double theta, const double tolerance, const double max_rel_mu_change, const int max_iter, const bool try_recov_w_optim) {

  // Init beta and mu
  VectorXd beta_hat = beta_out;
  VectorXd mu_hat = calculate_mu_mult<VectorXd>(model_matrix, beta_hat, exp_off);

  if (beta_hat.array().isNaN().any() || std::isnan(theta)) {
    beta_out.fill(NAN);
    dev_out = NAN;
    iters_out = 0;
    return;
  }

  // Init deviance
  double dev_old = conf.gpd_sum(counts, mu_hat, theta, beta_hat);
  int iter = 0;
  for (; iter < max_iter; iter++) {
    // Find good direction to optimize beta
    // but clamp at reasonable size to avoid huge betas
    // TODO: find out what's causing fisher steps to be on the scale of > 1e8
    const VectorXd step = conf.fs_step(model_matrix, counts, mu_hat, theta * mu_hat, beta_hat).cwiseMin(FS_CLAMP_HI);
    // Find step size that actually decreases the deviance
    const double dev = decrease_deviance(beta_hat, mu_hat, step, model_matrix, conf, exp_off, counts, theta, dev_old, tolerance, max_rel_mu_change);

    const double conv_test = std::abs(dev - dev_old) / (std::abs(dev) + 0.1);
    dev_old = dev;
    if (std::isnan(conv_test)) {
      // This should not happen
      beta_hat.fill(NAN);
      iter = max_iter;
      break;
    }
    if (conv_test < tolerance) {
      break;
    }
  }

  if (try_recov_w_optim && (iter == max_iter)) {
    beta_hat = beta_out; // re-copy since from reference since values have been overwritten
    fitBeta_FS_optim_step(beta_hat, dev_old, iter, model_matrix, counts, exp_off, conf, theta, max_iter);
  }

  beta_out = beta_hat;
  dev_out = dev_old;
  iters_out = iter;
}

template <class D1, class D2> inline double fitBeta_NR_optim_step(const EMB<D1> &counts, const EMB<D2> &off, const double theta) {
  return optimize_fmax([&counts = std::as_const(counts), &off = std::as_const(off), &theta = std::as_const(theta)](double beta_hat) -> double {
    double out = 0.0;
    // iterate through counts & offset simultaneously
    auto counts_begin = counts.cbegin();
    auto off_begin = off.cbegin();
    const auto counts_end = counts.cend();
    // we known that counts and off are of the same size, so we only check one vector
    for (; counts_begin != counts_end; ++counts_begin, ++off_begin) {
      out += dnbinom_impl(*counts_begin, beta_hat, *off_begin, theta);
    }
    return out;
  });
}

template <class D1, class D2>
inline void fitBeta_NR_internal_step(
    // in-out parameters
    double &beta_out, double &dev_out, int &iters_out,
    // const params
    const EMB<D1> &counts, const EMB<D2> &off,
    // other params
    const double theta, const double tolerance, const int max_iter) {

  double beta_hat = beta_out;
  if (std::isnan(beta_hat) || std::isnan(theta)) {
    beta_out = NAN;
    dev_out = NAN;
    iters_out = 0;
    return;
  }
  if (counts.isZero()) {
    beta_out = -INFINITY;
    dev_out = 0;
    iters_out = 0;
    return;
  }

  // Newton-Raphson
  int iter = 0;
  for (; iter < max_iter; iter++) {

    const ArrayXd mu = (off.array() + beta_hat).exp();
    const ArrayXd denom = 1.0 + (mu * theta).array();

    const double dl = ((counts.array() - mu) / denom).sum();
    const double ddl = (mu * (1.0 + (counts.array() * theta)) / denom / denom).sum();

    const double step = dl / ddl;
    beta_hat += step;
    if ((std::abs(step) < tolerance) || (std::isnan(beta_hat))) {
      break;
    }
  }

  if ((iter == max_iter) || std::isnan(beta_hat)) {
    beta_hat = fitBeta_NR_optim_step(counts, off, theta);
  }

  beta_out = beta_hat;
  dev_out = compute_gp_deviance_sum(counts, (off.array() + beta_hat).exp().matrix(), theta);
  iters_out = iter;
}

#endif