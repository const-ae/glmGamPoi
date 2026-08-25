#ifndef CALC_HELPERS_H
#define CALC_HELPERS_H

#include <cfloat>
#include <cmath>
const double SQRT_DBL_EPS = std::sqrt(DBL_EPSILON);

#include <Eigen/Dense>
template <class D> using EMB = Eigen::MatrixBase<D>;
const double C_MU_MULT_LO = 1e-50;
const double C_MU_MULT_HI = 1e50;

template <class D4, class D1, class D2, class D3>
inline D4 calculate_mu_mult(const EMB<D1> &model_matrix, const EMB<D2> &beta_hat, const EMB<D3> &exp_off) {
  return ((model_matrix * beta_hat).array().exp() * exp_off.array()).unaryExpr([](double exp_x) {
    // clamp to avoid large / zero values
    return std::clamp(exp_x, C_MU_MULT_LO, C_MU_MULT_HI);
  });
}

template <class D4, class D1, class D2, class D3>
inline D4 calculate_mu_add(const EMB<D1> &model_matrix, const EMB<D2> &beta_hat, const EMB<D3> &off) {
  return ((model_matrix * beta_hat).array() + off.array()).exp();
}

#include <Rmath.h>

inline double dnbinom_impl(const double y, const double beta_hat, const double off, const double theta) {
  return dnbinom_mu(y, 1.0 / theta, std::exp(beta_hat + off), 1);
}
inline double lgamma_impl(const double x) { return lgammafn(x); } // we do not use std::lgamma because it is explicitly NOT thread safe
inline double digamma_impl(const double x) { return digamma(x); }
inline double trigamma_impl(const double x) { return trigamma(x); }

#endif