#ifndef DEVIANCE_H
#define DEVIANCE_H

#include <Eigen/Dense>
template <class D> using EMB = Eigen::MatrixBase<D>;
using Eigen::MatrixX;
using Eigen::MatrixXd;

template <class T> inline int sgn(const T val) { return (T(0) < val) - (val < T(0)); }

inline double compute_gp_deviance(const double y, const double mu, const double theta) {
  if (theta < 1e-6) {
    // If theta is so small, calculate Poisson deviance
    if (y == 0) {
      return 2.0 * mu;
    } else {
      // the max is necessary because some combination of y and mu give negative results:
      // e.g. y = 1, mu = 0.99999999999994
      const auto d = 2.0 * (y * std::log(y / mu) - (y - mu));
      return (d > 0. ? d : 0.);
    }
  } else {
    // Otherwise calculate Gamma-Poisson deviance
    if (y == 0) {
      return 2.0 / theta * std::log1p(mu * theta);
    } else {
      const double s1 = y * std::log((mu + y * mu * theta) / (y + y * mu * theta));
      const double s2 = 1.0 / theta * std::log((1 + mu * theta) / (1 + y * theta));
      const auto d = -2.0 * (s1 - s2);
      return (d > 0. ? d : 0.);
    }
  }
}

template <class D1, class D2> inline double compute_gp_deviance_sum(const EMB<D1> &y, const EMB<D2> &mu, const double theta) {
  double dev = 0.0;
  auto y_begin = y.cbegin();
  auto mu_begin = mu.cbegin();
  const auto y_end = y.cend();
  for (; y_begin != y_end; ++y_begin, ++mu_begin) {
    dev += compute_gp_deviance(*y_begin, *mu_begin, theta);
  }
  return dev;
}

template <class D1, class D2, class D3>
inline MatrixXd compute_gp_deviance_residuals_matrix_impl(const Eigen::MatrixBase<D1> &Y, const Eigen::MatrixBase<D2> &Mu,
                                                          const Eigen::MatrixBase<D3> &thetas) {
  const auto nr = Y.rows();
  const auto nc = Y.cols();
  MatrixXd result(nr, nc);
  for (auto i = 0; i < nr; ++i) {
    const auto theta = thetas(i);
    for (auto j = 0; j < nc; ++j) {
      const auto y = Y(i, j);
      const auto mu = Mu(i, j);
      result(i, j) = sgn(y - mu) * std::sqrt(compute_gp_deviance(y, mu, theta));
    }
  }
  return result;
}

#endif
