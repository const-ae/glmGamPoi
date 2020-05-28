#ifndef DEVIANCE_H
#define DEVIANCE_H

// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;


template <typename T>
inline int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

inline double compute_gp_deviance (double y, double mu, double theta) {
  if(theta < 1e-6){
    // If theta is so small, calculate Poisson deviance
    if(y == 0){
      return 2.0 * mu;
    }else{
      // the max is necessary because some combination of y and mu give negative results:
      // e.g. y = 1, mu = 0.99999999999994
      return std::max(2.0 * (y * std::log(y/mu) - (y - mu)), 0.0);
    }
  }else{
    // Otherwise calculate Gamma-Poisson deviance
    if(y == 0){
      return 2.0/theta * std::log((1 + mu * theta));
    } else {
      double s1 = y * std::log((mu + y * mu * theta) / (y +  y * mu * theta));
      double s2 = 1.0/theta * std::log((1 + mu * theta) / (1 + y * theta));
      return std::max(-2.0 * (s1 - s2), 0.0);
    }
  }
}

template<class NumericType>
inline double compute_gp_deviance_sum(const arma::Mat<NumericType>& Y,
                               const arma::Mat<double>& Mu,
                               const NumericVector& thetas){
  double dev = 0.0;
  int nrows = Y.n_rows;
  for (int i = 0; i < Y.n_elem; i++) {
    dev += compute_gp_deviance(Y.at(i), Mu.at(i), thetas(i % nrows));
  }
  return dev;
}

template<class NumericType>
inline double compute_gp_deviance_sum(const arma::Mat<NumericType>& Y,
                               const arma::Mat<double>& Mu,
                               double theta){
  double dev = 0.0;
  for (int i = 0; i < Y.n_elem; i++) {
    dev += compute_gp_deviance(Y.at(i), Mu.at(i), theta);
  }
  return dev;
}


template<class NumericType>
inline arma::Mat<double> compute_gp_deviance_residuals_matrix_impl(const arma::Mat<NumericType> Y, const arma::Mat<double> Mu, NumericVector thetas) {
  arma::Mat<double> result(Y.n_rows, Y.n_cols);
  int nrows = Y.n_rows;
  for(int i = 0; i < Y.n_elem; i++){
    result(i) = sgn(Y.at(i) - Mu.at(i)) * sqrt(compute_gp_deviance(Y.at(i), Mu.at(i), thetas.at(i % nrows)));
  }
  return result;
}

inline arma::Mat<double> compute_gp_deviance_residuals_matrix(const SEXP Y_SEXP, const arma::Mat<double>& Mu, NumericVector thetas) {
  SEXP dims = Rf_getAttrib(Y_SEXP, R_DimSymbol);
  int nrow = INTEGER(dims)[0];
  int ncol = INTEGER(dims)[1];
  if(TYPEOF(Y_SEXP) == INTSXP){
    arma::Mat<int> Y(INTEGER(Y_SEXP), nrow, ncol, false);
    return compute_gp_deviance_residuals_matrix_impl<int>(Y, Mu, thetas);
  }else if(TYPEOF(Y_SEXP) == REALSXP){
    arma::Mat<double> Y(REAL(Y_SEXP), nrow, ncol, false);
    return compute_gp_deviance_residuals_matrix_impl<double>(Y, Mu, thetas);
  }else{
    stop("Cannot handle Y_SEXP of this type.");
  }
}


#endif
