#include <deviance.h>

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export(name = "compute_gp_deviance")]]
double compute_gp_deviance_mask(const double y, const double mu,const double theta) { return compute_gp_deviance(y, mu, theta); }

// [[Rcpp::export(name = "compute_gp_deviance_sum")]]
double compute_gp_deviance_sum_mask(const NumericVector y,const NumericVector mu,const double theta) {
  double dev = 0.0;
  const size_t n_elem = y.size();
  for (size_t i = 0; i < n_elem; i++) {
    dev += compute_gp_deviance(y[i], mu[i], theta);
  }
  return dev;
}

// [[Rcpp::export(name = "compute_gp_deviance_residuals_matrix")]]
Eigen::MatrixXd compute_gp_deviance_residuals_matrix_mask(const SEXP Y_SEXP, const Eigen::Map<Eigen::MatrixXd> &Mu,
                                                          const Eigen::Map<Eigen::MatrixXd> &thetas) {
  const SEXP dims = Rf_getAttrib(Y_SEXP, R_DimSymbol);
  const int nrow = INTEGER(dims)[0];
  const int ncol = INTEGER(dims)[1];
  if (TYPEOF(Y_SEXP) == INTSXP) {
    const Eigen::Map<const MatrixX<int>> Y(INTEGER(Y_SEXP), nrow, ncol);
    return compute_gp_deviance_residuals_matrix_impl(Y, Mu, thetas);
  } else if (TYPEOF(Y_SEXP) == REALSXP) {
    const Eigen::Map<const MatrixX<double>> Y(REAL(Y_SEXP), nrow, ncol);
    return compute_gp_deviance_residuals_matrix_impl(Y, Mu, thetas);
  } else {
    stop("Cannot handle Y_SEXP of this type.");
  }
}