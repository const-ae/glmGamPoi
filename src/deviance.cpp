// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <deviance.h>

using namespace Rcpp;


// [[Rcpp::export(name = "compute_gp_deviance")]]
double compute_gp_deviance_mask(double y, double mu, double theta) {
  return compute_gp_deviance(y, mu, theta);
}

// [[Rcpp::export(name = "compute_gp_deviance_sum")]]
double compute_gp_deviance_sum_mask(NumericVector y, NumericVector mu, double theta) {
  return compute_gp_deviance_sum(y, mu, theta);
}


// [[Rcpp::export(name = "compute_gp_deviance_residuals_matrix")]]
arma::Mat<double> compute_gp_deviance_residuals_matrix_mask(const SEXP Y_SEXP, const arma::Mat<double>& Mu, NumericVector thetas) {
  return compute_gp_deviance_residuals_matrix(Y_SEXP, Mu, thetas);
}

