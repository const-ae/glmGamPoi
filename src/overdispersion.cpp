// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]



// This correction factor is necessary to avoid estimates of
// theta that are basically +Inf. The problem is that for
// some combination of the y, mu, and X the term
// lgamma(1/theta) and the log(det(t(X) %*% W %*% X))
// with W = diag(1/(1/mu + theta)) canceled each other
// exactly out for large theta.
const double cr_correction_factor = 0.99;


// [[Rcpp::export]]
List make_table(NumericVector x){
  std::unordered_map<long, size_t> counts;
  for (double v : x){
    ++counts[(long) v];
  }
  NumericVector keys(counts.size());
  NumericVector values(counts.size());
  transform(counts.begin(), counts.end(), keys.begin(), [](std::pair<int, size_t> pair){return (double) pair.first;});
  transform(counts.begin(), counts.end(), values.begin(), [](std::pair<int, size_t> pair){return (double) pair.second;});
  return List::create(keys, values);
}


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

// [[Rcpp::export]]
double conventional_loglikelihood_fast(NumericVector y, NumericVector mu, double log_theta, const arma::mat& model_matrix, bool do_cr_adj,
                                       NumericVector unique_counts = NumericVector::create(),
                                       NumericVector count_frequencies = NumericVector::create()) {
  double theta = exp(log_theta);
  double cr_term = 0.0;
  if(do_cr_adj){
    arma::vec w_diag = 1.0 / (1.0 / mu + theta);
    arma::mat b = model_matrix.t() * (model_matrix.each_col() % w_diag);
    cr_term = -0.5 * log(det(b))  * cr_correction_factor;
  }
  double theta_neg1 = R_pow_di(theta, -1);
  double lgamma_term = 0;
  // If summarized counts are available use those to calculate sum(lgamma(y + theta_neg1))
  if(unique_counts.size() > 0 && unique_counts.size() == count_frequencies.size()){
    for(size_t iter = 0; iter < count_frequencies.size(); ++iter){
      lgamma_term += count_frequencies[iter] * lgamma(unique_counts[iter] + theta_neg1);
    }
  }else{
    lgamma_term = sum(lgamma(y + theta_neg1));
  }
  lgamma_term -=  y.size() * lgamma(theta_neg1);
  double ll_part = 0.0;
  for(size_t i = 0; i < y.size(); ++i){
    ll_part += (-y[i] - theta_neg1) * log(mu[i] + theta_neg1);
  }
  ll_part -= y.size() * theta_neg1 * log(theta);
  return lgamma_term + ll_part + cr_term;
}


// this function returns the derivative of the log posterior with respect to the log of the
// dispersion parameter alpha, given the same inputs as the previous function

// [[Rcpp::export]]
double conventional_score_function_fast(NumericVector y, NumericVector mu, double log_theta, const arma::mat& model_matrix, bool do_cr_adj,
                                        NumericVector unique_counts = NumericVector::create(),
                                        NumericVector count_frequencies = NumericVector::create()) {
  double theta = exp(log_theta);
  double theta_neg1 = 1.0 / theta;
  double theta_neg2 = 1.0 / R_pow_di(theta, 2);

  double cr_term = 0.0;
  if(do_cr_adj){
    arma::vec w_diag = 1.0 / (1.0 / mu + theta);
    arma::vec dw_diag = -1 * w_diag % w_diag;
    arma::mat b = model_matrix.t() * (model_matrix.each_col() % w_diag);
    arma::mat db = model_matrix.t() * (model_matrix.each_col() % dw_diag);
    cr_term = -0.5 * trace(b.i() * db) * cr_correction_factor;
  }


  double digamma_term = 0;
  // If summarized counts are available use those to calculate sum(digamma(y + theta_neg1))
  if(unique_counts.size() > 0 && unique_counts.size() == count_frequencies.size()){
    for(size_t iter = 0; iter < count_frequencies.size(); ++iter){
      digamma_term += count_frequencies[iter] * Rf_digamma(unique_counts[iter] + theta_neg1);
    }
    digamma_term -= y.size() * Rf_digamma(theta_neg1);
  }else{
    digamma_term = sum(digamma(y + theta_neg1));
    digamma_term -= y.size() * Rf_digamma(theta_neg1);
  }
  digamma_term *= theta_neg2;
  double ll_part = 0.0;
  for(size_t i = 0; i < y.size(); ++i){
    ll_part += log(1 + mu[i] * theta) + (y[i] - mu[i]) / (mu[i] + theta_neg1);
  }
  ll_part *= theta_neg2;
  return (ll_part - digamma_term + cr_term) * theta;
}


// this function returns the second derivative of the log posterior with respect to the log of the
// dispersion parameter alpha, given the same inputs as the previous function

// [[Rcpp::export]]
double conventional_deriv_score_function_fast(NumericVector y, NumericVector mu, double log_theta, const arma::mat& model_matrix, bool do_cr_adj,
                                              NumericVector unique_counts = NumericVector::create(),
                                              NumericVector count_frequencies = NumericVector::create()) {
  double theta = exp(log_theta);
  double cr_term = 0.0;
  double cr_term2 = 0.0;
  if(do_cr_adj){
    arma::vec w_diag = 1/(1/mu + theta);
    arma::vec dw_diag = -1 * w_diag % w_diag;
    arma::vec d2w_diag = -2 * dw_diag % w_diag;

    arma::mat b = model_matrix.t() * (model_matrix.each_col() % w_diag);
    arma::mat db = model_matrix.t() * (model_matrix.each_col() % dw_diag);
    arma::mat d2b = model_matrix.t() * (model_matrix.each_col() % d2w_diag);
    arma::mat b_i = b.i();
    arma::mat d_i_db = b_i * db;
    double ddetb = trace(d_i_db);
    double d2detb = ((R_pow_di(ddetb, 2) - trace(d_i_db * d_i_db) + trace(b_i * d2b)) );
    cr_term = (0.5 * R_pow_di(ddetb, 2) - 0.5 * d2detb)  * cr_correction_factor;
    cr_term2 = -0.5 * ddetb * cr_correction_factor;
  }

  double theta_neg1 = R_pow_di(theta, -1);
  double theta_neg2 = R_pow_di(theta, -2);
  double digamma_term = 0.0;
  double trigamma_term = 0.0;

  // If summarized counts are available use those to calculate sum(digamma()) and sum(trigamma())
  if(unique_counts.size() > 0 && unique_counts.size() == count_frequencies.size()){
    for(size_t iter = 0; iter < count_frequencies.size(); ++iter){
      digamma_term += count_frequencies[iter] * Rf_digamma(unique_counts[iter] + theta_neg1);
      trigamma_term += count_frequencies[iter] * Rf_trigamma(unique_counts[iter] + theta_neg1);
    }
    trigamma_term *= theta_neg2;

    digamma_term -= y.size() * Rf_digamma(theta_neg1);
    trigamma_term -=  theta_neg2 * y.size() * Rf_trigamma(theta_neg1);
  }else{
    digamma_term = sum(digamma(y + theta_neg1));
    digamma_term -= y.size() * Rf_digamma(theta_neg1);

    trigamma_term = theta_neg2 * sum(trigamma(y + theta_neg1));
    trigamma_term -=  theta_neg2 * y.size() * Rf_trigamma(theta_neg1);
  }

  double ll_part_1 = 0.0;
  double ll_part_2 = 0.0;
  for(size_t i = 0; i < y.size(); ++i){
    ll_part_1 += log(1 + mu[i] * theta) + (y[i] - mu[i]) / (mu[i] + theta_neg1);
    ll_part_2 += (mu[i] * mu[i] * theta + y[i]) / (1 + mu[i] * theta) / (1 + mu[i] * theta);
  }
  double ll_part = -2 * R_pow_di(theta, -3) * (ll_part_1 - digamma_term) + theta_neg2 * (ll_part_2 + trigamma_term);

  double res = ((ll_part + cr_term) * R_pow_di(theta, 2) +
                (ll_part_1 - digamma_term) * theta_neg1 + cr_term2 * theta);
  return res;
}


// ------------------------------------------------------------------------------------------------

template<class NumericType>
NumericVector estimate_overdispersions_fast_internal(RObject Y, RObject mean_matrix, NumericMatrix model_matrix, bool do_cox_reid_adjustment,
                                       double n_subsamples){
  auto Y_bm = beachmat::create_matrix<NumericType>(Y);
  auto mean_mat_bm = beachmat::create_numeric_matrix(mean_matrix);
  int n_samples = Y_bm->get_ncol();
  int n_genes = Y_bm->get_nrow();
  NumericVector result(n_genes);
  if(n_genes != mean_mat_bm->get_nrow() || n_samples != mean_mat_bm->get_ncol()){
    throw std::runtime_error("Dimensions of Y and mean_matrix do not match");
  }

  // This is calling back to R, which simplifies my code a lot
  Environment glmGamPoiEnv = Environment::namespace_env("glmGamPoi");
  Function gampoi_overdispersion_mle = glmGamPoiEnv["gampoi_overdispersion_mle"];

  for(int gene_idx = 0; gene_idx < n_genes; gene_idx++){
    if (gene_idx % 100 == 0) checkUserInterrupt();
    typename NumericType::vector counts(n_samples);
    Y_bm->get_row(gene_idx, counts.begin());
    NumericVector mu(n_samples);
    mean_mat_bm->get_row(gene_idx, mu.begin());

    SEXP dispRes = gampoi_overdispersion_mle(counts, mu, model_matrix, do_cox_reid_adjustment, n_subsamples);
    SEXP disp = Rcpp::as<List>(dispRes)["estimate"];
    result(gene_idx) = Rcpp::as<double>(disp);
  }

  return result;
}

// [[Rcpp::export]]
NumericVector estimate_overdispersions_fast(RObject Y, RObject mean_matrix, NumericMatrix model_matrix, bool do_cox_reid_adjustment,
                              double n_subsamples){
  auto mattype=beachmat::find_sexp_type(Y);
  if (mattype==INTSXP) {
    return estimate_overdispersions_fast_internal<beachmat::integer_matrix>(Y, mean_matrix, model_matrix, do_cox_reid_adjustment, n_subsamples);
  } else if (mattype==REALSXP) {
    return estimate_overdispersions_fast_internal<beachmat::numeric_matrix>(Y, mean_matrix, model_matrix, do_cox_reid_adjustment, n_subsamples);
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
}
