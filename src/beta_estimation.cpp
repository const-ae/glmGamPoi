// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"

#include <deviance.h>
#include <fisher_scoring_steps.h>

using namespace Rcpp;



template<class NumericType>
void clamp_inplace(arma::Mat<NumericType>& v, double min, double max){
  for(int i = 0; i < v.n_elem; i++){
    if(v.at(i) < min){
      v.at(i) = min;
    }else if(v.at(i) > max){
      v.at(i) = max;
    }
  }
}


arma::vec calculate_mu(const arma::mat& model_matrix, const arma::vec& beta_hat, const arma::vec& exp_off){
  arma::vec mu_hat = exp(model_matrix * beta_hat) % exp_off;
  clamp_inplace(mu_hat, 1e-50, 1e50);
  return mu_hat;
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
template<class NumericType>
double decrease_deviance(/*In-Out Parameter*/ arma::vec& beta_hat,
                         /*In-Out Parameter*/ arma::vec& mu_hat,
                         const arma::vec& step,
                         const arma::mat& model_matrix,
                         const arma::mat& exp_off,
                         const arma::Col<NumericType>& counts,
                         const double theta, const double dev_old, const double tolerance){
  double speeding_factor = 1.0;
  int line_iter = 0;
  double dev = 0;
  beta_hat = beta_hat + step;
  while(true){
    mu_hat = calculate_mu(model_matrix, beta_hat, exp_off);
    dev = compute_gp_deviance_sum(counts, mu_hat, theta);
    double conv_test = fabs(dev - dev_old)/(fabs(dev) + 0.1);
    if(dev < dev_old || conv_test < tolerance){
      break; // while loop
    }else if(line_iter >= 100){
      // speeding factor is very small, something is going wrong here
      dev = std::numeric_limits<double>::quiet_NaN();
      break; // while loop
    }else{
      // Halfing the speed
      speeding_factor = speeding_factor / 2.0;
      beta_hat = beta_hat - step * speeding_factor;
    }
    line_iter++;
  }
  return dev;
}




//--------------------------------------------------------------------------------------------------//
// The following code was originally copied from https://github.com/mikelove/DESeq2/blob/master/src/DESeq2.cpp
// I adapted it to the needs of this project by:
//  * remove lambda / ridge penality
//  * remove weights
//  * Calculate actual deviance (2 * (log(f_NB(y | mu, theta)) - log(f_NB(y | y, theta))))
//    instead of just 2 * log(f_NB(y | mu, theta)),
//  * Support DelayedArrays
//  * Remove unncessary outputs: beta_mat_var, hat_diagonals, deviance
//  * Remove beta divergence check if abs(beta) very large
//  * Add line search that ensures that deviance is decreasing at every step


// fit the Negative Binomial GLM with Fisher scoring
// note: the betas are on the natural log scale
//
template<class NumericType, class BMNumericType>
List fitBeta_fisher_scoring_impl(RObject Y, const arma::mat& model_matrix, RObject exp_offset_matrix,
                                 NumericVector thetas, SEXP beta_matSEXP, double ridge_penalty,
                                 double tolerance, int max_iter, bool use_diagonal_approx) {
  auto Y_bm = beachmat::create_matrix<BMNumericType>(Y);
  auto exp_offsets_bm = beachmat::create_numeric_matrix(exp_offset_matrix);
  int n_samples = Y_bm->get_ncol();
  int n_genes = Y_bm->get_nrow();
  // The result
  arma::mat beta_mat = as<arma::mat>(beta_matSEXP);
  // deviance, convergence and tolerance


  NumericVector iterations(n_genes);
  NumericVector deviance(n_genes);
  for (int gene_idx = 0; gene_idx < n_genes; gene_idx++) {
    if (gene_idx % 100 == 0) checkUserInterrupt();
    // Fill count and offset vector from beachmat matrix
    arma::Col<NumericType> counts(n_samples);
    Y_bm->get_row(gene_idx, counts.begin());
    arma::Col<double> exp_off(n_samples);
    exp_offsets_bm->get_row(gene_idx, exp_off.begin());
    // Init beta and mu
    arma::vec beta_hat = beta_mat.row(gene_idx).t();
    arma::vec mu_hat = calculate_mu(model_matrix, beta_hat, exp_off);
    // Init deviance
    double dev_old = compute_gp_deviance_sum(counts, mu_hat, thetas(gene_idx));
    for (int t = 0; t < max_iter; t++) {
      iterations(gene_idx)++;
      // Find good direction to optimize beta
      arma::vec step;
      if(use_diagonal_approx){
        step = fisher_scoring_diagonal_step(model_matrix, counts, mu_hat, thetas(gene_idx) * mu_hat);
      }else{
        if(ridge_penalty == 0){
          step = fisher_scoring_qr_step(model_matrix, counts, mu_hat, thetas(gene_idx) * mu_hat);
        }else{
          step = fisher_scoring_qr_ridge_step(model_matrix, counts, mu_hat, thetas(gene_idx) * mu_hat, ridge_penalty);
        }
      }
      // Find step size that actually decreases the deviance
      double dev = decrease_deviance(beta_hat, mu_hat, step, model_matrix, exp_off, counts,
                                     thetas(gene_idx), dev_old, tolerance);
      double conv_test = fabs(dev - dev_old)/(fabs(dev) + 0.1);
      if (std::isnan(conv_test)) {
        // This should not happen
        beta_hat.fill(NA_REAL);
        iterations(gene_idx) = max_iter;
        break;
      }
      if (conv_test < tolerance) {
        break;
      }
      dev_old = dev;
    }
    beta_mat.row(gene_idx) = beta_hat.t();
  }

  return List::create(
    Named("beta_mat", beta_mat),
    Named("iter", iterations));
}


// [[Rcpp::export]]
List fitBeta_fisher_scoring(RObject Y, const arma::mat& model_matrix, RObject exp_offset_matrix,
                                  NumericVector thetas, SEXP beta_matSEXP, double ridge_penalty,
                                  double tolerance, int max_iter) {
  auto mattype=beachmat::find_sexp_type(Y);
  if (mattype==INTSXP) {
    return fitBeta_fisher_scoring_impl<int, beachmat::integer_matrix>(Y, model_matrix, exp_offset_matrix,
                                                                      thetas,  beta_matSEXP,
                                                                      /*ridge_penalty=*/ ridge_penalty,
                                                                      tolerance, max_iter,
                                                                      /*use_diagonal_approx=*/ false);
  } else if (mattype==REALSXP) {
    return fitBeta_fisher_scoring_impl<double, beachmat::numeric_matrix>(Y, model_matrix, exp_offset_matrix,
                                                                         thetas,  beta_matSEXP,
                                                                         /*ridge_penalty=*/ ridge_penalty,
                                                                         tolerance, max_iter,
                                                                         /*use_diagonal_approx=*/ false);
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
}



// [[Rcpp::export]]
List fitBeta_diagonal_fisher_scoring(RObject Y, const arma::mat& model_matrix, RObject exp_offset_matrix,
                                     NumericVector thetas, SEXP beta_matSEXP,
                                     double tolerance, int max_iter) {
  auto mattype=beachmat::find_sexp_type(Y);
  if (mattype==INTSXP) {
    return fitBeta_fisher_scoring_impl<int, beachmat::integer_matrix>(Y, model_matrix, exp_offset_matrix,
                                                                      thetas,  beta_matSEXP,
                                                                      /*ridge_penalty=*/ 0,
                                                                      tolerance, max_iter,
                                                                      /*use_diagonal_approx=*/ true);
  } else if (mattype==REALSXP) {
    return fitBeta_fisher_scoring_impl<double, beachmat::numeric_matrix>(Y, model_matrix, exp_offset_matrix,
                                                                         thetas,  beta_matSEXP,
                                                                         /*ridge_penalty=*/ 0,
                                                                         tolerance, max_iter,
                                                                         /*use_diagonal_approx=*/ true);
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
}



// If there is only one group, there is no need to do the full Fisher-scoring
// Instead a simple Newton-Raphson algorithm will do
template<class NumericType>
List fitBeta_one_group_internal(SEXP Y_SEXP, SEXP offsets_SEXP,
                       NumericVector thetas, NumericVector beta_start_values,
                       double tolerance, int maxIter) {
  auto Y_bm = beachmat::create_matrix<NumericType>(Y_SEXP);

  auto offsets_bm = beachmat::create_numeric_matrix(offsets_SEXP);
  int n_samples = Y_bm->get_ncol();
  int n_genes = Y_bm->get_nrow();
  NumericVector result(n_genes);
  IntegerVector iterations(n_genes);

  for(int gene_idx = 0; gene_idx < n_genes; gene_idx++){
    if (gene_idx % 100 == 0) checkUserInterrupt();

    double beta = beta_start_values(gene_idx);
    const double& theta = thetas(gene_idx);

    typename NumericType::vector counts(n_samples);
    Y_bm->get_row(gene_idx, counts.begin());
    NumericVector off(n_samples);
    offsets_bm->get_row(gene_idx, off.begin());
    // Newton-Raphson
    int iter = 0;
    for(; iter < maxIter; iter++){
      double dl = 0.0;
      double ddl = 0.0;
      bool all_zero = true;
      for(int sample_iter = 0; sample_iter < n_samples; sample_iter++){
        const auto count = counts[sample_iter];
        all_zero = all_zero && count == 0;
        const double mu = std::exp(beta + off[sample_iter]);
        const double denom = 1.0 + mu * theta;
        dl += (count - mu) / denom;
        ddl += mu * (1.0 + count * theta) / denom / denom;
        // ddl += mu / denom;           // This is what edgeR is using
      }
      if(all_zero){
        beta = R_NegInf;
        break;
      }
      const double step = dl / ddl;
      beta += step;
      if(std::abs(step) < tolerance){
        break;
      }
    }
    result(gene_idx) = beta;
    iterations(gene_idx) = iter;
  }
  return List::create(
    Named("beta", result),
    Named("iter", iterations)
  );
}

// [[Rcpp::export(rng = false)]]
List fitBeta_one_group(RObject Y, RObject offset_matrix,
                        NumericVector thetas, NumericVector beta_start_values,
                        double tolerance, int maxIter) {
  auto mattype=beachmat::find_sexp_type(Y);
  if (mattype==INTSXP) {
    return fitBeta_one_group_internal<beachmat::integer_matrix>(Y, offset_matrix, thetas, beta_start_values, tolerance, maxIter);
  } else if (mattype==REALSXP) {
    return fitBeta_one_group_internal<beachmat::numeric_matrix>(Y, offset_matrix, thetas, beta_start_values, tolerance, maxIter);
  } else {
    throw std::runtime_error("unacceptable matrix type");
  }
}









