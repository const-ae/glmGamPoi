#include <fit_beta.h>
#include <par_helpers.h>
#include <tatami_helpers.h>

#include "Rtatami.h"

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::Map;
using Eigen::VectorXi;

// [[Rcpp::export(name = "fisher_scoring_qr_step")]]
Eigen::MatrixXd fisher_scoring_qr_step_mask(const Eigen::Map<Eigen::MatrixXd> &model_matrix, const Eigen::Map<Eigen::VectorXd> &counts,
                                            const Eigen::Map<Eigen::VectorXd> &mu, const Eigen::Map<Eigen::VectorXd> &theta_times_mu) {
  return fisher_scoring_qr_step(model_matrix, counts, mu, theta_times_mu);
}

// [[Rcpp::export(name = "fisher_scoring_qr_ridge_step")]]
Eigen::MatrixXd fisher_scoring_qr_ridge_step_mask(const Eigen::Map<Eigen::MatrixXd> &model_matrix, const Eigen::Map<Eigen::VectorXd> &counts,
                                                  const Eigen::Map<Eigen::VectorXd> &mu, const Eigen::Map<Eigen::VectorXd> &theta_times_mu,
                                                  const Eigen::Map<Eigen::MatrixXd> &ridge_penalty, const Eigen::Map<Eigen::VectorXd> &ridge_target,
                                                  const Eigen::Map<Eigen::VectorXd> &beta_hat) {
  return fisher_scoring_qr_ridge_step(model_matrix, counts, mu, theta_times_mu, ridge_penalty, ridge_target, beta_hat);
}

// [[Rcpp::export(name = "fisher_scoring_diagonal_step")]]
Eigen::MatrixXd fisher_scoring_diagonal_step_mask(const Eigen::Map<Eigen::MatrixXd> &model_matrix, const Eigen::Map<Eigen::VectorXd> &counts,
                                                  const Eigen::Map<Eigen::VectorXd> &mu, const Eigen::Map<Eigen::VectorXd> &theta_times_mu) {
  return fisher_scoring_diagonal_step(model_matrix, counts, mu, theta_times_mu);
}

// Check how many unique rows are in a matrix and if this number is less than or equal to n
// This is important to determine if the model can be solved by group averages
// (ie. the numer of unique rows == number of columns)
// [[Rcpp::export(rng = false)]]
bool lte_n_equal_rows(const NumericMatrix &matrix, const int n, const double tolerance = 1e-10) {
  NumericMatrix reference_matrix(n, matrix.ncol());
  size_t n_matches = 0;
  for (size_t row_idx = 0; row_idx < matrix.nrow(); row_idx++) {
    bool matched = false;
    NumericMatrix::ConstRow vec = matrix(row_idx, _);
    for (size_t ref_idx = 0; ref_idx < n_matches; ref_idx++) {
      NumericMatrix::Row ref_vec = reference_matrix(ref_idx, _);
      if (sum(abs(vec - ref_vec)) < tolerance) {
        matched = true;
        break;
      }
    }
    if (!matched) {
      ++n_matches;
      if (n_matches > n) {
        return false;
      }
      reference_matrix(n_matches - 1, _) = vec;
    }
  }
  return true;
}

// [[Rcpp::export(rng = false)]]
IntegerVector get_row_groups(const NumericMatrix &matrix, const int n_groups, const double tolerance = 1e-10) {
  NumericMatrix reference_matrix(n_groups, matrix.ncol());
  IntegerVector groups(matrix.nrow());
  size_t n_matches = 0;
  for (size_t row_idx = 0; row_idx < matrix.nrow(); row_idx++) {
    bool matched = false;
    NumericMatrix::ConstRow vec = matrix(row_idx, _);
    for (size_t ref_idx = 0; ref_idx < n_matches; ref_idx++) {
      NumericMatrix::Row ref_vec = reference_matrix(ref_idx, _);
      if (sum(abs(vec - ref_vec)) < tolerance) {
        groups(row_idx) = ref_idx;
        matched = true;
        break;
      }
    }
    if (!matched) {
      groups(row_idx) = n_matches;
      reference_matrix(n_matches, _) = vec;
      ++n_matches;
    }
  }
  return groups + 1;
}

//--------------------------------------------------------------------------------------------------//
// The following code was originally copied from https://github.com/mikelove/DESeq2/blob/master/src/DESeq2.cpp
// I adapted it to the needs of this project by:
//  * remove weights
//  * Calculate actual deviance (2 * (log(f_NB(y | mu, theta)) - log(f_NB(y | y, theta))))
//    instead of just 2 * log(f_NB(y | mu, theta)),
//  * Support DelayedArrays
//  * Remove unncessary outputs: beta_mat_var, hat_diagonals, deviance
//  * Remove beta divergence check if abs(beta) very large
//  * Add line search that ensures that deviance is decreasing at every step
//  * Add optional ridge penalty

// fit the Negative Binomial GLM with Fisher scoring
// note: the betas are on the natural log scale
//
template <class FSConf>
List fitBeta_fisher_scoring_impl(const RObject Y, const Eigen::Map<Eigen::MatrixXd> &model_matrix, const RObject exp_offset_matrix,
                                 const NumericVector thetas, const Eigen::Map<Eigen::MatrixXd> &beta_mat_v,
                                 // callbacks
                                 const FSConf &conf,
                                 // misc other parameters
                                 const double tolerance, const double max_rel_mu_change, const int max_iter, const bool try_recov_w_optim,
                                 const int do_parallel) {
  Rtatami::BoundNumericPointer Y_bm_ptr(Y);
  const auto &Y_bm = *(Y_bm_ptr->ptr);

  Rtatami::BoundNumericPointer exp_offsets_bm_ptr(exp_offset_matrix);
  const auto &exp_offsets_bm = *(exp_offsets_bm_ptr->ptr);

  const auto n_samples = Y_bm.ncol();
  const auto n_genes = Y_bm.nrow();

  // out matrix, initialized as copy from view of initial betas matrix
  MatrixXd beta_mat = beta_mat_v;

  // deviance, convergence and tolerance
  IntegerVector iterations(n_genes);
  NumericVector deviance(n_genes);

  Map<VectorXi> iterations_v(iterations.begin(), iterations.size());
  Map<VectorXd> deviance_v(deviance.begin(), deviance.size());

  const auto run = [&](const auto start, const auto length) -> void {
    auto Y_ext = tatami::consecutive_extractor<false>(Y_bm, true, start, length);
    std::unique_ptr<tatami::OracularDenseExtractor<double, tatami::NumericMatrix::index_type>> exp_offsets_ext;
    if (exp_offsets_bm.nrow() > 1) {
      exp_offsets_ext = tatami::consecutive_extractor<false>(exp_offsets_bm, true, start, length);
    } else {
      exp_offsets_ext = tatami::new_extractor<false, true>(exp_offsets_bm, true, std::make_shared<ConstIndexOracle<0>>(length));
    }
    // buffers needed in case pointers to data from extractors are not available (e.g. DelayedArray)
    VectorXd counts(n_samples), exp_off(n_samples);

    const auto grp_max_id = start + length;
    for (size_t gene_idx = start; gene_idx < grp_max_id; gene_idx++) {
      if (gene_idx % 100 == 0) {
        if (check_interrupt()) {
          return;
        }
      }

      // not using copy_n to avoid copies when not necessary, using Eigen::Map type instead.
      const auto cptr = Y_ext->fetch(counts.data());
      const auto eptr = exp_offsets_ext->fetch(exp_off.data());
      const Map<const VectorXd> counts_v(cptr, n_samples), exp_off_v(eptr, n_samples);

      auto beta_row = beta_mat.row(gene_idx); // using auto on purpose to get a block expression

      fitBeta_FS_internal_step(beta_row, deviance_v(gene_idx), iterations_v(gene_idx), model_matrix, counts_v, exp_off_v, conf, thetas(gene_idx),
                               tolerance, max_rel_mu_change, max_iter, try_recov_w_optim);
    }
  };

  if (do_parallel > 0) {
    run_par(run, n_genes, do_parallel);
  } else {
    run(0, n_genes);
  }
  // check if we got an interrupt, if yes, re-raise it
  if (check_interrupt()) {
    std::raise(SIGINT);
    Rcpp::checkUserInterrupt();
  }

  return List::create(_["beta_mat"] = beta_mat, _["iter"] = iterations, _["deviance"] = deviance);
}

// [[Rcpp::export(rng = false)]]
List fitBeta_fisher_scoring(const RObject Y, const Eigen::Map<Eigen::MatrixXd> model_matrix, const RObject exp_offset_matrix,
                            const NumericVector thetas, const Eigen::Map<Eigen::MatrixXd> &beta_mat_v, const Nullable<NumericMatrix> ridge_penalty_nl,
                            const double tolerance, const double max_rel_mu_change, const int max_iter, const bool try_recov_w_optim = false,
                            const int do_parallel = 0) {
  if (ridge_penalty_nl.isNotNull()) {
    const NumericMatrix tmp = ridge_penalty_nl.get();
    const auto ridge_penalty = Rcpp::as<Map<MatrixXd>>(tmp);
    const auto nc = ridge_penalty.cols();
    if (model_matrix.cols() != nc) {
      stop("Number of columns in model_matrix does not match the columns of the ridge_penalty");
    }
    const auto n_samples = static_cast<double>(model_matrix.rows());

    if (tmp.hasAttribute("target")) {
      const NumericVector ridge_target_buf = tmp.attr("target");
      const Map<const VectorXd> ridge_target(ridge_target_buf.begin(), ridge_target_buf.size());

      const FisherScoreQRwRidge conf(ridge_target, ridge_penalty, n_samples);
      return fitBeta_fisher_scoring_impl(Y, model_matrix, exp_offset_matrix, thetas, beta_mat_v, conf, tolerance, max_rel_mu_change, max_iter,
                                         try_recov_w_optim, do_parallel);
    } else {
      const VectorXd ridge_target = VectorXd::Constant(nc, 0.0);

      const FisherScoreQRwRidge conf(ridge_target, ridge_penalty, n_samples);
      return fitBeta_fisher_scoring_impl(Y, model_matrix, exp_offset_matrix, thetas, beta_mat_v, conf, tolerance, max_rel_mu_change, max_iter,
                                         try_recov_w_optim, do_parallel);
    }
  }

  const auto conf = FisherScoreQR();
  return fitBeta_fisher_scoring_impl(Y, model_matrix, exp_offset_matrix, thetas, beta_mat_v, conf, tolerance, max_rel_mu_change, max_iter,
                                     try_recov_w_optim, do_parallel);
}

// [[Rcpp::export(rng = false)]]
List fitBeta_diagonal_fisher_scoring(const RObject Y, const Eigen::Map<Eigen::MatrixXd> model_matrix, const RObject exp_offset_matrix,
                                     const NumericVector thetas, const Eigen::Map<Eigen::MatrixXd> &beta_mat_v, const double tolerance,
                                     const double max_rel_mu_change, const int max_iter, const bool try_recov_w_optim = false,
                                     const int do_parallel = 0) {
  const auto conf = FisherScoreDiagApprox();
  return fitBeta_fisher_scoring_impl(Y, model_matrix, exp_offset_matrix, thetas, beta_mat_v, conf, tolerance, max_rel_mu_change,
                                     max_iter, try_recov_w_optim, do_parallel);
}

template <class FSConf>
List fitBeta_optim_impl(const RObject Y, const Eigen::Map<Eigen::MatrixXd> &model_matrix, const RObject exp_offset_matrix, const NumericVector thetas,
                        const Eigen::Map<Eigen::MatrixXd> &beta_mat_v, const FSConf &conf, const int max_iter, const int do_parallel) {
  Rtatami::BoundNumericPointer Y_bm_ptr(Y);
  const auto &Y_bm = *(Y_bm_ptr->ptr);

  Rtatami::BoundNumericPointer exp_offsets_bm_ptr(exp_offset_matrix);
  const auto &exp_offsets_bm = *(exp_offsets_bm_ptr->ptr);

  const auto n_samples = Y_bm.ncol();
  const auto n_genes = Y_bm.nrow();

  MatrixXd beta_mat = beta_mat_v;

  IntegerVector iterations(n_genes);
  NumericVector deviance(n_genes);

  Map<VectorXi> iterations_v(iterations.begin(), iterations.size());
  Map<VectorXd> deviance_v(deviance.begin(), deviance.size());

  const auto run = [&](const auto start, const auto length) -> void {
    auto Y_ext = tatami::consecutive_extractor<false>(Y_bm, true, start, length);
    std::unique_ptr<tatami::OracularDenseExtractor<double, tatami::NumericMatrix::index_type>> exp_offsets_ext;
    if (exp_offsets_bm.nrow() > 1) {
      exp_offsets_ext = tatami::consecutive_extractor<false>(exp_offsets_bm, true, start, length);
    } else {
      exp_offsets_ext = tatami::new_extractor<false, true>(exp_offsets_bm, true, std::make_shared<ConstIndexOracle<0>>(length));
    }
    // buffers needed in case pointers to data from extractors are not available (e.g. DelayedArray)
    VectorXd counts(n_samples), exp_off(n_samples);

    const auto grp_max_id = start + length;
    for (size_t gene_idx = start; gene_idx < grp_max_id; gene_idx++) {
      if (gene_idx % 100 == 0) {
        if (check_interrupt()) {
          return;
        }
      }

      // not using copy_n to avoid copies when not necessary, using Eigen::Map type instead.
      const auto cptr = Y_ext->fetch(counts.data());
      const auto eptr = exp_offsets_ext->fetch(exp_off.data());
      const Map<const VectorXd> counts_v(cptr, n_samples), exp_off_v(eptr, n_samples);

      VectorXd beta_row = beta_mat.row(gene_idx);
      fitBeta_FS_optim_step(beta_row, deviance_v(gene_idx), iterations_v(gene_idx), model_matrix, counts_v, exp_off_v, conf, thetas(gene_idx),
                            max_iter);
      beta_mat.row(gene_idx) = beta_row;
    }
  };

  if (do_parallel > 0) {
    run_par(run, n_genes, do_parallel);
  } else {
    run(0, n_genes);
  }
  // check if we got an interrupt, if yes, re-raise it
  if (check_interrupt()) {
    std::raise(SIGINT);
    Rcpp::checkUserInterrupt();
  }

  return List::create(_["Beta"] = beta_mat, _["iterations"] = iterations, _["deviances"] = deviance);
}

//[[Rcpp::export(rng = false)]]
List fitBeta_optim(const RObject Y, const Eigen::Map<Eigen::MatrixXd> &model_matrix, const RObject exp_offset_matrix, const NumericVector thetas,
                   const Eigen::Map<Eigen::MatrixXd> &beta_mat_v, const Nullable<NumericMatrix> ridge_penalty_nl, const int max_iter,
                   const int do_parallel = 0) {
  if (ridge_penalty_nl.isNotNull()) {
    const NumericMatrix tmp = ridge_penalty_nl.get();
    const auto ridge_penalty = Rcpp::as<Map<MatrixXd>>(tmp);
    const auto nc = ridge_penalty.cols();
    if (model_matrix.cols() != nc) {
      stop("Number of columns in model_matrix does not match the columns of the ridge_penalty");
    }
    const auto n_samples = static_cast<double>(model_matrix.rows());

    if (tmp.hasAttribute("target")) {
      const NumericVector ridge_target_buf = tmp.attr("target");
      const Map<const VectorXd> ridge_target(ridge_target_buf.begin(), ridge_target_buf.size());

      const FisherScoreQRwRidge conf(ridge_target, ridge_penalty, n_samples);
      return fitBeta_optim_impl(Y, model_matrix, exp_offset_matrix, thetas, beta_mat_v, conf, max_iter, do_parallel);
    } else {
      const VectorXd ridge_target = VectorXd::Constant(nc, 0.0);

      const FisherScoreQRwRidge conf(ridge_target, ridge_penalty, n_samples);
      return fitBeta_optim_impl(Y, model_matrix, exp_offset_matrix, thetas, beta_mat_v, conf, max_iter, do_parallel);
    }
  }

  const auto conf = FisherScoreQR();
  return fitBeta_optim_impl(Y, model_matrix, exp_offset_matrix, thetas, beta_mat_v, conf, max_iter, do_parallel);
}

// If there is only one group, there is no need to do the full Fisher-scoring
// Instead a simple Newton-Raphson algorithm will do
//
//[[Rcpp::export(rng = false)]]
List fitBeta_one_group(const RObject Y, const RObject offset_matrix, const NumericVector thetas, const NumericVector beta_start_values,
                       const double tolerance, const int max_iter, const int do_parallel = 0) {
  Rtatami::BoundNumericPointer Y_bm_ptr(Y);
  const auto &Y_bm = *(Y_bm_ptr->ptr);

  Rtatami::BoundNumericPointer offsets_bm_ptr(offset_matrix);
  const auto &offsets_bm = *(offsets_bm_ptr->ptr);

  const auto n_samples = Y_bm.ncol();
  const auto n_genes = Y_bm.nrow();

  NumericVector result = Rcpp::clone(beta_start_values);
  NumericVector deviance(n_genes);
  IntegerVector iterations(n_genes);

  Map<VectorXd> result_v(result.begin(), result.size()), deviance_v(deviance.begin(), deviance.size());
  Map<VectorXi> iterations_v(iterations.begin(), iterations.size());

  const auto run = [&](const auto start, const auto length) -> void {
    auto Y_ext = tatami::consecutive_extractor<false>(Y_bm, true, start, length);
    std::unique_ptr<tatami::OracularDenseExtractor<double, tatami::NumericMatrix::index_type>> offset_ext;
    if (offsets_bm.nrow() > 1) {
      offset_ext = tatami::consecutive_extractor<false>(offsets_bm, true, start, length);
    } else {
      offset_ext = tatami::new_extractor<false, true>(offsets_bm, true, std::make_shared<ConstIndexOracle<0>>(length));
    }
    VectorXd counts_vec(n_samples), off_vec(n_samples);

    const auto grp_max_id = start + length;
    for (size_t gene_idx = start; gene_idx < grp_max_id; gene_idx++) {
      if (gene_idx % 100 == 0) {
        if (check_interrupt()) {
          return;
        }
      }
      const auto cptr = Y_ext->fetch(counts_vec.data());
      const auto optr = offset_ext->fetch(off_vec.data());
      // not using copy_n to avoid copies when not necessary, using Eigen::Map type instead.
      const Map<const VectorXd> counts_v(cptr, n_samples), off_v(optr, n_samples);

      fitBeta_NR_internal_step(result_v(gene_idx), deviance_v(gene_idx), iterations_v(gene_idx), counts_v, off_v, thetas(gene_idx), tolerance,
                               max_iter);
    }
  };
  if (do_parallel > 0) {
    run_par(run, n_genes, do_parallel);
  } else {
    run(0, n_genes);
  }
  // check if we got an interrupt, if yes, re-raise it
  if (check_interrupt()) {
    std::raise(SIGINT);
    Rcpp::checkUserInterrupt();
  }

  return List::create(_["beta"] = result, _["iter"] = iterations, _["deviance"] = deviance);
}