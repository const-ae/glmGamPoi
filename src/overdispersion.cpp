#include <calc_helpers.h>
#include <overdispersion.h>
#include <par_helpers.h>
#include <tatami_helpers.h>

#include <algorithm> // for std::shuffle
#include <numeric>   // for std::iota
#include <random>    // for std::default_random_engine

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::ArrayXi;
using Eigen::Map;
using Eigen::VectorXd;

#include "Rtatami.h"

// [[Rcpp::export(rng = false)]]
List make_table_if_small(const NumericVector x, int stop_if_larger) {
  VectorXd keys, values;

  make_map_if_small(keys, values, x, stop_if_larger);
  return List::create(keys, values);
}

// we take y as a NumericVector and not a Map<VectorXd> to allow for implicit conversions when y is given as integers if this is called directly from
// R all numerical methods on the C++ side work assuming floating point inputs, so these copies are theoretically unavoidable (TODO: find way of maybe
// avoiding this?)
// [[Rcpp::export(rng = false)]]
double conventional_loglikelihood_fast(const NumericVector y, const Eigen::Map<Eigen::VectorXd> &mu, double log_theta,
                                       const Eigen::Map<Eigen::MatrixXd> &model_matrix, bool do_cr_adj,
                                       NumericVector unique_counts = NumericVector::create(),
                                       NumericVector count_frequencies = NumericVector::create()) {
  const Map<const VectorXd> unique_counts_v(unique_counts.begin(), unique_counts.size()),
      count_frequencies_v(count_frequencies.begin(), count_frequencies.size()), y_v(y.begin(), y.size());
  return conventional_loglikelihood_fast_impl(y_v, mu, log_theta, model_matrix, do_cr_adj, unique_counts_v, count_frequencies_v);
}
// [[Rcpp::export(rng = false)]]
double conventional_score_function_fast(const NumericVector y, const Eigen::Map<Eigen::VectorXd> &mu, double log_theta,
                                        const Eigen::Map<Eigen::MatrixXd> &model_matrix, bool do_cr_adj,
                                        NumericVector unique_counts = NumericVector::create(),
                                        NumericVector count_frequencies = NumericVector::create()) {
  const Map<const VectorXd> unique_counts_v(unique_counts.begin(), unique_counts.size()),
      count_frequencies_v(count_frequencies.begin(), count_frequencies.size()), y_v(y.begin(), y.size());
  return conventional_score_function_fast_impl(y_v, mu, log_theta, model_matrix, do_cr_adj, unique_counts_v, count_frequencies_v);
}
// [[Rcpp::export(rng = false)]]
double conventional_deriv_score_function_fast(const NumericVector y, const Eigen::Map<Eigen::VectorXd> &mu, double log_theta,
                                              const Eigen::Map<Eigen::MatrixXd> &model_matrix, bool do_cr_adj,
                                              const NumericVector unique_counts = NumericVector::create(),
                                              const NumericVector count_frequencies = NumericVector::create()) {
  const Map<const VectorXd> unique_counts_v(unique_counts.begin(), unique_counts.size()),
      count_frequencies_v(count_frequencies.begin(), count_frequencies.size()), y_v(y.begin(), y.size());
  return conventional_deriv_score_function_fast_impl(y_v, mu, log_theta, model_matrix, do_cr_adj, unique_counts_v, count_frequencies_v);
}
// [[Rcpp::export(rng = false)]]
List NR_overdispersion_mle(const NumericVector y, const Eigen::Map<Eigen::VectorXd> &mu_vector, const Eigen::Map<Eigen::MatrixXd> &model_matrix,
                           const bool do_cox_reid_adjustment, const int max_iter, const double tolerance = 1e-8) {
  const Map<const VectorXd> y_v(y.begin(), y.size());

  double est = NAN;
  int iters = -1;
  std::string msg = "";
  overdispersion_mle_NR_impl(est, iters, msg, y_v, mu_vector, model_matrix, do_cox_reid_adjustment, max_iter, tolerance);

  return List::create(_["estimate"] = est, _["iterations"] = iters, _["message"] = msg);
}

// [[Rcpp::export(rng = false)]]
List estimate_overdispersions_fast(const RObject Y, const RObject mean_matrix, const NumericMatrix model_matrix, const bool do_cox_reid_adjustment,
                                   const double n_subsamples, const int max_iter) {
  Rtatami::BoundNumericPointer Y_bm_ptr(Y);
  const auto &Y_bm = *(Y_bm_ptr->ptr);
  Rtatami::BoundNumericPointer mean_mat_bm_ptr(mean_matrix);
  const auto &mean_mat_bm = *(mean_mat_bm_ptr->ptr);

  const int n_samples = Y_bm.ncol();
  const int n_genes = Y_bm.nrow();

  NumericVector estimates(n_genes);
  IntegerVector iterations(n_genes);
  CharacterVector messages(n_genes);

  if (n_genes != mean_mat_bm.nrow() || n_samples != mean_mat_bm.ncol()) {
    throw std::runtime_error("Dimensions of Y and mean_matrix do not match");
  }

  auto Y_ext = tatami::consecutive_extractor<false>(Y_bm, true, 0, n_genes);
  auto mean_mat_ext = tatami::consecutive_extractor<false>(mean_mat_bm, true, 0, n_genes);
  NumericVector counts(n_samples), mu(n_samples);

  // This is calling back to R, which simplifies my code a lot
  Environment glmGamPoiEnv = Environment::namespace_env("glmGamPoi2");
  Function overdispersion_mle_impl = glmGamPoiEnv["overdispersion_mle_impl"];
  for (int gene_idx = 0; gene_idx < n_genes; gene_idx++) {
    if (gene_idx % 100 == 0)
      checkUserInterrupt();

    // Using copy_n to ensure that the vectors are actually filled.
    const auto cptr = Y_ext->fetch(counts.begin());
    tatami::copy_n(cptr, n_samples, counts.begin());
    const auto mptr = mean_mat_ext->fetch(mu.begin());
    tatami::copy_n(mptr, n_samples, mu.begin());

    // Check if the first value is NA, if yes all of them will be
    if (n_samples > 0 && std::isnan(mu[0])) {
      estimates(gene_idx) = NAN;
      iterations(gene_idx) = max_iter;
      messages(gene_idx) = "Mean estimate was NA. Cannot estimate overdispersion";
    } else {
      List dispRes = Rcpp::as<List>(overdispersion_mle_impl(counts, mu, model_matrix, do_cox_reid_adjustment, n_subsamples, max_iter));
      estimates(gene_idx) = Rcpp::as<double>(dispRes["estimate"]);
      iterations(gene_idx) = Rcpp::as<double>(dispRes["iterations"]);
      messages(gene_idx) = Rcpp::as<String>(dispRes["message"]);
    }
  }
  return List::create(_["estimate"] = estimates, _["iterations"] = iterations, _["message"] = messages);
}

// [[Rcpp::export(rng = false)]]
List estimate_overdispersions_fast_delayed(const RObject Y, const NumericMatrix model_matrix, const RObject offset_matrix,
                                           const Eigen::Map<Eigen::MatrixXd> &beta_mat_v, const bool do_cox_reid_adjustment, const int n_subsamples,
                                           const int max_iter) {
  Rtatami::BoundNumericPointer Y_bm_ptr(Y);
  const auto &Y_bm = *(Y_bm_ptr->ptr);

  Rtatami::BoundNumericPointer offsets_bm_ptr(offset_matrix);
  const auto &offsets_bm = *(offsets_bm_ptr->ptr);

  const int n_samples = Y_bm.ncol();
  const int n_genes = Y_bm.nrow();

  NumericVector estimates(n_genes);
  IntegerVector iterations(n_genes);
  CharacterVector messages(n_genes);

  auto Y_ext = tatami::consecutive_extractor<false>(Y_bm, true, 0, n_genes);
  std::unique_ptr<tatami::OracularDenseExtractor<double, int>> offsets_ext;
  if (offsets_bm.nrow() > 1) {
    offsets_ext = tatami::consecutive_extractor<false>(offsets_bm, true, 0, n_genes);
  } else {
    offsets_ext = tatami::new_extractor<false, true>(offsets_bm, true, std::make_shared<ConstIndexOracle<0>>(n_genes));
  }
  NumericVector counts(n_samples), mu(n_samples);
  VectorXd off(n_samples);

  // This is calling back to R, which simplifies my code a lot
  Environment glmGamPoiEnv = Environment::namespace_env("glmGamPoi2");
  Function overdispersion_mle_impl = glmGamPoiEnv["overdispersion_mle_impl"];

  const Map<const MatrixXd> model_matrix_m(model_matrix.cbegin(), model_matrix.nrow(), model_matrix.ncol());

  for (int gene_idx = 0; gene_idx < n_genes; gene_idx++) {
    if (gene_idx % 100 == 0)
      Rcpp::checkUserInterrupt();

    // Using copy_n to ensure that the vectors are actually filled.
    const auto cptr = Y_ext->fetch(counts.begin());
    tatami::copy_n(cptr, n_samples, counts.begin());
    const auto optr = offsets_ext->fetch(off.data());

    const Map<const VectorXd> counts_v(cptr, n_samples), off_v(optr, n_samples);
    Map<VectorXd> mu_v(mu.begin(), n_samples);
    const auto beta_hat = beta_mat_v.row(gene_idx).transpose();
    mu_v = calculate_mu_add<VectorXd>(model_matrix_m, beta_hat, off_v);

    // Check if the first value is NA, if yes all of them will be
    if (n_samples > 0 && std::isnan(mu[0])) {
      estimates(gene_idx) = NAN;
      iterations(gene_idx) = max_iter;
      messages(gene_idx) = "Mean estimate was NA. Cannot estimate overdispersion";
    } else {
      List dispRes =
          Rcpp::as<List>(overdispersion_mle_impl(counts, mu, model_matrix, do_cox_reid_adjustment, n_subsamples, max_iter, gene_idx == 9520));
      estimates(gene_idx) = Rcpp::as<double>(dispRes["estimate"]);
      iterations(gene_idx) = Rcpp::as<double>(dispRes["iterations"]);
      messages(gene_idx) = Rcpp::as<String>(dispRes["message"]);
    }
  }
  return List::create(_["estimate"] = estimates, _["iterations"] = iterations, _["message"] = messages);
}

// [[Rcpp::export(rng = false)]]
NumericVector estimate_global_overdispersions_fast(const RObject Y, const RObject mean_matrix, const Eigen::Map<Eigen::MatrixXd> &model_matrix,
                                                   const bool do_cox_reid_adjustment, const NumericVector log_thetas, const int do_parallel = 0) {
  Rtatami::BoundNumericPointer Y_bm_ptr(Y);
  const auto &Y_bm = *(Y_bm_ptr->ptr);
  Rtatami::BoundNumericPointer mean_mat_bm_ptr(mean_matrix);
  const auto &mean_mat_bm = *(mean_mat_bm_ptr->ptr);

  const int n_samples = Y_bm.ncol();
  const int n_genes = Y_bm.nrow();
  const int n_spline_points = log_thetas.size();

  NumericVector log_likelihoods(n_spline_points);

  const auto run = [&](const int start, const int length) -> void {
    auto Y_ext = tatami::consecutive_extractor<false>(Y_bm, true, start, length);
    auto mean_mat_ext = tatami::consecutive_extractor<false>(mean_mat_bm, true, start, length);
    VectorXd counts(n_samples), mu(n_samples);

    const auto grp_max_id = start + length;
    for (int gene_idx = start; gene_idx < grp_max_id; gene_idx++) {
      if (gene_idx % 100 == 0) {
        if (check_interrupt()) {
          return;
        }
      }

      const auto cptr = Y_ext->fetch(counts.data());
      const auto mptr = mean_mat_ext->fetch(mu.data());
      // not using copy_n to avoid copies when not necessary, using Eigen::Map type instead.
      const Map<const VectorXd> counts_v(cptr, n_samples), mu_v(mptr, n_samples);

      VectorXd unique_counts, count_frequencies;
      make_map_if_small(unique_counts, count_frequencies, counts_v, /*stop_if_larger = */ n_samples / 2);

      for (int point_idx = 0; point_idx < n_spline_points; point_idx++) {
        log_likelihoods[point_idx] += conventional_loglikelihood_fast_impl(counts_v, mu_v, log_thetas[point_idx], model_matrix,
                                                                           do_cox_reid_adjustment, unique_counts, count_frequencies);
      }
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

  return log_likelihoods;
}

// [[Rcpp::export(rng = false)]]
NumericVector estimate_global_overdispersions_fast_delayed(const RObject Y, const Eigen::Map<Eigen::MatrixXd> &model_matrix,
                                                           const RObject offset_matrix, const Eigen::Map<Eigen::MatrixXd> &beta_mat_v,
                                                           const bool do_cox_reid_adjustment, const NumericVector log_thetas,
                                                           const int do_parallel = 0) {
  Rtatami::BoundNumericPointer Y_bm_ptr(Y);
  const auto &Y_bm = *(Y_bm_ptr->ptr);
  Rtatami::BoundNumericPointer offsets_bm_ptr(offset_matrix);
  const auto &offsets_bm = *(offsets_bm_ptr->ptr);

  const int n_samples = Y_bm.ncol();
  const int n_genes = Y_bm.nrow();
  const int n_spline_points = log_thetas.size();

  NumericVector log_likelihoods(n_spline_points);

  const auto run = [&](const int start, const int length) -> void {
    auto Y_ext = tatami::consecutive_extractor<false>(Y_bm, true, start, length);
    std::unique_ptr<tatami::OracularDenseExtractor<double, int>> offsets_ext;
    if (offsets_bm.nrow() > 1) {
      offsets_ext = tatami::consecutive_extractor<false>(offsets_bm, true, start, length);
    } else {
      offsets_ext = tatami::new_extractor<false, true>(offsets_bm, true, std::make_shared<ConstIndexOracle<0>>(n_genes));
    }
    VectorXd counts(n_samples), off(n_samples);

    const auto grp_max_id = start + length;
    for (int gene_idx = start; gene_idx < grp_max_id; gene_idx++) {
      if (gene_idx % 100 == 0) {
        if (check_interrupt()) {
          return;
        }
      }

      const auto cptr = Y_ext->fetch(counts.data());
      const auto optr = offsets_ext->fetch(off.data());

      // not using copy_n to avoid copies when not necessary, using Eigen::Map type instead.
      const Map<const VectorXd> counts_v(cptr, n_samples), off_v(optr, n_samples);

      // not using copy_n to avoid copies when not necessary, using Eigen::Map type instead.
      const auto beta_hat = beta_mat_v.row(gene_idx).transpose();
      const auto mu = calculate_mu_add<VectorXd>(model_matrix, beta_hat, off_v);

      VectorXd unique_counts, count_frequencies;
      make_map_if_small(unique_counts, count_frequencies, counts_v, /*stop_if_larger = */ n_samples / 2);

      for (int point_idx = 0; point_idx < n_spline_points; point_idx++) {
        log_likelihoods[point_idx] += conventional_loglikelihood_fast_impl(counts_v, mu, log_thetas[point_idx], model_matrix, do_cox_reid_adjustment,
                                                                           unique_counts, count_frequencies);
      }
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

  return log_likelihoods;
}

// [[Rcpp::export(rng = true)]]
List estimate_overdispersions_nr_fast(const RObject Y, const RObject mean_matrix, const Eigen::Map<Eigen::MatrixXd> &model_matrix,
                                      const bool do_cox_reid_adjustment, const int n_subsamples, const int max_iter, const double tolerance = 1e-8,
                                      const int do_parallel = 0) {
  Rtatami::BoundNumericPointer Y_bm_ptr(Y);
  const auto &Y_bm = *(Y_bm_ptr->ptr);
  Rtatami::BoundNumericPointer mean_mat_bm_ptr(mean_matrix);
  const auto &mean_mat_bm = *(mean_mat_bm_ptr->ptr);

  const int n_samples = Y_bm.ncol();
  const int n_genes = Y_bm.nrow();

  if (n_genes != mean_mat_bm.nrow() || n_samples != mean_mat_bm.ncol()) {
    throw std::runtime_error("Dimensions of Y and mean_matrix do not match");
  }

  const bool do_sub = n_subsamples < n_samples;

  const unsigned int seed = do_sub ? (int)Rcpp::runif(1, 0., (double)INT_MAX)[0] : 0;

  const auto run = [&](const int start, const int length, auto &estimates, auto &iterations, auto &messages) -> void {
    auto Y_ext = tatami::consecutive_extractor<false>(Y_bm, true, start, length);
    auto mean_mat_ext = tatami::consecutive_extractor<false>(mean_mat_bm, true, start, length);
    VectorXd counts(n_samples), mu(n_samples);

    // use pointers to avoid initializing an RNG when not subsampling
    std::unique_ptr<std::default_random_engine> gen;
    std::unique_ptr<std::vector<int>> idx;
    if (do_sub) {
      // seed rng w/ value from R RNG to ensure reproducibility w/ set.seed from R
      gen = std::unique_ptr<std::default_random_engine>(new std::default_random_engine(seed + ((unsigned int)(start))));

      idx = std::unique_ptr<std::vector<int>>(new std::vector<int>(n_samples));
      std::iota(idx->begin(), idx->end(), 0);
    }

    const auto grp_max_id = start + length;
    for (int gene_idx = start; gene_idx < grp_max_id; gene_idx++) {
      if (gene_idx % 100 == 0) {
        if (check_interrupt()) {
          return;
        }
      }

      // Using copy_n to ensure that the vectors are actually filled.
      const auto cptr = Y_ext->fetch(counts.data());
      const auto mptr = mean_mat_ext->fetch(mu.data());
      // not using copy_n to avoid copies when not necessary, using Eigen::Map type instead.
      const Map<const VectorXd> counts_v(cptr, n_samples), mu_v(mptr, n_samples);

      // important to use && to keep reference semantics for &string AND allow for Rcpp's string proxy type that is by-value
      auto &&msg_out = messages[gene_idx];

      // Check if the first value is NA, if yes all of them will be
      // std::isnan is valid here since the R NA value for a double is NaN
      if (n_samples > 0 && std::isnan(mu_v(0))) {
        estimates[gene_idx] = NAN;
        iterations[gene_idx] = max_iter;
        msg_out = "Mean estimate was NA. Cannot estimate overdispersion";
        continue;
      }

      if (do_sub) {
        // dereferencing gen & idx is safe when gated against do_sub
        std::shuffle(idx->begin(), idx->end(), *gen);
        const std::vector<int> idx_s(idx->begin(), idx->begin() + n_subsamples);

        const VectorXd counts_s = counts_v(idx_s);
        const VectorXd mu_s = mu_v(idx_s);
        const MatrixXd mm_s = model_matrix(idx_s, Eigen::all);
        overdispersion_mle_NR_impl(estimates[gene_idx], iterations[gene_idx], msg_out, counts_s, mu_s, mm_s, do_cox_reid_adjustment, max_iter,
                                   tolerance);
      } else {
        overdispersion_mle_NR_impl(estimates[gene_idx], iterations[gene_idx], msg_out, counts_v, mu_v, model_matrix, do_cox_reid_adjustment, max_iter,
                                   tolerance);
      }
    }
  };

  if (do_parallel > 0) {
    std::vector<double> estimates(n_genes);
    std::vector<int> iterations(n_genes);
    std::vector<std::string> messages(n_genes);

    const auto run_w = [&run, &estimates, &iterations, &messages](const int start, const int length) -> void {
      run(start, length, estimates, iterations, messages);
    };

    run_par(run_w, n_genes, do_parallel);
    if (check_interrupt()) {
      std::raise(SIGINT);
      Rcpp::checkUserInterrupt();
    }
    return List::create(_["estimate"] = estimates, _["iterations"] = iterations, _["message"] = messages);
  } else {
    NumericVector estimates(n_genes);
    IntegerVector iterations(n_genes);
    CharacterVector messages(n_genes);

    run(0, n_genes, estimates, iterations, messages);
    if (check_interrupt()) {
      std::raise(SIGINT);
      Rcpp::checkUserInterrupt();
    }
    return List::create(_["estimate"] = estimates, _["iterations"] = iterations, _["message"] = messages);
  }
}

// [[Rcpp::export(rng = true)]]
List estimate_overdispersions_nr_fast_delayed(const RObject Y, const Eigen::Map<Eigen::MatrixXd> &model_matrix, const RObject offset_matrix,
                                              const Eigen::Map<Eigen::MatrixXd> &beta_mat_v, const bool do_cox_reid_adjustment,
                                              const int n_subsamples, const int max_iter, const double tolerance = 1e-8, const int do_parallel = 0) {
  Rtatami::BoundNumericPointer Y_bm_ptr(Y);
  const auto &Y_bm = *(Y_bm_ptr->ptr);

  Rtatami::BoundNumericPointer offsets_bm_ptr(offset_matrix);
  const auto &offsets_bm = *(offsets_bm_ptr->ptr);

  const int n_samples = Y_bm.ncol();
  const int n_genes = Y_bm.nrow();

  std::vector<double> estimates(n_genes);
  std::vector<int> iterations(n_genes);
  std::vector<std::string> messages(n_genes);

  const bool do_sub = n_subsamples < n_samples;

  const unsigned int seed = do_sub ? (unsigned int)Rcpp::runif(1, 0., (double)UINT_MAX)[0] : 0;

  const auto run = [&](const int start, const int length, auto &estimates, auto &iterations, auto &messages) -> void {
    auto Y_ext = tatami::consecutive_extractor<false>(Y_bm, true, start, length);
    std::unique_ptr<tatami::OracularDenseExtractor<double, int>> offsets_ext;
    if (offsets_bm.nrow() > 1) {
      offsets_ext = tatami::consecutive_extractor<false>(offsets_bm, true, start, length);
    } else {
      offsets_ext = tatami::new_extractor<false, true>(offsets_bm, true, std::make_shared<ConstIndexOracle<0>>(length));
    }
    VectorXd counts(n_samples), off(n_samples);

    // use pointers to avoid initializing an RNG when not subsampling
    std::unique_ptr<std::default_random_engine> gen;
    std::unique_ptr<std::vector<int>> idx;
    if (do_sub) {
      // seed rng w/ value from R RNG to ensure reproducibility w/ set.seed from R
      gen = std::unique_ptr<std::default_random_engine>(new std::default_random_engine(seed + ((unsigned int)(start))));

      idx = std::unique_ptr<std::vector<int>>(new std::vector<int>(n_samples));
      std::iota(idx->begin(), idx->end(), 0);
    }

    const auto grp_max_id = start + length;
    for (int gene_idx = start; gene_idx < grp_max_id; gene_idx++) {
      if (gene_idx % 100 == 0) {
        if (check_interrupt()) {
          return;
        }
      }

      // Using copy_n to ensure that the vectors are actually filled.
      const auto cptr = Y_ext->fetch(counts.data());
      const auto optr = offsets_ext->fetch(off.data());

      // not using copy_n to avoid copies when not necessary, using Eigen::Map type instead.
      const Map<const VectorXd> counts_v(cptr, n_samples), off_v(optr, n_samples);

      // important to use && to keep reference semantics for &string AND allow for Rcpp's string proxy type that is by-value
      auto &&msg_out = messages[gene_idx];

      const auto beta_hat = beta_mat_v.row(gene_idx).transpose();
      if (do_sub) {
        // dereferencing gen & idx is safe when gated against do_sub
        std::shuffle(idx->begin(), idx->end(), *gen);
        const std::vector<int> idx_s(idx->begin(), idx->begin() + n_subsamples);

        const VectorXd counts_s = counts_v(idx_s);
        const VectorXd off_s = off_v(idx_s);
        const MatrixXd mm_s = model_matrix(idx_s, Eigen::all);

        const auto mu = calculate_mu_add<VectorXd>(mm_s, beta_hat, off_s);
        if (n_samples > 0 && std::isnan(mu(0))) {
          estimates[gene_idx] = NAN;
          iterations[gene_idx] = max_iter;
          msg_out = "Mean estimate was NA. Cannot estimate overdispersion";
          continue;
        }
        overdispersion_mle_NR_impl(estimates[gene_idx], iterations[gene_idx], msg_out, counts_s, mu, mm_s, do_cox_reid_adjustment, max_iter,
                                   tolerance);
      } else {
        const auto mu = calculate_mu_add<VectorXd>(model_matrix, beta_hat, off);
        if (n_samples > 0 && std::isnan(mu(0))) {
          estimates[gene_idx] = NAN;
          iterations[gene_idx] = max_iter;
          msg_out = "Mean estimate was NA. Cannot estimate overdispersion";
          continue;
        }
        overdispersion_mle_NR_impl(estimates[gene_idx], iterations[gene_idx], msg_out, counts_v, mu, model_matrix, do_cox_reid_adjustment, max_iter,
                                   tolerance);
      }
    }
  };

  if (do_parallel > 0) {
    std::vector<double> estimates(n_genes);
    std::vector<int> iterations(n_genes);
    std::vector<std::string> messages(n_genes);

    const auto run_w = [&run, &estimates, &iterations, &messages](const int start, const int length) -> void {
      run(start, length, estimates, iterations, messages);
    };

    run_par(run_w, n_genes, do_parallel);
    if (check_interrupt()) {
      std::raise(SIGINT);
      Rcpp::checkUserInterrupt();
    }
    return List::create(_["estimate"] = estimates, _["iterations"] = iterations, _["message"] = messages);
  } else {
    NumericVector estimates(n_genes);
    IntegerVector iterations(n_genes);
    CharacterVector messages(n_genes);

    run(0, n_genes, estimates, iterations, messages);
    if (check_interrupt()) {
      std::raise(SIGINT);
      Rcpp::checkUserInterrupt();
    }
    return List::create(_["estimate"] = estimates, _["iterations"] = iterations, _["message"] = messages);
  }
}
