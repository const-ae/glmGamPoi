local({
  for (seed in c(2, 8, 16)) {
    set.seed(seed)
    Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10)
    annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE), cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))

    for (cr_adj in c(FALSE, TRUE)) {
      logliks.fn <- function(mu, thetas, mm) {
        vapply(
          seq_len(nrow(Y)),
          function(i) {
            mu_i <- if (is.function(mu)) mu(i, ) else mu[i, ]
            conventional_loglikelihood_fast(Y[i, ], ifelse(mu_i > 1e-128, mu_i, 1e-6), log(thetas[[i]]), mm, cr_adj)
          },
          numeric(1)
        )
      }

      for (od in list(TRUE, "global")) {
        set.seed(1)
        ref_nr <- glm_gp(Y, col_data = annot, overdispersion = od, do_cox_reid_adjustment = cr_adj, design = ~group)
        set.seed(1)
        ref_fs <- glm_gp(Y, col_data = annot, overdispersion = od, do_cox_reid_adjustment = cr_adj, design = ~ group + cont1 + cont2)

        for (off_as_vec in c(FALSE, TRUE)) {
          for (delay_mu in c(FALSE, TRUE)) {
            for (p in c(1L, 32L)) {
              for (cast_to_dgR in c(FALSE, TRUE)) {
                conf_str <- sprintf(
                  "(seed=%s) perf_optim=list(offset_as_vec=%s, cast_dgC_Y_to_dgR=%s, delay_mu=%s, do_parallel=%s) w/ overdispersion=%s, do_cox_reid_adjustment=%s",
                  seed,
                  off_as_vec,
                  cast_to_dgR,
                  delay_mu,
                  p,
                  od,
                  cr_adj
                )

                test_w_ref <- function(ref, templ) {
                  set.seed(1)
                  res <- suppressWarnings(glm_gp(
                    Y,
                    design = ref[["design_formula"]],
                    col_data = annot,
                    overdispersion = od,
                    do_cox_reid_adjustment = cr_adj,
                    perf_optim = list(offset_as_vec = off_as_vec, cast_dgC_Y_to_dgR = cast_to_dgR, delay_mu = delay_mu, do_parallel = p)
                  ))

                  attr(res, "perf_optim") <- attr(ref, "perf_optim")
                  if (off_as_vec) {
                    res[["Offset"]] <- matrix(res[["Offset"]], nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
                  }
                  if (delay_mu) {
                    res[["Mu"]] <- res[["Mu"]]()
                  }

                  test_that(sprintf("%s : yields results consistent w/ reference", templ), expect_equal(ref, res, tolerance = if (delay_mu) 1e-7 else 0))

                  if (isTRUE(od)) {
                    res_alt <- suppressWarnings(glm_gp(
                      Y,
                      design = ref[["design_formula"]],
                      col_data = annot,
                      overdispersion = od,
                      do_cox_reid_adjustment = cr_adj,
                      perf_optim = list(offset_as_vec = off_as_vec, cast_dgC_Y_to_dgR = cast_to_dgR, delay_mu = delay_mu, do_parallel = p, use_nr_overdisp_impl = TRUE)
                    ))

                    test_that(sprintf("%s : perf_optim$use_nr_overdisp_impl generates at worst deviances larger by 0.01%%", templ), {
                      expect_all_true(((res_alt[["deviances"]] - res[["deviances"]]) <= (1e-4 * abs(res[["deviances"]]))))
                    })
                    test_that(sprintf("%s : perf_optim$use_nr_overdisp_impl generates at worst loglikehoods smaller by 0.01%%", templ), {
                      liks.ref <- logliks.fn(res[["Mu"]], res[["overdispersions"]], res[["model_matrix"]])
                      liks.alt <- logliks.fn(res_alt[["Mu"]], res_alt[["overdispersions"]], res_alt[["model_matrix"]])
                      mask <- !is.na(liks.ref)
                      liks.ref <- liks.ref[mask]
                      liks.alt <- liks.alt[mask]
                      expect_all_true(((liks.ref - liks.alt)) <= (1e-4 * abs(liks.ref)))
                    })
                    test_that(sprintf("%s : perf_optim$use_nr_overdisp_impl doesn't significantly change resulting Beta values", templ), {
                      skip_if(cr_adj, "known issue, default overdispersion estimator often ends up falling back to version w/o cox-reid adjustment, which is not true for C++ NR-based estimator")
                      expect_equal(res[["Beta"]], res_alt[["Beta"]], tolerance = 1e-8)
                    })
                  }
                }

                test_w_ref(ref_nr, sprintf("%s, NR case", conf_str))
                test_w_ref(ref_fs, sprintf("%s, FS case", conf_str))
              }
            }
          }
        }
      }
    }
  }
})

local({
  set.seed(1)
  y <- matrix(rnbinom(n = 1e4 * 4, mu = 5, size = 1 / 0.7), nrow = 4)
  annot <- data.frame(cont1 = runif(1e4))
  mm <- handle_design_parameter(~cont1, y, annot, NULL)$model_matrix
  off <- combine_size_factors_and_offset(0, "normed_sum", y)$offset_matrix
  beta <- estimate_betas_roughly(y, mm, off)

  test_w_mu <- function(mu, n) {
    r1 <- overdispersion_mle(y, mu, mm, use_nr_overdisp_impl = TRUE, do_parallel = n)

    set.seed(1)
    r2 <- overdispersion_mle(y, mu, mm, subsample = 1000, use_nr_overdisp_impl = TRUE, do_parallel = n)
    # check that shuffle respects seeding
    set.seed(1)
    r3 <- overdispersion_mle(y, mu, mm, subsample = 1000, use_nr_overdisp_impl = TRUE, do_parallel = n)

    expect_equal(r1$estimate, r2$estimate, tolerance = 0.1)
    expect_equal(r2$estimate, r3$estimate, tolerance = 0)
  }
  for (p in c(0L, 1L, 2L, 32L)) {
    test_that(sprintf("use_nr_overdisp_impl w/ n_subsamples parameter works (do_parallel=%s)", p), {
      test_w_mu(calculate_mu(beta, mm, off), p)
    })
    test_that(sprintf("use_nr_overdisp_impl w/ n_subsamples parameter works w/ delayed Mu (do_parallel=%s)", p), {
      test_w_mu(mk_delayed_mu(beta, mm, off), p)
    })
  }
})


local({
  set.seed(1)
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10)
  rownames(Y) <- sprintf("gene_%s", seq_len(nrow(Y)))
  colnames(Y) <- sprintf("cell_%s", seq_len(ncol(Y)))

  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE), cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- ~ group + cont1 + cont2

  res <- glm_gp(Y, design = design, col_data = annot)

  for (fn in c("predict", "residuals")) {
    fn.type <- get(sprintf("%s.%s", fn, class(res)))

    for (pred.type in as.list(eval(formals(fn.type)[["type"]]))) {
      test_that(sprintf("%s w/ type=%s (feat|obs).sub arguments works", fn, pred.type), {
        skip_if(pred.type == "randomized_quantile", "subsetting w/ randomized_quantile residual type is unsupported")

        fn.w <- function(...) fn.type(..., type = pred.type)
        ref <- fn.w(res)

        # using expect_all_equal instead of expect_equal because otherwise this generates
        # an obnoxious "number" of checks that pollutes reporter metrics
        expect_all_equal(
          vapply(
            list(
              fn.w(res, feat.sub = seq_len(nrow(Y))),
              fn.w(res, obs.sub = seq_len(ncol(Y))),
              fn.w(res, feat.sub = rownames(Y)),
              Reduce(rbind, lapply(seq_len(nrow(Y)), function(i) fn.w(res, feat.sub = i))),
              Reduce(cbind, lapply(seq_len(ncol(Y)), function(j) fn.w(res, obs.sub = j))),
              Reduce(rbind, lapply(rownames(Y), function(i) fn.w(res, feat.sub = i))),
              Reduce(cbind, lapply(colnames(Y), function(j) fn.w(res, obs.sub = j)))
            ),
            function(e) {
              cmp <- waldo::compare(ref, e, tolerance = 0)
              if (length(cmp) == 0) {
                ""
              } else {
                cmp
              }
            },
            character(1)
          ),
          ""
        )
      })
    }
  }
})
