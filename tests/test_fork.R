# jarl-ignore-file internal_function:>
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

# incredibly hacky way to source header files w/ generics
sourceH <- function(pth, monomorph = list(), export.all = TRUE, ...) {
  patch.list <- list(
    "#include <RcppEigen.h>" = "// [[Rcpp::depends(RcppEigen)]]",
    "#include <RcppArmadillo.h>" = "// [[Rcpp::depends(RcppArmadillo)]]"
  )

  lines <- purrr::reduce2(
    names(monomorph),
    monomorph,
    function(lines, trg, dst) {
      is_export <- export.all
      swp <- FALSE
      templ_pat <- stringr::regex(sprintf("(template ?<.*)( *class %s *,? *)(.*>)", trg))

      stringr::str_subset(
        purrr::map(lines, function(l) {
          if (stringr::str_detect(l, stringr::fixed("// [[Rcpp::export]]"))) {
            is_export <<- TRUE
            return(l)
          }
          if (is_export && stringr::str_detect(l, templ_pat)) {
            swp <<- TRUE
            return(stringr::str_replace(l, templ_pat, "\\1\\3"))
          }
          if (swp) {
            return(stringr::str_replace_all(l, stringr::fixed(trg), dst))
          }

          if (swp && (l == "}")) {
            is_export <<- export.all
            swp <<- FALSE
          }
          l
        }),
        stringr::regex("template ?<>"),
        negate = TRUE
      )
    },
    .init = readLines(pth)
  ) |>
    purrr::map(function(l) {
      # patch in header dependency comments
      for (inc in names(patch.list)) {
        if (stringr::str_detect(l, stringr::fixed(inc))) {
          return(paste0(c(patch.list[[inc]], l), collapse = "\n"))
        }
      }
      # incredibly hack way to find function declarations
      if (export.all && stringr::str_starts(l, stringr::regex("[^ ]+ [^ ]+\\("))) {
        return(paste0(c("// [[Rcpp::export]]", l), collapse = "\n"))
      }
      l
    })
  env <- new.env()
  Rcpp::sourceCpp(code = paste0(lines, collapse = "\n"), env = env, ...)

  env
}

obj.pl <- function(res, i, min_x = -10, max_x = 10, n = 1000, fn.ns = c("glmGamPoi", "glmGamPoi2"), cr.adjs = c(TRUE, FALSE)) {
  inp.range <- rev(seq(min_x, max_x, (max_x - min_x) / n))
  mu <- if (is.function(res[["Mu"]])) res[["Mu"]](i = i) else res[["Mu"]][i, ]
  mu[mu <= 1e-128] <- 1e-6
  y <- SummarizedExperiment::assay(res[["data"]])[i, ]

  df <- rbindlist(
    purrr::map(setNames(nm = fn.ns), function(ns) {
      rbindlist(
        purrr::map(
          list(
            "objective" = getFromNamespace("conventional_loglikelihood_fast", ns),
            "derivative" = getFromNamespace("conventional_score_function_fast", ns),
            "hessian" = getFromNamespace("conventional_deriv_score_function_fast", ns)
          ),
          function(f) {
            rbindlist(
              purrr::map(setNames(nm = cr.adjs), function(cr.adj) {
                data.table(
                  x = inp.range,
                  value = vapply(inp.range, function(e) f(y, mu, e, res[["model_matrix"]], cr.adj), numeric(1))
                )
              }),
              idcol = "cr.adj"
            )
          }
        ),
        idcol = "fn"
      )
    }),
    idcol = "ns"
  )
  max.obj <- df[fn == "objective", .SD[which.max(value)], by = .(ns, cr.adj)]
  ggplot(df) +
    aes(x = x, y = value, colour = ns, linetype = ns) +
    geom_line(alpha = 0.7) +
    geom_vline(
      inherit.aes = TRUE,
      data = rbind(
        max.obj,
        df[
          (fn == "derivative")
        ][
          max.obj[, .(x, ns, cr.adj)],
          on = .(ns, cr.adj)
        ][, .SD[which.min(abs(value) + 1e-03 * (x - i.x)^2)], by = .(ns, cr.adj), .SDcols = -c("i.x")]
      ),
      aes(xintercept = x),
      show.legend = FALSE
    ) +
    facet_grid(rows = vars(fn = forcats::fct_inorder(fn)), cols = vars(cr.adj), scales = "free_y", labeller = label_both) +
    labs(colour = NULL, linetype = NULL)
}


diff.dt <- function(ref, fork, drop.na = FALSE) {
  if (class(ref)[[1]] == "list") {
    stopifnot(class(fork) == "list")
    ref.de <- ref[[2]]
    fork.de <- fork[[2]]
    ref <- ref[[1]]
    fork <- fork[[1]]
  }

  DT <- melt(
    cbind(
      as.data.table(ref[["Beta"]], "gene"),
      overdispersion = ref[["overdispersions"]],
      deviance = ref[["deviances"]]
    ),
    value.name = "value.glmGamPoi",
    id.vars = "gene"
  )[
    melt(
      cbind(
        as.data.table(fork[["Beta"]], "gene"),
        overdispersion = fork[["overdispersions"]],
        deviance = fork[["deviances"]]
      ),
      id.vars = "gene",
      value.name = "value.glmGamPoi2"
    ),
    on = .(gene, variable)
  ]

  if (exists("ref.de")) {
    DT <- rbind(
      DT,
      melt(
        setnames(as.data.table(ref.de), old = "name", new = "gene"),
        id.vars = "gene",
        measure.vars = c("adj_pval", "f_statistic", "lfc", "pval"),
        value.name = "value.glmGamPoi"
      )[
        melt(
          setnames(as.data.table(fork.de), old = "name", new = "gene"),
          id.vars = "gene",
          measure.vars = c("adj_pval", "f_statistic", "lfc", "pval"),
          value.name = "value.glmGamPoi2"
        ),
        on = .(gene, variable)
      ]
    )
  }

  if (drop.na) {
    DT <- na.omit(DT)
  }

  split(DT, by = "variable", keep.by = FALSE)
}

stack.envs <- function(x, y) as.environment(append(as.list(x), as.list(y)))

fn.w.env <- function(fn, env) {
  env.c <- list2env(as.list(environment(fn)), parent = parent.env(environment(fn)))
  purrr::iwalk(as.list(env), function(x, nm) assign(nm, x, env.c))
  environment(fn) <- env.c

  fn
}

# using clang as default instead of gcc because clang's error messages are more readable
load.fork <- function(cc = "clang++", quiet = TRUE, compile = TRUE, attach = FALSE) {
  if (isTRUE(compile)) {
    pkgbuild::clean_dll()
  }
  withr::with_makevars(
    assignment = "=",
    c(
      CXX17 = sprintf("%s -fuse-ld=lld -flto -ffat-lto-objects -fopenmp", cc),
      CXX17STD = "-std=c++17",
      CXX17FLAGS = "-O3 -march=native",
      LDFLAGS = "-L/usr/lib64/R/lib -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,-z,pack-relative-relocs"
    ),
    pkgload::load_all(quiet = quiet, compile = compile, attach = attach, debug = FALSE)
  )
}
load.fork()

GLMGAMPOI.REF.ENV <- stack.envs(
  getNamespace("glmGamPoi"),
  # original package does not provide the fisher scoring step methods, so we compile/export them manually
  sourceH("https://raw.githubusercontent.com/const-ae/glmGamPoi/refs/heads/devel/inst/include/fisher_scoring_steps.h", list(NumericType = "double"))
)
GLMGAMPOI.FORK.ENV <- getNamespace("glmGamPoi2")
PERF_OPTIM_CONF <- list(delay_mu = TRUE, offset_as_vec = TRUE, use_nr_overdisp_impl = TRUE)


sce.diff <- function(sce, design, tol = 1e-8, ridge_vec = function(mm) c(5e-11 / nrow(mm), seq_len(ncol(mm) - 1)), ref.env = GLMGAMPOI.REF.ENV, fork.env = GLMGAMPOI.FORK.ENV) {
  Y <- SummarizedExperiment::assay(sce)
  col_data <- SummarizedExperiment::colData(sce)

  off <- glmGamPoi:::combine_size_factors_and_offset(0, "normed_sum", Y)[["offset_matrix"]]
  mm <- glmGamPoi:::handle_design_parameter(design, Y, col_data, NULL)[["model_matrix"]]
  betas <- glmGamPoi:::estimate_betas_roughly(Y, mm, off)
  disps <- glmGamPoi:::estimate_dispersions_roughly(Y, mm, off)
  ridge <- glmGamPoi:::handle_ridge_penalty_parameter(if (is.function(ridge_vec)) ridge_vec(mm) else ridge_vec, mm)
  disc.cols <- which(vapply(seq_len(ncol(mm)), function(j) length(unique(mm[, j])) < nrow(mm), logical(1)))

  ridge_t <- attr(ridge, "target")
  ridge <- diag(ridge, nrow = length(ridge))
  attr(ridge, "target") <- ridge_t

  mu <- exp(Matrix::tcrossprod(betas, mm) + off)

  fns <- list(
    compute_gp_deviance_residuals_matrix = function() compute_gp_deviance_residuals_matrix(as.matrix(Y), mu, disps),
    fitBeta_fisher_scoring = function() fitBeta_fisher_scoring(beachmat::initializeCpp(Y), mm, beachmat::initializeCpp(exp(off)), disps, betas, NULL, 1e-8, 1e5, 1000),
    fitBeta_fisher_scoring_ridge = function() fitBeta_fisher_scoring(beachmat::initializeCpp(Y), mm, beachmat::initializeCpp(exp(off)), disps, betas, ridge, 1e-8, 1e5, 1000),
    fitBeta_fisher_scoring_diagonal = function() fitBeta_diagonal_fisher_scoring(beachmat::initializeCpp(Y), mm, beachmat::initializeCpp(exp(off)), disps, betas, 1e-8, 1e5, 1000),
    estimate_betas_fisher_scoring = function() estimate_betas_fisher_scoring(Y, mm, off, disps, betas, NULL),
    estimate_betas_fisher_scoring_ridge = function() estimate_betas_fisher_scoring(Y, mm, off, disps, betas, ridge),
    estimate_betas_optim = function() estimate_betas_optim(Y, mm, off, disps, betas, NULL),
    estimate_betas_optim_ridge = function() estimate_betas_optim(Y, mm, off, disps, betas, ridge),
    overdispersion_mle = function() overdispersion_mle(Y, mu, mm, FALSE, FALSE),
    overdispersion_mle_wcr = function() overdispersion_mle(Y, mu, mm, TRUE, FALSE),
    overdispersion_mle_global = function() overdispersion_mle(Y, mu, mm, FALSE, TRUE),
    overdispersion_mle_global_wcr = function() overdispersion_mle(Y, mu, mm, TRUE, TRUE),
    overdispersion_mle_nr = function() {
      (if ("use_nr_overdisp_impl" %in% names(formals(overdispersion_mle))) {
        function(...) overdispersion_mle(..., use_nr_overdisp_impl = TRUE)
      } else {
        overdispersion_mle
      })(Y, mu, mm, FALSE, FALSE)[["estimate"]]
    },
    overdispersion_mle_nr_wcr = function() {
      (if ("use_nr_overdisp_impl" %in% names(formals(overdispersion_mle))) {
        function(...) overdispersion_mle(..., use_nr_overdisp_impl = TRUE)
      } else {
        overdispersion_mle
      })(Y, mu, mm, TRUE, FALSE)[["estimate"]]
    }
  )

  # if we have columns w/ discrete values that can be split into groups
  # also run diff for estimate_betas_group_wise
  grps <- glmGamPoi:::get_groups_for_model_matrix(mm[, disc.cols, drop = FALSE])
  if (!is.null(grps)) {
    fns <- append(
      fns,
      list(
        estimate_betas_group_wise = function() {
          estimate_betas_group_wise(
            Y,
            off,
            disps,
            beta_mat_init = betas[, disc.cols, drop = FALSE],
            groups = grps,
            model_matrix = mm[, disc.cols, drop = FALSE]
          )
        }
      )
    )
  } else {
    warning("diff for estimate_betas_group_wise not run because no component of model matrix can be split into discrete groups")
  }

  i.fns <- purrr::list_flatten(list(
    purrr::list_flatten(purrr::map(setNames(nm = c("conventional_loglikelihood_fast", "conventional_score_function_fast", "conventional_deriv_score_function_fast")), function(fn) {
      purrr::map(c("with_cr" = TRUE, "w/o_cr" = FALSE), function(cr.adj) function(i) get(fn)(Y[i, ], mu[i, ], log(disps[[i]]), mm, cr.adj))
    })),
    purrr::imap(
      list(
        "fisher_scoring_qr_step" = function(i) list(),
        "fisher_scoring_qr_ridge_step" = function(i) list(ridge, ridge_t, betas[i, ]),
        "fisher_scoring_diagonal_step" = function(i) list()
      ),
      function(extra.args, fn) {
        function(i) do.call(get(fn), append(list(mm, Y[i, ], c(mu[i, ]), c(mu[i, ] * disps[[i]])), extra.args(i)))
      }
    ),
    list(compute_gp_deviance_sum = function(i) compute_gp_deviance_sum(Y[i, ], mu[i, ], disps[[i]]))
  ))

  cmp.fn <- function(fn, nm) {
    if (!missing(nm)) {
      message(sprintf("running test for - %s\n", nm))
    }
    fn.ref <- fn.w.env(fn, ref.env)
    fn.new <- fn.w.env(fn, fork.env)
    if ((length(formals(fn)) != 0)) {
      fn.ref.o <- fn.ref
      fn.new.o <- fn.new
      fn.ref <- function() purrr::map(seq_len(nrow(Y)), function(i) fn.ref.o(i))
      fn.new <- function() purrr::map(seq_len(nrow(Y)), function(i) fn.new.o(i))
    }
    perf <- bench::mark(
      ref = {
        ref <- fn.ref()
      },
      new = {
        new <- fn.new()
      },
      check = FALSE,
      iterations = 1,
      filter_gc = FALSE
    )

    out <- waldo::compare(ref, new, tolerance = tol)
    attr(out, "ref") <- ref
    attr(out, "new") <- new
    attr(out, "bench") <- perf

    out
  }

  purrr::compact(append(
    # per gene functions
    purrr::imap(i.fns, cmp.fn),
    # whole data functions
    purrr::imap(fns, cmp.fn)
  ))
}
