# glmGamPoi(2/fork) explainer

code in this repo after commit [`8d9a005`](https://github.com/goeva-lab/glmGamPoi/compare/8d9a005..perf_optim_eigen) represents work done to optimize glmGamPoi for enhanced compatibility w/ large-scale sparse/scRNA-seq data.

an explicit aim of this project has been to allow for upstream-ing of the relevant code into the upstream ([`const-ae/glmGamPoi`](https://github.com/const-ae/glmGamPoi)) repo.

> [!IMPORTANT]
> to allow for easier testing / comparison against the original package to catch regressions, the package & main exported class have been renamed to `glmGamPoi2`.
> this change is intended to be temporary and to be reverted should we merge into upstream.
> moreover, the `R CMD CHECK` workflow for the repository has been disabled

## user-facing changes

the main user-facing consequence of these changes is the addition of the `perf_optim` list parameter to `glm_gp`, which controls various steps taken for optimization.

parameters for this option include:

- `do_parallel`: allows setting core count for parallelizing various steps

- `offset_as_vec`: storing the offsets as a vector in the case where offsets are set uniquely in a per-cell manner (which is true for the overwhelming majority of use-cases).
  
  **important**: this leads to a **_breaking_** change in output format, as the result's `Offset` field is no longer a matrix of shape `(nrow(Y), ncol(Y))` but instead a vector of length `ncol(Y)`.

- `delay_mu`: never explicitly form the Mu (prediction) matrix, instead returning a function which can return slices of the matrix as necessary

  **important**: this leads to a **_breaking_** change in output format, as the result's `Mu` field is no longer a matrix of shape `(nrow(Y), ncol(Y))` but instead a function which returns this matrix (or slice thereof) when called.

- `use_nr_overdisp_impl`: utilize an alternative implementation for overdispersion estimation, which is written in C++ and is thread safe, allowing for parallelization of this step

  **note**: due to being specificially aware of the possibility between inconsistencies between the objective and analytic gradient functions for overdispersion estimation, it often succeeds in finding cox-reid adjusted estimations of overdispersion which the original implementation fails to find, leading to non-insignificant changes in output values (including resulting model betas) when the cox-reid adjustment is enabled.

  **warning**: while initial testing shows that the method overall yields highly concordant results w/ the reference overdispersion estimation method, further testing/benchmarking is still ongoing and thus it should be considered to be in an unstable/beta status.

  **note**: switching to using the NR-based overdispersion estimator causes one (1) test to fail which previously passed:
  - [`test-overdispersion.R:149:3`](./tests/testthat/test-overdispersion.R#L149): yields a larger than 1e08 (1e10) value for "weird" input

by default, none of these options are enabled, as an explicit goal is that no changes should affect user outputs after having upgraded unless explicitly desired.

## internal changes

internally, these changes have required relatively significant changes to package internals (R and C++), with the main changes being:

- a re-write of C++ internals to utilize [Eigen](https://libeigen.gitlab.io/) (via [RcppEigen](https://cran.r-project.org/web/packages/RcppEigen/index.html)) instead of [Armadillo](https://arma.sourceforge.net/) (via [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html))for algebraic/matrix operations
  - the reasoning for this is that Eigen makes stronger guarantees around thread-safety of its data structures ([c.f.](https://libeigen.gitlab.io/eigen/docs-nightly/TopicMultiThreading.html)) compared to Armadillo, which was necessary for this project given the aim to enable cross-gene parallelization of beta/overdispersion estimation steps
  - while Armadillo offers some inherent parallelization (by using OPENMP for specific operations), it appears heavily dependent on system configuration (e.g. OPENMP presence, etc.), leading to difficulties w/ consistent benchmarking, further motivating this design decision, given our interest in direct control of parallelization levels
- a general refactoring of the C++ internals to move R-unaware logic into header files in [`inst/include`](./inst/include/)
- sets a strict clamp of +/- 1e8 to fisher step size when fitting model betas to avoid huge values (this might be reverted, true necessity unclear)
- adds [LBFGSpp](https://github.com/yixuan/LBFGSpp/) as a git submodule (kept in [`inst/include/LBFGSpp`](./inst/include/LBFGSpp/)), to allow for from-C++ fitting of model betas when fisher scoring based routine fails
- vendors an adapted version of the `Brent_fmin` routine ([c.f.](https://svn.r-project.org/R/trunk/src/library/stats/src/optimize.c)) from the R source code (kept in [`inst/include/opt_max.h`](./inst/include/opt_max.h)), to allow for from-C++ fitting of model beta when the newton-raphson based routine fails
- threading the `perf_optim` relevant variable(s) throughout package internals as necessary
- adding control flow to deal with alternate representations of Mu (prediction) and offset matrices

> [!NOTE]
> the latter two points non trivially increase (cyclomatic) complexity in the package (e.g. `is.vector(offset_matrix)`/`is.function(Mu)` checks in various function preludes, the branches in `overdispersion_mle`, etc.), worsening readability/parsability

a test suite for relevant changes has also been added at [`test-perf_optim.R`](./tests/testthat/test-perf_optim.R).

moreover, a file to drive benchmarking / regression analysis is provided at [`test_fork.R`](./tests/test_fork.R).
it will most likely be deleted and maybe moved to a separate repo prior to upstreaming.

## miscleanea

### dismissed directions

other notable experiments were undertaken before being dismissed as non worthwhile:

- `mu_dropout_thresh`: an attempt to improve memory usage by adding a dropout threshold to Mu matrix construction, allowing the matrix to be represented sparsely (reverted in `9b4cd99`)
  - this was reverted because it lead to a noticeable drop in accuracy while leading to performance regressions, and it significantly increased the internal complexity of `calculate_mu`
  - moreover, decreases in memory use appeared negligeable
- `cast_dgC_Y_to_dgR`: an attempt to improve memory access patterns in main functions by internally representing the counts matrix (when sparse) to row-major format, given that all internal accesses are by-row
  - does not appear to lead to any performance gains, instead seems to cause notable regressions in runtime

### singularity definition file

for easier benchmarking and evaluation, a singularity definition file [`./glmgp.def`](./glmgp.def) has been made.
it generates an environment w/ both forked (under the `glmGamPoi2` name) original (as fetched from bioconductor) packages available.

the build instructions are tuned to generate maximally-optimized builds (i.e. `-O3`) and leveraging all available CPU instructions available (i.e. `-march=native`).

> [!WARNING]
> this means that the resulting `sif` file should **_not_** be considered portable across CPUs (even of the same base architecture), and should be built on a per-node/machine basis
