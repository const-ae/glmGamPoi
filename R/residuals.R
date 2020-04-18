
#' Extract Residuals of Gamma Poisson Model
#'
#' @param object a fit of type `glmGamPoi`. It is usually produced with a call to
#'   `glm_gp()`.
#' @param Y any matrix-like object (e.g. `matrix()`, `DelayedArray()`, `HDF5Matrix()`) with
#'   one column per sample and row per gene.
#' @param type the type of residual that is calculated. See details for more information.
#'   Default: `"deviance"`.
#' @param ... currently ignored.
#'
#' @details
#' This method can calculate a range of different residuals:
#' \describe{
#'   \item{deviance}{The deviance for the Gamma-Poisson model is
#'     \deqn{dev = 2 * (1/theta * log((1 + m * theta) / (1 + y * theta)) - y log((m + y * theta) / (y + y * m * theta)))}
#'   and the residual accordingly is
#'     \deqn{res = sign(y - m) sqrt(dev).}
#'   }
#'   \item{pearson}{The Pearson residual is \eqn{res = (y - m) / sqrt(m + m^2 * theta)}}
#'   \item{randomized_quantile}{The randomized quantile residual was originally developed
#'   by Dunn & Smyth, 1995. Please see that publication or [statmod::qresiduals()] for more
#'   information.}
#'   \item{working}{The working residuals are \eqn{res = (y - m) / m}.}
#'   \item{response}{The response residuals are \eqn{res = y - m}}
#' }
#'
#' @return a matrix with the same size as `Y`. If `Y` is a `DelayedArray` than the
#'   result will be as well.
#'
#' @seealso [glm_gp()] and `stats::residuals.glm()
#' @export
residuals.glmGamPoi <- function(object, Y, type = c("deviance", "pearson", "randomized_quantile", "working", "response"), ...){
  type <- match.arg(type, c("deviance", "pearson", "randomized_quantile", "working", "response"))
  if(type == "deviance"){
    make_resid_hdf5_mat <- is(Y, "DelayedMatrix") && is(DelayedArray::seed(Y), "HDF5ArraySeed")
    if(! make_resid_hdf5_mat){
      compute_gp_deviance_residuals_matrix(Y, object$Mu, object$overdispersion)
    }else{
      delayed_matrix_apply_block(Y, object$Mu, object$overdispersion, compute_gp_deviance_residuals_matrix)
    }
  }else if(type == "pearson"){
    (Y - object$Mu) / sqrt(object$Mu + multiply_vector_to_each_column(object$Mu^2, object$overdispersion))
  }else if(type == "randomized_quantile"){
    make_resid_hdf5_mat <- is(Y, "DelayedMatrix") && is(DelayedArray::seed(Y), "HDF5ArraySeed")
    if(! make_resid_hdf5_mat){
      qres.gampoi(Y, object$Mu, object$overdispersion)
    }else{
      delayed_matrix_apply_block(Y, object$Mu, object$overdispersion, qres.gampoi)
    }
  }else if(type == "working"){
    (Y - object$Mu) / object$Mu
  }else if(type == "response"){
    Y - object$Mu
  }
}


qres.gampoi <- function(Y, Mu, dispersion){
  tmp <- vapply(seq_len(ncol(Y)), function(idx){
    a <- pnbinom(Y[, idx] - 1, mu = Mu[, idx], size = 1/dispersion)
    b <- pnbinom(Y[, idx], mu = Mu[, idx], size = 1/dispersion)
    rbind(a, b)
  }, FUN.VALUE = matrix(0, ncol = nrow(Y), nrow = 2))

  u <- runif(n = length(Y), min = pmin(tmp[1,,], tmp[2,,]), max = pmax(tmp[1,,], tmp[2,,]))
  res <- qnorm(u)

  # Get rid of infinite value by sampling on the log scale
  # This is technically wrong, but shouldn't matter.
  inf_res <- which(is.infinite(res))
  a_alt <- pnbinom(Y[inf_res] - 1, mu = Mu[inf_res], size = 1/dispersion[(inf_res - 1) %% nrow(Y) + 1],
                   log.p = TRUE, lower.tail = Mu[inf_res] > Y[inf_res])
  b_alt <- pnbinom(Y[inf_res], mu = Mu[inf_res], size = 1/dispersion[(inf_res - 1) %% nrow(Y) + 1],
                   log.p = TRUE, lower.tail = Mu[inf_res] > Y[inf_res])
  u_alt <- runif(n = length(inf_res), min = pmin(a_alt, b_alt), max = pmax(a_alt, b_alt))
  res[inf_res] <- qnorm(u_alt, log.p = TRUE, lower.tail = Mu[inf_res] > Y[inf_res])

  array(res, dim(Y))
}



