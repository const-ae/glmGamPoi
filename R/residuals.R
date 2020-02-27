
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
#'     \deqn{dev = 2 * (1/θ log((1 + µ θ) / (1 + y θ)) - y log((µ + y µ θ) / (y + y µ θ)))}
#'   and the residual accordingly is
#'     \deqn{res = sign(y - µ) sqrt(dev).}
#'   }
#'   \item{pearson}{The Pearson residual is \eqn{res = (y - µ) / sqrt(µ + µ^2 θ)}}
#'   \item{random_quantile}{The random quantile residual was originally developed
#'   by Dunn & Smyth, 1995. Please see that publicaton or [statmod::qresiduals()] for more
#'   information.}
#'   \item{working}{The working residuals are \eqn{res = (y - µ) / µ}.}
#'   \item{response}{The response residuals are \eqn{res = y - µ}}
#' }
#'
#' @return a matrix with the same size as `Y`. If `Y` is a `DelayedArray` than the
#'   result will be as well.
#'
#' @seealso [glm_gp()] and [stats::residuals.glm()]
#' @export
residuals.glmGamPoi <- function(object, Y, type = c("deviance", "pearson", "random_quantile", "working", "response"), ...){
  type <- match.arg(type, c("deviance", "pearson", "random_quantile", "working", "response"))
  if(type == "deviance"){
    make_resid_hdf5_mat <- is(Y, "DelayedMatrix") && is(DelayedArray::seed(Y), "HDF5ArraySeed")
    if(! make_resid_hdf5_mat){
      compute_gp_deviance_residuals_matrix(Y, object$Mu, object$overdispersion)
    }else{
      delayed_matrix_apply_block(Y, object$Mu, object$overdispersion, compute_gp_deviance_residuals_matrix)
    }
  }else if(type == "pearson"){
    (Y - object$Mu) / sqrt(object$Mu + multiply_vector_to_each_column(object$Mu^2, object$overdispersion))
  }else if(type == "random_quantile"){
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
  p <- 1 / (1 + multiply_vector_to_each_column(Mu, dispersion))
  tmp <- vapply(seq_len(ncol(Y)), function(idx){
    a <- pbeta(p[, idx], 1/dispersion, Y[, idx])
    b <- pbeta(p[, idx], 1/dispersion, Y[, idx] + 1)
    rbind(a, b)
  }, FUN.VALUE = matrix(0, ncol = nrow(Y), nrow = 2))

  u <- runif(n = length(Y), min = tmp[1,,], max = tmp[2,,])
  res <- qnorm(u)
  array(res, dim(Y))
}



