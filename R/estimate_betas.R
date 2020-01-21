
#' Make a quick first guess where reasonable beta would be
#'
#' @importFrom Matrix t
#'
#' @keywords internal
estimate_betas_roughly <- function(Y, model_matrix, offset_matrix, pseudo_count = 1){
  qrx <- qr(model_matrix)
  Q <- qr.Q(qrx)
  R <- qr.R(qrx)

  norm_log_count_mat <- t(log((Y / exp(offset_matrix) + pseudo_count)))
  t(solve(R, t(Q) %*% norm_log_count_mat))
}


#' Estimate the Betas for Fixed Dispersions
#'
#'
#' @keywords internal
estimate_betas <- function(Y, model_matrix, offset_matrix,
                           dispersions = NULL, beta_mat_init = NULL){

  if(is.null(beta_mat_init)){
    beta_mat_init<- estimate_betas_roughly(Y, model_matrix, offset_matrix = offset_matrix)
  }
  if(is.null(dispersions)){
    mu_mat <- exp(beta_mat_init %*% t(model_matrix) + offset_matrix)
    disps <- estimate_overdispersions(Y, mu_mat, model_matrix, do_cox_reid_adjustment = TRUE)
  }

  useWeights <- FALSE
  betaTol <- 1e-8
  maxit <- 100
  # minmu <- 0.5
  minmu <- 1e-6

  # if(engine == "DESeq2"){
  #   betaRes <- DESeq2:::fitBetaWrapper(ySEXP = as.matrix(Y), xSEXP = model_matrix,
  #                                      nfSEXP = exp(offset_matrix),
  #                                      alpha_hatSEXP = disps_fixed,
  #                                      beta_matSEXP = beta_mat_init,
  #                                      lambdaSEXP = lambdaNatLogScale,
  #                                      weightsSEXP = weights,
  #                                      useWeightsSEXP = useWeights,
  #                                      tolSEXP = betaTol, maxitSEXP = maxit,
  #                                      useQRSEXP=useQR, minmuSEXP=minmu)
  # }else{
  #   betaRes <- edgeR::glmFit.default(as.matrix(Y), design = model_matrix,
  #                                    dispersion = disps_fixed, offset = offset_matrix[1,],
  #                                    prior.count = 0, weights=NULL,
  #                                    start = beta_mat_init)
  #   betaRes$beta_mat <- coef(betaRes)
  # }
  # List fitBeta(SEXP ySEXP, SEXP xSEXP, SEXP nfSEXP, SEXP alpha_hatSEXP, SEXP beta_matSEXP,
  #              SEXP tolSEXP, SEXP maxitSEXP, SEXP minmuSEXP) {
  fitBeta(Y, model_matrix, exp(offset_matrix), dispersions, beta_mat_init, tolSEXP = betaTol,
          maxitSEXP = maxit, minmuSEXP = minmu)

  betaRes
}
