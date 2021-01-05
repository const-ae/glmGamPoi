
#' @export
predict.glmGamPoi <- function(object, newdata = NULL,
                              type = c("link", "response"),
                              se.fit = FALSE, ...){

  type <- match.arg(type, c("link", "response"))
  if(is.null(newdata)){
    # Easy, just return mu
    Mu <- object$Mu
  }else{
    # Do something with newdata

  }

  if(type == "response"){
    fit <- Mu
  }else if(type == "link"){
    fit <- log(Mu)
  }


  if(se.fit){
    design <- object$model_matrix
    p_idxs <- seq_len(ncol(design))
    # This could be adapted to the quasi-GamPoi value
    scale <- 1
    se_fit <- t(vapply(seq_len(nrow(fit)), function(gene_idx){
      disp <- object$overdispersions[gene_idx]
      mu <- Mu[gene_idx, ]
      weights <- mu / (1 + mu * disp)
      weighted_Design <-  object$model_matrix * sqrt(weights)
      R <- qr.R(qr(weighted_Design))[p_idxs, p_idxs]
      XRinv <- design %*% qr.solve(R)
      sqrt(matrixStats::rowSums2(XRinv^2))
    }, FUN.VALUE = rep(0.0, ncol(fit))))

    if(type == "response"){
      se_fit <- se_fit * Mu
    }
    list(fit = fit, se.fit = se_fit, residual.scale = scale)
  }else{
    fit
  }
}


