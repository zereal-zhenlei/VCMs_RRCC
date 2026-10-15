#' Fit analysis model to each imputed dataset and pool the results using Rubin's rules
#'
#' Fit analysis model to each imputed dataset and pool the results using Rubin's rules
#'
#' @param analysis_model imputation model formula (or coercible to formula) passed through to \code{survreg}, a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See \code{survreg} documentation for more details.
#' @param mult_imp a list (or list of lists) of imputed datasets returned from one of the \code{cmi_fp} functions
#'
#' @return A dataframe containing the pooled coefficients and standard errors for \code{analysis_model}

#'
#' @export
cmi_fp_pool_fit = function(analysis_model, mult_imp) {
  # Define the number of imputations
  B = length(mult_imp)

  # Fit analysis model to each imputed dataset
  mult_fit = do.call(what = rbind,
                     args = lapply(X = mult_imp,
                                   FUN = function(l) {
                                     summary(lm(formula = as.formula(analysis_model),
                                                data = l$imputed_data))$coefficients
                                   }
                     )
  )

  # Define the dimension of the coefficient vector
  p = nrow(summary(lm(formula = as.formula(analysis_model),
                      data = mult_imp[[1]]$imputed_data))$coefficients)

  # Separate coefficient and variances estimates
  beta_b = matrix(data = mult_fit[, 1],
                  nrow = B,
                  byrow = TRUE)
  vbeta_b = matrix(data = mult_fit[, 2],
                   nrow = B,
                   byrow = TRUE)

  # Pool coefficient and variance estimates
  beta_pooled = colMeans(beta_b)
  beta_pooled_rep = matrix(data = beta_pooled,
                           nrow = B,
                           ncol = p,
                           byrow = TRUE)
  vbeta_within = colMeans(vbeta_b)
  vbeta_between = colSums((beta_b - beta_pooled_rep) ^ 2) / (B - 1)
  vbeta_pooled = vbeta_within + vbeta_between

  # Create table of estimates
  res = data.frame(Est = beta_pooled,
                   SE = sqrt(vbeta_pooled))
  res$t = res$Est / res$SE ## t-statistic
  res$p = 2 * pnorm(q = abs(res$t), mean = 0, sd = 1, lower.tail = FALSE) ## two-sided p-value
  rownames(res) = rownames(summary(lm(formula = as.formula(analysis_model),
                                      data = mult_imp[[1]]$imputed_data))$coefficients) ## take row names from one of the fits

  # Return it
  return(res)
}
