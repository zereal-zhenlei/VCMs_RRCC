#' Single, fully parametric imputation for a right-censored covariate with conditional means following Equation (1)
#'
#' Single, fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then integrates over the estimated survival function from \code{W} to \code{Inf} to compute conditional means, as in Equation (1) of the manuscript.
#'
#' @param imputation_model imputation model formula (or coercible to formula) passed through to \code{survreg}, a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See \code{survreg} documentation for more details.
#' @param dist imputation model distribution passed through to \code{survreg}. See \code{survreg} documentation for more details.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and any other variables in \code{imputation_model}.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#'
#' @return A list containing:
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#' \item{aic}{Akaike information criterion (AIC) from the \code{imputation_model} fit.}
#' \item{bic}{Bayesian information criterion (BIC) from the \code{imputation_model} fit.}
#' \item{coefficients}{Vector of coefficients from the \code{imputation_model} fit.}
#' \item{scale}{Scale from the \code{imputation_model} fit.}
#'
#' @importFrom survival survreg
#' @importFrom survival Surv
#' @importFrom survival psurvreg

cmi_fp_original_single = function(imputation_model, dist, data, maxiter = 100) {
  # Extract variable names from imputation_model
  W = all.vars(imputation_model)[1] ## censored covariate
  Delta = all.vars(imputation_model)[2] ## corresponding event indicator

  # Initialize imputed values
  data$imp = data[, W] ## start with imp = W

  # Fit AFT imputation model for X ~ Z
  fit = survreg(formula = imputation_model,
                data = data,
                dist = dist,
                maxiter = maxiter)

  # Calculate linear predictor for AFT imputation model
  lp = fit$linear.predictors ## linear predictors from the survreg fit

  # Create an indicator variable for being uncensored
  uncens = data[, Delta] == 1

  # Use adaptive quadrature to estimate
  ## integral from X = 0 to X = Wi
  int_surv_W_to_Inf = sapply(
    X = which(!uncens),
    FUN = function(i) {
      integrate(f = function(x, mean, scale, distribution) {
        1 - psurvreg(q = x, mean = mean, scale = scale, distribution = distribution)
        },
                lower = data[i, W],
                upper = Inf,
                mean = lp[i],
                scale = fit$scale,
                distribution = dist)$value
    }
  )

  # Calculate S(W|Z)
  surv = 1 - psurvreg(q = data[which(!uncens), W],
                      mean = lp[which(!uncens)],
                      scale = fit$scale,
                      distribution = dist)

  # Calculate MRL(W) = int_surv / S(W|Z)
  est_mrl = int_surv_W_to_Inf / surv

  # Calculate imputed value E(X|X>W,Z) = W + MRL(W)
  data[which(!uncens), "imp"] = data[which(!uncens), W] + est_mrl

  # Return input dataset with appended column imp containing imputed values
  return_list = list(imputed_data = data,
                     code = !any(is.na(data$imp)),
                     aic = AIC(object = fit),
                     bic = BIC(object = fit),
                     coefficients = fit$coefficients,
                     scale = fit$scale)
  return(return_list)
}
