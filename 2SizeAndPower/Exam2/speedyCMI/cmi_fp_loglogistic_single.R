#' Single, fully parametric conditional mean imputation for a right-censored covariate (log-logistic distribution) with conditional means following Equation (6)
#'
#' Single, fully parametric conditional mean imputation for a right-censored covariate using a log-logistic model to estimate the conditional survival function and then uses an analytic solution to compute conditional means, as in Equation (5) of the manuscript.
#'
#' @param imputation_model imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param W character, column name for observed values of the censored covariate
#' @param Delta character, column name for censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and any other variables in \code{imputation_model}.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#'
#' @return
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#' \item{aic}{Akaike information criterion (AIC) from the \code{imputation_model} fit.}
#' \item{bic}{Bayesian information criterion (AIC) from the \code{imputation_model} fit.}
#' \item{coefficients}{Vector of coefficients from the \code{imputation_model} fit.}
#' \item{scale}{Scale from the \code{imputation_model} fit.}
#'
#' @importFrom survival survreg
#' @importFrom survival Surv
#' @importFrom survival psurvreg

cmi_fp_loglogistic_single = function(imputation_model, W, Delta, data, maxiter = 100) {
  # Initialize imputed values
  data$imp = data[, W] ## start with imp = W

  # Fit AFT imputation model for X ~ Z
  fit = survreg(formula = imputation_model,
                data = data,
                dist = "loglogistic",
                maxiter = maxiter)

  # Calculate linear predictor for AFT imputation model
  lp = fit$linear.predictors ## linear predictors

  # Transform parameters to agree with paper's parameterization
  alpha =  1 / fit$scale
  lambda = exp(lp)

  # Create an indicator variable for being uncensored
  uncens = data[, Delta] == 1

  # Use closed-form to compute the conditional means
  ## Save quantities for use in formula
  beta_shape1 = (alpha - 1) / alpha
  beta_shape2 = 1 / alpha
  z = (1 + (data[which(!uncens), W] / lambda[which(!uncens)]) ^ alpha) ^ (- 1)
  B = beta(a = beta_shape1,
           b = beta_shape2)
  incB = B * pbeta(q = z,
                   shape1 = beta_shape1,
                   shape2 = beta_shape2,
                   lower.tail = TRUE) ## CDF of a beta
  data[which(!uncens), "imp"] = data[which(!uncens), W] +
    lambda[which(!uncens)] / alpha * (1 / z) * incB

  # Return input dataset with appended column imp containing imputed values
  return_list = list(imputed_data = data,
                     code = !any(is.na(data$imp)),
                     aic = AIC(object = fit),
                     bic = BIC(object = fit),
                     coefficients = fit$coefficients,
                     scale = fit$scale)
  return(return_list)
}
