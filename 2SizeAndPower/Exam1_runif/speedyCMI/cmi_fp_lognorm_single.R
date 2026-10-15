#' Single, fully parametric conditional mean imputation for a right-censored covariate (log-normal distribution) with conditional means following Equation (10)
#'
#' Single, fully parametric conditional mean imputation for a right-censored covariate using a log-normal model to estimate the conditional survival function and then uses an analytic solution to compute conditional means, as in Equation (10) of the manuscript.
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
#' \item{bic}{Bayesian information criterion (BIC) from the \code{imputation_model} fit.}
#' \item{coefficients}{Vector of coefficients from the \code{imputation_model} fit.}
#' \item{scale}{Scale from the \code{imputation_model} fit.}
#'
#' @importFrom survival survreg
#' @importFrom survival Surv
#' @importFrom survival psurvreg

cmi_fp_lognorm_single = function(imputation_model, W, Delta, data, maxiter = 100) {
  # Extract variable names from imputation_model
  W = all.vars(imputation_model)[1] ## censored covariate
  Delta = all.vars(imputation_model)[2] ## corresponding event indicator

  # Initialize imputed values
  data$imp = data[, W] ## start with imp = W

  # Fit AFT imputation model for X ~ Z
  fit = survreg(formula = imputation_model,
                data = data,
                dist = "lognormal",
                maxiter = maxiter)

  # Calculate linear predictor for AFT imputation model
  mu = fit$linear.predictors ## linear predictors from the survreg fit = location
  sigma = fit$scale ## scale from the survreg fit = scale

  # Create an indicator variable for being uncensored
  uncens = data[, Delta] == 1

  # Use closed-form to compute the conditional means
  ## Save quantities for use in formula
  stand_W = (log(data[, W]) - mu) / sigma ### standardized W values (to go into standard normal CDF)
  mult = exp(mu + (sigma ^ 2) / 2) ### multiplier
  num = pnorm(q = (stand_W - sigma), mean = 0, sd = 1, lower.tail = FALSE) ### numerator
  denom = pnorm(q = stand_W, mean = 0, sd = 1, lower.tail = FALSE) ### denominator

  ## Put them together to compute conditional means
  condl_means = mult * num / denom

  # Replace censored values with their conditional means
  data[which(!uncens), "imp"] = condl_means[which(!uncens)]

  # Return input dataset with appended column imp containing imputed values
  return_list = list(imputed_data = data,
                     code = !any(is.na(data$imp)),
                     aic = AIC(object = fit),
                     bic = BIC(object = fit),
                     coefficients = fit$coefficients,
                     scale = fit$scale)
  return(return_list)
}
