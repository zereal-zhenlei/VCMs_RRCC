#' Single, fully parametric conditional mean imputation for a right-censored covariate (exponential distribution) with conditional means following Equation (6)
#'
#' Single, fully parametric conditional mean imputation for a right-censored covariate using an exponential model to estimate the conditional survival function and then uses an analytic solution to compute conditional means, as in Equation (5) of the manuscript.
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

cmi_fp_expo_single = function(imputation_model, W, Delta, data, maxiter = 100) {
  # Initialize imputed values
  data$imp = data[, W] ## start with imp = W

  # Fit AFT imputation model for X ~ Z
  fit = survreg(formula = imputation_model,
                data = data,
                dist = "exponential",
                maxiter = maxiter)

  # Calculate linear predictor for AFT imputation model
  lp = fit$linear.predictors ## linear predictors

  # Transform parameters to agree with R's weibull parameterization
  expo_scale = exp(lp)

  # Create an indicator variable for being uncensored
  uncens = data[, Delta] == 1

  # Use closed-form to compute the conditional mean
  ## Transform parameters to agree with paper's parameterization
  lambda = 1 / expo_scale

  # Use closed-form to compute the conditional means
  data[which(!uncens), "imp"] = data[which(!uncens), "imp"] + 1 / lambda[which(!uncens)]

  # Return input dataset with appended column imp containing imputed values
  return_list = list(imputed_data = data,
                     code = !any(is.na(data$imp)),
                     aic = AIC(object = fit),
                     bic = BIC(object = fit),
                     coefficients = fit$coefficients,
                     scale = fit$scale)
  return(return_list)
}
