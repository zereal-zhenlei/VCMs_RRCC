#' Single, fully parametric imputation for a right-censored covariate with conditional means following Equation (15)
#'
#' Single, fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then integrates over the estimated survival function from \code{0} to \code{W} to compute conditional means, as in Equation (15) of the manuscript.
#'
#' @param imputation_model imputation model formula (or coercible to formula) passed through to \code{survreg}, a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See \code{survreg} documentation for more details.
#' @param dist imputation model distribution passed through to \code{survreg}. Options of \code{dist =} \code{"exponential"}, \code{"extreme"}, \code{"gaussian"}, \code{"logistic"}, \code{"loglogistic"}, \code{"loggaussian"}, \code{"lognormal"}, \code{"rayleigh"}, \code{"t"}, and \code{"weibull"} currently accepted.
#' @param data dataframe or named matrix containing columns \code{W}, \code{Delta}, and any other variables in \code{imputation_model}.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#' @param use_cumulative_hazard logical, if \code{TRUE} the survival function is transformed to the cumulative hazard before integration. Default is \code{TRUE}.
#'
#' @return A list containing:
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

cmi_fp_stabilized_wo_mean_single = function(imputation_model, dist, data, maxiter = 100, use_cumulative_hazard = FALSE) {
  # Fit AFT imputation model for X ~ Z
  fit = survreg(formula = imputation_model,
                data = data,
                dist = dist,
                maxiter = maxiter)

  # Extract variable names from imputation_model
  W = all.vars(imputation_model)[1] ## censored covariate
  Delta = all.vars(imputation_model)[2] ## corresponding event indicator

  # Initialize imputed values
  data$imp = data[, W] ## start with imp = W

  # Calculate linear predictor for AFT imputation model
  lp = fit$linear.predictors ## linear predictors from the survreg fit

  # Create an indicator variable for being uncensored
  uncens = data[, Delta] == 1

  # Use adaptive quadrature (and clever transformation) to estimate
  ## integral from X = Wi to X = Infinity
  int_surv_W_to_Inf = sapply(
    X = which(!uncens),
    FUN = function(i) {
      integrate(f = eq15_integrand,
                lower = 0,
                upper = 1,
                Wi = data[i, W],
                lpi = lp[i],
                dist = dist,
                fit = fit,
                use_cumulative_hazard = use_cumulative_hazard)$value
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
