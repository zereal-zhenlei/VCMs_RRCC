#' Single, fully parametric conditional mean imputation for a right-censored covariate (piecewise exponential distribution) with conditional means following Equation (11)
#'
#' Single, fully parametric conditional mean imputation for a right-censored covariate using a piecewise exponential model to estimate the conditional survival function and then uses an analytic solution to compute conditional means, as in Equation (11) of the manuscript.
#'
#' @param imputation_model imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and any other variables in \code{imputation_model}.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#' @param nintervals integer, number of disjoint subintervals used to split the hazard function of \code{W}. Must specify either this or \code{breaks}.
#' @param breaks vector, fixed subinterval boundaries used to split the hazard function of \code{W}. Must specify either this or \code{nintervals}.
#' @return
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#' \item{aic}{Akaike information criterion (AIC) from the imputation model fit.}
#' \item{bic}{Bayesian information criterion (AIC) from the \code{imputation_model} fit.}
#' \item{coefficients}{Coefficients for covariates in the \code{imputation_model} fit.}
#' \item{blhaz}{Baseline hazards from the \code{imputation_model} fit.}
#'
#' @importFrom survival Surv
#' @importFrom eha pchreg

cmi_fp_pwe_single = function(imputation_model, data, maxiter = 100, nintervals = NULL, breaks = NULL) {
  ## Checks
  if (is.null(nintervals)) {
    if (is.null(breaks)) {
      stop('One of nintervals / breaks must be non-null')
    }
  }

  ## Perform checks and start setup
  setup = cmi_pwe_setup(imputation_model, data)

  # Initialize imputed values
  data$imp = data[, setup$Wname] ## start with imp = W

  ## Compute breaks based on nintervals quantiles (if specified)
  if (!(is.null(nintervals))) {
    if (is.null(breaks)) {
      wuncen = data[data[, setup$Deltaname] == 1, setup$Wname]
      breaks = quantile(x = unlist(wuncen),
                        probs = seq(1, nintervals-1) / nintervals)
      breaks = round(x = c(0, breaks, Inf),
                     digits = 3)
    } else {
      warning('Both nintervals and breaks were specified. Ignoring nintervals and using breaks...')
    }
  }

  # Fit imputation model for X ~ Z
  fit = pchreg(formula = imputation_model,
               cuts = breaks,
               data = data)

  ## Vectorized test
  wcen = data[data[, setup$Deltaname] == 0, setup$Wname]
  zcen = data[data[, setup$Deltaname] == 0, setup$Zname]
  impX = pwe_mean_imputation(centime = wcen,
                             breaks = breaks,
                             blhaz = fit$hazards,
                             X = zcen,
                             beta = fit$coefficients)
  data[data[, setup$Deltaname] == 0, "imp"] = impX

  # Return input dataset with appended column imp containing imputed values
  ll = fit$loglik[2] ## maximum of log-likelihood
  k = length(fit$hazards) + length(fit$coefficients) ## number of parameters
  return_list = list(imputed_data = data,
                     code = !any(is.na(data$imp)),
                     aic = (2 * k - 2 * ll),
                     bic = (k * log(nrow(data)) - 2 * ll),
                     coefficients = fit$coefficients,
                     blhaz = fit$hazards)
  return(return_list)
}
