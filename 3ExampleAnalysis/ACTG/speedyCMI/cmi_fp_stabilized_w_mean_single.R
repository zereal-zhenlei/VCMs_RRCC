#' Single, fully parametric imputation for a right-censored covariate with conditional means using the stabilized integral with mean
#'
#' Single, fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then integrates over the estimated survival function from \code{0} to \code{W} to compute conditional means using the stabilized integral.
#'
#' @param imputation_model imputation model formula (or coercible to formula) passed through to \code{survreg}, a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See \code{survreg} documentation for more details.
#' @param dist imputation model distribution passed through to \code{survreg}. Options of \code{dist =} \code{"exponential"}, \code{"loglogistic"}, \code{"lognormal"}, and \code{"weibull"} currently accepted.
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

cmi_fp_stabilized_w_mean_single = function(imputation_model, dist, data, maxiter = 100) {
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

  if (dist %in% c("weibull", "exponential")) {
    # Transform parameters to agree with R's weibull parameterization
    alpha = 1 / fit$scale
    lambda = exp(lp)

    # Calculate mean life = integral from 0 to \infty of S(t|Z)
    est_ml = lambda[which(!uncens)] * gamma(1 + fit$scale) ## using formula for Weibull distribution

    # Use adaptive quadrature to estimate
    ## integral from X = 0 to X = Wi
    int_surv_0_to_W = sapply(
      X = which(!uncens),
      FUN = function(i) {
        integrate(f = pweibull,
                  lower = 0,
                  upper = data[i, W],
                  shape = alpha,
                  scale = lambda[i],
                  lower.tail = FALSE)$value
      }
    )

    # Calculate S(W|Z)
    surv = pweibull(q = data[which(!uncens), W],
                    shape = alpha,
                    scale = lambda[which(!uncens)],
                    lower.tail = FALSE)
  } else if (dist == "lognormal") {
    mu = lp
    sigma = fit$scale

    # Calculate mean life = integral from 0 to \infty of S(t|Z)
    est_ml = exp(mu[which(!uncens)] + sigma ^ 2 / 2) ## using formula for log-normal distribution

    # Use adaptive quadrature to estimate
    ## integral from X = 0 to X = Wi
    int_surv_0_to_W = sapply(
      X = which(!uncens),
      FUN = function(i) {
        integrate(f = plnorm,
                  lower = 0,
                  upper = data[i, W],
                  meanlog = mu[i],
                  sdlog = sigma,
                  lower.tail = FALSE)$value
      }
    )

    # Calculate S(W|Z)
    surv = plnorm(q = data[which(!uncens), W],
                  meanlog = lp[which(!uncens)],
                  sdlog = sigma,
                  lower.tail = FALSE)
  } else if (dist == "loglogistic") {
    # Transform parameters to agree with R's weibull parameterization
    alpha = 1 / fit$scale
    lambda = exp(lp)

    if (alpha > 1) {
      # Calculate mean life = integral from 0 to \infty of S(t|Z)
      est_ml = (lambda[which(!uncens)] * pi / alpha) / sin(pi / alpha) ## using formula for log-logistic distribution

      # Use adaptive quadrature to estimate
      ## integral from X = 0 to X = Wi
      int_surv_0_to_W = sapply(
        X = which(!uncens),
        FUN = function(i) {
          integrate(f = function(x) 1 / (1 + (x / lambda[i]) ^ alpha),
                    lower = 0,
                    upper = data[i, W])$value
        }
      )

      # Calculate S(W|Z)
      surv = 1 / (1 + (data[which(!uncens), W] / lambda[which(!uncens)]) ^ alpha)
    } else {
      # Can't calculate mean life; expected value is undefined if shape alpha <= 1
      ## Make all imputed values NA
      data[which(!uncens), "imp"] = NA

      # Return input dataset with appended column imp containing imputed values
      return_list = list(imputed_data = data,
                         code = !any(is.na(data$imp)),
                         aic = AIC(object = fit),
                         bic = BIC(object = fit),
                         coefficients = fit$coefficients,
                         scale = fit$scale)
      return(return_list)
    }
  }

  # Subtract integral from mean life to get
  ## integral from X = Wi to X = Infinity
  int_surv_W_to_Inf = est_ml - int_surv_0_to_W

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
