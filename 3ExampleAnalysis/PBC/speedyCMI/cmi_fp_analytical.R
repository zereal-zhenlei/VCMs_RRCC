#' Fully parametric conditional mean imputation for a right-censored covariate with conditional means using the analytical solution
#'
#' Fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then uses an analytical solution to compute conditional means.
#'
#' @param imputation_model imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param dist imputation model distribution passed through to \code{survreg}. See \code{survreg} documentation for more details.
#' @param data dataframe or named matrix containing columns \code{W}, \code{Delta}, and any other variables in \code{imputation_model}.
#' @param nintervals (only if \code{dist = "PWE"}) integer, number of disjoint subintervals used to split the hazard function of \code{W}. Must specify either this or \code{breaks}.
#' @param breaks (only if \code{dist = "PWE"}) vector, fixed subinterval boundaries used to split the hazard function of \code{W}. Must specify either this or \code{nintervals}.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{max_iter = 100}.
#' @param boots (optional) numeric, for multiple imputation supply the desired number of imputations (obtained via bootstrapping) to \code{boots}. Default is \code{0}, which is single imputation.
#' @param seed (optional) numeric, for multiple imputation set the random seed for the bootstrapping with \code{seed}. Default is \code{NULL}, which does not reset \code{seed}.
#'
#' @return A list (or, in the case of multiple imputation, a list of \code{boots} lists) containing:
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#' \item{aic}{Akaike information criterion (AIC) from the \code{imputation_model} fit.}
#' \item{bic}{Bayesian information criterion (AIC) from the \code{imputation_model} fit.}
#' \item{coefficients}{Vector of coefficients from the \code{imputation_model} fit.}
#' \item{scale}{Scale from the \code{imputation_model} fit.}
#'
#' @export
#' @importFrom survival survreg
#' @importFrom survival Surv
#' @importFrom survival psurvreg

cmi_fp_analytical = function(imputation_model, dist, data, nintervals = NULL, breaks = NULL, maxiter = 100, boots = 0, seed = NULL) {
  # Extract variable names from imputation_model
  W = all.vars(imputation_model)[1] ## censored covariate
  Delta = all.vars(imputation_model)[2] ## corresponding event indicator

  # Convert data to dataframe (just in case)
  data = data.frame(data)

  # Single imputation
  if (boots == 0) {
    if (toupper(dist) == "WEIBULL") {
      ## If Weibull, use equation at the end of Section 2.4.1
      return_list = cmi_fp_weibull_single(imputation_model = imputation_model,
                                          W = W,
                                          Delta = Delta,
                                          data = data,
                                          maxiter = maxiter)
    } else if (toupper(dist) == "EXPONENTIAL") {
      ## If exponential, use equation from Section 2.4.2
      return_list = cmi_fp_expo_single(imputation_model = imputation_model,
                                       W = W,
                                       Delta = Delta,
                                       data = data,
                                       maxiter = maxiter)
    } else if (toupper(dist) %in% c("LOGNORMAL", "LOGGAUSSIAN")) {
      ## If log-normal, use equation at the end of Section 2.4.3
      return_list = cmi_fp_lognorm_single(imputation_model = imputation_model,
                                          W = W,
                                          Delta = Delta,
                                          data = data,
                                          maxiter = maxiter)
    } else if (toupper(dist) == "LOGLOGISTIC") {
      ## If log-logistic, use equation at the end of Section 2.4.4
      return_list = cmi_fp_loglogistic_single(imputation_model = imputation_model,
                                              W = W,
                                              Delta = Delta,
                                              data = data,
                                              maxiter = maxiter)
    } else if (toupper(dist) == "PWE") {
      ## If piecewise exponential, use equation at the end of Section 2.4.5
      return_list = cmi_fp_pwe_single(imputation_model = imputation_model,
                                      data = data,
                                      maxiter = maxiter,
                                      nintervals = nintervals,
                                      breaks = breaks)
    }
  } else { # Multiple imputation
    if (!is.null(seed)) {
      set.seed(seed)
    }
    return_list = list()
    for (b in 1:boots) {
      b_rows = sample(x = 1:nrow(data),
                      size = nrow(data),
                      replace = TRUE)
      b_data = data[b_rows, ]
      if (toupper(dist) == "WEIBULL") {
        ## If Weibull, use equation at the end of Section 2.4.1
        return_list[[b]] = cmi_fp_weibull_single(imputation_model = imputation_model,
                                                 W = W,
                                                 Delta = Delta,
                                                 data = b_data,
                                                 maxiter = maxiter)
      } else if (toupper(dist) == "EXPONENTIAL") {
        ## If exponential, use equation from Section 2.4.2
        return_list[[b]] = cmi_fp_expo_single(imputation_model = imputation_model,
                                              W = W,
                                              Delta = Delta,
                                              data = b_data,
                                              maxiter = maxiter)
      } else if (toupper(dist) %in% c("LOGNORMAL", "LOGGAUSSIAN")) {
        ## If log-normal, use equation at the end of Section 2.4.3
        return_list[[b]] = cmi_fp_lognorm_single(imputation_model = imputation_model,
                                                 W = W,
                                                 Delta = Delta,
                                                 data = b_data,
                                                 maxiter = maxiter)
      } else if (toupper(dist) == "LOGLOGISTIC") {
        ## If log-logistic, use equation at the end of Section 2.4.4
        return_list[[b]] = cmi_fp_loglogistic_single(imputation_model = imputation_model,
                                                     W = W,
                                                     Delta = Delta,
                                                     data = b_data,
                                                     maxiter = maxiter)
      } else if (toupper(dist) == "PWE") {
        ## If piecewise exponential, use equation at the end of Section 2.4.5
        return_list[[b]] = cmi_fp_pwe_single(imputation_model = imputation_model,
                                             data = b_data,
                                             maxiter = maxiter,
                                             nintervals = nintervals,
                                             breaks = breaks)
      }
    }
  }
  return(return_list)
}
