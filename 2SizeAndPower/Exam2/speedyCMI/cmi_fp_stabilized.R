#' Fully parametric imputation for a right-censored covariate with conditional means using the stabilized integral
#'
#' Fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then computes conditional means using the stabilized integral.
#'
#' @param imputation_model imputation model formula (or coercible to formula) passed through to \code{survreg}, a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See \code{survreg} documentation for more details.
#' @param dist imputation model distribution passed through to \code{survreg}. Options of \code{dist =} \code{"exponential"}, \code{"loglogistic"}, \code{"lognormal"}, and \code{"weibull"} currently accepted.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and any other variables in \code{imputation_model}.
#' @param with_mean logical, if \code{TRUE} the stabilized integral with the mean is used. Default is \code{FALSE}.
#' @param use_cumulative_hazard (optional) logical, if \code{use_cumulative_hazard = TRUE} and \code{with_mean = FALSE} the survival function is transformed to the cumulative hazard before integration. Default is \code{TRUE}.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#' @param boots (optional) numeric, for multiple imputation supply the desired number of imputations (obtained via bootstrapping) to \code{boots}. Default is \code{0}, which is single imputation.
#' @param seed (optional) numeric, for multiple imputation set the random seed for the bootstrapping with \code{seed}. Default is \code{NULL}, which does not reset \code{seed}.
#'
#' @return A list (or, in the case of multiple imputation, a list of lists) containing:
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

cmi_fp_stabilized = function(imputation_model, dist, data, with_mean = FALSE, use_cumulative_hazard = TRUE, maxiter = 100, boots = 0, seed = NULL) {
  if (with_mean) {
    if (boots == 0) {
      return_list = cmi_fp_stabilized_w_mean_single(imputation_model = imputation_model,
                                                    dist = dist,
                                                    data = data,
                                                    maxiter = maxiter)
    } else {
      if (!is.null(seed)) {
        set.seed(seed)
      }
      return_list = list()
      for (b in 1:boots) {
        re_data = data[sample(x = 1:nrow(data), size = nrow(data), replace = TRUE), ]
        return_list[[b]] = cmi_fp_stabilized_w_mean_single(imputation_model = imputation_model,
                                                           dist = dist,
                                                           data = re_data,
                                                           maxiter = maxiter)
      }
    }
  } else if (!with_mean) {
    if (boots == 0) {
      return_list = cmi_fp_stabilized_wo_mean_single(imputation_model = imputation_model,
                                                     dist = dist,
                                                     data = data,
                                                     maxiter = maxiter,
                                                     use_cumulative_hazard = use_cumulative_hazard)
    } else {
      if (!is.null(seed)) {
        set.seed(seed)
      }
      return_list = list()
      for (b in 1:boots) {
        re_data = data[sample(x = 1:nrow(data), size = nrow(data), replace = TRUE), ]
        return_list[[b]] = cmi_fp_stabilized_wo_mean_single(imputation_model = imputation_model,
                                                            dist = dist,
                                                            data = re_data,
                                                            maxiter = maxiter,
                                                            use_cumulative_hazard = use_cumulative_hazard)
      }
    }
  }
  return(return_list)
}
