#' Fully parametric imputation for a right-censored covariate with conditional means using the original integral
#'
#' Fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then integrates over the estimated survival function from \code{W} to \code{Inf} to compute conditional means.
#'
#' @param imputation_model imputation model formula (or coercible to formula) passed through to \code{survreg}, a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See \code{survreg} documentation for more details.
#' @param dist imputation model distribution passed through to \code{survreg}. See \code{survreg} documentation for more details.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and any other variables in \code{imputation_model}.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#' @param boots (optional) numeric, for multiple imputation supply the desired number of imputations (obtained via bootstrapping) to \code{boots}. Default is \code{0}, which is single imputation.
#' @param seed (optional) numeric, for multiple imputation set the random seed for the bootstrapping with \code{seed}. Default is \code{NULL}, which does not reset \code{seed}.
#'
#' @return A list (or, in the case of multiple imputation, a list of lists) containing:
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#' \item{aic}{Akaike information criterion (AIC) from the \code{imputation_model} fit.}
#' \item{bic}{Bayesian information criterion (BIC) from the \code{imputation_model} fit.}
#' \item{coefficients}{Vector of coefficients from the \code{imputation_model} fit.}
#' \item{scale}{Scale from the \code{imputation_model} fit.}
#'
#' @export
#' @importFrom survival survreg
#' @importFrom survival Surv
#' @importFrom survival psurvreg

cmi_fp_original = function(imputation_model, dist, data, maxiter = 100, boots = 0, seed = NULL) {
  if (boots == 0) {
    return_list = cmi_fp_original_single(imputation_model = imputation_model,
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
      return_list[[b]] = cmi_fp_original_single(imputation_model = imputation_model,
                                                dist = dist,
                                                data = re_data,
                                                maxiter = maxiter)
    }
  }
  return(return_list)
}
