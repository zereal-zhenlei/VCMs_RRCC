#' Multiple, parametric conditional mean imputation for a censored covariate
#'
#' Multiple, parametric conditional mean imputation for a censored covariate using parametric survival regression to estimate the conditional survival function.
#'
#' @param imputation_model Imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param dist Imputation model distribution passed through to \code{survreg}. See \code{survreg} documentation for more details.
#' @param analysis_model Analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be a continuous outcome for normal linear regression.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param nintervals (only if \code{dist = "PWE"}) integer, number of disjoint subintervals used to split the hazard function of \code{W}. Must specify either this or \code{breaks}.
#' @param breaks (only if \code{dist = "PWE"}) vector, fixed subinterval boundaries used to split the hazard function of \code{W}. Must specify either this or \code{nintervals}.
#' @param type String name of the type of imputation procedure to be used. Options include \code{"original"} (the default), \code{"analytical"}, \code{"stabilized (with mean)"}, and \code{"stabilized (without mean)"}.
#' @param B numeric, number of imputations. Default is \code{10}.
#'
#' @return A list containing:
#' \item{coeff}{dataframe containing the pooled coefficients and standard errors for \code{analysis_model}}
#' \item{aic}{Akaike information criterion (AIC) from the \code{imputation_model} fit to the original data.}
#' \item{bic}{Bayesian information criterion (AIC) from the \code{imputation_model} fit to the original data.}
#'
#'
#' @export

cmi_fp_bootstrap = function(imputation_model, dist, analysis_model, data, nintervals = NULL, breaks = NULL, type = "original", B = 10) {
  # Size of resample -----------------------------------------------------------
  n = nrow(data)

  # Make data a data frame (in case it's a tibble) -----------------------------
  data = data.frame(data)

  # Extract variable names from imputation_model -------------------------------
  W = all.vars(imputation_model)[1] ## censored covariate
  Delta = all.vars(imputation_model)[2] ## corresponding event indicator
  Z = all.vars(imputation_model)[-c(1:2)] ## additional covariates

  # Take names and dimension from naive fit
  data[, "imp"] = data[, W] ## start imputed value with observed
  naive_fit = lm(formula = analysis_model, data = data) ## fit naive model
  p = length(naive_fit$coefficients) ## dimension of beta vector

  # Initialize empty dataframe to hold results
  mult_fit = data.frame()

  # Loop through replicates
  for (b in 1:B) {
    # Sample with replacement from the original data ---------------------------
    re_rows = ceiling(runif(n = n, min = 0, max = 1) * nrow(data))
    re_data = data[re_rows, ]

    # Calculate proportion of censored observations in resampled data ----------
    re_prop_cens = 1 - mean(re_data[, Delta])

    # Check for censoring in resampled data
    if (0 < re_prop_cens) {
      # Impute censored x in re_data -------------------------------------------
      if (type == "analytical") {
        re_data_imp = cmi_fp_analytical(imputation_model = imputation_model,
                                        dist = dist,
                                        data = re_data,
                                        nintervals = nintervals,
                                        breaks = breaks,
                                        boots = 0)
      } else if (type == "stabilized (with mean)") {
        re_data_imp = cmi_fp_stabilized(imputation_model = imputation_model,
                                        dist = dist,
                                        data = re_data,
                                        with_mean = TRUE,
                                        boots = 0)
      } else if (type == "stabilized (without mean)") {
        re_data_imp = cmi_fp_stabilized(imputation_model = imputation_model,
                                        dist = dist,
                                        data = re_data,
                                        with_mean = FALSE,
                                        use_cumulative_hazard = TRUE,
                                        boots = 0)
      } else {
        re_data_imp = cmi_fp_original(imputation_model = imputation_model,
                                      dist = dist,
                                      data = re_data,
                                      boots = 0)
      }
    } else {
      # If no censored, just make imp = W and fit the usual model --------------
      re_data_imp = list(imputed_data = re_data,
                         code = TRUE)
      re_data_imp$imputed_data$imp = re_data[, W]
    }

    # Check for imputation model not converging --------------------------------
    while (!re_data_imp$code) {
      # Sample with replacement from the original data -------------------------
      re_rows = ceiling(runif(n = n, min = 0, max = 1) * nrow(data))
      re_data = data[re_rows, ]

      # Calculate proportion of censored observations in resampled data --------
      re_prop_cens = 1 - mean(re_data[, Delta])

      # Check for censoring in resampled data
      if (0 < re_prop_cens) {
        # Impute censored x in re_data -----------------------------------------
        if (type == "analytical") {
          re_data_imp = cmi_fp_analytical(imputation_model = imputation_model,
                                          dist = dist,
                                          data = re_data,
                                          nintervals = nintervals,
                                          breaks = breaks,
                                          boots = 0)
        } else if (type == "stabilized (with mean)") {
          re_data_imp = cmi_fp_stabilized(imputation_model = imputation_model,
                                          dist = dist,
                                          data = re_data,
                                          with_mean = TRUE,
                                          boots = 0)
        } else if (type == "stabilized (without mean)") {
          re_data_imp = cmi_fp_stabilized(imputation_model = imputation_model,
                                          dist = dist,
                                          data = re_data,
                                          with_mean = FALSE,
                                          use_cumulative_hazard = TRUE,
                                          boots = 0)
        } else {
          re_data_imp = cmi_fp_original(imputation_model = imputation_model,
                                        dist = dist,
                                        data = re_data,
                                        boots = 0)
        }
      } else {
        # If no censored, just make imp = W and fit the usual model ------------
        re_data_imp = list(imputed_data = re_data,
                           code = TRUE)
        re_data_imp$imputed_data$imp = re_data[, W]
      }
    }

    # Once imputation was successful, fit the analysis model -------------------
    re_fit = lm(formula = analysis_model,
                data = re_data_imp$imputed_data)

    # Save coefficients to results matrix --------------------------------------
    mult_fit = rbind(mult_fit,
                     summary(re_fit)$coefficients)
  }

  ## Pool estimates
  pooled_est = rowsum(x = mult_fit[, 1],
                      group = rep(x = 1:p, times = B)) / B

  ## Calculate within-imputation variance
  within_var = rowsum(x = mult_fit[, 2] ^ 2,
                      group = rep(x = 1:p, times = B)) / B

  ## Calculate between-imputation variance
  between_var = rowsum(x = (mult_fit[, 1] - rep(pooled_est, times = B)) ^ 2,
                       group = rep(x = 1:p, times = B)) / (B - 1)

  ## Pool variance
  pooled_var = within_var + between_var + (between_var / B)

  # Return table of pooled estimates
  tab = data.frame(Coefficient = names(naive_fit$coefficients),
                   Est = pooled_est,
                   SE = sqrt(pooled_var))

  # Fit imputation model to original data to get model diagnostics
  # Impute censored x in re_data -------------------------------------------
  if (type == "analytical") {
    data_imp = cmi_fp_analytical(imputation_model = imputation_model,
                                 dist = dist,
                                 data = data,
                                 boots = 0,
                                 nintervals = nintervals,
                                 breaks = breaks)
  } else if (type == "stabilized (with mean)") {
    data_imp = cmi_fp_stabilized(imputation_model = imputation_model,
                                 dist = dist,
                                 data = data,
                                 with_mean = TRUE,
                                 boots = 0)
  } else if (type == "stabilized (without mean)") {
    data_imp = cmi_fp_stabilized(imputation_model = imputation_model,
                                 dist = dist,
                                 data = data,
                                 with_mean = FALSE,
                                 use_cumulative_hazard = TRUE,
                                 boots = 0)
  } else {
    data_imp = cmi_fp_original(imputation_model = imputation_model,
                               dist = dist,
                               data = data,
                               boots = 0)
  }

  return(list(coeff = tab,
              aic = data_imp$aic,
              bic = data_imp$bic))
}
