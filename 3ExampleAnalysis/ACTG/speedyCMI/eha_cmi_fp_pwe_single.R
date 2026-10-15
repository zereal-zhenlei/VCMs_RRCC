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
#' \item{logblhaz}{Log baseline hazards from the \code{imputation_model} fit.}
#'
#' @importFrom eha pchreg

eha_cmi_fp_pwe_single = function(imputation_model, data, maxiter = 100, nintervals = NULL, breaks = NULL) {
  ## Checks
  if (is.null(nintervals)) {
    if (is.null(breaks)) {
      stop('One of nintervals / breaks must be non-null')
    }
  } else if (nintervals == 1) {
    stop('For one interval please use the exponential distribution option in cmi_fp_analytical()')
  }

  # Extract variable names from imputation_model
  Wname = all.vars(imputation_model)[1] ## censored covariate
  Deltaname = all.vars(imputation_model)[2] ## corresponding event indicator

  # Initialize imputed values
  data$imp = data[, Wname] ## start with imp = W

  # Compute breaks based on nintervals quantiles (if specified)
  if (!(is.null(nintervals))) {
    if (is.null(breaks)) {
      Wobs <- data[data[, Deltaname] == 1, Wname]
      breaks <- quantile(x = unlist(Wobs),
                         probs = (1:(nintervals - 1)) / nintervals)
      breaks <- c(0, breaks, Inf)
      breaks <- round(breaks, 3)
    } else {
      warning('Both nintervals and breaks were specified. Ignoring nintervals and using breaks...')
    }
  }

  # Fit AFT imputation model for X ~ Z
  fit = pchreg(formula = imputation_model,
               cuts = breaks,
               data = data)

  ## Perform checks and start setup
  W  = data[[Wname]]
  Delta = data[[Deltaname]]

  ## -----------------------
  ## Impute censored times
  ## -----------------------
  ## Obtain censored values
  cenindx = data[[Deltaname]] == 0
  lower = data[[Wname]][cenindx]
  idcen = (1:nrow(data))[cenindx]

  ## Filter pdata so that only censored individuals are included
  data2 = data
  data2$id = 1:nrow(data2)
  data2 = data2[which(data2$id %in% idcen), ]
  C = data2[[Wname]] ## Censored times

  ## Extract design matrix (ignoring intercept terms / basline hazards)
  X = model.matrix(imputation_model, data2)
  if ( '(Intercept)' %in% colnames(X) ) {
    intindx = which(colnames(X) == '(Intercept)')
    X = X[, -intindx]
    if (is.null(ncol(X))) {
      X = data.matrix(frame = X)
     }
  }
  ## Extract coefficients (ignoring intercept terms)
  beta = coef(fit)
  #blhazindx = which(names(beta) %in% paste0('interval', intvlnames))
  logblhaz = fit$hazards #beta[blhazindx]
  #beta = beta[-blhazindx]
  if (length(beta) == 0) {
    lp = rep(0, nrow(X))
  } else {
    lp = (X %*% beta)[, 1]
  }
  J = length(breaks) - 1

  ## Compute survival function at each endpoint
  ## and at censored time
  cumHaz = matrix(0, nrow = nrow(X), ncol = J+1)   ## cumulative hazard at each cut point
  cumHazCen = matrix(0, nrow = nrow(X), ncol = J+1)   ## cumulative hazard at each point until censoring
  survProbs = matrix(1, nrow = nrow(X), ncol = J+1)   ## survival probabilities at each cut point
  survProbsCen = matrix(1, nrow = nrow(X), ncol = J+1)   ## survival probabilities at each cut point
  survInt = numeric(nrow(X))  ## \int_{0}^{\infty} S(t) dt
  survIntCen = numeric(nrow(X))  ## \int_{0}^{C} S(t) dt
  logSurvCen = numeric(nrow(X))  ## log S(C)
  for (j in 1:J) {
    ## Compute cumulative hazard and survival at interval j
    cumHaz[, j+1]    = cumHaz[, j] + exp( logblhaz[j] + lp + log(breaks[j+1] - breaks[j]) )
    survProbs[, j+1] = exp(-cumHaz[, j+1])
    survInt          = survInt + exp(log(survProbs[, j] - survProbs[, j+1]) - logblhaz[j])
    ## See if subject was at risk in interval j;
    ##  compute time at risk and cumulative hazard in that interval
    indx_j                    = (C > breaks[j])
    risktime_j                = pmin(C[indx_j], breaks[j+1]) - breaks[j]
    cumHazCen[indx_j, j+1]    = cumHazCen[indx_j, j] + exp( logblhaz[j] + lp[indx_j] + log(risktime_j) )
    survProbsCen[indx_j, j+1] = exp(-cumHazCen[indx_j, j+1])
    survIntCen[indx_j]        = survIntCen[indx_j] +
      exp(log(survProbsCen[indx_j, j] - survProbsCen[indx_j, j+1]) - logblhaz[j])
    logSurvCen[indx_j] = logSurvCen[indx_j] - exp( logblhaz[j] + lp[indx_j] + log(risktime_j) )
  }
  ## Impute censoring time
  imp = C + exp( log( survInt - survIntCen ) - logSurvCen )

  ## Return data set with imputed time
  data$imp  = data[[Wname]]
  data$imp[data[[Deltaname]] == 0] = imp

  # Return input dataset with appended column imp containing imputed values
  k = length(beta) + length(logblhaz)
  return_list = list(imputed_data = data,
                     code = !any(is.na(data$imp)),
                     aic = (2 * k - 2 * fit$loglik)[2],
                     bic = (k * log(nrow(X)) - 2 * fit$loglik)[2],
                     coefficients = fit$coefficients,
                     logblhaz = fit$hazards)
  return(return_list)
}
