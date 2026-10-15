eq15_integrand = function(t, Wi, lpi, dist, fit, use_cumulative_hazard) {
  if (use_cumulative_hazard) {
    ## H{Wi + (1 - t) / t|Zi}
    num = -log(1 - psurvreg(q = Wi + (1 - t) / t,
                            mean = lpi,
                            scale = fit$scale,
                            dist = dist))
    return(exp(- (num + 2 * log(t))))
  } else {
    ## S{Wi + (1 - t) / t|Zi}
    num = 1 - psurvreg(q = Wi + (1 - t) / t,
                       mean = lpi,
                       scale = fit$scale,
                       dist = dist)
    denom = t ^ 2
    return(num / denom)
  }
}

cmi_pwe_setup = function(imputation_formula, data) {
  Terms = terms(imputation_formula, data = data)
  if ( attr(Terms, 'response') == 0 )
    stop('formula must have a Surv response')
  mf = model.frame(imputation_formula, data = data)
  Y  = model.extract(mf, 'response')
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")
  ## Extract names of variables
  Wname = all.vars(imputation_formula)[1]
  Deltaname = all.vars(imputation_formula)[2]   ## assumes right-censoring
  Zname = all.vars(imputation_formula)[-c(1:2)]
  ## Return list
  list('Wname' = Wname,
       'Deltaname' = Deltaname,
       'Zname' = Zname)
}

loglik_pweph <- function(y, event, X, breaks, logblhaz, beta) {

  ## Compute hazard for each interval
  lp     <- (as.matrix(X) %*% beta)[, 1]

  ## Get interval start and end points
  J     <- length(breaks) - 1
  start <- breaks[1:J]
  end   <- breaks[2:(J+1)]

  ## Initialize log likelihood contribution
  loglik <- numeric(length(y))
  for ( j in 1:J ) {
    loghaz_j <- logblhaz[j] + lp
    atrisk_j <- ( y >= start[j] )
    action_j <- atrisk_j & ( y < end[j] ) & ( event == 1 )

    ## Compute survival likelihood contribution
    logrisktime_j    <- log( pmin(y[atrisk_j], end[j]) - start[j] )
    loglik[atrisk_j] <- loglik[atrisk_j] - exp(logrisktime_j + loghaz_j[atrisk_j] )

    ## Compute hazard likelihood contribution
    loglik[action_j] <- loglik[action_j] + loghaz_j[action_j]
  }

  return(sum(loglik))
}

survInt_pwe <- function(centime, breaks, blhaz, X, beta) {
  X     <- as.matrix(X)
  lp    <- (X %*% beta)[, 1]
  J     <- length(breaks) - 1
  start <- breaks[1:J]
  end   <- breaks[2:(J+1)]
  diff  <- end - start

  fullInt  <- numeric(nrow(X))   ## \int_{0}^{\infty} S(t) dt
  lowerInt <- numeric(nrow(X))   ## \int_{0}^{c} S(t) dt
  cumHaz   <- numeric(nrow(X))

  for ( j in 1:J ) {

    ## Compute hazard for interval j
    haz_j   <- exp( log(blhaz[j]) + lp )

    ## Update full integral
    timeatrisk_j <- end[j] - start[j]
    fullInt      <- fullInt + ( 1 - exp( -haz_j * timeatrisk_j ) ) / haz_j * exp( -cumHaz )

    ## Compute lower integral
    indx_atrisk   <- centime > start[j]
    if ( sum(indx_atrisk) > 0 ) {
      timeatrisk_j  <- pmin( centime[indx_atrisk], end[j] ) - start[j]
      lowerInt[indx_atrisk] <- (
        lowerInt[indx_atrisk]
        + ( 1 - exp( -haz_j[indx_atrisk] * timeatrisk_j ) ) / haz_j[indx_atrisk] * exp( -cumHaz[indx_atrisk] )
      )
    }

    ## Compute cumulative hazard for interval j
    cumHaz <- cumHaz + haz_j * (end[j] - start[j])
  }
  return(fullInt - lowerInt)
}

Spweph <- function(x, breaks, hazards) {
  J     <- length(breaks) - 1
  start <- breaks[1:J]
  end   <- breaks[2:(J+1)]
  diff  <- end - start

  ## Cumulative hazard computation
  n      <- length(x)
  cumHaz <- numeric(n)
  for ( j in 1:J ) {
    indx_surv           <- x >= end[j]
    cumHaz[indx_surv]   <- cumHaz[indx_surv] + hazards[j] * diff[j]

    indx_action         <- (x >= start[j]) & (x < end[j])
    cumHaz[indx_action] <- cumHaz[indx_action] + hazards[j] * (x[indx_action] - start[j])
  }
  exp(-cumHaz)
}

pwe_mean_imputation <- function(centime, breaks, blhaz, X, beta) {
  X     <- as.matrix(X)
  lp    <- (X %*% beta)[, 1]
  J     <- length(breaks) - 1
  start <- breaks[1:J]
  end   <- breaks[2:(J+1)]
  diff  <- end - start

  fullInt  <- numeric(nrow(X))   ## \int_{0}^{\infty} S(t) dt
  lowerInt <- numeric(nrow(X))   ## \int_{0}^{c} S(t) dt
  cumHaz   <- numeric(nrow(X))
  negLogSurvtime <- numeric(nrow(X))  ## -log S(c)

  for ( j in 1:J ) {

    ## Compute hazard for interval j
    haz_j   <- exp( log(blhaz[j]) + lp )

    ## Update full integral
    timeatrisk_j <- end[j] - start[j]
    fullInt      <- fullInt + ( 1 - exp( -haz_j * timeatrisk_j ) ) / haz_j * exp( -cumHaz )

    ## Compute lower integral and update hazard
    indx_atrisk   <- centime > start[j]
    if ( sum(indx_atrisk) > 0 ) {
      timeatrisk_j                <- pmin( centime[indx_atrisk], end[j] ) - start[j]
      cumhazcontrib_j             <- haz_j[indx_atrisk] * timeatrisk_j
      lowerInt[indx_atrisk] <- (
        lowerInt[indx_atrisk]
        + ( 1 - exp( -cumhazcontrib_j ) ) / haz_j[indx_atrisk] * exp( -cumHaz[indx_atrisk] )
      )

      negLogSurvtime[indx_atrisk] <- negLogSurvtime[indx_atrisk] + cumhazcontrib_j

    }

    ## Compute cumulative hazard for interval j
    cumHaz <- cumHaz + haz_j * (end[j] - start[j])
  }
  Int <- fullInt - lowerInt

  ## Return L + [\int_{L}^{\infty} S(t) dt] / S(L)
  survTime <- exp(-negLogSurvtime)
  res      <- centime + Int / survTime
  return(res)
}
