#' Conditional mean single imputation
#'
#' Imputes censored covariates with their conditional mean given censored value and additional covariates (where supplied).
#'
#' @param fit A \code{coxph} or \code{survfit} imputation model object.
#' @param obs String column name for the censored covariate.
#' @param delta String column name for the censoring indicator of the covariate.
#' @param addl_covar (Optional) string or vector of strings for the additional fully-observed covariates. Default is \code{NULL}.
#' @param data Datafrane containing columns \code{obs}, \code{delta}, and (if provided) \code{addl_covar}.
#' @param approx_beyond Choice of approximation used to extrapolate the survival function beyond the last observed covariate value. Default is \code{"expo"} for the exponential approximation. Other choices include \code{"zero"} or \code{"carryforward"}.
#'
#' @return A dataframe augmenting \code{data} with a column of imputed covariate values called \code{imp}.
#'
#' @export
condl_mean_impute <- function(fit, obs, delta, addl_covar = NULL, data, approx_beyond = "expo") {
  if (is.null(addl_covar)) {
    # Estimate baseline survival from Kaplan-Meier estimator
    surv_df <- with(fit, data.frame(t = time, surv = surv))
  } else {
    # Estimate baseline survival from Cox model fit
    data$hr <- exp(fit$coefficients * data[, addl_covar])
    cox_surv <- breslow_estimator(time = obs, delta = delta, hr = "hr", data = data)
    surv_df <- with(cox_surv, data.frame(t = times, surv = basesurv))
  }
  colnames(surv_df)[1] <- obs
  # Merge survival estimates into data
  data <- merge(x = data, y = surv_df, all.x = TRUE, sort = FALSE)
  # Order data by obs
  data <- data[order(data[, obs]), ]
  # Create an indicator variable for being uncensored
  uncens <- data[, delta] == 1

  if (any(is.na(data[, "surv"]))) {
    # For censored subjects, survival is average of times right before/after
    suppressWarnings(
      data[is.na(data[, "surv"]), "surv"] <- sapply(X = data[is.na(data[, "surv"]), obs], FUN = impute_censored_surv, time = obs, delta = delta, surv = "surv", data = data)
    )
  }

  # Extrapolate survival beyond last observed covariate
  if (any(data[, obs] > max(data[uncens, obs]))) {
    cens_after <- which(data[, obs] > max(data[uncens, obs]))
    t_cens_after <- data[cens_after, obs]
    last_event_surv <- data[max(which(uncens)), "surv"]
    # Efron (1967) immediately goes to zero
    if (approx_beyond == "zero") { data[cens_after, "surv"] <- 0 }
    # Gill (1980) carry forward survival at last event
    if (approx_beyond == "carryforward") { data[cens_after, "surv"] <- last_event_surv }
    # Brown, Hollander, and Kowar (1974) exponential approx
    if (approx_beyond == "expo") {
      max_event <- max(which(data[, delta] == 1))
      t_max_event <- data[max_event, obs]
      surv_max_event <- data[max_event, "surv"]
      data[cens_after, "surv"] <- exp(t_cens_after * log(surv_max_event) / t_max_event)
    }
  }
  # Distinct rows (in case of non-unique obs values)
  data_dist <- unique(data[, c(obs, delta, addl_covar, "surv")])
  # [T_{(i+1)} - T_{(i)}]
  t_diff <- data_dist[-1, obs] - data_dist[-nrow(data_dist), obs]
  # Censored subject values (to impute)
  t_cens <- data[data[, delta] == 0, obs]
  # For people with events, obs = X
  data$imp <- data[, obs]

  if (is.null(addl_covar)) {
    for (x in which(!uncens)) {
      Zj <- data[x, addl_covar]
      Cj <- data[x, obs]
      Sj <- data_dist[-1, "surv"] + data_dist[-nrow(data_dist), "surv"]
      num <- sum((data_dist[-nrow(data_dist), obs] > Cj) * Sj * t_diff)
      denom <- data[x, "surv"]
      data$imp[x] <- (1 / 2) * (num / denom) + Cj
    }
  } else {
    for (x in which(!uncens)) {
      Zj <- data[x, addl_covar]
      Cj <- data[x, obs]
      Sj <- data_dist[-1, "surv"] ^ (exp(fit$coefficients * Zj)) + data_dist[-nrow(data_dist), "surv"] ^ (exp(fit$coefficients * Zj))
      num <- sum((data_dist[-nrow(data_dist), obs] > Cj) * Sj * t_diff)
      denom <- data[x, "surv"] ^ (exp(fit$coefficients * Zj))
      data$imp[x] <- (1 / 2) * (num / denom) + Cj
    }
  }
  return(data)
}