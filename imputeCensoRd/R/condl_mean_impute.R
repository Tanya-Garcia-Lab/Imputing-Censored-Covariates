#' Conditional mean single imputation
#'
#' Imputes censored covariates with their conditional mean given censored value and additional covariates (where supplied).
#'
#' @param fit A \code{coxph} or \code{survfit} imputation model object.
#' @param obs String column name for the censored covariate.
#' @param event String column name for the censoring indicator of the covariate.
#' @param addl_covar (Optional) string or vector of strings for the additional fully-observed covariates. Default is \code{NULL}.
#' @param data Datafrane containing columns \code{obs}, \code{event}, and (if provided) \code{addl_covar}.
#' @param approx_beyond Choice of approximation used to extrapolate the survival function beyond the last observed covariate value. Default is \code{"expo"} for the exponential approximation. Other choices include \code{"zero"} or \code{"carryforward"}.
#' @param sample_lambda Logical; if true, lambda will be randomly drawn from asymptotic normal distribution of lambda based on Cox model fit
#'
#' @return A dataframe augmenting \code{data} with a column of imputed covariate values called \code{imp}.
#'
#' @export
#' @importFrom MASS mvrnorm
condl_mean_impute <- function(fit, obs, event, addl_covar = NULL, data, approx_beyond = "expo", sample_lambda = F) {
  # test for bad input
  if (!is.character(obs)) { stop("argument obs must be a character") }
  if (!is.character(event)) { stop("argument event must be a character") }
  if (!is.null(addl_covar) & !is.character(addl_covar)) { stop("when supplied, argument addl_covar must be a character") }
  if (!is.data.frame(data) & !is.matrix(data)) { stop("argument data must be a data frame or a matrix") }
  # test that data contains columns with specified names
  if (!(obs %in% colnames(data))) { stop(paste("data does not have column with name", obs)) }
  if (!(event %in% colnames(data))) { stop(paste("data does not have column with name", event)) }
  if (!is.null(addl_covar) & !all(addl_covar %in% colnames(data))) { stop(paste("data does not have columns with names:", paste(addl_covar, collapse = ", ")))}
  # test that approx_beyond is one of three acceptable options
  if (!(approx_beyond %in% c("expo", "zero", "carryforward"))) { stop("argument approx_beyond must be expo, zero, or carryforward") }
  # test for improper entries in columns of data
  #### Still deciding whether the following conditions should produce errors (stop()) or warnings
  if (any(data[, obs] < 0)) { warning(paste("elements of column", obs, "must be positive")) }
  if (!all(data[, event] == 0 | data[, event] == 1)) { warning(paste("elements of column", event, "must be either 0 or 1")) }
  
  if (sample_lambda) {
    # sample lambda from its asymptotic distribution
    lambda = matrix(mvrnorm(1, fit$coefficients, vcov(fit)), ncol = 1)
  } else { 
    # use point estimate of lambda
    lambda = matrix(fit$coefficients, ncol = 1)
  }

  if (is.null(addl_covar)) {
    # Estimate baseline survival from Kaplan-Meier estimator
    surv_df <- with(fit, data.frame(t = time, surv = surv))
  } else {
    # Calculate linear predictor \lambda %*% addl_covar for Cox model
    lp <- data.matrix(data[, addl_covar]) %*% lambda
    data$hr <- exp(lp)
    # Estimate baseline survival from Cox model fit using Breslow estimator
    cox_surv <- breslow_estimator(time = obs, event = event, hr = "hr", data = data)
    surv_df <- with(cox_surv, data.frame(t = times, surv = basesurv))
  }
  
  colnames(surv_df)[1] <- obs
  # Merge survival estimates into data
  data <- merge(x = data, y = surv_df, all.x = TRUE, sort = FALSE)
  # Order data by obs
  data <- data[order(data[, obs]), ]
  # Create an indicator variable for being uncensored
  uncens <- data[, event] == 1

  if (any(is.na(data[, "surv"]))) {
    # For censored subjects, survival is average of times right before/after
    suppressWarnings(
      data[is.na(data[, "surv"]), "surv"] <- sapply(X = data[is.na(data[, "surv"]), obs], FUN = impute_censored_surv, time = obs, event = event, surv = "surv", data = data)
    )
  }

  # Assume survival = 1 before first observed covariate value
  if (any(data[, obs] < min(data[uncens, obs]))) {
    cens_before <- which(data[, obs] < min(data[uncens, obs]))
    data[cens_before, "surv"] <- 1
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
      max_event <- max(which(data[, event] == 1))
      t_max_event <- data[max_event, obs]
      surv_max_event <- data[max_event, "surv"]
      data[cens_after, "surv"] <- exp(t_cens_after * log(surv_max_event) / t_max_event)
    }
  }
  # Distinct rows (in case of non-unique obs values)
  data_dist <- unique(data[, c(obs, event, addl_covar, "surv")])
  # [T_{(i+1)} - T_{(i)}]
  t_diff <- data_dist[-1, obs] - data_dist[-nrow(data_dist), obs]
  # Censored subject values (to impute)
  t_cens <- data[data[, event] == 0, obs]
  # For people with events, obs = X
  data$imp <- data[, obs]

  if (is.null(addl_covar)) {
    # Follow formula without additional covariates Z
    for (x in which(!uncens)) {
      Cj <- data[x, obs]
      Sj <- data_dist[-1, "surv"] + data_dist[-nrow(data_dist), "surv"]
      num <- sum((data_dist[-nrow(data_dist), obs] >= Cj) * Sj * t_diff)
      denom <- data[x, "surv"]
      data$imp[x] <- (1 / 2) * (num / denom) + Cj
    }
  } else {
    # Follow formula assuming Cox model with additional covariates Z
    for (x in which(!uncens)) {
      Zj <- data[x, addl_covar]
      lp <- as.numeric(data.matrix(Zj) %*% lambda)
      Cj <- data[x, obs]
      Sj <- data_dist[-1, "surv"] ^ (exp(lp)) + data_dist[-nrow(data_dist), "surv"] ^ (exp(lp))
      num <- sum((data_dist[-nrow(data_dist), obs] >= Cj) * Sj * t_diff)
      denom <- data[x, "surv"] ^ (exp(lp))
      data$imp[x] <- (1 / 2) * (num / denom) + Cj
    }
  }
  
  # Return input dataset with appended column imp containing imputed values 
  return(data)
}
