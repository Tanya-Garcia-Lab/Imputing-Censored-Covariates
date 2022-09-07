#' Conditional mean imputation with bootstrap
#'
#' Imputes censored covariates with their conditional mean given censored value and additional covariates (where supplied).
#'
#' @param obs String column name for the censored covariate.
#' @param event String column name for the censoring indicator of the covariate.
#' @param addl_covar (Optional) string or vector of strings for the additional fully-observed covariates. Default is \code{NULL}.
#' @param data Datafrane containing columns \code{obs}, \code{event}, and (if provided) \code{addl_covar}.
#' @param approx_beyond Choice of approximation used to extrapolate the survival function beyond the last observed covariate value. Default is \code{"expo"} for the exponential approximation. Other choices include \code{"zero"} or \code{"carryforward"}.
#' @param M number of bootstrap samples to be taken from \code{data}
#'
#' @return A list of \code{M} dataframes, each of which is sampled with replacement from \code{data} and then augmented with a column of imputed covariate values called \code{imp}.
#'
#' @export
condl_mean_impute_sample_lambda <- function(obs, event, addl_covar = NULL, data, approx_beyond = "expo", M) {
  # test for bad input
  if (!is.character(obs)) { stop("argument obs must be a character") }
  if (!is.character(event)) { stop("argument event must be a character") }
  if (!is.null(addl_covar) & !is.character(addl_covar)) { stop("when supplied, argument addl_covar must be a character") }
  if (!is.data.frame(data) & !is.matrix(data)) { stop("argument data must be a data frame or a matrix") }
  # if (!(M >= 1 & is.integer(M))) { stop("argument M must be a positive integer")}
  # test that data contains columns with specified names
  if (!(obs %in% colnames(data))) { stop(paste("data does not have column with name", obs)) }
  if (!(event %in% colnames(data))) { stop(paste("data does not have column with name", event)) }
  if (!is.null(addl_covar) & !all(addl_covar %in% colnames(data))) { stop(paste("data does not have columns with names:", paste(addl_covar, collapse = ", ")))}
  # test that approx_beyond is one of three acceptable options
  if (!(approx_beyond %in% c("expo", "zero", "carryforward"))) { stop("argument approx_beyond must be expo, zero, or carryforward") }
  # test for improper entries in columns of data
  #### Still deciding whether the following conditions should produce errors (stop()) or warnings
  if (any(data[, obs] < 0)) { warning(paste("elements of column", obs, "must be positive")) }
  if (!all(data[, event] %in% c(0, 1))) { warning(paste("elements of column", event, "must be either 0 or 1")) }
  
  # fit imputation model
  if (!is.null(addl_covar)) {
    fit.i <- survival::coxph(formula = as.formula(paste0("survival::Surv(time = ", obs, ", event = ", event, ") ~", paste0(addl_covar, collapse = "+"))),
                             data = data)
  } else {
    fit.i <- survival::survfit(formula = as.formula(paste0("survival::Surv(time = ", obs, ", event = ", event, ") ~ 1")),
                               data = data)
  }
  
  # list of imputed datasets to be returned
  imputed.datasets = list()
  n = nrow(data)
  
  # perform conditional mean imputation M times, sampling lambda from its asymptotic distribution
  for (i in 1:M) {
    # call single imputation function
    imputed.datasets[[i]] <- condl_mean_impute(fit = fit.i, obs = obs, event = event, addl_covar = addl_covar,
                                               data = data, approx_beyond = approx_beyond, sample_lambda = T)
  }
  
  # return list of imputed datasets
  return(imputed.datasets)
}
