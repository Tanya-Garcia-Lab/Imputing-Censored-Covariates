#' Impute survival at censored time
#'
#' Imputes survival at censored covariate with mean of uncensored covariates before and after.
#'
#' @param at_time A scalar time
#' @param time String column name for the censored covariate.
#' @param event String column name for the censoring indicator of the covariate.
#' @param surv String column name for the survival estimate.
#' @param data Datafrane containing columns \code{obs}, \code{event}, and \code{surv}.
#'
#' @return Scalar survival estimate for value \code{at_time}.
#'
#' @export
impute_censored_surv <- function(at_time, time, event, surv, data) {
  # test for bad input
  if (!is.character(time)) { stop("argument time must be a character") }
  if (!is.character(event)) { stop("argument event must be a character") }
  if (!is.character(surv)) { stop("argument surv must be a character") }
  if (!is.data.frame(data) & !is.matrix(data)) { stop("argument data must be a data frame or a matrix") }
  # test that data contains columns with specified names
  if (!(time %in% colnames(data))) { stop(paste("data does not have column with name", time)) }
  if (!(event %in% colnames(data))) { stop(paste("data does not have column with name", event)) }
  if (!(surv %in% colnames(data))) { stop(paste("data does not have column with name", surv)) }
  # test for improper entries in columns of data
  #### Still deciding whether the following conditions should produce errors (stop()) or warnings
  if (any(data[, time] < 0)) { warning(paste("elements of column", time, "must be positive")) }
  if (!all(data[, event] %in% c(0, 1))) { warning(paste("elements of column", event, "must be either 0 or 1")) }
  if (any(data[, surv][!is.na(data[, surv])] < 0 | data[, surv][!is.na(data[, surv])] > 1)) { warning(paste("elements of column", surv, "must be inclusively between 0 and 1"))}
  
  # which (if any) event times are equal to at_time?
  same_time <- which(round(data[, time] - at_time, 8) == 0 & data[, event] == 1 & !is.na(data[, surv]))
  
  # if no event times are equal to at_time, impute with the mean of values immediately before/after
  if (length(same_time) == 0) {
    # index of greatest event time less than at_time
    before <- which(data[, time] <= at_time & data[, event] == 1)
    # corresponding survival estimate
    surv_before <- data[max(before), surv]
   
    # index of smallest event time greater than at_time
    after <- which(data[, time] > at_time & data[, event] == 1)
    # corresponding survival estimate
    surv_after <- data[min(after), surv]
    
    # average the above survival estimates
    return((surv_before + surv_after) / 2)
  } else {
    return(data[max(same_time), surv])
  }
}
