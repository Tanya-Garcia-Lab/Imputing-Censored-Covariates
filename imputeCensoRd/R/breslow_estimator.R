#' Breslow estimator of baseline survival
#'
#' Estimates the baseline survival function from a Cox proportional hazards model following Breslow's estimator.
#'
#' @param time String column name for the observed times.
#' @param event String column name for the censoring indicator of the covariate.
#' @param hr String column name for the hazard ratios calculated from a \code{coxph} model fit.
#' @param data Dataframe containing columns \code{time}, \code{event}, and \code{hr}.
#'
#' @return A list with the following three elements:
#' \item{times}{a vector of the unique observed failure times}
#' \item{basesurv}{baseline survival estimates}
#' \item{basehaz}{baseline cumulative hazard estimates}
#'
#' @export

breslow_estimator <- function(time, event, hr, data) {
  # test for bad input
  if (!is.character(time)) { stop("argument time must be a character") }
  if (!is.character(event)) { stop("argument event must be a character") }
  if (!is.character(hr)) { stop("argument hr must be a character") }
  if (!is.data.frame(data) & !is.matrix(data)) { stop("argument data must be a data frame or a matrix") }
  # test that data contains columns with specified names
  if (!(time %in% colnames(data))) { stop(paste("data does not have column with name", time)) }
  if (!(event %in% colnames(data))) { stop(paste("data does not have column with name", event)) }
  if (!(hr %in% colnames(data))) { stop(paste("data does not have column with name", hr)) }
  # test for improper entries in columns of data
  #### Still deciding whether the following conditions should produce errors (stop()) or warnings
  if (any(data[, time] < 0)) { warning(paste("elements of column", time, "must be positive")) }
  if (!all(data[, event] == 0 | data[, event] == 1)) { warning(paste("elements of column", event, "must be either 0 or 1")) }
  if (any(data[, hr] < 0)) { warning(paste("elements of column", hr, "must be inclusively between 0 and 1"))}
  
  tj <- data[, time] # save vector of observed times 
  dj <- aggregate(x = data[, event], by = list(data[, time]), FUN = sum) # tabulate number of deaths per tj 
  colnames(dj) = c(time, event)
  dj <- dj[dj[, event] > 0, ] # subset to times with at least one death (i.e., dj > 0)
  tauj <- dj[, time] # observed failure times
  # create a dataframe for the riskset containing number of people still alive at or just before each tauj 
  riskset <- data.frame(tauj = rep(tauj, times = length(tj)),
                        tj = rep(tj, each = length(tauj)),
                        hr = rep(data[, hr], each = length(tauj)))
  riskset$atrisk <- with(riskset, tauj <= tj)
  riskset$atrisk_hr <- with(riskset, atrisk * hr)
  # denominator for baseline hazard estimate
  haz0_denom <- aggregate(x = riskset[, "atrisk_hr"], by = list(riskset[, "tauj"]), FUN = sum) # tabulate number of deaths per tj 
  colnames(haz0_denom) = c("tauj", "atrisk_hr")
  
  # baseline hazard estimate
  haz0 <- dj[, event] / haz0_denom$atrisk_hr
  # baseline cumulative hazard estimate
  cumhaz0 <- cumsum(haz0)
  # baseline survival estimate
  surv0 <- exp(- cumhaz0)
  
  return(list(times = tauj, basesurv = surv0, basecumhaz = cumhaz0))
}
