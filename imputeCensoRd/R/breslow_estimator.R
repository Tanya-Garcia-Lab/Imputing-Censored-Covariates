#' Breslow estimator of baseline survival
#'
#' Estimates the baseline survival function from a Cox proportional hazards model following Breslow's estimator.
#'
#' @param x (Optional) Value at which the baseline survival function is to be calculated. If \code{x = NULL} (the default), values at all unique event times is returned.
#' @param time Column name of observed predictor values (including censored opens). 
#' @param event Column name of censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param hr Column name of hazard ratios (HR) from a Cox proportional hazards model. 
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{HR}.
#'
#' @return If \code{x} is not \code{NULL}, then a scalar value for the baseline survival function at value \code{x}. Otherwise, a list with the following three elements is returned:
#' \item{times}{a vector of the unique observed failure times}
#' \item{basesurv}{baseline survival estimates}
#' \item{basehaz}{baseline cumulative hazard estimates}
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by summarize

breslow_estimator <- function (x = NULL, time, event, hr, data)
{
  tj <- data[, time] # all times (censored or event)
  data %>%
    dplyr::group_by(get(time)) %>%
    dplyr::summarize(d = sum(get(event))) %>%
    data.frame() -> dj
  colnames(dj) <- c(time, event)
  dj <- dj[dj[, event] > 0, ] # filter to event times
  
  if (!is.null(x)) {
    # For times prior to first event, survival = 1
    if ((min(dj[, time]) - x) > 1E-10) {
      return(1)
    } else {
      dj <- dj[(dj[, time] - x) <= 1E-10, ] # filter to event times up to at_x
    }
  }
  tauj <- dj[, time] # event times
  atrisk_hr <- sapply(X = tauj, FUN = at_risk, t = data[, time], hr = data[, hr])
  haz0 <- dj[, event]/atrisk_hr
  cumhaz0 <- cumsum(haz0)
  surv0 <- exp(- cumhaz0)
  if (!is.null(x)) {
    return(surv0[length(surv0)])
  } else {
    return(list(times = tauj, basesurv = surv0, basecumhaz = cumhaz0))
  }
}