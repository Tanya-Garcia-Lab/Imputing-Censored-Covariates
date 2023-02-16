#' Calculate size of risk set at \code{x} based on a vector of times \code{t}
#'
#' Calculates the size of the risk set at specified times. Alternatively, if \code{hr} is not \code{NULL} then it sums up HR of at-risk people (for use with \code{breslow_estimator()}).
#'
#' @param x Value at which the size of the risk set is to be calculated. 
#' @param t Numeric vector of observed predictor values (including censored opens). 
#' @param hr (Optional) Numeric vector of hazard ratios (HR) corresponding to values in \code{t}. 
#'
#' @return A scalar
#'
#' @export
#' 
at_risk <- function(x, t, hr = NULL) {
  if (!is.null(hr)) {
    sum((t >= x) * hr)
  } else {
    sum(t >= x)
  }
}