#'
#' Semi-analytical model of Lee et al JGR 2005.
#' WARNING: This is a temporary model that doesn't considered the SZA...
#'
#'
#' @param a is a vector of absorption coefficient
#' @param b os a vector of backscattering coefficient (length is the same as a)
#'
#'
#'@author Simon Bélanger
#'@export

Kd.Lee  <- function(a,bb) {

  coeff <- c(1.320, 4.12, 0.504, 10.304)

  Kd <- coeff[1]*a + coeff[2] ^ (1 - coeff[3] * exp(-coeff[4] * a)) * bb
  return(Kd)
}
