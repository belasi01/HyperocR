#'
#' @export
Kd.Lee2005  <- function(a,bb, theta0, LAYER="10p") {

  if (LAYER == "10p") {
    # From Eq 11 in Lee et al JGR 2005
    m1 = 1.0 + 0.005*theta0
    coeff <- c(m1, 4.18, 0.52, 10.8)

    Kd <- coeff[1]*a + coeff[2] * (1 - coeff[3] * exp(-coeff[4] * a)) * bb

  }

  if (LAYER == "1m") {

    # Use Eq 10 with coefficient from table 2
    coeff10 <- c(1.06,4.307,0.675, 1.484)
    coeff30 <- c(1.118,4.373,0.657, 1.498)
    coeff60 <- c(1.311,4.461,0.587, 2.98)

    Kd10 <- coeff10[1]*a + coeff10[2] * (1 - coeff10[3] * exp(-coeff10[4] * a)) * bb
    Kd30 <- coeff30[1]*a + coeff30[2] * (1 - coeff30[3] * exp(-coeff30[4] * a)) * bb
    Kd60 <- coeff60[1]*a + coeff60[2] * (1 - coeff60[3] * exp(-coeff60[4] * a)) * bb

    Ks = rbind(Kd10,Kd30,Kd60)
    Ks[is.na(Ks)] <-  -999
    # interpolate to the good thetas
    Kd <- apply(Ks, 2, function(x) {spline (c(10,30,60),x, xout=theta0, method="natural")$y})

    Kd[Kd < 0] <- NA

  }
  return(Kd)
}
