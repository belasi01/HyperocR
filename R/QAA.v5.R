#'
#' Quasi-Analytical Algorithm to retreive spectral absorption and backscattering coefficients.
#' The current paramterisation is version 5
#'
#' @param waves is a vector of wavelenghts
#' @param Rrs os a vector of Rrs (length is the same as waves)
#'
#'
#'@author Simon Bélanger
#'@export

QAA.v5 <- function(waves, Rrs)
{
	#wl  <- c(412,443,490,510,555,670)
	bbw <- spectral.bw(waves) * 0.5
	aw <- spectral.aw(waves)

  ix443 = which.min(abs(waves-443))
  ix490 = which.min(abs(waves-490))
  ix510 = which.min(abs(waves-510))
  ix555 = which.min(abs(waves-555))
  #ix640 = which.min(abs(waves-640))
  ix667 = which.min(abs(waves-667))

  nwaves = length(waves)
	a <- rep(0,nwaves)
	bb <- rep(0,nwaves)
	a2 <- rep(0,nwaves)
	bb2 <- rep(0,nwaves)
  g0 = 0.0895
  g1 = 0.1247

	# Check if 667 chanel is available
	if (min(abs(waves-667)) > 5 ) {
	  print("No 667 nm channel")
	  Rrs667 = 1.27*Rrs[ix555]^1.47 + 0.00018*(Rrs[ix490]/Rrs[ix555])^-3.19
	  rrs667 = Rrs667 / (0.52 + 1.7*Rrs667)
	} else {
	  up.lim <- 20 * Rrs[ix555]^1.5
	  low.lim <- 0.9 * Rrs[ix555]^1.7
	  if (Rrs[ix667] > up.lim | Rrs[ix667] < low.lim ) {
	    Rrs667 = 1.27*Rrs[ix555]^1.47 + 0.00018*(Rrs[ix490]/Rrs[ix555])^-3.19
	    rrs667 = Rrs667 / (0.52 + 1.7*Rrs667)
	  } else {
	    Rrs667 = Rrs[ix667]
	    rrs667 = Rrs[ix667] / (0.52 + 1.7*Rrs[ix667])
	  }
	}


	#  STEP 0 - compute RRS from below Sea Surface
	rrs <-  Rrs / (0.52 + 1.7*Rrs)

	# Step 1 - Compute  bb/a+bb ratio
	X <-  (-g0 + sqrt(g0^2 + 4*g1*rrs)) / (2*g1)

	# Step 2 - Empirical estimation of a(555) from v5
	xsi <- log10((rrs[ix443]+rrs[ix490])/
               (rrs[ix555]+5*rrs667/rrs[ix490]*rrs667))
  a.ref <- aw[ix555] + 10^(-1.146 - 1.366*xsi -0.469*xsi^2)

	 # Step 3  - Analytical estimation of bb at 555 nm
	 bbp.ref <- (X[ix555]*a.ref)/(1-X[ix555]) - bbw[ix555]

	 # Step 4 - Empirical estimation of Nu
	 nu <- 2.0*(1 - 1.2*exp(-0.9*rrs[ix443]/rrs[ix555]))

	 # Step 5 - Spectral bbp
  bb <- bbw + bbp.ref*(waves[ix555]/waves)^nu

  # Step 6 - Spectral a
	a <-  ((1 - X)*bb)/X

	return(list(a=a,bb=bb))
	}

Kd_QAA  <- function(a,bb) {

	coeff <- c(1.320, 4.12, 0.504, 10.304)

	Kd <- coeff[1]*a + coeff[2] ^ (1 - coeff[3] * exp(-coeff[4] * a)) * bb
	return(Kd)
	}
