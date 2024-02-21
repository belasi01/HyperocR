#'@title Compute the remote sensing reflectance from HyperSAS data
#'
#'@description
#'This is the main function of the package. It processes the SAS data
#'for a given data acquisition file.
#'
#'
#'@param SAS is a list return by \code{\link{read.horc.L2.SAS}}.
#'@param tilt.max is the maximum tilt tolerance for the SAS.
#'A default value is 3 degrees is used otherwise.
#'@param quantile.prob is a value (betwwen 0.25 to 1) for the maximum quantile probability
#'for which the surface radiance (Lt) values will be discarded.
#'The quantile probability is the value at which the probability
#'of the random variable being less than or equal to this value.
#'For example, a value of 0.5 means that every Lt value above the
#'50\% quantile probability will be discarded, while a value of 0.8 means
#'that only the values above the 80\% quantile probability will be discarded.
#'The rational is to eliminate outliers resulting from sun glitters, foam, sea ice, etc.
#'The default value is 0.5 to eliminate Lt value above the 50\% quantile probalibity.
#'@param windspeed is the wind speed in m/s.
#'A default value of 5 m/s is used otherwise.
#'@param thetaV is the viewing sensor zenith angle in degree.
#'A default value of 40 degrees is used otherwise.
#'@param Dphi is the azimuthal difference between the sun and the sensor view angle.
#'A default value of 90 degrees is used otherwise.
#'@param NIR.CORRECTION is the method used to eliminate residual "white reflectance".
#'Such residual is frequent under windy condisitons at sea.
#'Those may be due to sun glint, foam, ocean spray, etc.
#'The "NULL" correction subtracts the residual Rrs value at 800 nm.
#'The "SIMILARITY" correction is prefer in turbid waters. It is based on
#'the similarity spectrum (see Ruddick et al. L&O 2005) in the NIR. Here it
#'assumes a constant Rrs(780)/Rrs(720) ratio of 2.35.
#'@param use.COMPASS is a logical parameter indicating whether or not the azimuth difference
#'between the sun and the view geometry is calculated using the SAS compass data.
#'When FALSE, the azimuth angle difference is taken
#'from the cast.info.dat file. NOTE: The Default is FALSE because
#'compass usually doesn't work on ship.
#'
#'@details
#'This functions is the main part of the HyperSAS data processing.
#'
#'First it computes the tilt from the Pitch and Roll data recorded by the SAS
#'and stored in SAS$Anc data frame. Then tilt values areassociated to
#'each irradiance (Ed), sky radiance (Li) and surface radiance (Lt)
#'frame found whitin the file (each of them have their
#'specific time integration).
#'Second, the Lt frames are discarded if
#'1) the tilt is greater than the specified tilt.max;
#'2) the Lt value around 490 nm is greater than the specified quantile.prob
#'(to remove upper outliers); and
#'3) the Lt value around 490 nm is below the 10\% quantile probalility
#'(to remove lower outliers).
#'
#'Second, for each valid Lt frame, a corresponding value of Ed and Li is interpolated
#'using a spline function. Mean and standard devivation of valid Lt
#'and interpolated Ed and Li spectra are computed.
#'
#'Third, a spectral interpolation is performed on the average
#'Ed, Li and Lt spectra from 380 to 800 nm every 5 nm.
#'This step is done using a \code{\link{loess}} function with a span of 0.05.
#'
#'Forth, the sky reflectance (Li/Ed) and surface (Lt/Ed) are computed and smoothed using a
#'\code{\link{loess}} function with a span of 0.75 and 0.1, respectively. This
#'step is necessary to avoid spectral artefact resulting from the individual spectral
#'response of each sensor.
#'
#'Fift, the specular reflectance of the air-sea interface is determined
#'(i.e., rho_{sky}).
#'This is based on Ruddick et al. L&0 2005 and Mobley, App. Opt. 1999.
#'Under clear sky, i.e. when the sky reflectance at 750 nm is below 5\%,
#'rho_{sky} is interpolated from a Look-Up-Table provided by Mobley.
#'The rho_{sky} LUT has 4 dimensions for thetaS, thetaV, Dphi
#'and Windspeed. The current table is not spectrally dependent
#'(so strictly valid for 550 nm). Spectral dependency is therefore
#'not taken into account, which may be significant
#'(see Lee et al, Opt. Exp. 2010)
#'
#'
#'@seealso \code{\link{process.HyperSAS}} and \code{\link{read.hocr.SAS}}
#'@author Simon Bélanger
#'@export
#'@name compute.Rrs.SAS
compute.Rrs.SAS <- function(SAS,
                            tilt.max= 3,
                            quantile.prob = 0.5,
                            windspeed = 5,
                            thetaV = 40,
                            Dphi=90,
                            NIR.CORRECTION="NULL",
                            Good=1,
                            use.COMPASS=FALSE) {

  #### Compute Tilt from Roll and Pitch for each sensor
  d2r <- pi / 180

  # tilt
  if (is.na(SAS$Anc)) {

    ###### I AM HERE NEED TO DECIDE WHAT WE DO WHEN TILT IS ABSENT
    tilt=NA
    tilt.Lt = NA
    tilt.Ed = NA
    tilt.Li = NA
    phiv.mag.Lt = NA
    ix.Lt.good = which(SAS$Lt[,35] < quantile(SAS$Lt[,35], probs = quantile.prob) &
                       SAS$Lt[,35] > quantile(SAS$Lt[,35], probs = 0.1))
    phiv.mag = NA
    declination = NA
    phiv.true = NA
    delta.phi =NA

  } else {
    Roll <- SAS$Anc$ROLL
    Pitch <- SAS$Anc$PITCH
    tilt <- atan(sqrt(tan(Roll*d2r)^2+tan(Pitch*d2r)^2))/d2r

    tilt.Lt = spline(SAS$Anc.Time, tilt, xout=SAS$Lt.Time)$y
    tilt.Ed = spline(SAS$Anc.Time, tilt, xout=SAS$Ed.Time)$y
    tilt.Li = spline(SAS$Anc.Time, tilt, xout=SAS$Li.Time)$y

    #####Compute sensor azimuth for each Lt measurements
    phiv.mag.Lt = spline(SAS$Anc.Time, SAS$Anc$COMP, xout=SAS$Lt.Time)$y

    # Trim data to remove high tilt and the last quantile
    ix.Lt.good = which(tilt.Lt < tilt.max &
                         SAS$Lt[,35] < quantile(SAS$Lt[,35], probs = quantile.prob) &
                         SAS$Lt[,35] > quantile(SAS$Lt[,35], probs = 0.1))
  }






  Lt = SAS$Lt[ix.Lt.good,]

  Lt.time.mean = mean(SAS$Lt.Time[ix.Lt.good])

  Ed=apply(SAS$Ed,2, function(x){tmp <- smooth.spline(as.numeric(SAS$Ed.Time),x)
                                   predict(tmp,as.numeric(SAS$Lt.Time[ix.Lt.good]))$y})

  Li=apply(SAS$Li,2, function(x){tmp <- smooth.spline(as.numeric(SAS$Li.Time),x)
                                 predict(tmp,as.numeric(SAS$Lt.Time[ix.Lt.good]))$y})


  ######### Average data
  Ed.mean = apply(Ed, 2, mean)
  Li.mean = apply(Li, 2, mean)
  Lt.mean = apply(Lt, 2, mean)

  Ed.sd = apply(Ed, 2, sd)
  Li.sd = apply(Li, 2, sd)
  Lt.sd = apply(Lt, 2, sd)


  ######### Interpolated wavelengths
  waves=seq(350,810,5)


  # Wavelength interpolation
  df = as.data.frame(cbind(SAS$Ed.wl, Ed.mean))
  names(df) = c("wl","Ed")
  mod = loess(Ed ~ wl, data = df, span=0.05,
              control = loess.control(surface = "direct"))
  Ed.int = predict(mod,  waves)

  df = as.data.frame(cbind(SAS$Li.wl, Li.mean))
  names(df) = c("wl","Li")
  mod = loess(Li ~ wl, data = df, span=0.05,
              control = loess.control(surface = "direct"))
  Li.int = predict(mod,  waves)

  df = as.data.frame(cbind(SAS$Lt.wl, Lt.mean))
  names(df) = c("wl","Lt")
  mod = loess(Lt ~ wl, data = df, span=0.05,
              control = loess.control(surface = "direct"))
  Lt.int = predict(mod,  waves)

  ######## Compute Sky reflectance and rho
  sky = Li.int/Ed.int
  mod= loess(sky~waves, data=data.frame(waves=waves, sky=sky), span=0.75)
  sky.smooth = predict(mod, waves)

  sea =Lt.int/Ed.int
  sea[77:78] = NA # remove the oxygen band
  mod= loess(sea~waves, data=data.frame(waves=waves, sea=sea), span=0.10)
  sea.smooth = predict(mod, waves)




  ######### Get rho from MOBLEY LUT or use a constant rho of 0.0256 if cloudy
  # Find wavelength indices
  ix350 = which.min(abs(waves - 350))
  ix380 = which.min(abs(waves - 380))
  ix490 = which.min(abs(waves - 490))
  ix720 = which.min(abs(waves - 720))
  ix750 = which.min(abs(waves - 750))
  ix780 = which.min(abs(waves - 780))
  ix810 = which.min(abs(waves - 810))



  ##### Compute sun-viweing geometry
  thetaS = SAS$dd$sunzen
  phiS = SAS$dd$sunazim

  if (sky.smooth[ix750] >= 0.05){
    #Then  CLOUDY SKY (Ruddick et al L&O 2006, eq. 23, 24)
    print("Cloudy sky")
    rho = 0.0256
    CLEARSKY = FALSE

  }  else {
    CLEARSKY = TRUE

    # Then interpolate within the LUT provided by MOBLEY (for 550 nm)
    if (use.COMPASS) {
      if (is.na(delta.phi)) {
        print("WARNING: Calculated delta phi taken in cast.info.dat")
        use.COMPASS = FALSE
      } else {
        if (abs(Dphi - delta.phi) > 10) {
          print("WARNING: Calculated delta phi does not match the logged Dphi in cast.info.dat")
          use.COMPASS = FALSE
        }

        if (delta.phi < 80) {
          print("WARNING: delta phi is below 80 degrees")
          use.COMPASS = FALSE
        }

        if (delta.phi > 150) {
          print("WARNING: delta phi is greater than 150 degrees")
          use.COMPASS = FALSE
        }
      }

    }

    #

    if (use.COMPASS) {
      rho = get.rho550(thetaV, delta.phi, windspeed,thetaS)
    } else {
      rho = get.rho550(thetaV, Dphi, windspeed,thetaS)
    }

  }


  ###################
  Rrs = sea.smooth - (rho*sky.smooth) ### Mobley LUT standard method

  # Apply a correction

  #####
  ###### Method 1.  A standard NULL correction

  offset = Rrs[ix810]
  Rrs.NULL = Rrs - offset # Apply NIR correction

  #####
  ###### Methods 2. Estimation of the NIR Rrs offset correction based on
  #      Ruddick et al L&O 2006, SPIE 2005
  offset = 2.35*Rrs[ix780] - Rrs[ix720]/(2.35-1)
  Rrs.SIMILARITY = Rrs - offset # Apply NIR correction

  #####
  ###### Method 3. Estimation of the rho.fresnel assuming BLACK Pixel assumption in the NIR
  rho.sky.NIR =  (mean(sea.smooth[ix780:ix810], na.rm = T) /
                    mean(sky.smooth[ix780:ix810], na.rm = T))
  Rrs.BP = sea.smooth - (rho.sky.NIR*sky.smooth)

  #####
  ###### Method 4. Estimation of the rho.fresnel assuming BLACK Pixel assumption in the UV
  rho.sky.UV =  (mean(sea.smooth[ix350:ix380], na.rm = T) /
                   mean(sky.smooth[ix350:ix380], na.rm = T))
  Rrs.UV = sea.smooth - (rho.sky.UV*sky.smooth)

  #####
  ###### Method 5. Estimation of the rho.fresnel assuming BLACK Pixel assumption in both UV and NIR (spectrally dependent)
  rho.sky.UV.NIR = spline(c(mean(waves[ix350:ix380], na.rm = T),
                            mean(waves[ix780:ix810], na.rm = T)),
                          c(rho.sky.UV, rho.sky.NIR),
                          xout = waves)$y
  Rrs.UV.NIR = sea.smooth - (rho.sky.UV.NIR*sky.smooth)

  #####
  ##### Method 6: Implementation of Kutser et al. 2013 for removal of sky glint
  #if (VERBOSE) print("Begining the Kutser correction for glint")
  UVdata <- sea.smooth[ix350:ix380]
  NIRdata <- sea.smooth[ix780:ix810]
  UVwave <- waves[ix350:ix380]
  NIRwave <- waves[ix780:ix810]
  #if (VERBOSE) print("Wavelength and data at UV and NIR binned")
  UV.NIR.data <- c(UVdata,NIRdata)
  UV.NIR.wave <- c(UVwave,NIRwave)
  Kutserdata <- data.frame(UV.NIR.wave, UV.NIR.data)
  names(Kutserdata) <- c("waves", "urhow")
  #if (VERBOSE) print("Starting the NLS")
  glint.fit <- nls(urhow ~ b*waves^z,start = list(b = 1, z = -1),data=Kutserdata,
                   control = nls.control(maxiter = 300))
  p <- coef(glint.fit)
  print(p)
  Kutserestimate <- p[1]*(waves)^p[2]
  Rrs.Kutser <- sea.smooth - Kutserestimate


  return(list(
    Rrs.wl = waves,
    Rrs = Rrs,
    Rrs.NULL = Rrs.NULL,
    Rrs.SIMILARITY1 = Rrs.SIMILARITY,
    Rrs.NIR = Rrs.BP,
    Rrs.UV = Rrs.UV,
    Rrs.UV.NIR = Rrs.UV.NIR,
    Rrs.Kutser = Rrs.Kutser,
    rho.sky = rho,
    rho.sky.NIR = rho.sky.NIR,
    rho.sky.UV = rho.sky.UV,
    rho.sky.UV.NIR = rho.sky.UV.NIR,
    offset = offset,
    Ed=Ed,
    Lt=Lt,
    Li=Li,
    tilt=tilt,
    tilt.Ed=tilt.Ed,
    tilt.Li=tilt.Li,
    tilt.Lt=tilt.Lt,
    Ed.mean = Ed.mean,
    Ed.sd = Ed.sd,
    Lt.mean = Lt.mean,
    Lt.sd = Lt.sd,
    Li.mean = Li.mean,
    Li.sd = Li.sd,
    ix.Lt.good=ix.Lt.good,
    DateTime = Lt.time.mean,
    ThetaV = thetaV,
    ThetaS = thetaS,
    PhiS = phiS,
    PhiV.true = phiv.true,
    PhiV.mag = phiv.mag,
    Dphi.comp = delta.phi,
    Dphi.log = Dphi,
    Windspeed = windspeed,
    CLEARSKY = CLEARSKY,
    Good = Good,
    use.COMPASS = use.COMPASS,
    dd = SAS$dd
    ))

}
