#' Processed in-water HyperOCR data
#'
#' This routine compute AOPs from in-water hyperspectral
#' measurements performed using HyperOCR radiometers deployed
#' at the water surface with an arm and (optionaly) with a
#' tripod. Four set-ups have been implemented so far:
#'
#'
#' 1. Upwelling radiance near the air-water interface (i.e. Lu0) and
#' downwelling irradiance at depth Z (EdZ). Tripod set up.
#' 2. Upwelling radiance near the air-water interface (i.e. Lu0) and
#' upwelling radiance at depth Z (LuZ).
#' 3. Upwelling radiance near the air-water interface (i.e. Lu0) only
#' 4. Upwelling irradiance near the air-water interface (i.e. Eu0) only
#'
#' @param filen is the name of a L2 file in *.dat format procuded by
#' Prosoft
#' @param Station is the station ID
#' @param Lu0.Depth is the approximate depth (in meters) of the Lu0 sensor
#' (Default is NA assuming there is no Lu0 sensor)
#' @param Eu0.Depth is the approximate depth (in meters) of the Eu0 sensor
#' (Default is NA assuming there is no Eu0 sensor)
#' @param EdZ.Depth is the depth (in meters) of the EdZ sensor
#' (Default is NA assuming there is no EdZ sensor)
#' @param Delta.LuZ.Depth is the distance between Lu0 and LuZ sensors (in meters)
#' (Default is NA assuming there is no LuZ sensor)
#' @param Ag.file is the full file name and path of the RData file containing CDOM absorption
#' @param Ap.file is the full file name and path of the RData file containing particles absorption
#' @param Bbp.file is the full file name and path of the RData file containing bbp
#'
#' @details
#'    The HOCR set up needs at least on sensor looking down just beneath
#'    the air-water interface. It could be either a radiance sensor (Lu0) or
#'    an irradiance sensor (Eu0). If the set up uses a Lu0 sensor, you must
#'    provide the depth at which the sensor was submerged using the parameter
#'    Lu0.Depth. For example, if the Lu0 sensor was lowered at 10cm below the
#'    surface, Lu0.Depth = 0.1. If the set up uses a Eu0 sensor, then
#'    Eu0.Depth must be provided.
#'    NOTE: non-NA is therefore mandatory for Lu0.Depth or Eu0.Depth. This
#'    depth will be used to apply a correction to retrieve Lu0- or Eu0-. This
#'    correction needs an estimation of the diffuse attenuation coefficient
#'    for the upwelling radiance Lu or irradiance Eu (K_Lu or K_u).
#'    In the set up #1 ( Lu0 and EdZ sensors available), Kd is estimated
#'    using Ed0+ and EdZ sensors, and KLu=Kd (only true in opticaly deep water).
#'    In the set up #2 (Lu0 and LuZ), K_Lu is computed from the measurements
#'
#'    In the set up #3 and #4, K_Lu or K_u can only be approximated using the
#'    measured IOPs. No correction will be apply is IOPs are not provided.
#'
#'    If Lu0 is provided, then the code will compute Rrs=Lw/Ed.
#'    If Eu0 is provided, then the code will compute R=Eu/Ed.
#'
#'
#'    IMPORTANT: The instrument self-shadow correction has not been
#'    implemented yet.
#'@return a long list including raw and processed data.
#'
#' @author Simon Bélanger
#' @export
#' @name process.HOCR
process.HOCR <- function (filen=filen,
                          Station=Station,
                          Lu0.Depth = NA,
                          Eu0.Depth = NA,
                          EdZ.Depth = NA,
                          Delta.LuZ.Depth = NA,
                          Ag.file = NA,
                          Ap.file = NA,
                          Bbp.file = NA) {

  if (is.na(Lu0.Depth)&&
      is.na(Eu0.Depth)) {
    print("PROBLEM: you need at least the depth of
          a Lu or Eu sensor just below the sea surface.
          Lu0.Depth or Eu0.Depth should have a value.")
    return(0)
  }

  if (!is.na(Lu0.Depth)&&
      !is.na(EdZ.Depth)) {

    RADIOMETERS = "Es_Lu0_EdZ"
    print(paste("processing",filen))
    hocr <- read.hocr.L2(filen=filen, RADIOMETERS = RADIOMETERS)

    # average spectra
    Es <- apply(hocr$Meas[[1]], 2, mean)
    Lu0 <- apply(hocr$Meas[[2]], 2, mean)
    EdZ <- apply(hocr$Meas[[3]], 2, mean)

    ######## Interpolated wavelengths
    waves=seq(380,800,3)

    # Wavelenght interpolation
    df <- data.frame(waves=hocr$waves[,1], Es)
    mod = loess(Es ~ waves, data = df, span=0.04)
    Es.int = predict(mod,  waves)

    df <- data.frame(waves=hocr$waves[,2], Lu0)
    mod = loess(Lu0 ~ waves, data = df, span=0.04)
    Lu0.int = predict(mod,  waves)

    df <- data.frame(waves=hocr$waves[,3], EdZ)
    mod = loess(EdZ ~ waves, data = df, span=0.04)
    EdZ.int = predict(mod,  waves)

    # Calculate AOPs
    # Calculate Kd
    Ed0m <- Es.int *0.97
    Kd <- (log(Ed0m)-log(EdZ.int))/EdZ.Depth
    KLu <- NA


    # Uncorrected Rrs
    uRrs <- 0.54 * Lu0.int / Es.int

    # Set irradiance-derived quantity to NA
    R     <- NA
    uR    <- NA
    Eu.0m <- NA

    # Modeled Kd from measured IOPs
    if (!is.na(Ag.file)&&!is.na(Ap.file)&&!is.na(Bbp.file)) {
      if (file.exists(Ag.file)) {
        load(Ag.file)
        df.Ag <- data.frame(waves=Ag$Lambda, Ag=Ag$Ag.offset)
        Ag.int <- spline(df.Ag, xout=waves)$y
      }  else {
        print(paste("CDOM File not found", Ag.file))
        return(0)
      }

      if (file.exists(Ap.file)) {
        load(Ap.file)
        df.Ap <- data.frame(waves=A$Ap$Lambda, Ap=A$Ap$Ap.Stramski.mean)
        Ap.int <- spline(df.Ap, xout=waves)$y
      }  else {
        print(paste("Ap File not found", Ap.file))
        return(0)
      }

      if (file.exists(Bbp.file)) {
        load(Bbp.file)
        df.Bbp <-data.frame(waves=bb$waves, bb=bb$bb)
        Bbp.int <- df.Bbp$bb[2]*(df.Bbp$waves[2]/waves)^0.5
      }  else {
        print(paste("Bb File not found", Bbp.file))
        print("Ignoring bbp")
        Bbp.int <- rep(0, length(waves))
      }

      aw <- spectral.aw(waves)
      bbw <- spectral.bw(waves)
      iop <- data.frame(a=aw+Ag.int+Ap.int, bb=bbw+Bbp.int)

      Kd.mod <- Kd.Lee2005(iop$a, iop$bb, 45, LAYER="1m")

    } else {
      print("No measured IOPs or one iop file missing")
      print(paste(Ag.file, Ap.file, Bbp.file))
      Kd.mod <- rep(0, length(waves))
      iop <-NA
    }


    # Corrected Rrs assuming K_Lu == Kd
    Lw <- 0.54 * Lu0.int / exp(-Lu0.Depth*Kd)
    Rrs <- Lw / Es.int
    K.CORRECTION =  "Kd"
  }

  if (!is.na(Lu0.Depth)&&
      !is.na(Delta.LuZ.Depth)) {
    RADIOMETERS = "Es_Lu0_LuZ"
    print(paste("processing",filen))
    hocr <- read.hocr.L2(filen=filen, RADIOMETERS = RADIOMETERS)

    # average spectra
    Es <- apply(hocr$Meas[[1]], 2, mean)
    Lu0 <- apply(hocr$Meas[[2]], 2, mean)
    LuZ <- apply(hocr$Meas[[3]], 2, mean)

    ######## Interpolated wavelengths
    waves=seq(380,800,3)

    # Wavelenght interpolation
    df <- data.frame(waves=hocr$waves[,1], Es)
    mod = loess(Es ~ waves, data = df, span=0.04)
    Es.int = predict(mod,  waves)

    df <- data.frame(waves=hocr$waves[,2], Lu0)
    mod = loess(Lu0 ~ waves, data = df, span=0.04)
    Lu0.int = predict(mod,  waves)

    df <- data.frame(waves=hocr$waves[,3], LuZ)
    mod = loess(LuZ ~ waves, data = df, span=0.04)
    LuZ.int = predict(mod,  waves)

    # Calculate AOPs
    # Calculate KLu
    KLu <- (log(Lu0.int)-log(LuZ.int))/Delta.LuZ.Depth
    Kd <- NA


    # Uncorrected Rrs
    uRrs <- 0.54 * Lu0.int / Es.int


    # Set irradiance-derived quantity to NA
    R     <- NA
    uR    <- NA
    Eu.0m <- NA

    # Modeled Kd from measured IOPs
    if (!is.na(Ag.file)&&!is.na(Ap.file)&&!is.na(Bbp.file)) {
      if (file.exists(Ag.file)) {
        load(Ag.file)
        df.Ag <- data.frame(waves=Ag$Lambda, Ag=Ag$Ag.offset)
        Ag.int <- spline(df.Ag, xout=waves)$y
      }  else {
        print(paste("CDOM File not found", Ag.file))
        return(0)
      }

      if (file.exists(Ap.file)) {
        load(Ap.file)
        df.Ap <- data.frame(waves=A$Ap$Lambda, Ap=A$Ap$Ap.Stramski.mean)
        Ap.int <- spline(df.Ap, xout=waves)$y
      }  else {
        print(paste("Ap File not found", Ap.file))
        return(0)
      }

      if (file.exists(Bbp.file)) {
        load(Bbp.file)
        df.Bbp <-data.frame(waves=bb$waves, bb=bb$bb)
        Bbp.int <- df.Bbp$bb[2]*(df.Bbp$waves[2]/waves)^0.5
      }  else {
        print(paste("Bb File not found", Bbp.file))
        print("Ignoring bbp")
        Bbp.int <- rep(0, length(waves))
      }

      aw <- spectral.aw(waves)
      bbw <- spectral.bw(waves)
      iop <- data.frame(a=aw+Ag.int+Ap.int, bb=bbw+Bbp.int)

      Kd.mod <- Kd.Lee2005(iop$a, iop$bb, 45, LAYER="1m")

    } else {
      print("No measured IOPs or one iop file missing")
      print(paste(Ag.file, Ap.file, Bbp.file))
      Kd.mod <- rep(0, length(waves))
      iop <-NA
    }

    # Corrected Rrs using measured K_Lu
    Lw <- 0.54 * Lu0.int / exp(-Lu0.Depth*KLu)
    Rrs <- Lw / Es.int
    K.CORRECTION = "KLu"

  }

  if (!is.na(Lu0.Depth)&&
      is.na(Delta.LuZ.Depth)&&
      is.na(EdZ.Depth)) {
    RADIOMETERS = "Es_Lu0"
    print(paste("processing",filen))
    hocr <- read.hocr.L2(filen=filen, RADIOMETERS = RADIOMETERS)

    # average spectra
    Es <- apply(hocr$Meas[[1]], 2, mean)
    Lu0 <- apply(hocr$Meas[[2]], 2, mean)

    ######## Interpolated wavelengths
    waves=seq(380,800,3)

    # Wavelenght interpolation
    df <- data.frame(waves=hocr$waves[,1], Es)
    mod = loess(Es ~ waves, data = df, span=0.04)
    Es.int = predict(mod,  waves)

    df <- data.frame(waves=hocr$waves[,2], Lu0)
    mod = loess(Lu0 ~ waves, data = df, span=0.04)
    Lu0.int = predict(mod,  waves)


    # Calculate AOPs
    # Calculate KLu
    KLu <- NA
    Kd <- NA

    # Uncorrected Rrs
    uRrs <- 0.54 * Lu0.int / Es.int

    # Set irradiance-derived quantity to NA
    R     <- NA
    uR    <- NA
    Eu.0m <- NA

    # estimate IOPs using QAA and Kd
    # Modeled Kd from measured IOPs
    if (!is.na(Ag.file)&&!is.na(Ap.file)&&!is.na(Bbp.file)) {
      if (file.exists(Ag.file)) {
        load(Ag.file)
        df.Ag <- data.frame(waves=Ag$Lambda, Ag=Ag$Ag.offset)
        Ag.int <- spline(df.Ag, xout=waves)$y
      }  else {
        print(paste("CDOM File not found", Ag.file))
        return(0)
      }

      if (file.exists(Ap.file)) {
        load(Ap.file)
        df.Ap <- data.frame(waves=A$Ap$Lambda, Ap=A$Ap$Ap.Stramski.mean)
        Ap.int <- spline(df.Ap, xout=waves)$y
      }  else {
        print(paste("Ap File not found", Ap.file))
        return(0)
      }

      if (file.exists(Bbp.file)) {
        load(Bbp.file)
        df.Bbp <-data.frame(waves=bb$waves, bb=bb$bb)
        Bbp.int <- df.Bbp$bb[2]*(df.Bbp$waves[2]/waves)^0.5
      }  else {
        print(paste("Bb File not found", Bbp.file))
        print("Ignoring bbp")
        Bbp.int <- rep(0, length(waves))
      }

      aw <- spectral.aw(waves)
      bbw <- spectral.bw(waves)
      iop <- data.frame(a=aw+Ag.int+Ap.int, bb=bbw+Bbp.int)

      Kd.mod <- Kd.Lee2005(iop$a, iop$bb, 45, LAYER="1m")
      K.CORRECTION = "Kd.modeled"
    } else {
      print("No measured IOPs or one iop file missing")
      print(paste(Ag.file, Ap.file, Bbp.file))
      Kd.mod <- rep(0, length(waves))
      print("Kd.mod set to 0 so no correction applied")
      iop <-NA
      K.CORRECTION = NA
    }

    # Corrected Rrs using modeled Kd
    Lw <- 0.54 * Lu0.int / exp(-Lu0.Depth*Kd.mod)
    Rrs <- Lw / Es.int


  }

  if (!is.na(Eu0.Depth)) {
    RADIOMETERS = "Es_Eu0"
    print(paste("processing",filen))
    hocr <- read.hocr.L2(filen=filen, RADIOMETERS = RADIOMETERS)

    # average spectra
    Es <- apply(hocr$Meas[[1]], 2, mean)
    Eu0 <- apply(hocr$Meas[[2]], 2, mean)

    ######## Interpolated wavelengths
    waves=seq(380,800,3)

    # Wavelenght interpolation
    df <- data.frame(waves=hocr$waves[,1], Es)
    mod = loess(Es ~ waves, data = df, span=0.04)
    Es.int = predict(mod,  waves)

    df <- data.frame(waves=hocr$waves[,2], Eu0)
    mod = loess(Eu0 ~ waves, data = df, span=0.04)
    Eu0.int = predict(mod,  waves)


    # Calculate AOPs
    # Calculate KLu
    KLu <- NA
    Kd <- NA

    # Uncorrected R above water
    Rrs <- NA
    uRrs <-NA
    Lw <- NA
    uR <- Eu0.int / Es.int / 0.96

    # estimate IOPs using QAA and Kd
    # Modeled Kd from measured IOPs
    if (!is.na(Ag.file)&&!is.na(Ap.file)&&!is.na(Bbp.file)) {
      if (file.exists(Ag.file)) {
        load(Ag.file)
        df.Ag <- data.frame(waves=Ag$Lambda, Ag=Ag$Ag.offset)
        Ag.int <- spline(df.Ag, xout=waves)$y
      }  else {
        print(paste("CDOM File not found", Ag.file))
        return(0)
      }

      if (file.exists(Ap.file)) {
        load(Ap.file)
        df.Ap <- data.frame(waves=A$Ap$Lambda, Ap=A$Ap$Ap.Stramski.mean)
        Ap.int <- spline(df.Ap, xout=waves)$y
      }  else {
        print(paste("Ap File not found", Ap.file))
        return(0)
      }

      if (file.exists(Bbp.file)) {
        load(Bbp.file)
        df.Bbp <-data.frame(waves=bb$waves, bb=bb$bb)
        Bbp.int <- df.Bbp$bb[2]*(df.Bbp$waves[2]/waves)^0.5
      }  else {
        print(paste("Bb File not found", Bbp.file))
        print("Ignoring bbp")
        Bbp.int <- rep(0, length(waves))
      }

      aw <- spectral.aw(waves)
      bbw <- spectral.bw(waves)
      iop <- data.frame(a=aw+Ag.int+Ap.int, bb=bbw+Bbp.int)

      Kd.mod <- Kd.Lee2005(iop$a, iop$bb, 45, LAYER="1m")
      K.CORRECTION = "Kd.modeled"
    } else {
      print("No measured IOPs or one iop file missing")
      print(paste(Ag.file, Ap.file, Bbp.file))
      Kd.mod <- rep(0, length(waves))
      print("Kd.mod set to 0 so no correction applied")
      iop <-NA
      K.CORRECTION = NA
    }
    # Corrected Rrs using modeled Kd
    Eu.0m <- Eu0.int / exp(-Eu0.Depth*Kd.mod)
    R <- Eu.0m / Es.int / 0.96
  }

  return(list(hocr = hocr,
              waves= waves,
              Lw   = Lw,
              Eu.0m= Eu.0m,
              Rrs  = Rrs,
              uRrs = uRrs,
              R    = R,
              uR   = uR,
              Kd   = Kd,
              KLu  = KLu,
              iop  = iop,
              Kd.mod = Kd.mod,
              K.CORRECTION = K.CORRECTION,
              Es   = Es.int,
              DateTime = mean(mean(hocr$Meas.Time[[1]])),
              Station=Station,
              Lu0.Depth = Lu0.Depth,
              Eu0.Depth = Eu0.Depth,
              EdZ.Depth = EdZ.Depth,
              Delta.LuZ.Depth = Delta.LuZ.Depth))
}



