#' Processed in-water HyperOCR data
#'
#' This routine compute AOPs from in-water hyperspectral
#' measurements performed using HyperOCR radiometers deployed
#' at the water surface with an arm and (optionaly) with a
#' tripod. Three set-ups have been implemented so far:
#'
#' 1. Upwelling radiance near the water surface (i.e. Lu0) only
#' 2. Upwelling radiance near the water surface (i.e. Lu0) and
#' downwelling irradiance at depth Z (EdZ)
#' 3. Upwelling radiance near the water surface (i.e. Lu0) and
#' upwelling radiance at depth Z (LuZ).
#'
#'
#' @param filen is the name of a L2 file in *.dat format procuded by
#' Prosoft
#' @param Station is the station ID
#' @param Lu0.Depth is the approximate depth of the Lu0 sensor
#' (Default is 0.1 m)
#' @param EdZ.Depth is the depth of the EdZ sensor
#' (Default is NA assuming there is no EdZ sensor)
#' @param Delta.LuZ.Depth is the distance between Lu0 and LuZ sensors
#' (Default is NA assuming there is no LuZ sensor)
#' @param Ag.file is the full file name and path of the RData file containing CDOM absorption
#' @param Ap.file is the full file name and path of the RData file containing particles absorption
#' @param Bbp.file is the full file name and path of the RData file containing bbp
#'
#' @author Simon Bélanger
#' @export
#' @name process.HOCR
process.HOCR <- function (filen=filen,
                          Station=Station,
                          Lu0.Depth = 0.1,
                          EdZ.Depth = NA,
                          Delta.LuZ.Depth = NA,
                          Ag.file = NA,
                          Ap.file = NA,
                          Bbp.file = NA) {

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

  return(list(hocr = hocr,
              waves = waves,
              Lw   = Lw,
              Rrs  = Rrs,
              uRrs = uRrs,
              Kd   = Kd,
              KLu  = KLu,
              iop  = iop,
              Kd.mod = Kd.mod,
              K.CORRECTION = K.CORRECTION,
              Es   = Es.int,
              DateTime = mean(mean(hocr$Meas.Time[[1]])),
              Station=Station,
              Lu0.Depth = Lu0.Depth,
              EdZ.Depth = EdZ.Depth,
              Delta.LuZ.Depth = Delta.LuZ.Depth))
}



