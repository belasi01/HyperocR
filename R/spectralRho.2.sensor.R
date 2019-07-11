#'
#' Compute reflectance (or radiance) at a sensor wavebands
#' from an hyperspectral or multispectral measurement using
#' relative spectral responses (RSR)
#'
#' @param waves vector of wavelenghts for the input spectral reflectance spectra
#' @param rho vector of the input spectral reflectance spectra (same lenght as waves)
#' @param SENSOR is a string for the SENSOR names. Supported sensors are
#' OLI (default),
#' MERIS,
#' S2A,
#' S2B,
#' MODISA,
#' MODIST,
#' VIIRS,
#' SeaWIFS,
#' RapidEye,
#' PlanetScope_0C0D,
#' PlanetScope_0E,
#' PlanetScope_0F_10
#'
#' @author Simon Bélanger
#' @export
#' @name spectralRho.2.sensor
spectralRho.2.sensor <- function(waves,rho,SENSOR="OLI", PLOT=FALSE) {

  if (SENSOR == "OLI") {
    nbands = ncol(OLI)-1
    rho.sensor = rep(NA, nbands)
    waves.sensor = rep(NA, nbands)

    # interpolate the input spectra to RSR
    rho.int = spline(waves, rho, xout=OLI$waves, method = "natural")$y
    rho.int[rho.int < 0] = 0
    rho.int[OLI$waves < min(waves)] = 0
    rho.int[OLI$waves > max(waves)] = 0

    # loop on wavebands
    for (i in 2:ncol(OLI)) {
      print(paste("Band : ",names(OLI)[i]))

      # integrate on non-zero indices
      ix = which(OLI[,i] > 0) # where RSR is not zero
      fx.linear <- approxfun(OLI$waves[ix], OLI[ix,i]*rho.int[ix])
      X = integrate(fx.linear, min(OLI$waves[ix]), max(OLI$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value

      fx.linear <- approxfun(OLI$waves[ix], OLI[ix,i])
      X = integrate(fx.linear, min(OLI$waves[ix]), max(OLI$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      denom = X$value

      rho.sensor[(i-1)] = num/denom

      fx.linear <- approxfun(OLI$waves[ix], OLI[ix,i]*OLI$waves[ix])
      X = integrate(fx.linear, min(OLI$waves[ix]), max(OLI$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value


      waves.sensor[(i-1)] = num/denom


    }

    # plot rho on top of RSRs
    if (PLOT) {
      p <- plot_ly(x=~OLI$waves,y=~OLI$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
        add_trace(x=~OLI$waves,y=~OLI$B2, name = 'B2', color = I("blue"), fill = 'tozeroy') %>%
        add_trace(x=~OLI$waves,y=~OLI$B3, name = 'B3', color = I("green"), fill = 'tozeroy') %>%
        add_trace(x=~OLI$waves,y=~OLI$B4, name = 'B4', color = I("red"), fill = 'tozeroy') %>%
        add_trace(x=~OLI$waves,y=~OLI$B5, name = 'B5', color = I("darkred"), fill = 'tozeroy') %>%
        add_trace(x=~OLI$waves,y=~OLI$B6, name = 'B6', color = I("grey"), fill = 'tozeroy') %>%
        add_trace(x=~OLI$waves,y=~OLI$B7, name = 'B7', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~OLI$waves,y=~OLI$B8, name = 'B8', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~OLI$waves,y=~rho.int/max(rho.int, na.rm = T), name = 'Rho', color = I("black")) %>%
        layout(xaxis = list(title = 'Wavelenght', range=c(400,2500)),
               yaxis = list(title = 'RSR'))
      print(p)
    }

    return(data.frame(waves=waves.sensor,rho=rho.sensor))
  }

  if (SENSOR == "MERIS") {
    nbands = ncol(MERIS)-1
    rho.sensor = rep(NA, nbands)
    waves.sensor = rep(NA, nbands)

    # interpolate the input spectra to RSR
    rho.int = spline(waves, rho, xout=MERIS$waves, method = "natural")$y
    rho.int[rho.int < 0] = 0
    rho.int[MERIS$waves < min(waves)] = 0
    rho.int[MERIS$waves > max(waves)] = 0

    # loop on wavebands
    for (i in 2:ncol(MERIS)) {
      print(paste("Band : ",names(MERIS)[i]))

      # integrate on non-zero indices
      ix = which(MERIS[,i] > 0) # where RSR is not zero
      fx.linear <- approxfun(MERIS$waves[ix], MERIS[ix,i]*rho.int[ix])
      X = integrate(fx.linear, min(MERIS$waves[ix]), max(MERIS$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value

      fx.linear <- approxfun(MERIS$waves[ix], MERIS[ix,i])
      X = integrate(fx.linear, min(MERIS$waves[ix]), max(MERIS$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      denom = X$value

      rho.sensor[(i-1)] = num/denom

      fx.linear <- approxfun(MERIS$waves[ix], MERIS[ix,i]*MERIS$waves[ix])
      X = integrate(fx.linear, min(MERIS$waves[ix]), max(MERIS$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value


      waves.sensor[(i-1)] = num/denom


    }

    # plot rho on top of RSRs
    if (PLOT) {
      p <- plot_ly(x=~MERIS$waves,y=~MERIS$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B2, name = 'B2', color = I("blue"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B3, name = 'B3', color = I("cyan"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B4, name = 'B4', color = I("green"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B5, name = 'B5', color = I("darkgreen"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B6, name = 'B6', color = I("orange"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B7, name = 'B7', color = I("darkorange"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B8, name = 'B8', color = I("red"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B9, name = 'B9', color = I("darkred"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B10, name = 'B10', color = I("grey"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B11, name = 'B11', color = I("grey"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B12, name = 'B12', color = I("grey"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B13, name = 'B13', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B14, name = 'B14', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~MERIS$B15, name = 'B15', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~MERIS$waves,y=~rho.int/max(rho.int, na.rm = T), name = 'Rho', color = I("black")) %>%
        layout(xaxis = list(title = 'Wavelenght', range=c(400,920)),
               yaxis = list(title = 'RSR'))
      print(p)
    }

    return(data.frame(waves=waves.sensor,rho=rho.sensor))
  }

  if (SENSOR == "S2A") {
    nbands = ncol(S2A)-1
    rho.sensor = rep(NA, nbands)
    waves.sensor = rep(NA, nbands)

    # interpolate the input spectra to RSR
    rho.int = spline(waves, rho, xout=S2A$waves, method = "natural")$y
    rho.int[rho.int < 0] = 0
    rho.int[S2A$waves < min(waves)] = 0
    rho.int[S2A$waves > max(waves)] = 0

    # loop on wavebands
    for (i in 2:ncol(S2A)) {
      print(paste("Band : ",names(S2A)[i]))

      # integrate on non-zero indices
      ix = which(S2A[,i] > 0) # where RSR is not zero
      fx.linear <- approxfun(S2A$waves[ix], S2A[ix,i]*rho.int[ix])
      X = integrate(fx.linear, min(S2A$waves[ix]), max(S2A$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value

      fx.linear <- approxfun(S2A$waves[ix], S2A[ix,i])
      X = integrate(fx.linear, min(S2A$waves[ix]), max(S2A$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      denom = X$value

      rho.sensor[(i-1)] = num/denom

      fx.linear <- approxfun(S2A$waves[ix], S2A[ix,i]*S2A$waves[ix])
      X = integrate(fx.linear, min(S2A$waves[ix]), max(S2A$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value


      waves.sensor[(i-1)] = num/denom


    }

    # plot rho on top of RSRs
    if (PLOT) {
      p <- plot_ly(x=~S2A$waves,y=~S2A$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
        add_trace(x=~S2A$waves,y=~S2A$B2, name = 'B2', color = I("blue"), fill = 'tozeroy') %>%
        add_trace(x=~S2A$waves,y=~S2A$B3, name = 'B3', color = I("green"), fill = 'tozeroy') %>%
        add_trace(x=~S2A$waves,y=~S2A$B4, name = 'B4', color = I("orange"), fill = 'tozeroy') %>%
        add_trace(x=~S2A$waves,y=~S2A$B5, name = 'B5', color = I("red"), fill = 'tozeroy') %>%
        add_trace(x=~S2A$waves,y=~S2A$B6, name = 'B6', color = I("darkred"), fill = 'tozeroy') %>%
        add_trace(x=~S2A$waves,y=~S2A$B7, name = 'B7', color = I("grey"), fill = 'tozeroy') %>%
        add_trace(x=~S2A$waves,y=~S2A$B8, name = 'B8', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~S2A$waves,y=~S2A$B9, name = 'B9', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~S2A$waves,y=~S2A$B10, name = 'B10', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~S2A$waves,y=~S2A$B11, name = 'B11', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~S2A$waves,y=~S2A$B12, name = 'B12', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~S2A$waves,y=~S2A$B13, name = 'B13', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~S2A$waves,y=~rho.int/max(rho.int, na.rm = T), name = 'Rho', color = I("black")) %>%
        layout(xaxis = list(title = 'Wavelenght', range=c(400,2350)),
               yaxis = list(title = 'RSR'))
      print(p)
    }

    return(data.frame(waves=waves.sensor,rho=rho.sensor))
  }

  if (SENSOR == "S2B") {
    nbands = ncol(S2B)-1
    rho.sensor = rep(NA, nbands)
    waves.sensor = rep(NA, nbands)

    # interpolate the input spectra to RSR
    rho.int = spline(waves, rho, xout=S2B$waves, method = "natural")$y
    rho.int[rho.int < 0] = 0
    rho.int[S2B$waves < min(waves)] = 0
    rho.int[S2B$waves > max(waves)] = 0

    # loop on wavebands
    for (i in 2:ncol(S2B)) {
      print(paste("Band : ",names(S2B)[i]))

      # integrate on non-zero indices
      ix = which(S2B[,i] > 0) # where RSR is not zero
      fx.linear <- approxfun(S2B$waves[ix], S2B[ix,i]*rho.int[ix])
      X = integrate(fx.linear, min(S2B$waves[ix]), max(S2B$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value

      fx.linear <- approxfun(S2B$waves[ix], S2B[ix,i])
      X = integrate(fx.linear, min(S2B$waves[ix]), max(S2B$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      denom = X$value

      rho.sensor[(i-1)] = num/denom

      fx.linear <- approxfun(S2B$waves[ix], S2B[ix,i]*S2B$waves[ix])
      X = integrate(fx.linear, min(S2B$waves[ix]), max(S2B$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value


      waves.sensor[(i-1)] = num/denom


    }

    # plot rho on top of RSRs
    if (PLOT) {
      p <- plot_ly(x=~S2B$waves,y=~S2B$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
        add_trace(x=~S2B$waves,y=~S2B$B2, name = 'B2', color = I("blue"), fill = 'tozeroy') %>%
        add_trace(x=~S2B$waves,y=~S2B$B3, name = 'B3', color = I("green"), fill = 'tozeroy') %>%
        add_trace(x=~S2B$waves,y=~S2B$B4, name = 'B4', color = I("orange"), fill = 'tozeroy') %>%
        add_trace(x=~S2B$waves,y=~S2B$B5, name = 'B5', color = I("red"), fill = 'tozeroy') %>%
        add_trace(x=~S2B$waves,y=~S2B$B6, name = 'B6', color = I("darkred"), fill = 'tozeroy') %>%
        add_trace(x=~S2B$waves,y=~S2B$B7, name = 'B7', color = I("grey"), fill = 'tozeroy') %>%
        add_trace(x=~S2B$waves,y=~S2B$B8, name = 'B8', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~S2B$waves,y=~S2B$B9, name = 'B9', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~S2B$waves,y=~S2B$B10, name = 'B10', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~S2B$waves,y=~S2B$B11, name = 'B11', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~S2B$waves,y=~S2B$B12, name = 'B12', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~S2B$waves,y=~S2B$B13, name = 'B13', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~S2B$waves,y=~rho.int/max(rho.int, na.rm = T), name = 'Rho', color = I("black")) %>%
        layout(xaxis = list(title = 'Wavelenght', range=c(400,2350)),
               yaxis = list(title = 'RSR'))
      print(p)
    }

    return(data.frame(waves=waves.sensor,rho=rho.sensor))
  }

  if (SENSOR == "MODISA") {
    nbands = ncol(MODISA)-1
    rho.sensor = rep(NA, nbands)
    waves.sensor = rep(NA, nbands)

    # interpolate the input spectra to RSR
    rho.int = spline(waves, rho, xout=MODISA$waves, method = "natural")$y
    rho.int[rho.int < 0] = 0
    rho.int[MODISA$waves < min(waves)] = 0
    rho.int[MODISA$waves > max(waves)] = 0

    # loop on wavebands
    for (i in 2:ncol(MODISA)) {
      print(paste("Band : ",names(MODISA)[i]))

      # integrate on non-zero indices
      ix = which(MODISA[,i] > 0) # where RSR is not zero
      fx.linear <- approxfun(MODISA$waves[ix], MODISA[ix,i]*rho.int[ix])
      X = integrate(fx.linear, min(MODISA$waves[ix]), max(MODISA$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value

      fx.linear <- approxfun(MODISA$waves[ix], MODISA[ix,i])
      X = integrate(fx.linear, min(MODISA$waves[ix]), max(MODISA$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      denom = X$value

      rho.sensor[(i-1)] = num/denom

      fx.linear <- approxfun(MODISA$waves[ix], MODISA[ix,i]*MODISA$waves[ix])
      X = integrate(fx.linear, min(MODISA$waves[ix]), max(MODISA$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value


      waves.sensor[(i-1)] = num/denom


    }

    # plot rho on top of RSRs
    if (PLOT) {
      p <- plot_ly(x=~MODISA$waves,y=~MODISA$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B2, name = 'B2', color = I("darkblue"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B3, name = 'B3', color = I("blue"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B4, name = 'B4', color = I("cyan"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B5, name = 'B5', color = I("turquoise"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B6, name = 'B6', color = I("green"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B7, name = 'B7', color = I("darkgreen"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B8, name = 'B8', color = I("darkorange"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B9, name = 'B9', color = I("red"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B10, name = 'B10', color = I("darkred"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B11, name = 'B11', color = I("grey"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B12, name = 'B12', color = I("grey"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B13, name = 'B13', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B14, name = 'B14', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B15, name = 'B15', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~MODISA$B16, name = 'B16', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~MODISA$waves,y=~rho.int/max(rho.int, na.rm = T), name = 'Rho', color = I("black")) %>%
        layout(xaxis = list(title = 'Wavelenght', range=c(400,920)),
               yaxis = list(title = 'RSR'))
      print(p)
    }

    return(data.frame(waves=waves.sensor,rho=rho.sensor))
  }

  if (SENSOR == "MODIST") {
    nbands = ncol(MODIST)-1
    rho.sensor = rep(NA, nbands)
    waves.sensor = rep(NA, nbands)

    # interpolate the input spectra to RSR
    rho.int = spline(waves, rho, xout=MODIST$waves, method = "natural")$y
    rho.int[rho.int < 0] = 0
    rho.int[MODIST$waves < min(waves)] = 0
    rho.int[MODIST$waves > max(waves)] = 0

    # loop on wavebands
    for (i in 2:ncol(MODIST)) {
      print(paste("Band : ",names(MODIST)[i]))

      # integrate on non-zero indices
      ix = which(MODIST[,i] > 0) # where RSR is not zero
      fx.linear <- approxfun(MODIST$waves[ix], MODIST[ix,i]*rho.int[ix])
      X = integrate(fx.linear, min(MODIST$waves[ix]), max(MODIST$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value

      fx.linear <- approxfun(MODIST$waves[ix], MODIST[ix,i])
      X = integrate(fx.linear, min(MODIST$waves[ix]), max(MODIST$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      denom = X$value

      rho.sensor[(i-1)] = num/denom

      fx.linear <- approxfun(MODIST$waves[ix], MODIST[ix,i]*MODIST$waves[ix])
      X = integrate(fx.linear, min(MODIST$waves[ix]), max(MODIST$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value


      waves.sensor[(i-1)] = num/denom


    }

    # plot rho on top of RSRs
    if (PLOT) {
      p <- plot_ly(x=~MODIST$waves,y=~MODIST$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B2, name = 'B2', color = I("darkblue"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B3, name = 'B3', color = I("blue"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B4, name = 'B4', color = I("cyan"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B5, name = 'B5', color = I("turquoise"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B6, name = 'B6', color = I("green"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B7, name = 'B7', color = I("darkgreen"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B8, name = 'B8', color = I("darkorange"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B9, name = 'B9', color = I("red"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B10, name = 'B10', color = I("darkred"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B11, name = 'B11', color = I("grey"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B12, name = 'B12', color = I("grey"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B13, name = 'B13', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B14, name = 'B14', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B15, name = 'B15', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~MODIST$B16, name = 'B16', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~MODIST$waves,y=~rho.int/max(rho.int, na.rm = T), name = 'Rho', color = I("black")) %>%
        layout(xaxis = list(title = 'Wavelenght', range=c(400,920)),
               yaxis = list(title = 'RSR'))
      print(p)
    }

    return(data.frame(waves=waves.sensor,rho=rho.sensor))
  }

  if (SENSOR == "VIIRS") {
    nbands = ncol(VIIRS)-1
    rho.sensor = rep(NA, nbands)
    waves.sensor = rep(NA, nbands)

    # interpolate the input spectra to RSR
    rho.int = spline(waves, rho, xout=VIIRS$waves, method = "natural")$y
    rho.int[rho.int < 0] = 0
    rho.int[VIIRS$waves < min(waves)] = 0
    rho.int[VIIRS$waves > max(waves)] = 0

    # loop on wavebands
    for (i in 2:ncol(VIIRS)) {
      print(paste("Band : ",names(VIIRS)[i]))

      # integrate on non-zero indices
      ix = which(VIIRS[,i] > 0) # where RSR is not zero
      fx.linear <- approxfun(VIIRS$waves[ix], VIIRS[ix,i]*rho.int[ix])
      X = integrate(fx.linear, min(VIIRS$waves[ix]), max(VIIRS$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value

      fx.linear <- approxfun(VIIRS$waves[ix], VIIRS[ix,i])
      X = integrate(fx.linear, min(VIIRS$waves[ix]), max(VIIRS$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      denom = X$value

      rho.sensor[(i-1)] = num/denom

      fx.linear <- approxfun(VIIRS$waves[ix], VIIRS[ix,i]*VIIRS$waves[ix])
      X = integrate(fx.linear, min(VIIRS$waves[ix]), max(VIIRS$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value


      waves.sensor[(i-1)] = num/denom


    }

    # plot rho on top of RSRs
    if (PLOT) {
      p <- plot_ly(x=~VIIRS$waves,y=~VIIRS$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
        add_trace(x=~VIIRS$waves,y=~VIIRS$B2, name = 'B2', color = I("blue"), fill = 'tozeroy') %>%
        add_trace(x=~VIIRS$waves,y=~VIIRS$B3, name = 'B3', color = I("cyan"), fill = 'tozeroy') %>%
        add_trace(x=~VIIRS$waves,y=~VIIRS$B4, name = 'B4', color = I("green"), fill = 'tozeroy') %>%
        add_trace(x=~VIIRS$waves,y=~VIIRS$B5, name = 'B5', color = I("darkred"), fill = 'tozeroy') %>%
        add_trace(x=~VIIRS$waves,y=~VIIRS$B6, name = 'B6', color = I("grey"), fill = 'tozeroy') %>%
        add_trace(x=~VIIRS$waves,y=~VIIRS$B7, name = 'B7', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~VIIRS$waves,y=~rho.int/max(rho.int, na.rm = T), name = 'Rho', color = I("black")) %>%
        layout(xaxis = list(title = 'Wavelenght', range=c(380,900)),
               yaxis = list(title = 'RSR'))
      print(p)
    }

    return(data.frame(waves=waves.sensor,rho=rho.sensor))
  }

  if (SENSOR == "SeaWiFS") {
    nbands = ncol(SeaWiFS)-1
    rho.sensor = rep(NA, nbands)
    waves.sensor = rep(NA, nbands)

    # interpolate the input spectra to RSR
    rho.int = spline(waves, rho, xout=SeaWiFS$waves, method = "natural")$y
    rho.int[rho.int < 0] = 0
    rho.int[SeaWiFS$waves < min(waves)] = 0
    rho.int[SeaWiFS$waves > max(waves)] = 0

    # loop on wavebands
    for (i in 2:ncol(SeaWiFS)) {
      print(paste("Band : ",names(SeaWiFS)[i]))

      # integrate on non-zero indices
      ix = which(SeaWiFS[,i] > 0) # where RSR is not zero
      fx.linear <- approxfun(SeaWiFS$waves[ix], SeaWiFS[ix,i]*rho.int[ix])
      X = integrate(fx.linear, min(SeaWiFS$waves[ix]), max(SeaWiFS$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value

      fx.linear <- approxfun(SeaWiFS$waves[ix], SeaWiFS[ix,i])
      X = integrate(fx.linear, min(SeaWiFS$waves[ix]), max(SeaWiFS$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      denom = X$value

      rho.sensor[(i-1)] = num/denom

      fx.linear <- approxfun(SeaWiFS$waves[ix], SeaWiFS[ix,i]*SeaWiFS$waves[ix])
      X = integrate(fx.linear, min(SeaWiFS$waves[ix]), max(SeaWiFS$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value


      waves.sensor[(i-1)] = num/denom


    }

    # plot rho on top of RSRs
    if (PLOT) {
      p <- plot_ly(x=~SeaWiFS$waves,y=~SeaWiFS$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
        add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B2, name = 'B2', color = I("blue"), fill = 'tozeroy') %>%
        add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B3, name = 'B3', color = I("cyan"), fill = 'tozeroy') %>%
        add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B4, name = 'B4', color = I("turquoise"), fill = 'tozeroy') %>%
        add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B5, name = 'B5', color = I("green"), fill = 'tozeroy') %>%
        add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B6, name = 'B6', color = I("red"), fill = 'tozeroy') %>%
        add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B7, name = 'B7', color = I("darkgrey"), fill = 'tozeroy') %>%
        add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B8, name = 'B8', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~SeaWiFS$waves,y=~rho.int/max(rho.int, na.rm = T), name = 'Rho', color = I("black")) %>%
        layout(xaxis = list(title = 'Wavelenght', range=c(389,950)),
               yaxis = list(title = 'RSR'))
      print(p)
    }

    return(data.frame(waves=waves.sensor,rho=rho.sensor))
  }

  if (SENSOR == "RapidEye") {
    nbands = ncol(RapidEye)-1
    rho.sensor = rep(NA, nbands)
    waves.sensor = rep(NA, nbands)

    # interpolate the input spectra to RSR
    rho.int = spline(waves, rho, xout=RapidEye$waves, method = "natural")$y
    rho.int[rho.int < 0] = 0
    rho.int[RapidEye$waves < min(waves)] = 0
    rho.int[RapidEye$waves > max(waves)] = 0

    # loop on wavebands
    for (i in 2:ncol(RapidEye)) {
      print(paste("Band : ",names(RapidEye)[i]))

      # integrate on non-zero indices
      ix = which(RapidEye[,i] > 0) # where RSR is not zero
      fx.linear <- approxfun(RapidEye$waves[ix], RapidEye[ix,i]*rho.int[ix])
      X = integrate(fx.linear, min(RapidEye$waves[ix]), max(RapidEye$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value

      fx.linear <- approxfun(RapidEye$waves[ix], RapidEye[ix,i])
      X = integrate(fx.linear, min(RapidEye$waves[ix]), max(RapidEye$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      denom = X$value

      rho.sensor[(i-1)] = num/denom

      fx.linear <- approxfun(RapidEye$waves[ix], RapidEye[ix,i]*RapidEye$waves[ix])
      X = integrate(fx.linear, min(RapidEye$waves[ix]), max(RapidEye$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value


      waves.sensor[(i-1)] = num/denom


    }

    # plot rho on top of RSRs
    if (PLOT) {
      p <- plot_ly(x=~RapidEye$waves,y=~RapidEye$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("blue"), name ="B1") %>%
        add_trace(x=~RapidEye$waves,y=~RapidEye$B2, name = 'B2', color = I("green"), fill = 'tozeroy') %>%
        add_trace(x=~RapidEye$waves,y=~RapidEye$B3, name = 'B3', color = I("red"), fill = 'tozeroy') %>%
        add_trace(x=~RapidEye$waves,y=~RapidEye$B4, name = 'B4', color = I("darkred"), fill = 'tozeroy') %>%
        add_trace(x=~RapidEye$waves,y=~RapidEye$B5, name = 'B5', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~RapidEye$waves,y=~rho.int/max(rho.int, na.rm = T), name = 'Rho', color = I("black")) %>%
        layout(xaxis = list(title = 'Wavelenght', range=c(400,900)),
               yaxis = list(title = 'RSR'))
      print(p)
    }

    return(data.frame(waves=waves.sensor,rho=rho.sensor))
  }

  if (SENSOR == "PlanetScope_0C0D") {
    nbands = ncol(PlanetScope_0C0D)-1
    rho.sensor = rep(NA, nbands)
    waves.sensor = rep(NA, nbands)

    # interpolate the input spectra to RSR
    rho.int = spline(waves, rho, xout=PlanetScope_0C0D$waves, method = "natural")$y
    rho.int[rho.int < 0] = 0
    rho.int[PlanetScope_0C0D$waves < min(waves)] = 0
    rho.int[PlanetScope_0C0D$waves > max(waves)] = 0

    # loop on wavebands
    for (i in 2:ncol(PlanetScope_0C0D)) {
      print(paste("Band : ",names(PlanetScope_0C0D)[i]))

      # integrate on non-zero indices
      ix = which(PlanetScope_0C0D[,i] > 0) # where RSR is not zero
      fx.linear <- approxfun(PlanetScope_0C0D$waves[ix], PlanetScope_0C0D[ix,i]*rho.int[ix])
      X = integrate(fx.linear, min(PlanetScope_0C0D$waves[ix]), max(PlanetScope_0C0D$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value

      fx.linear <- approxfun(PlanetScope_0C0D$waves[ix], PlanetScope_0C0D[ix,i])
      X = integrate(fx.linear, min(PlanetScope_0C0D$waves[ix]), max(PlanetScope_0C0D$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      denom = X$value

      rho.sensor[(i-1)] = num/denom

      fx.linear <- approxfun(PlanetScope_0C0D$waves[ix], PlanetScope_0C0D[ix,i]*PlanetScope_0C0D$waves[ix])
      X = integrate(fx.linear, min(PlanetScope_0C0D$waves[ix]), max(PlanetScope_0C0D$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value


      waves.sensor[(i-1)] = num/denom


    }

    # plot rho on top of RSRs
    if (PLOT) {
      p <- plot_ly(x=~PlanetScope_0C0D$waves,y=~PlanetScope_0C0D$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("blue"), name ="B1") %>%
        add_trace(x=~PlanetScope_0C0D$waves,y=~PlanetScope_0C0D$B2, name = 'B2', color = I("green"), fill = 'tozeroy') %>%
        add_trace(x=~PlanetScope_0C0D$waves,y=~PlanetScope_0C0D$B3, name = 'B3', color = I("red"), fill = 'tozeroy') %>%
        add_trace(x=~PlanetScope_0C0D$waves,y=~PlanetScope_0C0D$B4, name = 'B4', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~PlanetScope_0C0D$waves,y=~rho.int/max(rho.int, na.rm = T), name = 'Rho', color = I("black")) %>%
        layout(xaxis = list(title = 'Wavelenght', range=c(400,900)),
               yaxis = list(title = 'RSR'))
      print(p)
    }

    return(data.frame(waves=waves.sensor,rho=rho.sensor))
  }

  if (SENSOR == "PlanetScope_0E") {
    nbands = ncol(PlanetScope_0E)-1
    rho.sensor = rep(NA, nbands)
    waves.sensor = rep(NA, nbands)

    # interpolate the input spectra to RSR
    rho.int = spline(waves, rho, xout=PlanetScope_0E$waves, method = "natural")$y
    rho.int[rho.int < 0] = 0
    rho.int[PlanetScope_0E$waves < min(waves)] = 0
    rho.int[PlanetScope_0E$waves > max(waves)] = 0

    # loop on wavebands
    for (i in 2:ncol(PlanetScope_0E)) {
      print(paste("Band : ",names(PlanetScope_0E)[i]))

      # integrate on non-zero indices
      ix = which(PlanetScope_0E[,i] > 0) # where RSR is not zero
      fx.linear <- approxfun(PlanetScope_0E$waves[ix], PlanetScope_0E[ix,i]*rho.int[ix])
      X = integrate(fx.linear, min(PlanetScope_0E$waves[ix]), max(PlanetScope_0E$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value

      fx.linear <- approxfun(PlanetScope_0E$waves[ix], PlanetScope_0E[ix,i])
      X = integrate(fx.linear, min(PlanetScope_0E$waves[ix]), max(PlanetScope_0E$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      denom = X$value

      rho.sensor[(i-1)] = num/denom

      fx.linear <- approxfun(PlanetScope_0E$waves[ix], PlanetScope_0E[ix,i]*PlanetScope_0E$waves[ix])
      X = integrate(fx.linear, min(PlanetScope_0E$waves[ix]), max(PlanetScope_0E$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value


      waves.sensor[(i-1)] = num/denom


    }

    # plot rho on top of RSRs
    if (PLOT) {
      p <- plot_ly(x=~PlanetScope_0E$waves,y=~PlanetScope_0E$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("blue"), name ="B1") %>%
        add_trace(x=~PlanetScope_0E$waves,y=~PlanetScope_0E$B2, name = 'B2', color = I("green"), fill = 'tozeroy') %>%
        add_trace(x=~PlanetScope_0E$waves,y=~PlanetScope_0E$B3, name = 'B3', color = I("red"), fill = 'tozeroy') %>%
        add_trace(x=~PlanetScope_0E$waves,y=~PlanetScope_0E$B4, name = 'B4', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~PlanetScope_0E$waves,y=~rho.int/max(rho.int, na.rm = T), name = 'Rho', color = I("black")) %>%
        layout(xaxis = list(title = 'Wavelenght', range=c(400,900)),
               yaxis = list(title = 'RSR'))
      print(p)
    }

    return(data.frame(waves=waves.sensor,rho=rho.sensor))
  }

  if (SENSOR == "PlanetScope_0F_10") {
    nbands = ncol(PlanetScope_0F_10)-1
    rho.sensor = rep(NA, nbands)
    waves.sensor = rep(NA, nbands)

    # interpolate the input spectra to RSR
    rho.int = spline(waves, rho, xout=PlanetScope_0F_10$waves, method = "natural")$y
    rho.int[rho.int < 0] = 0
    rho.int[PlanetScope_0F_10$waves < min(waves)] = 0
    rho.int[PlanetScope_0F_10$waves > max(waves)] = 0

    # loop on wavebands
    for (i in 2:ncol(PlanetScope_0F_10)) {
      print(paste("Band : ",names(PlanetScope_0F_10)[i]))

      # integrate on non-zero indices
      ix = which(PlanetScope_0F_10[,i] > 0) # where RSR is not zero
      fx.linear <- approxfun(PlanetScope_0F_10$waves[ix], PlanetScope_0F_10[ix,i]*rho.int[ix])
      X = integrate(fx.linear, min(PlanetScope_0F_10$waves[ix]), max(PlanetScope_0F_10$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value

      fx.linear <- approxfun(PlanetScope_0F_10$waves[ix], PlanetScope_0F_10[ix,i])
      X = integrate(fx.linear, min(PlanetScope_0F_10$waves[ix]), max(PlanetScope_0F_10$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      denom = X$value

      rho.sensor[(i-1)] = num/denom

      fx.linear <- approxfun(PlanetScope_0F_10$waves[ix], PlanetScope_0F_10[ix,i]*PlanetScope_0F_10$waves[ix])
      X = integrate(fx.linear, min(PlanetScope_0F_10$waves[ix]), max(PlanetScope_0F_10$waves[ix]), subdivisions=1000, stop.on.error = FALSE)[1]
      num = X$value


      waves.sensor[(i-1)] = num/denom


    }

    # plot rho on top of RSRs
    if (PLOT) {
      p <- plot_ly(x=~PlanetScope_0F_10$waves,y=~PlanetScope_0F_10$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("blue"), name ="B1") %>%
        add_trace(x=~PlanetScope_0F_10$waves,y=~PlanetScope_0F_10$B2, name = 'B2', color = I("green"), fill = 'tozeroy') %>%
        add_trace(x=~PlanetScope_0F_10$waves,y=~PlanetScope_0F_10$B3, name = 'B3', color = I("red"), fill = 'tozeroy') %>%
        add_trace(x=~PlanetScope_0F_10$waves,y=~PlanetScope_0F_10$B4, name = 'B4', color = I("black"), fill = 'tozeroy') %>%
        add_trace(x=~PlanetScope_0F_10$waves,y=~rho.int/max(rho.int, na.rm = T), name = 'Rho', color = I("black")) %>%
        layout(xaxis = list(title = 'Wavelenght', range=c(400,900)),
               yaxis = list(title = 'RSR'))
      print(p)
    }

    return(data.frame(waves=waves.sensor,rho=rho.sensor))
  }

}
