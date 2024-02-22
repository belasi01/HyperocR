#'
#'@title Process a HyperSAS data folder
#'
#'@description
#'This function is called by the higher level function \code{\link{HyperSAS.go}}. .
#'It can be also called in the command line to process a given data directory.
#'
#'
#'@param dirdat is the current directory path contaning the files to process.
#'@param TYPE is either equal to "STATION" or "TRANSIT" depending
#'whether the folder gather replicate Rrs of the same station, or the data
#'collected during the course of a ship transit. If the function is
#'called from \code{\link{HyperSAS.go}}, the TYPE will be read in the file named
#'directories.for.HyperSAS.dat. Otherwise the default is "STATION".
#'@param use.SAS.RData is a logical parameter indicating whether or not a SAS.raw.RData file
#'already exists in the data directory and do not need to be generated again. Note that
#'SAS.raw.RData is an output of the function \code{\link{process.HyperSAS}}.  Default is FALSE.
#'@param use.COMPASS is a logical parameter indicating whether or not the azimuth difference
#'between the sun and the view geometry is calculated using the SAS compass data.
#'When FALSE, the azimuth angle difference is taken
#'from the cast.info.dat file. NOTE: The Default is FALSE because
#'compass usually doesn't work on ship.
#'
#'@return The function returns a list containing the Rrs data.
#'The same list (RRS) is saved into a RRS.RData file in dirdat.
#'If use.SAS.RData is set to FALSE, a list with raw data (SAS)
#'will be saved in file SAS.raw.RData in dirdar.
#'
#'@details The function first looks into the cast.info.dat file found
#'in the dirdat and check whether all the mandatory fields are present.
#'If not, the processing may be stop or a default value is taken.
#'Next the HyperSAS files to process are red and stored in a list of lists
#'named SAS. Finally, the function calls \code{\link{compute.Rrs.SAS}} to
#'compute an Rrs spectrum for each data file.
#'
#'@seealso
#'\code{\link{compute.Rrs.SAS}} and \code{\link{HyperSAS.go}}
#'
#'@author Simon Bélanger
#'@export
#'@name process.HyperSAS
process.HyperSAS<- function(dirdat,
                            TYPE="STATION",
                            use.SAS.RData=FALSE,
                            use.COMPASS=FALSE) {

  # Cast info file
  default.cast.info.file <- paste( Sys.getenv("R_HOCR_DATA_DIR"), "cast.info.dat", sep = "/")

  cast.info.file <- paste(dirdat, "cast.info.dat", sep = "/")
  if(!file.exists(cast.info.file)) {
    file.copy(from = default.cast.info.file, to = cast.info.file)
    cat("EDIT file", cast.info.file, "and CUSTOMIZE IT\n")
    cat("This file contains the information on viewing geometry,\n")
    cat("wind speed, wind speed units, as well as processing parameters\n")
    cat("such as the quantile probility, maximum tilt tolerance,\n")
    cat("NIR correction method. See help for more details.")
    stop("Abort processing...")
  }

  cast.info <- read.table(cast.info.file, header=T)


  # Check if cast.info contains all required information
  if (is.null(cast.info$Year)) {
    print("Year not found in cast.info.dat")
    stop()
  }
  if (is.null(cast.info$Day)) {
    print("Day (of year) not found in cast.info.dat")
    stop()
  } else {
    nreplicates = length(cast.info$Day)
  }
  if (is.null(cast.info$Time)) {
    print("Time (HHMMSS) not found in cast.info.dat")
    stop()
  }
  if (is.null(cast.info$Dphi)) {
    print("Dphi (delta azimuth betwen sun and sensor viewing aximuth) not found in cast.info.dat")
    print("Mandatory if USE.COMPASS=FALSE")
    stop()
  }
  if (is.null(cast.info$ThetaV)) {
    print("ThetaV not found in cast.info.dat")
    stop()
  }
  if (is.null(cast.info$Windspeed)) {
    print("Windspeed not found in cast.info.dat")
    print("Assumes wind = 4 m/s")
    cast.info$Windspeed = rep(4,nreplicates)
  } else {
    if (is.null(cast.info$Wind.units)) {
      print("Windspeed Units (i.e., Kts or m.s or km.h) not found in cast.info.dat")
      stop()
    } else
      {
        for (i in 1:nreplicates) {
          if (cast.info$Wind.units[i] == "Kts") {
            print("Convert wind speed from Knots to m/s" )
            cast.info$Windspeed[i] = cast.info$Windspeed[i]/1.9426
          }
          if (cast.info$Wind.units[i] == "Km.h" | cast.info$Wind.units[i] == "Km/h") {
            print("Convert wind speed from Km/h to m/s" )
            cast.info$Windspeed[i] = cast.info$Windspeed[i]/3.6
          }
        }
    }
  }

  if (is.null(cast.info$quantile.prob)) {
    print("quantile.prob not found in cast.info.dat")
    print("Assumes 0.5")
    cast.info$quantile.prob = rep(0.5,nreplicates)
  }
  if (is.null(cast.info$tilt.max)) {
    print("tilt.max not found in cast.info.dat")
    print("Assumes tilt.max of 3 degrees")
    cast.info$tilt.max = rep(3,nreplicates)
  }
  if (is.null(cast.info$NIR.CORRECTION)) {
    print("NIR.CORRECTION method not found in cast.info.dat")
    print("Assumes NULL")
  cast.info$NIR.CORRECTION = rep("NULL",nreplicates)
  }

  if (is.null(cast.info$good)) {
    print("Good spectra not found in cast.info.dat")
    print("Assumes that all replicates are good (i.e. =1)")
    cast.info$good = rep(1,nreplicates)
  }

  if (TYPE == "STATION") {

    # construct a file name vector from cast info
    filen = paste(cast.info$Year,"-", cast.info$Day,"-", cast.info$Time, "_L2.dat", sep="")
    nb.rec = length(filen)

    if (nb.rec > 2) {
      print("More than two Ed0+ file")

      if (use.SAS.RData) {
        print(paste("Load: ",dirdat, "/SAS.raw.RData", sep="" ))
        load(paste(dirdat, "SAS.raw.RData", sep="/"))
        nb.rec = length(SAS)
        } else {
          SAS = lapply(filen, read.hocr.L2.SAS, VERBOSE=TRUE)
          save(file = paste(dirdat, "SAS.raw.RData", sep="/"), SAS)
      }

      RRS = list()
      all.rrs = matrix(NA, nrow=length(filen), ncol=93)
      method.selected <- cast.info$NIR.CORRECTION
      AVW.selected <- rep(NA,nb.rec)
      NDI.selected <- rep(NA,nb.rec)
      QWIP.selected <- rep(NA,nb.rec)
      QWIP.score.selected <- rep(NA,nb.rec)
      FU.selected <- rep(NA,nb.rec)
      cast.datetime<-rep(NA,nb.rec)

      for (i in 1:nb.rec) {
        print(paste("Process data collected at:", SAS[[i]]$dd$date))
        tmp = compute.Rrs.SAS(SAS[[i]],
                      tilt.max= cast.info$tilt.max[i],
                      quantile.prob=cast.info$quantile.prob[i],
                      windspeed=cast.info$Windspeed[i],
                      NIR.CORRECTION=cast.info$NIR.CORRECTION[i],
                      thetaV=cast.info$ThetaV[i],
                      Dphi=cast.info$Dphi[i],
                      Good=cast.info$good[i],
                      use.COMPASS)

        plot.SAS.Rrs(tmp, PNG=TRUE)
        RRS[[i]] <- tmp # store the list in a list of list
        # store the Rrs in a matrix for further mean calculation using the selected method
        ix.method <- which(tmp$methods == method.selected[i])
        all.rrs[i,]           <- tmp$Rrs[ix.method,]
        AVW.selected[i]       <- tmp$AVW[ix.method]
        NDI.selected[i]       <- tmp$NDI[ix.method]
        QWIP.selected[i]      <- tmp$QWIP[ix.method]
        QWIP.score.selected[i]<- tmp$QWIP.score[ix.method]
        FU.selected[i]        <- tmp$FU[ix.method]
        cast.datetime[i]      <- RRS[[i]]$dd$date
      }

  # Average the good RRS
      ix.good = which(cast.info$good == 1)
      Rrs.mean = apply(all.rrs[ix.good,], 2, mean, na.rm=T)
      Rrs.sd = apply(all.rrs[ix.good,], 2, sd, na.rm=T)
      RRS$Rrs.mean =Rrs.mean
      RRS$Rrs.sd =Rrs.sd
      RRS$Rrs.wl =tmp$Rrs.wl
      RRS$QWIP <- QWIP.Rrs(tmp$Rrs.wl, Rrs.mean)
      RRS$FU   <- Rrs2FU(tmp$Rrs.wl, Rrs.mean)$FU


    # Plot selected Rrs
      rrs.df <- as.data.frame(t(all.rrs))
      cast.names <- paste(as_datetime(cast.datetime, origin = "1970-01-01"), method.selected)
      names(rrs.df) <- cast.names
      Df.stat = as.data.frame(cbind(wavelength=RRS$Rrs.wl,
                               Rrs.mean =Rrs.mean,
                               Rrs.sd =Rrs.sd))

      Df.cast = as.data.frame(cbind(wavelength=RRS$Rrs.wl,rrs.df))
      Dfm.cast = melt(Df.cast, id.vars = c("wavelength"))
      names(Dfm.cast) = c("wavelength", "Methods", "value" )


      plot.rrs <- ggplot() +
        geom_line(data = Dfm.cast, aes(x = wavelength, y = value, group = Methods, color=Methods)) +
        scale_color_viridis_d() +
        geom_line(data = Df.stat, aes(x = wavelength, y = Rrs.mean), color="black")+
        geom_ribbon(data = Df.stat, aes(x = wavelength, y = Rrs.mean, ymax = Rrs.mean + Rrs.sd, ymin = Rrs.mean - Rrs.sd), fill = "grey", alpha = 0.3) +
        scale_x_continuous(limits = c(350, 800))+
        labs(x=expression(lambda), y=expression(paste(rho[w])), colour="Methods", title = "Selected Rrs casts comparison")


      # ##### generate the QWIP plot
      # QWIP coefficients
      p1 <- -8.399885e-9
      p2 <- 1.715532e-5
      p3 <- -1.301670e-2
      p4 <- 4.357838e0
      p5 <- -5.449532e2

      predicted.AVW <- 440:600
      predicted.NDI <- p1*(predicted.AVW^4) +
        p2*(predicted.AVW^3) +
        p3*(predicted.AVW^2) +
        p4*predicted.AVW   + p5

      # My line data frame
      df <- data.frame(AVW = predicted.AVW,
                       NDI = predicted.NDI,
                       NDI.minus.0.1 = predicted.NDI-0.1,
                       NDI.plus.0.1 = predicted.NDI+0.1,
                       NDI.minus.0.2 = predicted.NDI-0.2,
                       NDI.plus.0.2 = predicted.NDI+0.2)

      # Reshaping
      dfm <- melt(df,id.vars = "AVW")
      names(dfm) <- c("AVW", "Predicted", "NDI")


      # my point data frame
      df.qwip <- data.frame(AVW = AVW.selected,
                                     NDI = NDI.selected,
                                     FU  = FU.selected,
                                     Methods = cast.names)

      # New row to be added
      new_row <- data.frame(
        AVW = RRS$QWIP$AVW,
        NDI = RRS$QWIP$NDI,
        FU = RRS$FU,
        Methods = "Averaged Spectrum",
        stringsAsFactors = FALSE
      )

      # Add the new row to the existing data frame
      df.qwip <- rbind(df.qwip, new_row)

      # Define meaningful colors for the points and match them to the levels of the Methods variable
      method_colors <- c(viridis(nb.rec), "black")
      names(method_colors) <- c(cast.names, "Averaged Spectrum")



      # Plotting
      plot.QWIP <- ggplot() +
        geom_line(data = dfm, aes(x = AVW, y = NDI, color = Predicted, linetype = Predicted)) +
        geom_point(data = df.qwip, aes(x = AVW, y = NDI, fill = Methods), shape = 21, size = 2.5, color = "black") +
        geom_text(data = df.qwip, aes(x = AVW, y = NDI, label = FU), vjust = -0.5) + # Add labels
        scale_color_manual(name = "Lines",
                           labels = c("Predicted", "-0.1", "+0.1", "-0.2", "+0.2"),
                           values = c("black", "orange", "orange", "red", "red")) +
        scale_fill_manual(name = "Methods",
                          values = method_colors,
                          labels = unique(df.qwip$Methods)) +
        scale_linetype_manual(name = "Lines",
                              labels = c("Predicted", "-0.1", "+0.1", "-0.2", "+0.2"),
                              values = c("solid", "dashed", "dashed", "dotted", "dotted"))

      fullplot <- plot.rrs / plot.QWIP +
        plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5))) +
        plot_layout(guides = "collect")
      suppressMessages(plot(fullplot))


      ggsave(paste("PNG/Rrs_casts_comparison.png",sep=""), units = "in",
             width = 8, height = 7)

    } else {
      SAS = read.hocr.L2.SAS(filen, VERBOSE=FALSE)
      save(file = paste(dirdat, "SAS.raw.RData", sep="/"), SAS)

      RRS = compute.Rrs.SAS(SAS, tilt.max= cast.info$tilt.max,
                            quantile.prob=cast.info$quantile.prob,
                            windspeed=cast.info$Windspeed,
                            NIR.CORRECTION=cast.info$NIR.CORRECTION,
                            thetaV=cast.info$ThetaV,
                            Dphi=cast.info$Dphi,
                            Good=cast.info$good[i],
                            use.COMPASS)
      RRS$Rrs.mean = RRS$Rrs
      RRS$Rrs.sd =   rep(NA,length(RRS$Rrs.mean))
      RRS$Rrs.wl =   RRS$Rrs.wl
      plot.SAS.Rrs(RRS, PNG=TRUE)
    }
  }

  if (TYPE == "TRANSIT") {

    # construct a file name vector from cast info
    if (use.SAS.RData) {
      print(paste("Load: ",dirdat, "/SAS.raw.RData", sep="" ))
      load(paste(dirdat, "SAS.raw.RData", sep="/"))
      nb.rec = length(SAS)
      print(paste(nb.rec, "records along the transit...."))
    } else {
      # construct a file name vector from cast info
      filen = paste(cast.info$Year,"-", cast.info$Day,"-", cast.info$Time, "_L2.dat", sep="")
      SAS = lapply(filen, read.hocr.L2.SAS, VERBOSE=FALSE)
      save(file = paste(dirdat, "SAS.raw.RData", sep="/"), SAS)
      nb.rec = length(filen)
    }

    RRS = list()
    all.rrs = matrix(NA, nrow=nb.rec, ncol=85)
    all.lat.lon =  matrix(NA, nrow=nb.rec, ncol=2)
    Anc =  matrix(NA,nrow=nb.rec, ncol=10)
    DateTime = rep(NA, nb.rec)


      for (i in 1:nb.rec) {
        print(paste("Process data collected at:", SAS[[i]]$dd$date))
        tmp = compute.Rrs.SAS(SAS[[i]],
                              tilt.max= cast.info$tilt.max[i],
                              quantile.prob=cast.info$quantile.prob[i],
                              windspeed=cast.info$Windspeed[i],
                              NIR.CORRECTION=cast.info$NIR.CORRECTION[i],
                              thetaV=cast.info$ThetaV[i],
                              Dphi=cast.info$Dphi[i],
                              Good=cast.info$good[i],
                              use.COMPASS)
        RRS[[i]] <- tmp
        all.rrs[i,] <- tmp$Rrs
        all.lat.lon[i,] <- c(SAS[[i]]$dd$latitude, SAS[[i]]$dd$longitude)
        Anc[i,1] = tmp$ThetaS
        Anc[i,3] = tmp$ThetaV
        Anc[i,2] = tmp$PhiS
        Anc[i,4] = tmp$Dphi.log
        Anc[i,5] = tmp$Dphi.comp
        Anc[i,6] = tmp$Windspeed
        Anc[i,7] = tmp$CLEARSKY
        Anc[i,8] = tmp$rho.sky
        Anc[i,9] = tmp$offset
        DateTime[i] = tmp$DateTime
      }
    RRS$rrs <- all.rrs
    RRS$lat.lon <- all.lat.lon
    RRS$Anc <- Anc
    RRS$DateTime <- DateTime

  }

  save(file = paste(dirdat,"RRS.RData", sep="/"), RRS)

  return(RRS)

}

