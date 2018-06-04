#'
#'@title HyperSAS data processing launcher
#'
#'@description This function launches the HyperSAS data processing,
#'which constists in calculating the hyperspectral remote sensing reflectance (Rrs)
#'from radiometric measurements performed using a Satlantic HyperSAS system. It processes
#'data files found in each directories specify in the file named directories.for.HyperSAS.dat.
#'The TYPE "STATION" or "TRANSIT" need to be specify for each directory in the file
#'directories.for.HyperSAS.dat. It depends whether the directory gather replicates of Rrs
#'of the same station, or the data collected along the course of a ship transit.
#'
#'For each directory, the logging information must be stored in an
#'ASCII file named cast.info.dat. Note that the HyperSAS file were previously
#'processed from raw to L2 using ProSoft and exported to ASCII format (*L2.dat).
#'
#'
#'@param report is a logical parameter indicating whether or not a
#'PDF report is produced using knitr.  Default is FALSE.
#'@param use.SAS.RData is a logical parameter indicating whether or not a SAS.raw.RData file
#'already exists in the data directory and do not need to be generated again. Note that
#'SAS.raw.RData is an output of the function \code{\link{process.HyperSAS}}.  Default is FALSE.
#'@param use.COMPASS is a logical parameter indicating whether or not the azimuth difference
#'between the sun and the view geometry is calculated using the SAS compass data. Default is FALSE
#'because compass usually don't work on ship. When FALSE, the azimuth angle difference
#'is taken from the cast.info.dat file.
#'@param  DIAGNOSTIC.PLOTS is a logical parameter indicating whether or not diagnostic plots are saved in PNG format.
#'It is used only when TYPE="TRANSIT".
#'The plots shows the temporal variability in Ed and tilt,
#'as well as the spectra for each Ed, Li and Lt recorded during the data acquisition.
#'Default is FALSE.
#'
#'@details First when HyperSAS.go is executed, it reads a file named directories.for.HyperSAS.dat
#'in the working directory from which \code{\link{HyperSAS.go}} was launched.
#'
#'Second, in each folder found in directories.for.HyperSAS.dat, the programm will look for a
#'file named cast.info.dat. This file contains the logging information need to process each *L2.dat
#'The cast.info.dat file contains the information on viewing geometry, windspeed,
#'dark current file, as well as processing parameters such as the quantile probility,
#'the maximum tilt tolerance, the "white correction method" to eliminate sun glint, foam,
#'ocean spray etc. In details, the following 11 fields should be found in cast.info.dat:
#'Year,
#'Day,
#'Time,
#'Dphi,
#'ThetaV,
#'Windspeed,
#'Wind.units,
#'quantile.prob,
#'tilt.max,
#'NIR.CORRECTION,
#'good.
#'
#'Finally, the function \code{\link{process.HyperSAS}} will be called to process each data folder.
#'
#'
#'@seealso See \code{\link{process.HyperSAS}}  and \code{\link{compute.Rrs.SAS}} for more details
#'about the processing parameters of the cast.info.dat file listed above.
#'
#'@author Simon Bélanger
#'
#'@export

HyperSAS.go <- function(report=FALSE,
                        use.SAS.RData=FALSE,
                        use.COMPASS=FALSE,
                        DIAGNOSTIC.PLOTS=FALSE) {

  #data(rho550)

  if(!file.exists("directories.for.HyperSAS.dat")) {
    cat("CREATE a file named directories.for.HyperSAS.dat in current directory (where R is launched)\n")
    cat("  and put in it the names of the directories where data files can be found (one by line)\n")
    stop()
  } else {
    dirdats <- read.table("directories.for.HyperSAS.dat", sep=" ",
                          comment.char = "#", colClasses = "character")
    for(i in 1:length(dirdats$V1)) {
      dirdat  = dirdats$V1[i]
      if(!file.exists(dirdat)) {
        cat(dirdat, "does not exist")
        stop()
      } else setwd(dirdat)

      TYPE = dirdats$V2[i]
      if (TYPE != "STATION" & TYPE != "TRANSIT"){
        cat("The TYPE:", TYPE, " does not exist.  ")
        cat("Add either STATION or TRANSIT after each directory\n")
        cat("path in directories.for.HyperSAS.dat\n")
        stop()
      }

      mymessage(paste("PROCESSING DIRECTORY", dirdat), head = "@", tail = "@")
      Rrs = process.HyperSAS(dirdat, TYPE, use.SAS.RData, use.COMPASS)
      if (report) {
        if (is.list(Rrs)) {
            create.HyperSAS.report(dirdat, TYPE, DIAGNOSTIC.PLOTS)
        } else print("No Rrs to report!")

      } else print("No report requested; Add report=TRUE to generate one")

    }
  }
}

