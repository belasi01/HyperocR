#' @title Generate a data base for rhow (pi*Rrs) from the HyperSAS
#'
#'
#' @param path is the path where the file directories.for.HyperSAS.dat containing the data folders
#' to merge in the data base.
#' @param mission is a string for the name of the mission. It will be used for the file names of the output.
#' @param wave.range is a vector of two integers indicating the wavelength range to output in the database.
#' The default is wave.range=c(350,810) to put the data from 350 nm to 810 nm.
#'
#' @return It returns a list object named SAS.DB containing a matrix of rhow (rhow.m) and vectors for
#' StationID, data, lat, lon sunzen, windspeed, rhow.Method.
#'
#' The object SAS.DB is saved in RData format. The data are also saved in ASCII (.dat with comma separator)
#' and a figure showing the measured rho_w spectra of the data base is produced.
#'
#' @author Simon Bélanger
#' @export
#' @name generate.SAS.rhow.DB

generate.SAS.rhow.DB <- function(path="./",mission="XXX", wave.range=c(350,810)) {

  setwd(path)
  path<- getwd()
  dirdats <- read.table("directories.for.HyperSAS.dat", sep=" ",
                        comment.char = "#", colClasses = "character")
  dirs <- dirdats$V1


  ndirs = length(dirs)

  # Read one file to get the waves configuration of the SAS
  setwd(as.character(dirs[1]))
  file.name = "RRS.RData"
  load(file.name)

  # find the indices corresponding to the wavelength range to output
  ix.min <- which.min(abs(RRS$Rrs.wl -wave.range[1]))
  ix.max <- which.min(abs(RRS$Rrs.wl -wave.range[2]))
  waves = RRS$Rrs.wl[ix.min:ix.max]
  nwaves = length(waves)

  #####
  setwd(path)
  print(paste("The number of station found is : ", ndirs))

  rhow.m = matrix(ncol=nwaves, nrow = ndirs)
  ID = rep("ID", ndirs)
 # Replicate = rep("A", ncasts)
  date = rep(NA, ndirs)
  sunzen = rep(NA, ndirs)
  viewzen = rep(NA, ndirs)
  dphi = rep(NA, ndirs)
  lat = rep(NA, ndirs)
  lon = rep(NA, ndirs)
  windspeed = rep(NA, ndirs)
  rhow.Method = rep(NA, ndirs)
  cast=1
  for (i in 1:ndirs) {
      setwd(as.character(dirs[i]))
      cast.info <- read.table(file="cast.info.dat", header=T, comment.char = "#")
      file.name = "RRS.RData"
      load(file.name)

      ix.good = which(cast.info$good == 1)
      rhow.Method[cast] <- paste(cast.info$NIR.CORRECTION[ix.good], collapse = " ") # retrieve the method from the first cast
      rhow.m[cast,] <- RRS$Rrs.mean[ix.min:ix.max] * pi # to convert into rhow
      ID[cast] <- RRS$ID
      date[cast] <- RRS$DateTime
      sunzen[cast] <- RRS$sunzen
      viewzen[cast]<- RRS$viewzen
      dphi[cast]<- RRS$dphi
      lat[cast] <- RRS$lat
      lon[cast] <- RRS$lon
      windspeed[cast] <- RRS$windspeed


      rec.info = data.frame(ID[cast],
                            date[cast],
                            lat[cast],
                            lon[cast],
                            sunzen[cast],
                            viewzen[cast],
                            dphi[cast],
                            windspeed[cast],
                            rhow.Method[cast])


      if (cast == 1) {
        all = data.frame(rec.info,t(rhow.m[cast,]))

        col.names = c(paste("rhow_", waves,sep=""))

        names(all) <- c("StationID","DateTime",  "latitude", "longitude", "sunzen", "viewzen", "dphi",  "WindSpeed", "rhow.Method", col.names)
      } else
      {
        rec = data.frame(rec.info,t(rhow.m[cast,]))
        names(rec) <-  c("StationID","DateTime",  "latitude", "longitude", "sunzen",  "viewzen", "dphi", "WindSpeed", "rhow.Method", col.names)
        all = rbind(all,rec)
      }


      cast = cast + 1




    }



  SAS.BD <- list(ID=ID,
                 #Replicate=Replicate,
                 waves=waves,
                 rhow.m=rhow.m,
                 date=as.POSIXct(date, origin="1970-01-01"),
                 lat=lat,
                 lon=lon,
                 sunzen=sunzen,
                 viewzen=viewzen,
                 dphi=dphi,
                 windspeed=windspeed,
                 rhow.Method=rhow.Method)

  # Save the data
  setwd(path)
  save(SAS.BD, file=paste("SAS.DB.PackageVersion.",packageVersion("HyperocR"),".", mission,".RDATA",sep=""))
  all$DateTime=as.POSIXct(all$DateTime, origin="1970-01-01")
  write.table(all, file = paste("SAS.DB.PackageVersion.",packageVersion("HyperocR"),".", mission,".dat",sep=""), sep=",", quote=F, row.names=F)



  # plot the data

  png(paste("SAS.DB.PackageVersion.",packageVersion("HyperocR"),".",  mission,".png",sep=""), res=300, height = 6, width = 8, units = "in")

  Df = as.data.frame(cbind(wavelength=waves, t(rhow.m)))
  colnames(Df) <- c("wavelength", ID)
  Dfm = reshape2::melt(Df, id.vars = c("wavelength"))
  names(Dfm) = c("wavelength", "rho_w", "value" )

  p1 <- ggplot(data=Dfm, aes(x=wavelength, y=value, colour=rho_w)) + geom_line()
  p1 <- p1 + labs(x=expression(lambda), y=expression(paste(rho[w])), colour="Station")
  p1 <- p1 + ggtitle(paste(mission))
  print(p1)

  dev.off()

  return(SAS.BD)


}
