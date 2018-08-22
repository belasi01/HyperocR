#'
#' Run batch processing for HyperOCR
#'
#'  @param log.file is the name of the ASCII file containing the
#' list of L2 data to processed.
#' @param data.path is the full path where the L2 data are stored (*.dat)
#'
#' @author Simon Bélanger
#' @export



run.process.HOCR.batch <- function(log.file="log.txt",
                                   data.path="./",
                                   Ag.path="./",
                                   Ap.path="./",
                                   Bbp.path="./") {
  if (file.exists(data.path)) {
    print("Data path exists")
    print(data.path)
  } else {
    print("data.path does not exist!")
    print("STOP processing")
    return(0)
  }

  # Check the output directories and Create them if they are not available.
  path.png = file.path(data.path,"png/") #paste(data.path,"/png/", sep="")
  path.out = file.path(data.path,"RData/") #paste(data.path,"/RData/", sep="")
  path.raw = file.path(data.path,"dat/")
  if (!file.exists(path.raw)) {
    print("No folder named raw found")
    print("STOP processing")
    print("Create a folder:")
    print(path.raw)
    print("and put your data files in it")
    return(0)
  } else print(paste("Raw data in :",path.raw))

  if (!file.exists(path.png)) dir.create(path.png)
  if (!file.exists(path.out)) dir.create(path.out)

  # Lecture des informations dans un fichier texte
  ext = unlist(strsplit(log.file, "[.]"))[2]
  if (ext == "csv") {
    log = read.table(file=log.file, header=T, sep=",")
  } else {
    log = read.table(file=log.file, header=T)
  }

  nfiles = length(log$HOCR.filename)
  print(paste(nfiles, "files to process."))
  for (i in 1:nfiles) {
    #print(paste("Processing ", log$HOCR.filename[i]))

    raw.file <- paste(path.raw,log$HOCR.filename[i], sep="/")
    Ag.file <- paste(Ag.path,log$Ag.file[i], sep="/")
    if (!file.exists(Ag.file)) {
      print("WARNING no CDOM file found")
      Ag.file = NA
    }
    Ap.file <- paste(Ap.path,log$Ap.file[i], sep="/")
    if (!file.exists(Ap.file)) {
      print("WARNING no Ap file found")
      Ap.file = NA
    }

    Bbp.file <- paste(Bbp.path,log$bbp.file[i], sep="/")
    if (!file.exists(Bbp.file)) {
      print("WARNING no Bbp file found")
      Bbp.file = NA
    }
    if (file.exists(raw.file)) {
      print(paste("Processing ", raw.file))
      HOCR <- process.HOCR(filen=raw.file,
                           Station=log$Station[i],
                           Lu0.Depth = log$Lu0.Depth[i],
                           EdZ.Depth = log$Ed.Depth[i],
                           Delta.LuZ.Depth = log$Delta.Depth[i],
                           Ag.file =  Ag.file,
                           Ap.file =  Ap.file,
                           Bbp.file = Bbp.file)

      # Plot data
      time <- paste(unlist(str_split(log$UTCTime[i],":"))[1],
                    unlist(str_split(log$UTCTime[i],":"))[2],sep="")
      png(file=paste(path.png,log$Station[i],"_",log$Date[i],"_",time,".png", sep=""),
          width=6, height=6, units="in", res=300)

      par(mfrow=c(2,1))
      par(mar=c(4,5,1,1))

      plot(HOCR$waves, HOCR$Rrs,
           type="l", lwd=2,
           main=paste(HOCR$Station, HOCR$DateTime),
           ylab=expression(R[rs](sr^-1)),
           xlab=expression(lambda))
      lines(HOCR$waves,HOCR$uRrs, col=2, lwd=1)
      legend("topleft", c("uncorrected", "corrected"), lwd=c(1,2), col=c(2,1))

      if (!is.na(HOCR$K.CORRECTION)) {
        if (HOCR$K.CORRECTION == "Kd") {
          plot(HOCR$waves, HOCR$Kd,
               type="l", lwd=2,
               ylab=expression(K[d](m^-1)),
               xlab=expression(lambda))
        }
        if (HOCR$K.CORRECTION == "KLu") {
          plot(HOCR$waves, HOCR$KLu,
               type="l", lwd=2,
               ylab=expression(K[Lu](m^-1)),
               xlab=expression(lambda))
        }
        if (HOCR$K.CORRECTION == "Kd.modeled") {
          plot(HOCR$waves, HOCR$Kd.mod,
               type="l", lwd=2,
               ylab=expression(paste("Modeled",K[d](m^-1))),
               xlab=expression(lambda))
        }
      }

      dev.off()
    } else {
      print("WARNING no data file")
      print(paste("Skipping", raw.file))
    }
  }

}
