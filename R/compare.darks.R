

compare.darks <- function(dirdat) {

  filen = list.files(pattern = "*DARK*")

  nb.file = length(filen)

  Ed.dark.mean = matrix(NA, ncol=nb.file, nrow=137)
  Li.dark.mean = matrix(NA, ncol=nb.file, nrow=137)
  Lt.dark.mean = matrix(NA, ncol=nb.file, nrow=137)

  for (i in 1:nb.file) {
    SAS.dark = read.hocr.L2.SAS(filen[i], VERBOSE=FALSE)
    Ed.dark.mean[,i] = apply(SAS.dark$Ed, 2, mean)
    Li.dark.mean[,i] = apply(SAS.dark$Li, 2, mean)
    Lt.dark.mean[,i] = apply(SAS.dark$Lt, 2, mean)
  }

  # plot the data

  png(file=paste(dirdat, "Darks.comparison.png", sep="/"),
      res=300, height = 5, width = 6, units = "in")

  par(mfrow=c(3,1))
  par(mar=c(0,5,1,1))
  matplot(SAS.dark$Ed.wl, Ed.dark.mean, pch=19,
          col=rainbow.modified(nb.file), xlab="", xaxt = "n")
  grid(col="darkgrey")
  par(mar=c(1,5,1,1))
  matplot(SAS.dark$Li.wl, Li.dark.mean, pch=19,
          col=rainbow.modified(nb.file), xlab="", xaxt = "n")
  grid(col="darkgrey")
  par(mar=c(2,5,0,1))
  matplot(SAS.dark$Li.wl, Li.dark.mean, pch=19,
          col=rainbow.modified(nb.file), xlab="Wavelength")
  grid(col="darkgrey")
  dev.off()

  print(paste("See PNG file at :", dirdat, "/Darks.comparison.png", sep=""))


}
