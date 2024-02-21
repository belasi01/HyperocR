
plot.SAS.Rrs <- function (SAS, PNG=FALSE, RADIANCES=FALSE) {

  ix.wl = which(SAS$Rrs.wl > 350 & SAS$Rrs.wl <900)
  if (PNG & !dir.exists("PNG")) dir.create("PNG")

  if (RADIANCES) {

      # first plot of the raw radiances with error bars
      Df = as.data.frame(cbind(wavelength=SAS$Rrs.wl[ix.wl],
                               Li=SAS$Li.mean[ix.wl], Li.sd=SAS$Li.sd[ix.wl],
                               Lt=SAS$Lt.mean[ix.wl], Lt.sd=SAS$Lt.sd[ix.wl],
                               Ed=SAS$Ed.mean[ix.wl], Ed.sd=SAS$Ed.sd[ix.wl],
                               rhosky=pi*SAS$rho.sky*(SAS$Li.mean[ix.wl]/SAS$Ed.mean[ix.wl]),
                               rhosurf=pi*SAS$Lt.mean[ix.wl]/SAS$Ed.mean[ix.wl]))

      p1=ggplot(data=Df, aes(x=wavelength, y=Ed))  + geom_line()  + scale_x_continuous(limits = c(350, 800))  + geom_ribbon(aes(ymin=Ed-Ed.sd, ymax=Ed+Ed.sd, x=wavelength), alpha = 0.5)
      p2=ggplot(data=Df, aes(x=wavelength, y=Li))  + geom_line()  + scale_x_continuous(limits = c(350, 800))+ geom_ribbon(aes(ymin=Li-Li.sd, ymax=Li+Li.sd, x=wavelength), alpha = 0.5)
      p3=ggplot(data=Df, aes(x=wavelength, y=Lt))  + geom_line()  + scale_x_continuous(limits = c(350, 800))+ geom_ribbon(aes(ymin=Lt-Lt.sd, ymax=Lt+Lt.sd, x=wavelength), alpha = 0.5)
      p4=ggplot(data=Df, aes(x=wavelength, y=rhosky))  + geom_line()
      p4 = p4  + geom_line(aes(x=wavelength, y=rhosurf), linetype=2) + labs(x=expression(lambda),
                                                                            y=expression(paste(rho[sky],rho[surf])))

      p1 / p2 / p3 / p4

      if (PNG) ggsave(paste("PNG/Radiances_",SAS$DateTime,".png",sep=""), units = "in",
             width = 5, height = 4)

  } else {
      Df = as.data.frame(cbind(wavelength=SAS$Rrs.wl[ix.wl],
                               None=SAS$Rrs[ix.wl],
                               Null_900=SAS$Rrs.NULL[ix.wl],
                               Similarity_720_780=SAS$Rrs.SIMILARITY[ix.wl],
                               NIR = SAS$Rrs.NIR[ix.wl],
                               UV = SAS$Rrs.UV[ix.wl],
                               UV.NIR = SAS$Rrs.UV.NIR[ix.wl],
                               Kutser13 = SAS$Rrs.Kutser[ix.wl]))

    Dfm = melt(Df, id.vars = c("wavelength"))
    names(Dfm) = c("wavelength", "rho_w", "value" )

      p1 <- ggplot(data=Dfm, aes(x=wavelength, y=value, colour=rho_w)) + geom_line()
      p1 <- p1 + scale_x_continuous(limits = c(350, 800))
      p1 <- p1 + labs(x=expression(lambda), y=expression(paste(rho[w])), colour="Correction method")
      p1 <- p1 + ggtitle(paste(SAS$anc$StationID, SAS$DateTime, "Lat:", signif(SAS$dd$lat,5), "Lon:", signif(SAS$dd$lon,5)),
                         subtitle = bquote(rho[Fresnel]^Mobley2015 == .(signif(SAS$rho.sky,3)) ~
                                             "   "~ rho[Fresnel]^NIR == .(signif(SAS$rho.sky.NIR, 3))~
                                             "   "~ rho[Fresnel]^UV == .(signif(SAS$rho.sky.UV, 3))))
      print(p1)
      if (PNG) ggsave(paste("PNG/Rrs_",SAS$DateTime,".png",sep=""), units = "in",
                      width = 5, height = 4)

  }
}
