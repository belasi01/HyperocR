
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
      rrs.df <- as.data.frame(t(SAS$Rrs[,ix.wl]))
      names(rrs.df) <- SAS$methods
      Df = as.data.frame(cbind(wavelength=SAS$Rrs.wl[ix.wl],rrs.df))

    Dfm = melt(Df, id.vars = c("wavelength"))
    names(Dfm) = c("wavelength", "Methods", "value" )


    # Define meaningful colors for the points and match them to the levels of the Methods variable
    method_colors <- c("grey", "black", "orange", "darkred", "violet", "green", "yellow", "blue")
    names(method_colors) <- unique(SAS$methods)

      plot.rrs <- ggplot(data=Dfm, aes(x=wavelength, y=value, colour=Methods)) +
        geom_line() + scale_x_continuous(limits = c(350, 800)) +
        labs(x=expression(lambda), y=expression(paste(rho[w])), colour="Methods") +
        scale_color_manual(name = "Methods",
                          values = method_colors,
                          labels = unique(rrs.df$Methods)) +
        ggtitle(paste( "Lat:", signif(SAS$dd$lat,5), "Lon:", signif(SAS$dd$lon,6)),
                         subtitle = bquote(rho[Fresnel]^Mobley2015 == .(signif(SAS$rho.sky,3)) ~
                                             "   "~ rho[Fresnel]^NIR == .(signif(SAS$rho.sky.NIR, 3))~
                                             "   "~ rho[Fresnel]^UV == .(signif(SAS$rho.sky.UV, 3))))
      #print(p1)

      ##### generate the QWIP plot
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
      df.qwip <- data.frame(AVW = SAS$AVW,
                           NDI = SAS$NDI,
                           Methods = SAS$methods,
                           FU = SAS$FU)



      # Plotting
      plot.QWIP <- ggplot() +
        geom_line(data = dfm, aes(x = AVW, y = NDI, color = Predicted, linetype = Predicted)) +
        geom_point(data = df.qwip, aes(x = AVW, y = NDI, fill = Methods), shape = 21, size = 3, color = "black") +
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
        plot_annotation(title = SAS$DateTime,
                        theme = theme(plot.title = element_text(hjust = 0.5))) +
        plot_layout(guides = "collect")
      suppressMessages(plot(fullplot))


      if (PNG) ggsave(paste("PNG/Rrs_",SAS$DateTime,".png",sep=""), units = "in",
                      width = 8, height = 7)

  }
}
