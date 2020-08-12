library(ncdf4)
library(plotly)
nc_data <- nc_open('extra/S3A_OL_SRF_20160713_mean_rsr.nc4')


#sink('extra/S3A_OL_SRF_0_20180109_mean_rsr.txt')
#print(nc_data)
#sink()

SRF=ncvar_get(nc_data, "mean_spectral_response_function")
SRF.waves=ncvar_get(nc_data, "mean_spectral_response_function_wavelength")

# create a continous matrix similat to S3A and S2B
waves <- seq(387,1044,0.1)
nwaves<- length(waves)

S3A <- as.data.frame(matrix(0, nrow = nwaves, ncol=21))

S3A$waves <- waves

names(S3A) <- c("waves", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
                "B9", "B10", "B11", "B12", "B13", "B14", "B15", "B16", "B17",
                "B18", "B19", "B20", "B21")

for (i in 1:21) {
  tmp <- spline(SRF.waves[,i], SRF[,i] , xout=waves)$y
  tmp[tmp<0] = 0
  S3A[,i+1] <- tmp
}

S3A$waves <- waves
p <- plot_ly(x=~S3A$waves,y=~S3A$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
  add_trace(x=~S3A$waves,y=~S3A$B2, name = 'B2', color = I("darkviolet"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B3, name = 'B3', color = I("blue"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B4, name = 'B4', color = I("cyan"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B5, name = 'B5', color = I("green"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B6, name = 'B6', color = I("darkgreen"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B7, name = 'B7', color = I("orange"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B8, name = 'B8', color = I("red"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B9, name = 'B9', color = I("darkred"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B10, name = 'B10', color = I("darkred"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B11, name = 'B11', color = I("darkred"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B12, name = 'B12', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B13, name = 'B13', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B14, name = 'B14', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B15, name = 'B15', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B16, name = 'B16', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B17, name = 'B17', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B18, name = 'B18', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B19, name = 'B19', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B20, name = 'B20', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~S3A$waves,y=~S3A$B21, name = 'B21', color = I("black"), fill = 'tozeroy') %>%
  layout(xaxis = list(title = 'Wavelenght', range=c(380,1050)),
         yaxis = list(title = 'RSR'))
print(p)

save(S3A, file = "data/S3A_RSRs.RData")




