library(plotly)

##### convert RSR to RData

VIIRS=read.table("extra/VIIRSN_IDPSv3_RSRs.txt", skip=5)
VIIRS<-VIIRS[,1:8]
names(VIIRS) <- c("waves", "B1", "B2", "B3", "B4", "B5", "B6", "B7")

p <- plot_ly(x=~VIIRS$waves,y=~VIIRS$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
  add_trace(x=~VIIRS$waves,y=~VIIRS$B2, name = 'B2', color = I("blue"), fill = 'tozeroy') %>%
  add_trace(x=~VIIRS$waves,y=~VIIRS$B3, name = 'B3', color = I("green"), fill = 'tozeroy') %>%
  add_trace(x=~VIIRS$waves,y=~VIIRS$B4, name = 'B4', color = I("red"), fill = 'tozeroy') %>%
  add_trace(x=~VIIRS$waves,y=~VIIRS$B5, name = 'B5', color = I("darkred"), fill = 'tozeroy') %>%
  add_trace(x=~VIIRS$waves,y=~VIIRS$B6, name = 'B6', color = I("grey"), fill = 'tozeroy') %>%
  add_trace(x=~VIIRS$waves,y=~VIIRS$B7, name = 'B7', color = I("darkgrey"), fill = 'tozeroy') %>%
  layout(xaxis = list(title = 'Wavelenght', range=c(380,1000)),
         yaxis = list(title = 'RSR'))
print(p)

save(VIIRS, file = "data/VIIRS_RSRs.RData")
