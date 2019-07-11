library(plotly)

##### convert RSR to RData

MERIS=read.table("extra/MERIS_RSRs_avg.txt", skip=5)
names(MERIS) <- c("waves", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9",
                "B10", "B11", "B12", "B13", "B14", "B15")
MERIS[MERIS == -999] <- 0

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

  layout(xaxis = list(title = 'Wavelenght', range=c(400,920)),
         yaxis = list(title = 'RSR'))
print(p)

save(MERIS, file = "data/MERIS_RSRs.RData")
