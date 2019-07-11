library(plotly)

##### convert RSR to RData

MODIST=read.table("extra/HMODIST_RSRs.txt", skip=8)
names(MODIST) <- c("waves", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8" , "B9",
                   "B10", "B11", "B12", "B13", "B14", "B15", "B16")

p <- plot_ly(x=~MODIST$waves,y=~MODIST$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B2, name = 'B2', color = I("darkblue"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B3, name = 'B3', color = I("blue"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B4, name = 'B4', color = I("cyan"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B5, name = 'B5', color = I("turquoise"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B6, name = 'B6', color = I("green"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B7, name = 'B7', color = I("darkgreen"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B8, name = 'B8', color = I("darkorange"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B9, name = 'B9', color = I("red"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B10, name = 'B10', color = I("darkred"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B11, name = 'B11', color = I("grey"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B12, name = 'B12', color = I("grey"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B13, name = 'B13', color = I("darkgrey"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B14, name = 'B14', color = I("darkgrey"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B15, name = 'B15', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~MODIST$waves,y=~MODIST$B16, name = 'B16', color = I("black"), fill = 'tozeroy') %>%
  layout(xaxis = list(title = 'Wavelenght', range=c(380,2200)),
         yaxis = list(title = 'RSR'))
print(p)

save(MODIST, file = "data/MODIST_RSRs.RData")
