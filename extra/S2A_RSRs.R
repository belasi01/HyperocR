library(plotly)

##### convert RSR to RData

S2A=read.table("extra/S2A_RSRs.txt", skip=5)
names(S2A) <- c("waves", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
                "B9", "B10", "B11", "B12", "B13")

p <- plot_ly(x=~S2A$waves,y=~S2A$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
  add_trace(x=~S2A$waves,y=~S2A$B2, name = 'B2', color = I("blue"), fill = 'tozeroy') %>%
  add_trace(x=~S2A$waves,y=~S2A$B3, name = 'B3', color = I("green"), fill = 'tozeroy') %>%
  add_trace(x=~S2A$waves,y=~S2A$B4, name = 'B4', color = I("orange"), fill = 'tozeroy') %>%
  add_trace(x=~S2A$waves,y=~S2A$B5, name = 'B5', color = I("red"), fill = 'tozeroy') %>%
  add_trace(x=~S2A$waves,y=~S2A$B6, name = 'B6', color = I("darkred"), fill = 'tozeroy') %>%
  add_trace(x=~S2A$waves,y=~S2A$B7, name = 'B7', color = I("grey"), fill = 'tozeroy') %>%
  add_trace(x=~S2A$waves,y=~S2A$B8, name = 'B8', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~S2A$waves,y=~S2A$B9, name = 'B9', color = I("darkgrey"), fill = 'tozeroy') %>%
  add_trace(x=~S2A$waves,y=~S2A$B10, name = 'B10', color = I("darkgrey"), fill = 'tozeroy') %>%
  add_trace(x=~S2A$waves,y=~S2A$B11, name = 'B11', color = I("darkgrey"), fill = 'tozeroy') %>%
  add_trace(x=~S2A$waves,y=~S2A$B12, name = 'B12', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~S2A$waves,y=~S2A$B13, name = 'B13', color = I("black"), fill = 'tozeroy') %>%
  layout(xaxis = list(title = 'Wavelenght', range=c(400,2350)),
         yaxis = list(title = 'RSR'))
print(p)

save(S2A, file = "data/S2A_RSRs.RData")
