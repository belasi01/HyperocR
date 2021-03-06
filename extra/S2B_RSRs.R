library(plotly)

##### convert RSR to RData

S2B=read.table("extra/S2B_RSRs.txt", skip=5)
names(S2B) <- c("waves", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
                "B9", "B10", "B11", "B12", "B13")

p <- plot_ly(x=~S2B$waves,y=~S2B$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
  add_trace(x=~S2B$waves,y=~S2B$B2, name = 'B2', color = I("blue"), fill = 'tozeroy') %>%
  add_trace(x=~S2B$waves,y=~S2B$B3, name = 'B3', color = I("green"), fill = 'tozeroy') %>%
  add_trace(x=~S2B$waves,y=~S2B$B4, name = 'B4', color = I("orange"), fill = 'tozeroy') %>%
  add_trace(x=~S2B$waves,y=~S2B$B5, name = 'B5', color = I("red"), fill = 'tozeroy') %>%
  add_trace(x=~S2B$waves,y=~S2B$B6, name = 'B6', color = I("darkred"), fill = 'tozeroy') %>%
  add_trace(x=~S2B$waves,y=~S2B$B7, name = 'B7', color = I("grey"), fill = 'tozeroy') %>%
  add_trace(x=~S2B$waves,y=~S2B$B8, name = 'B8', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~S2B$waves,y=~S2B$B9, name = 'B9', color = I("darkgrey"), fill = 'tozeroy') %>%
  add_trace(x=~S2B$waves,y=~S2B$B10, name = 'B10', color = I("darkgrey"), fill = 'tozeroy') %>%
  add_trace(x=~S2B$waves,y=~S2B$B11, name = 'B11', color = I("darkgrey"), fill = 'tozeroy') %>%
  add_trace(x=~S2B$waves,y=~S2B$B12, name = 'B12', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~S2B$waves,y=~S2B$B13, name = 'B13', color = I("black"), fill = 'tozeroy') %>%
  layout(xaxis = list(title = 'Wavelenght', range=c(400,2350)),
         yaxis = list(title = 'RSR'))
print(p)

save(S2B, file = "data/S2B_RSRs.RData")
