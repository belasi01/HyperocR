
library(plotly)

##### convert RSR to RData

OLI=read.table("extra/OLI_RSRs.txt", skip=5)
names(OLI) <- c("waves", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8")

p <- plot_ly(x=~OLI$waves,y=~OLI$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
  add_trace(x=~OLI$waves,y=~OLI$B2, name = 'B2', color = I("blue"), fill = 'tozeroy') %>%
  add_trace(x=~OLI$waves,y=~OLI$B3, name = 'B3', color = I("green"), fill = 'tozeroy') %>%
  add_trace(x=~OLI$waves,y=~OLI$B4, name = 'B4', color = I("red"), fill = 'tozeroy') %>%
  add_trace(x=~OLI$waves,y=~OLI$B5, name = 'B5', color = I("darkred"), fill = 'tozeroy') %>%
  add_trace(x=~OLI$waves,y=~OLI$B6, name = 'B6', color = I("grey"), fill = 'tozeroy') %>%
  add_trace(x=~OLI$waves,y=~OLI$B7, name = 'B7', color = I("darkgrey"), fill = 'tozeroy') %>%
  add_trace(x=~OLI$waves,y=~OLI$B8, name = 'B8', color = I("black"), fill = 'tozeroy') %>%
  layout(xaxis = list(title = 'Wavelenght', range=c(400,2500)),
         yaxis = list(title = 'RSR'))
print(p)

save(OLI, file = "data/OLI_RSRs.RData")
