library(plotly)

##### convert RSR to RData

PlanetScope_0E=read.table("extra/PlanetScope_RSR_SatID_0e.csv", skip=1, sep=",")
names(PlanetScope_0E) <- c("waves", "B1", "B2", "B3", "B4")
PlanetScope_0E$waves <- PlanetScope_0E$waves*1000

p <- plot_ly(x=~PlanetScope_0E$waves,y=~PlanetScope_0E$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("blue"), name ="B1") %>%
  add_trace(x=~PlanetScope_0E$waves,y=~PlanetScope_0E$B2, name = 'B2', color = I("green"), fill = 'tozeroy') %>%
  add_trace(x=~PlanetScope_0E$waves,y=~PlanetScope_0E$B3, name = 'B3', color = I("red"), fill = 'tozeroy') %>%
  add_trace(x=~PlanetScope_0E$waves,y=~PlanetScope_0E$B4, name = 'B4', color = I("black"), fill = 'tozeroy') %>%
  layout(xaxis = list(title = 'Wavelenght', range=c(400,900)),
         yaxis = list(title = 'RSR'))
print(p)

save(PlanetScope_0E, file = "data/PlanetScope_0E_RSRs.RData")
