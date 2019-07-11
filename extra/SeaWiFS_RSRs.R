library(plotly)

##### convert RSR to RData

SeaWiFS=read.table("extra/SeaWiFS_RSRs.txt", skip=10)
names(SeaWiFS) <- c("waves", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8")

# Normalise each RSR

for (i in 2:9) {
  SeaWiFS[,i] <- SeaWiFS[,i]/max(SeaWiFS[,i], na.rm = T)
}

p <- plot_ly(x=~SeaWiFS$waves,y=~SeaWiFS$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
  add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B2, name = 'B2', color = I("blue"), fill = 'tozeroy') %>%
  add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B3, name = 'B3', color = I("cyan"), fill = 'tozeroy') %>%
  add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B4, name = 'B4', color = I("turquoise"), fill = 'tozeroy') %>%
  add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B5, name = 'B5', color = I("green"), fill = 'tozeroy') %>%
  add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B6, name = 'B6', color = I("red"), fill = 'tozeroy') %>%
  add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B7, name = 'B7', color = I("darkgrey"), fill = 'tozeroy') %>%
  add_trace(x=~SeaWiFS$waves,y=~SeaWiFS$B8, name = 'B8', color = I("black"), fill = 'tozeroy') %>%
  layout(xaxis = list(title = 'Wavelenght', range=c(380,950)),
         yaxis = list(title = 'RSR'))
print(p)

save(SeaWiFS, file = "data/SeaWiFS_RSRs.RData")
