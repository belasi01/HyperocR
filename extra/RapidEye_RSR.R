library(plotly)

##### convert RSR to RData

RapidEye=read.table("extra/RapidEye_RSR.csv", skip=1, sep=",")
names(RapidEye) <- c("waves", "B1", "B2", "B3", "B4", "B5")
RapidEye$waves <- RapidEye$waves*1000

p <- plot_ly(x=~RapidEye$waves,y=~RapidEye$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("blue"), name ="B1") %>%
  add_trace(x=~RapidEye$waves,y=~RapidEye$B2, name = 'B2', color = I("green"), fill = 'tozeroy') %>%
  add_trace(x=~RapidEye$waves,y=~RapidEye$B3, name = 'B3', color = I("red"), fill = 'tozeroy') %>%
  add_trace(x=~RapidEye$waves,y=~RapidEye$B4, name = 'B4', color = I("darkred"), fill = 'tozeroy') %>%
  add_trace(x=~RapidEye$waves,y=~RapidEye$B5, name = 'B5', color = I("black"), fill = 'tozeroy') %>%
  layout(xaxis = list(title = 'Wavelenght', range=c(400,900)),
         yaxis = list(title = 'RSR'))
print(p)

save(RapidEye, file = "data/RapidEye_RSRs.RData")
