library(plotly)

##### convert RSR to RData

MODISA=read.table("extra/HMODISA_RSRs.txt", skip=8)
names(MODISA) <- c("waves", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8" , "B9",
                   "B10", "B11", "B12", "B13", "B14", "B15", "B16")

p <- plot_ly(x=~MODISA$waves,y=~MODISA$B1, type = 'scatter', mode = 'lines', fill = 'tozeroy', color = I("violet"), name ="B1") %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B2, name = 'B2', color = I("darkblue"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B3, name = 'B3', color = I("blue"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B4, name = 'B4', color = I("cyan"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B5, name = 'B5', color = I("turquoise"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B6, name = 'B6', color = I("green"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B7, name = 'B7', color = I("darkgreen"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B8, name = 'B8', color = I("darkorange"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B9, name = 'B9', color = I("red"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B10, name = 'B10', color = I("darkred"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B11, name = 'B11', color = I("grey"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B12, name = 'B12', color = I("grey"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B13, name = 'B13', color = I("darkgrey"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B14, name = 'B14', color = I("darkgrey"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B15, name = 'B15', color = I("black"), fill = 'tozeroy') %>%
  add_trace(x=~MODISA$waves,y=~MODISA$B16, name = 'B16', color = I("black"), fill = 'tozeroy') %>%
  layout(xaxis = list(title = 'Wavelenght', range=c(380,2200)),
         yaxis = list(title = 'RSR'))
print(p)

save(MODISA, file = "data/MODISA_RSRs.RData")
