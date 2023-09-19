library(lubridate)
library(ggplot2)
library(IDetect)

## Access data from https://finance.yahoo.com/quote/EURGBP%3DX/history?period1=915408000&period2=1695081600&interval=1wk&filter=history&frequency=1wk&includeAdjustedClose=true
data <- read.csv('path_to_data')

data <- read.csv("/Users/sophia.loizidou/Library/Mobile Documents/com~apple~CloudDocs/Sophia/Research Project/New code/Two sided only/Real data slope/EURGBP=X August.csv")

## Keep date and close prices
data <- data[,c(1,5)]
cpt <- DAIS(data$Close, contrast = 'slope')

## Define date format
data$Date <- as.Date(data$Date, format = "%Y-%m-%d")

ggplot(data, aes(x=Date, y=Close)) +
  geom_line() +
  xlab("Time") +
  ylab("Number of Crimes reported") +
  geom_line(aes(Date, IDetect:::est_signal(Close, cpt, type = "slope"), colour = "DAIS"), lwd=1) +
  theme_classic()

