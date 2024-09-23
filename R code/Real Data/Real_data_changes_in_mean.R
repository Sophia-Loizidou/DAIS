library(dplyr)
library(lubridate)
library(ggplot2)

## Access data from: https://catalog.data.gov/dataset/crime

## Read data 
data <- read.csv('path_to_data')

## Process the data
split_date <- function(data){ #remove time from date
  stringr::str_split(data, ' ')[[1]][1]
}
data$date <- sapply(data$Start_Date_Time, split_date)
crime <- data %>% count(date)
crime <- crime %>% arrange(mdy(crime$date))

## Only use data starting from 01/22/2020 until 08/31/2024
crime_small <- crime[1301:2984,]

## Transform data using Anscombe transformation as they are integer valued
cases_transformed <- 2*sqrt(crime_small$n + 3/8)
cpt <- DAIS(cases_transformed, contrast = 'mean')

## Define format of date
crime_small$date <- as.Date(crime_small$date, format = "%m/%d/%Y")

ggplot(crime_small, aes(x=date, y=n)) +
  geom_line() + 
  xlab("Time") +
  ylab("Number of Crimes reported") +
  geom_line(aes(date, IDetect:::est_signal(n, cpt, type = "mean")), col='red', lwd=1) +
  theme_classic()
