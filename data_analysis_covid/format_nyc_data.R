################################################################################
# NYC data formatting
################################################################################


library(lubridate)

### Read in data from GitHub
nyc <-read.csv("https://raw.githubusercontent.com/nychealth/coronavirus-data/master/trends/cases-by-day.csv") 

nyc <- nyc[,c('date_of_interest', 'CASE_COUNT')]
colnames(nyc) <- c('date', 'dailyCases')

# format dates
nyc$date <- as.Date(nyc$date, format = '%m/%d/%Y')

# population
nyc$Population <- 8804190


# peak identifier
nyc$peak <- NA

# peak 1 - Feb 29 - Jun 1, 2020
startDate <- as.Date('2020-02-29')
endDate <- as.Date('2020-06-16')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 1

# peak 2 - Oct 01, 2020 - Jun 1, 2021
startDate <- as.Date('2020-10-01')
endDate <- as.Date('2021-05-16')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 2

# peak 3 - Jul 1 - Nov 1, 2021
startDate <- as.Date('2021-07-01')
endDate <- as.Date('2021-11-01')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 3

# peak 4 - Dec 1, 2021 - Feb 15, 2022
startDate <- as.Date('2021-12-01')
endDate <- as.Date('2022-02-16')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 4


write.csv(nyc, './data/nycClean.csv', quote = F, row.names = F)

