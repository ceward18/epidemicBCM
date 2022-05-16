################################################################################
# NYC data formatting
################################################################################


library(lubridate)

# calculate moving average for smoothing
movingAverage <- function(x, bw) {
    
    n <- length(x)
    bw <- floor(bw)
    
    out <- rep(0, n)
    for (i in 1:n) {
        
        if (i < bw) {
            t1 = 1
            t2 = i
        } else {
            t1 = i - bw + 1
            t2 = i
        }
        
        out[i] <- mean(x[t1:t2])
    }
    
    return(out)
}

### Read in data from GitHub
nyc <-read.csv("https://raw.githubusercontent.com/nychealth/coronavirus-data/master/trends/cases-by-day.csv") 

nyc <- nyc[,c('date_of_interest', 'CASE_COUNT')]
colnames(nyc) <- c('date', 'dailyCases')

# format dates
nyc$date <- as.Date(nyc$date, format = '%m/%d/%Y')

# 7-day moving average of cases to account for reporting delays
nyc$smoothedCases <- round(movingAverage(nyc$dailyCases, 7))

# cumulative cases
nyc$cumulativeCases <- cumsum(nyc$smoothedCases)

# population
nyc$Population <- 8.419*1e6


# peak identifier
nyc$peak <- NA

# peak 1 - Feb 29 - Jun 1, 2020
startDate <- as.Date('2020-02-29')
endDate <- as.Date('2020-06-15')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 1

# peak 2 - Oct 15, 2020 - Jun 1, 2021
startDate <- as.Date('2020-10-01')
endDate <- as.Date('2021-06-01')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 2

# peak 3 - Jul 1 - Nov 1, 2021
startDate <- as.Date('2021-07-01')
endDate <- as.Date('2021-11-01')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 3

# peak 4 - Nov 1, 2021 - Mar 15, 2022
startDate <- as.Date('2021-11-01')
endDate <- as.Date('2022-03-15')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 4

#plot(nyc$date, nyc$smoothedCases, col = rainbow(4)[nyc$peak], pch = 16)


write.csv(nyc, './Data/nycClean.csv', quote = F, row.names = F)




