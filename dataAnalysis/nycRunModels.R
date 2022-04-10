################################################################################
# Run models by peak for NYC
# for each peak run 7 different models
#   thresh, hill, power, gp, spline, betat, basic
################################################################################

### read data
nyc <- read.csv('./Data/nycClean.csv')

peak <- c('full', '1', '2', '3', '4')
alarmFit <- c( 'thresh', 'hill', 'power', 'gp', 'spline', 'betat', 'basic')


################################################################################
### Peak 1

peak <- 2

# format data and initial values

incData <- nyc[nyc$peak == peak,]

N <- nyc$Population[1]

# use first 5 days as those initially infectious
I0 <- incData[1:5]
incData <- incData[-c(1:5)]

# use cumulative infectious before peak to determine those initially removed
if (peak %in% c('full', '1')) {
    R0 <- 0 
} else {
    idxStart <- min(which(nyc$peak == peak))
    R0 <- nyc$cumulativeCases[idxStart - 1]
}







