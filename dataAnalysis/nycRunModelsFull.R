################################################################################
# Run models by peak for NYC
# for each peak run 7 different models
#   thresh, hill, power, gp, spline, betat, basic
################################################################################

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
idx <- as.numeric(task_id)

### load libraries
library(parallel)

### read data
nyc <- read.csv('./Data/nycClean.csv')

peak <- c('full')
alarmFit <- c( 'thresh', 'hill', 'power', 'gp', 'spline', 'betat', 'basic')
infPeriod <- 'fixed'

# 7 possibilities (7 alarmFits each to full data)
allModels <- expand.grid(peak = peak,
                         alarmFit = alarmFit,
                         infPeriod = infPeriod)

# constants for all models
N <- nyc$Population[1]
lengthI <- 7
smoothWindow <- 30

peak_i <- allModels$peak[idx]
alarmFit_i <- allModels$alarmFit[idx]
infPeriod_i <- allModels$infPeriod[idx]

print(paste('Running alarm:', alarmFit_i,
            ', peak:', peak_i, 
            ', infPeriod:', infPeriod_i))

# get data 
incData <- nyc$smoothedCases

# initialize current number of infectious and removed individuals
idxStart <- 5
incData <- incData[-c(1:idxStart)]

# currently infectious
I0 <- sum(nyc$smoothedCases[max(1, (idxStart - lengthI + 1)):(idxStart)])
R0 <- nyc$cumulativeCases[idxStart] - I0

# run three chains in parallel
cl <- makeCluster(3)
clusterExport(cl, list('incData',  'infPeriod_i', 'alarmFit_i',
                       'N', 'I0', 'R0', 'lengthI', 'smoothWindow'))

resThree <- parLapplyLB(cl, 1:3, function(x) {
    
    library(nimble)
    
    # source relevant scripts
    source('./scripts/modelFits.R')
    
    fitAlarmModel(incData = incData, N = N, I0 = I0, R0 = R0, 
                  lengthI = lengthI, infPeriod = infPeriod_i, 
                  alarmFit = alarmFit_i, smoothWindow = smoothWindow,
                  seed = x)
    
})
stopCluster(cl)


source('./scripts/summarizePost.R')

# debugonce(summarizePost)
postSummaries <- summarizePost(resThree = resThree, incData = incData,
                               N = N, I0 = I0, R0 = R0, lengthI = lengthI,
                               alarmFit = alarmFit_i, infPeriod = infPeriod_i, 
                               smoothWindow = smoothWindow)

# save results in separate files
modelInfo <- data.frame(alarmFit = alarmFit_i,
                        infPeriod = infPeriod_i,
                        peak = peak_i)

gr <- cbind.data.frame(postSummaries$gdiag, modelInfo)
paramsPost <- cbind.data.frame(postSummaries$postParams, modelInfo)
alarmPost <- cbind.data.frame(postSummaries$postAlarm, modelInfo)
epiPredPost <- cbind.data.frame(postSummaries$postEpiPred, modelInfo)
betaPost <- cbind.data.frame(postSummaries$postBeta, modelInfo)
R0AlarmPost <- cbind.data.frame(postSummaries$postR0Alarm, modelInfo)
R0Post <- cbind.data.frame(postSummaries$postR0, modelInfo)
waicPost <- cbind.data.frame(postSummaries$waic, modelInfo)

# save output in RDS form
saveRDS(gr, paste0('./Output/grFull', idx, '.rds'))
saveRDS(paramsPost, paste0('./Output/paramsPostFull', idx, '.rds'))
saveRDS(alarmPost, paste0('./Output/alarmPostFull', idx, '.rds'))
saveRDS(epiPredPost, paste0('./Output/epiPredPostFull', idx, '.rds'))
saveRDS(betaPost, paste0('./Output/betaPostFull', idx, '.rds'))
saveRDS(R0AlarmPost, paste0('./Output/R0AlarmPostFull', idx, '.rds'))
saveRDS(R0Post, paste0('./Output/R0PostFull', idx, '.rds'))
saveRDS(waicPost, paste0('./Output/waicPostFull', idx, '.rds'))


