################################################################################
# Run models by peak for Montreal
# for each peak run 7 different models
#   thresh, hill, power, gp, spline, betatSpline, basic
################################################################################

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
idx <- as.numeric(task_id)

### load libraries
library(parallel)
library(nimble)

### source scripts (for movingAverage function)
source('../scripts/modelCodes.R')

### read data
montreal <- read.csv('./Data/montrealClean.csv')
montreal$smoothedCases <- round(movingAverage(montreal$dailyCases, 3))
montreal$cumulativeCases <- cumsum(montreal$smoothedCases)


peak <- c('1', '2', '3', '4', '5')
alarmFit <- c( 'thresh', 'hill', 'power', 'gp', 'spline', 'betatSpline', 'basic')
smoothWindow <- c(14, 30)

# 70 possibilities (7 alarmFits, 5 peaks, 2 infPeriods)
allModelsFixed <- expand.grid(peak = peak,
                              alarmFit = alarmFit,
                              smoothWindow = smoothWindow,
                              infPeriod = 'fixed')
allModelsExp <- expand.grid(peak = peak,
                            alarmFit = alarmFit,
                            smoothWindow = smoothWindow,
                            infPeriod = 'exp')

allModels <- rbind.data.frame(allModelsFixed, allModelsExp)

# constants for all models
N <- montreal$Population[1]
lengthI <- 5
smoothWindow <- 30

# batches by alarmFit (7 batches total)
batchSize <- 5
batchIdx <- batchSize * (idx - 1) + 1:batchSize


for (i in batchIdx) {
    
    peak_i <- allModels$peak[i]
    alarmFit_i <- allModels$alarmFit[i]
    infPeriod_i <- allModels$infPeriod[i]
    smoothWindow_i <- allModels$smoothWindow[i]
    
    print(paste('Running alarm:', alarmFit_i,
                ', peak:', peak_i, 
                ', infPeriod:', infPeriod_i,
                ', smoothWindow:', smoothWindow_i))
    
    # smoothed incidence to inform alarm function 
    # (shifted so alarm is informed only by data up to time t-1)
    montreal$smoothI <- head(movingAverage(c(0, montreal$dailyCases), smoothWindow_i), -1)
    
    # get data for the specified peak
    if (peak_i == 'full') {
        incData <- montreal$smoothedCases
        smoothI <- montreal$smoothI
    } else {
        incData <- montreal$smoothedCases[which(montreal$peak == peak_i)]
        smoothI <- montreal$smoothI[which(montreal$peak == peak_i)]
    }
    
    # initialize current number of infectious and removed individuals
    
    if (peak_i %in% c('full', '1')) {
        idxStart <- 15
        incData <- incData[-c(1:idxStart)]
        smoothI <- smoothI[-c(1:idxStart)]
    } else {
        idxStart <- min(which(montreal$peak == peak_i))
        incData <- incData[-1]
        smoothI <- smoothI[-1]
    }
    
    
    # currently infectious
    I0 <- sum(montreal$smoothedCases[max(1, (idxStart - lengthI + 1)):(idxStart)])
    R0 <- montreal$cumulativeCases[idxStart] - I0
    
    Rstar0 <- montreal$smoothedCases[max(1, (idxStart - lengthI + 1)):(idxStart)]
    
    # run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('incData', 'smoothI', 'infPeriod_i', 'alarmFit_i',
                           'N', 'I0', 'R0', 'Rstar0', 'lengthI', 'smoothWindow_i'))
    
    resThree <- parLapplyLB(cl, 1:3, function(x) {
        
        library(nimble)
        
        # source relevant scripts
        source('../scripts/modelFits.R')
        
        fitAlarmModel(incData = incData, smoothI = smoothI,
                      N = N, I0 = I0, R0 = R0, Rstar0 = Rstar0,
                      lengthI = lengthI, infPeriod = infPeriod_i, 
                      alarmFit = alarmFit_i, seed = x)
        
    })
    stopCluster(cl)
    
    source('../scripts/summarizePost.R')
    
    # debugonce(summarizePost)
    postSummaries <- summarizePost(resThree = resThree, incData = incData,
                                   smoothI = smoothI,
                                   N = N, I0 = I0, R0 = R0, Rstar0 = Rstar0,
                                   lengthI = lengthI, alarmFit = alarmFit_i, 
                                   infPeriod = infPeriod_i)
    
    # save results in separate files
    modelInfo <- data.frame(alarmFit = alarmFit_i,
                            infPeriod = infPeriod_i,
                            peak = peak_i,
                            smoothWindow = smoothWindow_i)
    
    if (i == batchIdx[1]) {
        gr <- cbind.data.frame(postSummaries$gdiag, modelInfo)
        paramsPost <- cbind.data.frame(postSummaries$postParams, modelInfo)
        alarmPost <- cbind.data.frame(postSummaries$postAlarm, modelInfo)
        epiPredPost <- cbind.data.frame(postSummaries$postEpiPred, modelInfo)
        betaPost <- cbind.data.frame(postSummaries$postBeta, modelInfo)
        R0Post <- cbind.data.frame(postSummaries$postR0, modelInfo)
        waicPost <- cbind.data.frame(postSummaries$waic, modelInfo)
        
    } else {
        gr <- rbind.data.frame(gr, 
                               cbind.data.frame(postSummaries$gdiag, modelInfo))
        paramsPost <- rbind.data.frame(paramsPost, 
                                       cbind.data.frame(postSummaries$postParams, modelInfo))
        alarmPost <- rbind.data.frame(alarmPost, 
                                      cbind.data.frame(postSummaries$postAlarm, modelInfo))
        epiPredPost <- rbind.data.frame(epiPredPost, 
                                        cbind.data.frame(postSummaries$postEpiPred, modelInfo))
        betaPost <- rbind.data.frame(betaPost, 
                                     cbind.data.frame(postSummaries$postBeta, modelInfo))
        R0Post <- rbind.data.frame(R0Post, 
                                   cbind.data.frame(postSummaries$postR0, modelInfo))
        waicPost <- rbind.data.frame(waicPost, 
                                     cbind.data.frame(postSummaries$waic, modelInfo))
    }
    
} # end loop

# save output in RDS form
saveRDS(gr, paste0('./Output/grBatch', idx, '.rds'))
saveRDS(paramsPost, paste0('./Output/paramsPostBatch', idx, '.rds'))
saveRDS(alarmPost, paste0('./Output/alarmPostBatch', idx, '.rds'))
saveRDS(epiPredPost, paste0('./Output/epiPredPostBatch', idx, '.rds'))
saveRDS(betaPost, paste0('./Output/betaPostBatch', idx, '.rds'))
saveRDS(R0Post, paste0('./Output/R0PostBatch', idx, '.rds'))
saveRDS(waicPost, paste0('./Output/waicPostBatch', idx, '.rds'))


