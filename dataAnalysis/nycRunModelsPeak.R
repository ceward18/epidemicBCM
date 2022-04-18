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

peak <- c('1', '2', '3', '4')
alarmFit <- c( 'thresh', 'hill', 'power', 'gp', 'spline', 'betat', 'basic')
infPeriod <- 'fixed'

# 28 possibilities (7 alarmFits, 4 peaks)
allModels <- expand.grid(peak = peak,
                         alarmFit = alarmFit,
                         infPeriod = infPeriod)

# constants for all models
N <- nyc$Population[1]
lengthI <- 7
smoothWindow <- 30

# batches by alarmFit (7 batches total)
batchSize <- 4
batchIdx <- batchSize * (idx - 1) + 1:batchSize


for (i in 1:batchIdx) {
    
    peak_i <- allModels$peak[i]
    alarmFit_i <- allModels$alarmFit[i]
    infPeriod_i <- allModels$infPeriod[i]
    
    print(paste('Running alarm:', alarmFit_i,
                ', peak:', peak_i, 
                ', infPeriod:', infPeriod_i))
    
    # get data for the specified peak
    if (peak_i == 'full') {
        incData <- nyc$smoothedCases
    } else {
        incData <- nyc$smoothedCases[which(nyc$peak == peak_i)]
    }
    
    # initialize current number of infectious and removed individuals
    
    if (peak_i %in% c('full', '1')) {
        idxStart <- 5
    } else {
        idxStart <- min(which(nyc$peak == peak_i))
    }
    
    # currently infectious
    I0 <- sum(nyc$smoothedCases[max(1, (idxStart - lengthI + 1)):(idxStart)])
    R0 <- nyc$cumulativeCases[idxStart] - I0
    
    # first time point is included in initial values
    incData <- incData[-1]
    
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
    
    if (i == batchIdx[1]) {
        gr <- cbind.data.frame(postSummaries$gdiag, modelInfo)
        paramsPost <- cbind.data.frame(postSummaries$postParams, modelInfo)
        alarmPost <- cbind.data.frame(postSummaries$postAlarm, modelInfo)
        epiPredPost <- cbind.data.frame(postSummaries$postEpiPred, modelInfo)
        betaPost <- cbind.data.frame(postSummaries$postBeta, modelInfo)
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
saveRDS(waicPost, paste0('./Output/waicPostBatch', idx, '.rds'))


