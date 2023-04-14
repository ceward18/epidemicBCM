################################################################################
# Run models by peak for NYC
# for each peak run 7 different models
#   thresh, hill, power, gp, spline, betat, basic
################################################################################

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
idx <- as.numeric(task_id)

### load libraries
library(parallel)
library(nimble)

### source scripts (for movingAverage function)
source('./scripts/modelCodes.R')

### read data
dat <- read.csv('./data/ebolaClean.csv')
dat$smoothedCases <- round(movingAverage(dat$dailyCases, 7))
dat$cumulativeCases <- cumsum(dat$smoothedCases)

alarmFit <- c('power', 'thresh', 'hill', 
              'spline', 'splineFixKnot', 'gp', 
              'basic', 'betatSpline')
smoothWindow <- c(30, 60)

# 12 possibilities (6 alarmFits, 2 smoothWindows)
allModelsAlarm <- expand.grid(alarmFit = alarmFit[c(1:6)],
                              smoothWindow = smoothWindow)

# 2 possibilities (2 alarmFits, 1 smoothWindow)
allModelsNoAlarm <- expand.grid(alarmFit = alarmFit[7:8],
                                smoothWindow = 1)

# 14
allModels <- rbind.data.frame(allModelsAlarm, allModelsNoAlarm)
allModels <- allModels[order(allModels$alarmFit,
                             allModels$smoothWindow),]
rownames(allModels) <- NULL

# list of batches
tmp <- allModels[seq(1, nrow(allModels), 2),]
rownames(tmp) <- NULL

# constants for all models
N <- 5363500
E0 <- 1
I0 <- 0
R0 <- 0

# batches by alarmFit (56 batches total)
batchSize <- 4
batchIdx <- batchSize * (idx - 1) + 1:batchSize


for (i in batchIdx) {
    
    alarmFit_i <- allModels$alarmFit[i]
    smoothWindow_i <- allModels$smoothWindow[i]
    
    print(paste('Running alarm:', alarmFit_i,
                ', smoothing:', smoothWindow_i))
    
    # get data 
    incData <- dat$onset
    deathData <- dat$death
    
    # smoothed incidence to inform alarm function 
    # (shifted so alarm is informed only by data up to time t-1)
    smoothI <-  head(movingAverage(c(0, incData), smoothWindow_i), -1)
    
    # run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('incData', 'deathData', 'smoothI', 'alarmFit_i', 
                           'smoothWindow_i', 'N', 'E0', 'I0', 'R0'))
    
    resThree <- parLapplyLB(cl, 1:3, function(x) {
        
        library(nimble)
        
        # source relevant scripts
        source('./scripts/modelFits.R')
        
        # debugonce(fitAlarmModel)
        fitAlarmModel(incData = incData, deathData = deathData,
                      smoothI = smoothI,
                      N = N, E0 = E0, I0 = I0, R0 = R0, 
                      smoothWindow = smoothWindow_i,
                      alarmFit = alarmFit_i, seed = x)
        
    })
    stopCluster(cl)
   
    source('./scripts/summarizePost.R')
    # debugonce(summarizePost)
    # debugonce(postPredFit)
    postSummaries <- summarizePost(resThree = resThree, incData = incData,
                                   deathData = deathData,
                                   smoothI = smoothI, smoothWindow = smoothWindow_i,
                                   N = N, E0 = E0, I0 = I0, R0 = R0, 
                                   alarmFit = alarmFit_i)
    
    # if the model did not converge save the chains so these can be examined later
    if (!all(postSummaries$gdiag$gr < 1.1)) {
        
        # create thinned version
        resThree[[1]] <- resThree[[1]][seq(1,nrow(resThree[[1]]), 10),]
        resThree[[2]] <- resThree[[2]][seq(1,nrow(resThree[[2]]), 10),]
        resThree[[3]] <- resThree[[3]][seq(1,nrow(resThree[[3]]), 10),]
        
        saveRDS(resThree, 
                paste0('./output/chains_', alarmFit_i, '_', 
                       smoothWindow_i, '.rds'))
    }
    
    
    # save results in separate files
    modelInfo <- data.frame(alarmFit = alarmFit_i,
                            smoothWindow = smoothWindow_i)
    
    if (i == batchIdx[1]) {
        gr <- cbind.data.frame(postSummaries$gdiag, modelInfo)
        paramsPost <- cbind.data.frame(postSummaries$postParams, modelInfo)
        alarmPost <- cbind.data.frame(postSummaries$postAlarm, modelInfo)
        postPredFits <- cbind.data.frame(postSummaries$postPredFit, modelInfo)
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
        postPredFits <- rbind.data.frame(postPredFits, 
                                         cbind.data.frame(postSummaries$postPredFit, modelInfo))
        betaPost <- rbind.data.frame(betaPost, 
                                     cbind.data.frame(postSummaries$postBeta, modelInfo))
        R0Post <- rbind.data.frame(R0Post, 
                                   cbind.data.frame(postSummaries$postR0, modelInfo))
        waicPost <- rbind.data.frame(waicPost, 
                                     cbind.data.frame(postSummaries$waic, modelInfo))
    }
    
} # end loop

idxPrint <- sprintf("%02d",idx)

# save output in RDS form
saveRDS(gr, paste0('./output/grBatch', idxPrint, '.rds'))
saveRDS(paramsPost, paste0('./output/paramsPostBatch', idxPrint, '.rds'))
saveRDS(alarmPost, paste0('./output/alarmPostBatch', idxPrint, '.rds'))
saveRDS(postPredFits, paste0('./output/postPredFitBatch', idxPrint, '.rds'))
saveRDS(betaPost, paste0('./output/betaPostBatch', idxPrint, '.rds'))
saveRDS(R0Post, paste0('./output/R0PostBatch', idxPrint, '.rds'))
saveRDS(waicPost, paste0('./output/waicPostBatch', idxPrint, '.rds'))


