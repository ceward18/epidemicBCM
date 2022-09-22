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
dat <- read.csv('./data/nycClean.csv')
dat$smoothedCases <- round(movingAverage(dat$dailyCases, 7))
dat$cumulativeCases <- cumsum(dat$smoothedCases)

peak <- c('1', '2', '3', '4')
alarmFit <- c('power', 'thresh', 'hill', 
              'spline', 'splineFixKnot', 'gp', 
              'basic', 'betatSpline')
smoothWindow <- c(30, 60)
prior <- 1:4

# 96 possibilities (6 alarmFits, 4 peaks, 4 priors, 2 smoothWindows)
allModelsAlarm <- expand.grid(peak = peak,
                              alarmFit = alarmFit[c(1:6)],
                              smoothWindow = smoothWindow,
                              prior = prior)

# 32 possibilities (2 alarmFits, 4 peaks, 4 priors, 1 smoothWindow)
allModelsNoAlarm <- expand.grid(peak = peak,
                                alarmFit = alarmFit[7:8],
                                smoothWindow = 1,
                                prior = prior)

# 112
allModels <- rbind.data.frame(allModelsAlarm, allModelsNoAlarm)
allModels <- allModels[order(allModels$alarmFit,
                             allModels$smoothWindow, allModels$peak),]
rownames(allModels) <- NULL

# list of batches
tmp <- allModels[seq(1, nrow(allModels), 4),]
rownames(tmp) <- NULL

# constants for all models
N <- dat$Population[1]

# used to initialize I0 and R0
lengthI <- 5

# batches by alarmFit (56 batches total)
batchSize <- 4
batchIdx <- batchSize * (idx - 1) + 1:batchSize


for (i in batchIdx) {
    
    peak_i <- allModels$peak[i]
    alarmFit_i <- allModels$alarmFit[i]
    smoothWindow_i <- allModels$smoothWindow[i]
    prior_i <- allModels$prior[i]
    
    print(paste('Running alarm:', alarmFit_i,
                ', peak:', peak_i, 
                ', smoothing:', smoothWindow_i, 
                ', prior:', prior_i))
    
    # smoothed incidence to inform alarm function 
    # (shifted so alarm is informed only by data up to time t-1)
    dat$smoothI <- head(movingAverage(c(0, dat$dailyCases), smoothWindow_i), -1)
    
    # get data for the specified peak
    incData <- dat$smoothedCases[which(dat$peak == peak_i)]
    smoothI <- dat$smoothI[which(dat$peak == peak_i)]
    
    # initialize current number of infectious and removed individuals
    
    if (peak_i == 1) {
        idxStart <- 5
        incData <- incData[-c(1:idxStart)]
        smoothI <- smoothI[-c(1:idxStart)]
    } else {
        idxStart <- min(which(dat$peak == peak_i))
        incData <- incData[-1]
        smoothI <- smoothI[-1]
    }
    
    # currently infectious
    I0 <- sum(dat$smoothedCases[max(1, (idxStart - lengthI + 1)):(idxStart)])
    R0 <- dat$cumulativeCases[idxStart] - I0 
    
    # used to initialize removal vector
    Rstar0 <- dat$smoothedCases[max(1, idxStart - lengthI + 1):(idxStart)]
    
    
    # run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('incData', 'smoothI', 'alarmFit_i', 'prior_i', 
                           'peak_i', 'smoothWindow_i',
                           'N', 'I0', 'R0', 'Rstar0', 'lengthI'))
    
    resThree <- parLapplyLB(cl, 1:3, function(x) {
        
        library(nimble)
        
        # source relevant scripts
        source('./scripts/modelFits.R')
        
        fitAlarmModel(incData = incData, smoothI = smoothI,
                      N = N, I0 = I0, R0 = R0, 
                      Rstar0 = Rstar0, lengthI = lengthI, 
                      prior = prior_i, peak = peak_i, 
                      smoothWindow = smoothWindow_i,
                      alarmFit = alarmFit_i, seed = x)
        
    })
    stopCluster(cl)
   
    
    source('./scripts/summarizePost.R')
    # debugonce(summarizePost)
    postSummaries <- summarizePost(resThree = resThree, incData = incData,
                                   smoothI = smoothI, smoothWindow = smoothWindow_i,
                                   N = N, I0 = I0, R0 = R0, 
                                   Rstar0 = Rstar0, lengthI = lengthI, 
                                   alarmFit = alarmFit_i, prior = prior_i,
                                   peak = peak_i)
    
    # if the model did not converge save the chains so these can be examined later
    if (!all(postSummaries$gdiag$gr < 1.1)) {
        
        # create thinned version
        resThree[[1]] <- resThree[[1]][seq(1,nrow(resThree[[1]]), 10),]
        resThree[[2]] <- resThree[[2]][seq(1,nrow(resThree[[2]]), 10),]
        resThree[[3]] <- resThree[[3]][seq(1,nrow(resThree[[3]]), 10),]
        
        saveRDS(resThree, 
                paste0('./output/chains_', alarmFit_i, '_peak', 
                       peak_i, '_', smoothWindow_i, '_', prior_i, '.rds'))
    }
    
    
    # save results in separate files
    modelInfo <- data.frame(alarmFit = alarmFit_i,
                            peak = peak_i,
                            smoothWindow = smoothWindow_i,
                            prior = prior_i)
    
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


