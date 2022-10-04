################################################################################
# Model fitting script
################################################################################

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
idx <- as.numeric(task_id)

# load libraries
library(nimble)
library(parallel)

# set up grid of models to fit
nSim <- 50
alarmGen <- c('power', 'thresh', 'hill')
alarmFit <- c('power', 'thresh', 'hill', 'spline', 'gp', 'basic', 'betatSpline')
smoothWindow <- c(14, 30)

allFits <- expand.grid(simNumber = 1:nSim,
                       alarmGen = alarmGen,
                       alarmFit = alarmFit,
                       smoothWindow = smoothWindow,
                       stringsAsFactors = FALSE)

# 1500 rows 5 alarmFits * 3 alarmGens * 50 epidemics 
allFits <- allFits[-which(allFits$alarmFit %in% alarmGen &
                              allFits$alarmFit != allFits$alarmGen),]
rownames(allFits) <- NULL

# tmp <- allFits[seq(1,nrow(allFits), 25),]
# rownames(tmp) <- NULL

# fit models in batches of 25 (60 batches total)
batchSize <- 25
batchIdx <- batchSize * (idx - 1) + 1:batchSize

for (i in batchIdx) {
    
    simNumber_i <- allFits$simNumber[i]
    alarmGen_i <- allFits$alarmGen[i]
    alarmFit_i <- allFits$alarmFit[i]
    smoothWindow_i <- allFits$smoothWindow[i]
    
    print(paste0('alarm Gen: ', alarmGen_i,
                 ', alarm fit: ', alarmFit_i, 
                 ', smoothing window: ', smoothWindow_i, 
                 ', simulation: ', simNumber_i))
    
    # load data
    incData <- readRDS(paste0('./data/', alarmGen_i, '_', 
                              smoothWindow_i, '.rds'))
    
    # only use data on incidence (ignore removal times)
    incData <- incData[,grep('Istar', colnames(incData))]
    
    # subset row corresponding to simulation number specified
    incData <- incData[simNumber_i,]
    
    # only use the first 50 time points for model fitting
    incDataFit <- incData[1:50]
    
    # run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('incDataFit', 'alarmFit_i',
                           'smoothWindow_i', 'simNumber_i'))
    
    resThree <- parLapplyLB(cl, 1:3, function(x) {
        
        library(nimble)
        
        # source relevant scripts
        source('./scripts/modelFits.R')
        
        fitAlarmModel(incData = incDataFit, alarmFit = alarmFit_i, 
                      smoothWindow = smoothWindow_i,  simNumber = simNumber_i,
                      seed = x)
    })
    stopCluster(cl)
    
    source('./scripts/summarizePost.R')
    
    # debugonce(summarizePost)
    # debugonce(postPredFit)
    postSummaries <- summarizePost(resThree = resThree, incData = incData,
                                   incDataFit = incDataFit, 
                                   alarmFit = alarmFit_i, 
                                   smoothWindow = smoothWindow_i)
    
    # if the model did not converge save the chains so these can be examined later
    if (!all(postSummaries$gdiag$gr < 1.1)) {
        
        # create thinned version
        resThree[[1]] <- resThree[[1]][seq(1,nrow(resThree[[1]]), 10),]
        resThree[[2]] <- resThree[[2]][seq(1,nrow(resThree[[2]]), 10),]
        resThree[[3]] <- resThree[[3]][seq(1,nrow(resThree[[3]]), 10),]
        
        saveRDS(resThree, 
                paste0('./output/chains_', alarmGen, '_', alarmFit_i,
                       '_', smoothWindow_i, '_', simNumber_i, '.rds'))
    }
    
    
    # save results in separate files
    modelInfo <- data.frame(alarmGen = alarmGen_i,
                            alarmFit = alarmFit_i,
                            smoothWindow = smoothWindow_i,
                            simNumber = simNumber_i)
    
    # posterior summaries
    if (i == batchIdx[1]) {
        gr <- cbind.data.frame(postSummaries$gdiag, modelInfo)
        paramsPost <- cbind.data.frame(postSummaries$postParams, modelInfo)
        alarmPost <- cbind.data.frame(postSummaries$postAlarm, modelInfo)
        epiPredPost <- cbind.data.frame(postSummaries$postEpiPred, modelInfo)
        predFitPost <- cbind.data.frame(postSummaries$postPredFit, modelInfo)
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
        predFitPost <- rbind.data.frame(predFitPost, 
                                        cbind.data.frame(postSummaries$postPredFit, modelInfo))
        betaPost <- rbind.data.frame(betaPost, 
                                     cbind.data.frame(postSummaries$postBeta, modelInfo))
        waicPost <- rbind.data.frame(waicPost, 
                                     cbind.data.frame(postSummaries$waic, modelInfo))
        
    }
    
} # end loop

idxPrint <- sprintf("%02d",idx)

# save output in RDS form
saveRDS(gr, paste0('./output/grBatch', idxPrint, '.rds'))
saveRDS(paramsPost, paste0('./output/paramsPostBatch', idxPrint, '.rds'))
saveRDS(alarmPost, paste0('./output/alarmPostBatch', idxPrint, '.rds'))
saveRDS(epiPredPost, paste0('./output/epiPredPostBatch', idxPrint, '.rds'))
saveRDS(predFitPost, paste0('./output/predFitPostBatch', idxPrint, '.rds'))
saveRDS(betaPost, paste0('./output/betaPostBatch', idxPrint, '.rds'))
saveRDS(waicPost, paste0('./output/waicPostBatch', idxPrint, '.rds'))




