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
alarmGen <- c('thresh', 'hill', 'power')
alarmFit <- c('thresh', 'hill', 'power', 'spline', 'gp', 'betatSpline', 'basic')
smoothWindow <- c(14, 30)
priors <- 1:4
epiSize <- c('small', 'large')

allFits <- expand.grid(simNumber = 1:nSim,
                       alarmGen = alarmGen,
                       alarmFit = alarmFit,
                       smoothWindow = smoothWindow,
                       prior = priors,
                       epiSize = epiSize,
                       stringsAsFactors = FALSE)


# 6000 rows
allFits <- allFits[-which(allFits$alarmFit %in% alarmGen &
                              allFits$alarmFit != allFits$alarmGen),]
rownames(allFits) <- 1:nrow(allFits)

tmp <- allFits[seq(1,nrow(allFits), 25),]
rownames(tmp) <- 1:nrow(tmp)

# fit models in batches of 25 (240 batches total)
batchSize <- 25
batchIdx <- batchSize * (idx - 1) + 1:batchSize

for (i in batchIdx) {
    
    simNumber_i <- allFits$simNumber[i]
    alarmGen_i <- allFits$alarmGen[i]
    alarmFit_i <- allFits$alarmFit[i]
    smoothWindow_i <- allFits$smoothWindow[i]
    prior_i <- allFits$prior[i]
    epiSize_i <- allFits$epiSize[i]
    
    print(paste0('alarm Gen: ', alarmGen_i,
                 ', alarm fit: ', alarmFit_i, 
                 ', smoothing window: ', smoothWindow_i, 
                 ', prior: ', prior_i, 
                 ', size: ', epiSize_i, 
                 ', simulation: ', simNumber_i))
    
    # load data
    if (epiSize_i == 'small') {
        print(paste0('./Data/', alarmGen_i, '_exp_', smoothWindow_i, '.rds'))
        incData <- readRDS(paste0('./Data/', alarmGen_i, '_exp_', smoothWindow_i, '.rds'))
        
    } else if (epiSize_i == 'large') {
        
        incData <- readRDS(paste0('./Data/', alarmGen_i, '_exp_', smoothWindow_i, '_large.rds'))
    }
    incData <- incData[,grep('Istar', colnames(incData))]
    
    # subset row corresponding to simulation number specified
    incData <- incData[simNumber_i,]
    
    # only use the first 50 time points for model fitting
    if (epiSize_i == 'small') {
        incDataFit <- incData[1:50]
    } else if (epiSize_i == 'large') {
        incDataFit <- incData[1:100]
    }
    
    
    # run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('incDataFit', 'alarmFit_i', 'smoothWindow_i',
                           'prior_i', 'simNumber_i'))
    
    resThree <- parLapplyLB(cl, 1:3, function(x) {
        
        library(nimble)
        
        # source relevant scripts
        source('./scripts/modelFits.R')
        
        fitAlarmModel(incData = incDataFit, alarmFit = alarmFit_i,
                      smoothWindow = smoothWindow_i, prior = prior_i,
                      simNumber = simNumber_i, seed = x)
    })
    stopCluster(cl)
    
    source('./scripts/summarizePost.R')
    
    # debugonce(summarizePost)
    postSummaries <- summarizePost(resThree = resThree, incData = incData,
                                   alarmFit = alarmFit_i, smoothWindow = smoothWindow_i,
                                   prior = prior_i, epiSize = epiSize_i)
    
    # save results in separate files
    modelInfo <- data.frame(alarmGen = alarmGen_i,
                            alarmFit = alarmFit_i,
                            smoothWindow = smoothWindow_i,
                            prior = prior_i,
                            epiSize = epiSize_i,
                            simNumber = simNumber_i)
    
    # gelman-rubin
    # posterior parameters
    # posterior alarm
    # posterior predictions
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
saveRDS(waicPost, paste0('./Output/waicPostBatch', idx, '.rds'))




