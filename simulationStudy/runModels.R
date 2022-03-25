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

infPeriod <- c('fixed')
alarmGen <- c('thresh', 'hill', 'power')
alarmFit <- c('thresh', 'hill', 'power', 'spline', 'gp', 'betat', 'basic')
smoothWindow <- c(14, 30)


allFits <- expand.grid(simNumber = 1:nSim,
                       infPeriod = infPeriod,
                       alarmGen = alarmGen,
                       alarmFit = alarmFit,
                       smoothWindow = smoothWindow,
                       stringsAsFactors = FALSE)

# 1500 rows
allFits <- allFits[-which(allFits$alarmFit %in% alarmGen &
                            allFits$alarmFit != allFits$alarmGen),]
rownames(allFits) <- 1:nrow(allFits)


# fit models in batches of 50 (30 batches total)
batchSize <- 50
batchIdx <- batchSize * (idx - 1) + 1:batchSize

for (i in batchIdx) {
  
  print(paste0('Now starting row: ', i))
  
  simNumber_i <- allFits$simNumber[i]
  infPeriod_i <- allFits$infPeriod[i]
  alarmGen_i <- allFits$alarmGen[i]
  alarmFit_i <- allFits$alarmFit[i]
  smoothWindow_i <- allFits$smoothWindow[i]
  
  # load data
  incData <- readRDS(paste0('./Data/', alarmGen_i, '_', 
                            infPeriod_i, '_', smoothWindow_i, '.rds'))
  
  # subset row corresponding to simulation number specified
  incData <- incData[simNumber_i,]
  
  # only use the first 50 time points for model fitting
  incDataFit <- incData[1:50]
  
  # run three chains in parallel
  cl <- makeCluster(3)
  clusterExport(cl, list('incDataFit',  'infPeriod_i', 'alarmFit_i',
                         'smoothWindow_i'))
  
  resThree <- parLapplyLB(cl, 1:3, function(x) {
    
    library(nimble)
    
    # source relevant scripts
    source('./scripts/modelFits.R')
    
    fitAlarmModel(incData = incDataFit, infPeriod = infPeriod_i, 
                  alarmFit = alarmFit_i, smoothWindow = smoothWindow_i,
                  seed = x)
  })
  stopCluster(cl)
  
  source('./scripts/summarizePost.R')
  
  # debugonce(summarizePost)
  postSummaries <- summarizePost(resThree = resThree, incData = incData,
                                 alarmFit = alarmFit_i, infPeriod = infPeriod_i, 
                                 smoothWindow = smoothWindow_i)
  
  # save results in separate files
  modelInfo <- data.frame(alarmGen = alarmGen_i,
                          alarmFit = alarmFit_i,
                          infPeriod = infPeriod_i,
                          smoothWindow = smoothWindow_i,
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
saveRDS(alarmPost, paste0('./Output/betaPostBatch', idx, '.rds'))
saveRDS(alarmPost, paste0('./Output/waicPostBatch', idx, '.rds'))




