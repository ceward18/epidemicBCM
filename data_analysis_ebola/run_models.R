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
# delete first row as this is used in initial values
dat <- dat[-1,]

alarmFit <- c('power', 'thresh', 'hill', 
              'spline', 'gp', 
              'basic', 'betatSpline', 'basicInt')

# 7
allModels <- data.frame(alarmFit = alarmFit)


# list of batches
tmp <- allModels[seq(1, nrow(allModels), 1),]
rownames(tmp) <- NULL

# constants for all models
N <- 5363500
E0 <- 2             # two infections at start - assume already exposed
I0 <- 1
R0 <- 2             # two deaths already occurred

# intervention time for intervention model
intTime <- which(dat$date == as.Date('1995-05-09'))


alarmFit_i <- allModels$alarmFit[idx]

print(paste('Running alarm:', alarmFit_i))

# get data 
incData <- dat$onset
deathData <- dat$death

# cumulative incidence to inform alarm function 
smoothI <- c(0, head(cumsum(incData), -1))

# run three chains in parallel
cl <- makeCluster(3)
clusterExport(cl, list('incData', 'deathData', 'smoothI', 'alarmFit_i', 
                       'N', 'E0', 'I0', 'R0', 'intTime'))

resThree <- parLapplyLB(cl, 1:3, function(x) {
    
    library(nimble)
    
    # source relevant scripts
    source('./scripts/modelFits.R')
    
    # debugonce(fitAlarmModel)
    fitAlarmModel(incData = incData, deathData = deathData,  smoothI = smoothI,
                  N = N, E0 = E0, I0 = I0, R0 = R0, intTime = intTime,
                  alarmFit = alarmFit_i, seed = x)
    
})
stopCluster(cl)

source('./scripts/summarizePost.R')
# debugonce(summarizePost)
# debugonce(postPredFit)
postSummaries <- summarizePost(resThree = resThree, incData = incData,
                               deathData = deathData,
                               smoothI = smoothI, 
                               N = N, E0 = E0, I0 = I0, R0 = R0, 
                               alarmFit = alarmFit_i)

# if the model did not converge save the chains so these can be examined later
if (!all(postSummaries$gdiag$gr < 1.1)) {
    
    # create thinned version
    resThree[[1]] <- resThree[[1]][seq(1,nrow(resThree[[1]]), 10),]
    resThree[[2]] <- resThree[[2]][seq(1,nrow(resThree[[2]]), 10),]
    resThree[[3]] <- resThree[[3]][seq(1,nrow(resThree[[3]]), 10),]
    
    saveRDS(resThree, 
            paste0('./output/chains_', alarmFit_i, '.rds'))
}


# save results in separate files
modelInfo <- data.frame(alarmFit = alarmFit_i)

gr <- cbind.data.frame(postSummaries$gdiag, modelInfo)
paramsPost <- cbind.data.frame(postSummaries$postParams, modelInfo)
alarmPost <- cbind.data.frame(postSummaries$postAlarm, modelInfo)
R0Post <- cbind.data.frame(postSummaries$postR0, modelInfo)
compsPost <- cbind.data.frame(postSummaries$postComps, modelInfo)
postPredFits <- cbind.data.frame(postSummaries$postPredFit, modelInfo)
betaPost <- cbind.data.frame(postSummaries$postBeta, modelInfo)
waicPost <- cbind.data.frame(postSummaries$waic, modelInfo)


# save output in RDS form
saveRDS(gr, paste0('./output/grBatch', idx, '.rds'))
saveRDS(paramsPost, paste0('./output/paramsPostBatch', idx, '.rds'))
saveRDS(alarmPost, paste0('./output/alarmPostBatch', idx, '.rds'))
saveRDS(R0Post, paste0('./output/R0PostBatch', idx, '.rds'))
saveRDS(compsPost, paste0('./output/compsPostBatch', idx, '.rds'))
saveRDS(postPredFits, paste0('./output/postPredFitBatch', idx, '.rds'))
saveRDS(betaPost, paste0('./output/betaPostBatch', idx, '.rds'))
saveRDS(waicPost, paste0('./output/waicPostBatch', idx, '.rds'))


