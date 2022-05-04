################################################################################
# function to do posterior prediction on remainder of the epidemic 
################################################################################


postPred <- function(incData, alarmFit, infPeriod, smoothWindow, 
                     paramsPost, alarmSamples, RstarPost = NULL) {
  
  source('./scripts/modelCodesSim.R')
  
  # constants that are true for all alarm types
  obsTime <- 50
  fullTime <- length(incData)
  predTime <- (obsTime + 1):fullTime
  dataObs <- incData[1:obsTime]
  
  # model code
  modelCode <- get(paste0('SIR_', alarmFit, '_', infPeriod))
  
  # model-specific constants, data, and inits
  modelInputs <- getModelInput(alarmFit, incData, infPeriod, smoothWindow)
  
  modelInputs$constantsList$bw <- smoothWindow
  
  # compile model and simulator
  if (alarmFit == 'spline') {
    # spline model needs initial values to avoid warning from NA knots
    myModelPred <- nimbleModel(modelCode, 
                               data = modelInputs$dataList, 
                               constants = modelInputs$constantsList,
                               inits = modelInputs$initsList)
  } else {
    myModelPred <- nimbleModel(modelCode, 
                               data = modelInputs$dataList, 
                               constants = modelInputs$constantsList)
  }
  
  compiledPred  <- compileNimble(myModelPred) 
  dataNodes <- paste0('Istar[', predTime, ']')
  
  if (infPeriod == 'exp') {
    dataNodes <- c(dataNodes, paste0('Rstar[', predTime, ']'))
  }
  
  sim_R <- simulator(myModelPred, dataNodes)
  sim_C <- compileNimble(sim_R)
  
  # get order of parameters
  parentNodes <- myModelPred$getParents(dataNodes, stochOnly = TRUE)
  parentNodes <- parentNodes[-which(parentNodes %in% dataNodes)]
  parentNodes <- myModelPred$expandNodeNames(parentNodes, returnScalarComponents = TRUE)
  
  nPost <- 10000
  postPredInc <- matrix(NA, nrow = 50, ncol = nPost)
  set.seed(1)
  for (j in 1:nPost) {
    
    postIdx <- sample(1:nrow(paramsPost), 1)
    
    betaPost <- paramsPost[postIdx, 'beta']
    
    if (infPeriod == 'exp') {
      rateIPost <- paramsPost[postIdx, 'rateI']
      RstarPostSamp <- RstarPost[postIdx, ]
    }
    
    # model specific parameters
    if (alarmFit == 'thresh') {
      
      alarmParamPost <- paramsPost[postIdx, c('delta', 'H')]
      trueVals <- c(betaPost, alarmParamPost, dataObs)
      
      
    } else if (alarmFit == 'hill') {
      
      alarmParamPost <- paramsPost[postIdx, c('delta', 'nu', 'x0')]
      trueVals <- c(betaPost, alarmParamPost, dataObs)
      
      
    } else if (alarmFit == 'power') {
      
      alarmParamPost <- paramsPost[postIdx, 'k']
      trueVals <- c(betaPost, alarmParamPost, dataObs)
      
      
    } else if (alarmFit == 'spline') {
      
      bPost <- paramsPost[postIdx, grep('b\\[', colnames(paramsPost))]
      knotsPost <- paramsPost[postIdx, grep('knots\\[', colnames(paramsPost))]
      trueVals <- c(betaPost, bPost, knotsPost, dataObs)
      
      
    } else if (alarmFit == 'gp') {
      
      logitAlarmPost <- logit(alarmSamples[,postIdx])[-1]
      names(logitAlarmPost) <- paste0('logit_', names(logitAlarmPost))
      trueVals <- c(betaPost, logitAlarmPost, dataObs)
      
    } else if (alarmFit == 'basic') {
      
      trueVals <- c(betaPost, dataObs)
      
    }
    
    # for exponential infectious period
    if (infPeriod == 'exp') {
      trueVals <- c(trueVals, rateIPost, RstarPostSamp)
    }
    
    trueVals <- trueVals[parentNodes]
    
    
    postPredInc[,j] <- apply(sim_C$run(trueVals, 10), 2, median)[grep('Istar', dataNodes)]
  }
  
  postPredInc
}
