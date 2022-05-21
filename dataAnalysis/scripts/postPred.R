################################################################################
# function to do posterior prediction on remainder of the epidemic 
################################################################################


postPred <- function(incData, N, I0, R0, Rstar0, lengthI,
                     alarmFit, infPeriod, smoothWindow, 
                     paramsPost, alarmSamples) {
  
  source('../scripts/modelCodesSim.R')
  
  # want to predict out 100 days from last observed time point
  nDaysPred <- 100
  obsTime <- length(incData)
  fullTime <- obsTime + nDaysPred
  predTime <- (obsTime + 1):fullTime
  dataObs <- incData

  # model code
  modelCode <- get(paste0('SIR_', alarmFit, '_', infPeriod))
  
  modelInputs <- getModelInput(alarmFit, 
                               c(incData, rep(1, nDaysPred)), 
                               infPeriod, smoothWindow, 
                               N, I0, R0, Rstar0, lengthI)
  
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
  
  names(dataObs) <- paste0('Istar[', 1:obsTime, ']')
  
  # get order of parameters
  parentNodes <- myModelPred$getParents(dataNodes, stochOnly = TRUE)
  parentNodes <- parentNodes[-which(parentNodes %in% dataNodes)]
  parentNodes <- myModelPred$expandNodeNames(parentNodes, returnScalarComponents = TRUE)
  
  nPost <- 10000
  postPredInc <- matrix(NA, nrow = nDaysPred, ncol = nPost)
  set.seed(1)
  for (j in 1:nPost) {
    
    postIdx <- sample(1:nrow(paramsPost), 1)
    
    betaPost <- paramsPost[postIdx,'beta']
    RstarPost <- paramsPost[postIdx,grep('Rstar', colnames(paramsPost))]
    
    if (infPeriod == 'exp') {
      rateIPost <- paramsPost[postIdx, 'rateI']
    }
    
    
    # model specific parameters
    if (alarmFit == 'thresh') {
      
      alarmParamPost <- paramsPost[postIdx, c('delta', 'H')]
      trueVals <- c(betaPost, alarmParamPost, dataObs, RstarPost)
      
      
    } else if (alarmFit == 'hill') {
      
      alarmParamPost <- paramsPost[postIdx, c('delta', 'nu', 'x0')]
      trueVals <- c(betaPost, alarmParamPost, dataObs, RstarPost)
      
      
    } else if (alarmFit == 'power') {
      
      alarmParamPost <- paramsPost[postIdx, 'k']
      trueVals <- c(betaPost, alarmParamPost, dataObs, RstarPost)
      
      
    } else if (alarmFit == 'spline') {
      
      bPost <- paramsPost[postIdx, grep('b\\[', colnames(paramsPost))]
      knotsPost <- paramsPost[postIdx, grep('knots\\[', colnames(paramsPost))]
      trueVals <- c(betaPost, bPost, knotsPost, dataObs, RstarPost)
      
      
    } else if (alarmFit == 'gp') {
      
      logitAlarmPost <- logit(alarmSamples[,postIdx])[-1]
      names(logitAlarmPost) <- paste0('logit_', names(logitAlarmPost))
      trueVals <- c(betaPost, logitAlarmPost, dataObs, RstarPost)
      
    } else if (alarmFit == 'basic') {
      
      trueVals <- c(betaPost, dataObs, RstarPost)
      
    }
    
    # for exponential infectious period
    if (infPeriod == 'exp') {
      trueVals <- c(trueVals, rateIPost)
    }
    
    trueVals <- trueVals[parentNodes]
    
    
    postPredInc[,j] <- apply(sim_C$run(trueVals, 10), 2, median)[grep('Istar', dataNodes)]
  }
  
  postPredInc
}
