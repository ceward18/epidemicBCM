################################################################################
# function to do posterior prediction on remainder of the epidemic 
################################################################################


postPred <- function(incData, N, I0, R0, lengthI,
                     alarmFit, infPeriod, smoothWindow, 
                     paramsPost, alarmSamples) {
  
  # want to predict out 100 days from last observed time point
  nDaysPred <- 100
  obsTime <- length(incData)
  fullTime <- obsTime + nDaysPred
  predTime <- (obsTime + 1):fullTime
  dataObs <- incData
  
  # constants that are the same for all models
  S0 <- N - I0 - R0
  tau <- length(incData) + nDaysPred
  
  # model code
  modelCode <- get(paste0('SIR_', alarmFit, '_', infPeriod))
  
  if (alarmFit == 'thresh') {
    
    n <- 50
    maxI <- ceiling(max(movingAverage(incData, smoothWindow)))
    xAlarm <- seq(0, maxI, length.out = n)
    
    constantsList <- list(tau = tau,
                          N = N,
                          S0 = S0,
                          I0 = I0,
                          probRstar = rep(1/lengthI, lengthI),
                          bw = smoothWindow,
                          lengthI = lengthI,
                          n = n,
                          xAlarm = xAlarm)
    
  } else if (alarmFit == 'hill') {
    
    n <- 50
    maxI <- ceiling(max(movingAverage(incData, smoothWindow)))
    xAlarm <- seq(0, maxI, length.out = n)
    
    constantsList <- list(tau = tau,
                          N = N,
                          S0 = S0,
                          I0 = I0,
                          probRstar = rep(1/lengthI, lengthI),
                          bw = smoothWindow,
                          lengthI = lengthI,
                          n = n,
                          xAlarm = xAlarm,
                          maxI = maxI)
    
  } else if (alarmFit == 'power') {
    
    n <- 50
    maxI <- ceiling(max(movingAverage(incData, smoothWindow)))
    xAlarm <- seq(0, maxI, length.out = n)
    
    constantsList <- list(tau = tau,
                          N = N,
                          S0 = S0,
                          I0 = I0,
                          probRstar = rep(1/lengthI, lengthI),
                          bw = smoothWindow,
                          lengthI = lengthI,
                          n = n,
                          xAlarm = xAlarm)
    
  } else if (alarmFit == 'spline') {
    
    n <- 50
    maxI <- ceiling(max(movingAverage(incData, smoothWindow)))
    xAlarm <- seq(0, maxI, length.out = n)
    nb <- 3
    
    constantsList <- list(tau = tau,
                          N = N,
                          S0 = S0,
                          I0 = I0,
                          probRstar = rep(1/lengthI, lengthI),
                          bw = smoothWindow,
                          xAlarm = xAlarm,
                          n = n,
                          maxI = maxI,
                          lengthI = lengthI,
                          nb = nb)
    
    initsList <- list(b = c(1, 1, 1), 
                      knots = quantile(1:maxI, probs = c(0.3, 0.6)), 
                      beta = 0.5)
    
  } else if (alarmFit == 'gp') {
    
    n <- 10
    maxI <- ceiling(max(movingAverage(incData, smoothWindow)))
    xAlarm <- seq(0, maxI, length.out = n)
    distMat <- as.matrix(dist(matrix(xAlarm)))
    
    vals <- c(1, 1) # doesn't actually matter here
    constantsList <- list(tau = tau,
                          N = N,
                          S0 = S0,
                          I0 = I0,
                          probRstar = rep(1/lengthI, lengthI),
                          bw = smoothWindow,
                          lengthI = lengthI,
                          dists = distMat,
                          mu0 = 1,
                          ones = logit(seq(0.0001, 0.9999, length.out= n)),
                          n = n,
                          xAlarm = xAlarm,
                          c = vals[1],
                          d = vals[2])
    
  } else if (alarmFit == 'basic') {
    
    constantsList <- list(tau = tau,
                          N = N,
                          S0 = S0,
                          I0 = I0,
                          probRstar = rep(1/lengthI, lengthI),
                          bw = smoothWindow,
                          lengthI = lengthI)
    
  }
  
  # true data and fake data during prediction interval just for initialization
  dataList <- list(Istar = c(incData, 
                             rep(1, nDaysPred)))
  
  
  # compile model and simulator
  if (alarmFit == 'spline') {
    # spline model needs initial values to avoid warning from NA knots
    myModelPred <- nimbleModel(modelCode, 
                               data = dataList, 
                               constants = constantsList,
                               inits = initsList)
  } else {
    myModelPred <- nimbleModel(modelCode, 
                               data = dataList, 
                               constants = constantsList)
  }
  
  compiledPred  <- compileNimble(myModelPred) 
  dataNodes <- paste0('Istar[', predTime, ']')
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
      trueVals <- c(betaPost, logitAlarmPost, dataObs, RstarPost)
      
    } else if (alarmFit == 'basic') {
      
      trueVals <- c(betaPost, dataObs, RstarPost)
      
    }
    
    trueVals <- trueVals[parentNodes]
    
    
    postPredInc[,j] <- apply(sim_C$run(trueVals, 10), 2, median)
  }
  
  postPredInc
}

# 
# parentNodes <- myModelPred$getParents(dataNodes, stochOnly = TRUE)
# # exclude data from parent nodes
# parentNodes <- parentNodes[-which(parentNodes %in% dataNodes)]
# parentNodes <- myModelPred$expandNodeNames(parentNodes, returnScalarComponents = TRUE)
# cat("Stochastic parents of data are: ", paste(parentNodes, sep = ','), ".\n")
# simNodes <- myModelPred$getDependencies(parentNodes, self = FALSE,
#                                   downstream = T)
# 
# values(myModelPred, parentNodes) <- trueVals
# myModelPred$simulate(simNodes, includeData = TRUE)
# values(myModelPred, dataNodes)
# 
