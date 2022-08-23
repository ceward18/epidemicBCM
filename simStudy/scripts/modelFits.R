################################################################################
# function to fit models 
# inputs: 
#   incData - observed incidence data
#   alarmFit - type of function to use to describe the alarm
#   smoothWindow - width of smoothing window 
#   prior - prior for rate parameter (1, 2, 3, or 4)
################################################################################


fitAlarmModel <- function(incData,  alarmFit, smoothWindow, prior, simNumber, seed) {
  
  source('./scripts/modelCodes.R')
  source('./scripts/getModelInputs.R')
  
  # get appropriate model code
  modelCode <- get(paste0('SIR_', alarmFit, '_exp'))
  
  # for reproducibility so inits are always the same
  set.seed(seed + simNumber)

  # model-specific constants, data, inits, and MCMC specs
  modelInputs <- getModelInput(alarmFit, incData, smoothWindow, prior)

  ### MCMC specifications
  niter <- modelInputs$niter
  nburn <- modelInputs$nburn
  nthin <- modelInputs$nthin
  
  ### create nimble model
  myModel <- nimbleModel(modelCode, 
                         data = modelInputs$dataList, 
                         constants = modelInputs$constantsList,
                         inits = modelInputs$initsList)
  myConfig <- configureMCMC(myModel)
  
  # need to ensure all stochastic nodes are monitored for WAIC calculation
  
  if (!alarmFit %in% c('betat', 'betatSpline', 'basic')) {
    myConfig$addMonitors(c('yAlarm', 'alarm'))
  }
  if (alarmFit %in% c('betat', 'betatSpline')) {
    myConfig$addMonitors(c('beta'))
  }
  
  # for exponential infectious period, need to use special proposal
  myConfig$removeSamplers('Rstar') # Nodes will be expanded
  myConfig$addSampler(target = c('Rstar'),
                      type = "RstarUpdate")
  myConfig$addMonitors(c('Rstar', 'R0'))
  
  myConfig$removeSampler('rateI')
  myConfig$addSampler(target = 'rateI', type = "slice")
  
  if (alarmFit == 'thresh') {
      
      # block sampler for transmission parameters
      paramsForBlock <- c('beta', 'delta', 'rateI')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock, type = "AF_slice")
      
  } else if (alarmFit == 'hill') {
      
      # block sampler for transmission parameters
      paramsForBlock <- c('beta', 'delta', 'nu', 'x0', 'rateI')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock, type = "AF_slice")
      
  } else if (alarmFit == 'power') {
      
      # block sampler for transmission parameters
      paramsForBlock <- c('beta', 'k', 'rateI')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock, type = "AF_slice")
      
  } else if (alarmFit == 'spline') {
      
      # block sampler for transmission parameters
      paramsForBlock <- c('b', 'knots')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock[1], type = "AF_slice")
      myConfig$addSampler(target = paramsForBlock[2], type = "AF_slice")
      
      # block sampler for transmission parameters
      paramsForBlock <- c('beta', 'rateI')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock, type = "AF_slice")
      
  } else if (alarmFit == 'gp') {
      
      # if gaussian process model, use slice sampling
      paramsForSlice <- c('l', 'sigma')
      myConfig$removeSampler(paramsForSlice)
      myConfig$addSampler(target = paramsForSlice[1], type = "slice")
      myConfig$addSampler(target = paramsForSlice[2], type = "slice")
      
      # block sampler for transmission parameters
      paramsForBlock <- c('beta', 'rateI')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock, type = "AF_slice")
      
  }
  
  myMCMC <- buildMCMC(myConfig)
  compiled <- compileNimble(myModel, myMCMC) 

  runMCMC(compiled$myMCMC, 
          niter = niter, 
          nburnin = nburn,
          thin = nthin,
          setSeed  = seed)

}
