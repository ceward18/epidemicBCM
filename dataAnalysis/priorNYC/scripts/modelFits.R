################################################################################
# function to fit models 
# inputs: 
#   incData - observed incidence data
#   infPeriod - type of infectious period (fixed or exponential)
#   alarmFit - type of function to use to describe the alarm
################################################################################

fitAlarmModel <- function(incData, smoothI, N, I0, R0, Rstar0, lengthI,
                          infPeriod, prior, alarmFit, seed) {
  
  source('./scripts/modelCodes.R')
  source('./scripts/getModelInputs.R')
  
  # get appropriate model code
  modelCode <- get(paste0('SIR_', alarmFit, '_', infPeriod))
  
  # for reproducibility so inits are always the same
  set.seed(seed + 3)

  # model-specific constants, data, and inits
  modelInputs <- getModelInput(alarmFit = alarmFit, 
                               incData = incData, smoothI = smoothI,
                               infPeriod = infPeriod, prior = prior,
                               N = N, I0 = I0, R0 = R0, Rstar0 = Rstar0, 
                               lengthI = lengthI)
  
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
  
  if (!alarmFit %in% c('betatSpline', 'basic')) {
    myConfig$addMonitors(c('yAlarm', 'alarm'))
  }
  if (alarmFit == 'betatSpline') {
    myConfig$addMonitors(c('beta'))
  }
  
  # if exponential infectious period, need to use special proposal
  if (infPeriod == 'exp') {
    
    myConfig$removeSamplers('Rstar') # Nodes will be expanded
    myConfig$addSampler(target = 'Rstar',
                        type = "RstarUpdate",
                        control = list(ignoreIdx = 1:length(Rstar0)))
    
    myConfig$removeSampler('rateI')
    myConfig$addSampler(target = 'rateI', type = "slice")
  }
  
  # customize samplers depending on model being fitted
  if (alarmFit == 'thresh') {
    
    # block sampler for transmission parameters
    paramsForBlock <- c('beta', 'delta')
    myConfig$removeSampler(paramsForBlock)
    myConfig$addSampler(target = paramsForBlock, type = "AF_slice")
    
  } else if (alarmFit == 'hill') {
    
    # block sampler for transmission parameters
    paramsForBlock <- c('beta', 'delta', 'nu', 'x0')
    myConfig$removeSampler(paramsForBlock)
    myConfig$addSampler(target = paramsForBlock, type = "AF_slice")
    
  } else if (alarmFit == 'power') {
      
      # block sampler for transmission parameters
      paramsForBlock <- c('beta', 'k')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock, type = "AF_slice")
      
  } else if (alarmFit == 'spline') {
      
      # block sampler for transmission parameters
      paramsForBlock <- c('beta', 'b', 'knots')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock[1:2], type = "AF_slice")
      myConfig$addSampler(target = paramsForBlock[3], type = "AF_slice")
     
  } else if (alarmFit == 'gp') {
      
      # if gaussian process model, use slice sampling
      paramsForSlice <- c('beta', 'l', 'sigma')
      myConfig$removeSampler(paramsForSlice)
      myConfig$addSampler(target = paramsForSlice[1], type = "slice")
      myConfig$addSampler(target = paramsForSlice[2], type = "slice")
      myConfig$addSampler(target = paramsForSlice[3], type = "slice")
      
  }  else if (alarmFit == 'betatSpline') {
      
      # block sampler for transmission parameters
      paramsForBlock <- c('b', 'knots')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock[1], type = "AF_slice")
      myConfig$addSampler(target = paramsForBlock[2], type = "AF_slice")
      
  } 
  
  myConfig$addMonitors(c('Rstar', 'R0'))
  myMCMC <- buildMCMC(myConfig)
  compiled <- compileNimble(myModel, myMCMC) 
  
  runMCMC(compiled$myMCMC, 
          niter = niter, 
          nburnin = nburn,
          thin = nthin,
          setSeed  = seed)
  
}
