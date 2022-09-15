################################################################################
# function to fit models 
# inputs: 
#   incData - observed incidence data
#   alarmFit - type of function to use to describe the alarm
################################################################################

fitAlarmModel <- function(incData, smoothI, N, I0, R0, Rstar0, lengthI,
                          prior, peak, alarmFit, seed) {
  
  source('./scripts/modelCodes.R')
  source('./scripts/getModelInputs.R')
  
  # get appropriate model code
  modelCode <- get(paste0('SIR_', alarmFit))
  
  # for reproducibility so inits are always the same
  set.seed(seed + 3)

  # model-specific constants, data, and inits
  modelInputs <- getModelInput(alarmFit = alarmFit, 
                               incData = incData, smoothI = smoothI,
                               prior = prior, peak = peak,
                               N = N, I0 = I0, R0 = R0,
                               Rstar0 = Rstar0, lengthI = lengthI)
  
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
  myConfig$removeSamplers('Rstar') # Nodes will be expanded
  myConfig$addSampler(target = 'Rstar',
                      type = "RstarUpdate")
  
  # customize samplers depending on model being fitted
  if (alarmFit == 'power') {
      
      # block sampler for beta, rateI parameters
      paramsForBlock <- c('beta', 'rateI')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock, 
                          type = "RW_block",
                          propCov = diag(0.4, length(paramsForBlock)))
      
      # slice sampler for k
      paramsForSlice <- c('k')
      myConfig$removeSampler(paramsForSlice)
      myConfig$addSampler(target = paramsForSlice[1], type = "slice")
      
  } else if (alarmFit == 'thresh') {
    
    # block sampler for transmission parameters
    paramsForBlock <- c('beta', 'delta', 'rateI')
    myConfig$removeSampler(paramsForBlock)
    myConfig$addSampler(target = paramsForBlock, 
                        type = "RW_block",
                        propCov = diag(0.4, length(paramsForBlock)))
    
    
  } else if (alarmFit == 'hill') {
    
    # block sampler for transmission parameters
    paramsForBlock <- c('beta', 'delta', 'rateI')
    myConfig$removeSampler(paramsForBlock)
    myConfig$addSampler(target = paramsForBlock, 
                        type = "RW_block",
                        propCov = diag(0.4, length(paramsForBlock)))
    
    # slice samplers for nu and x0
    paramsForSlice <- c('nu', 'x0')
    myConfig$removeSampler(paramsForSlice)
    myConfig$addSampler(target = paramsForSlice[1], type = "slice")
    myConfig$addSampler(target = paramsForSlice[2], type = "slice")
    
  } else if (alarmFit == 'spline') {
      
      # block sampler for transmission parameters
      paramsForBlock <- c('beta', 'b', 'rateI')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock, 
                          type = "RW_block",
                          propCov = diag(0.4, 
                                         length(myModel$expandNodeNames(paramsForBlock))))
      
      # block slice sampler for knots
      paramsForBlock <- c('knots')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock, 
                          type = "AF_slice")
     
  } else if (alarmFit == 'splineFixKnot') {
      
      # block sampler for transmission parameters
      paramsForBlock <- c('beta', 'b', 'rateI')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock, 
                          type = "RW_block",
                          propCov = diag(0.4, 
                                         length(myModel$expandNodeNames(paramsForBlock))))
      
      
  } else if (alarmFit == 'gp') {
      
      # if gaussian process model, use slice sampling for GP parameters
      paramsForSlice <- c('l', 'sigma')
      myConfig$removeSampler(paramsForSlice)
      myConfig$addSampler(target = paramsForSlice[1], type = "slice")
      myConfig$addSampler(target = paramsForSlice[2], type = "slice")
      
      # sample beta and rateI jointly
      paramsForBlock <- c('beta', 'rateI')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock, 
                          type = "RW_block",
                          propCov = diag(0.4, length(paramsForBlock)))
      
      
  } else if (alarmFit == 'betatSpline') {
      
      # block sampler for transmission parameters
      paramsForBlock <- c('b', 'rateI')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock, 
                          type = "RW_block",
                          propCov = diag(0.4, 
                                         length(myModel$expandNodeNames(paramsForBlock))))
      
      # block slice sampler for knots
      paramsForBlock <- c('knots')
      myConfig$removeSampler(paramsForBlock)
      myConfig$addSampler(target = paramsForBlock, 
                          type = "AF_slice")
      
  } 
  
  myConfig$addMonitors(c('Rstar', 'R0'))
  print(myConfig)
  myMCMC <- buildMCMC(myConfig)
  compiled <- compileNimble(myModel, myMCMC) 
  
  runMCMC(compiled$myMCMC, 
          niter = niter, 
          nburnin = nburn,
          thin = nthin,
          setSeed  = seed)
  
}
