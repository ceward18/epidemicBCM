################################################################################
# function to fit models 
# inputs: 
#   incData - observed incidence data
#   infPeriod - type of infectious period (fixed or exponential)
#   alarmFit - type of function to use to describe the alarm
#   smoothWindow - width of smoothing window 
################################################################################

fitAlarmModel <- function(incData, N, I0, R0, Rstar0, lengthI,
                          infPeriod, alarmFit, smoothWindow, seed) {
  
  source('../scripts/modelCodes.R')
  source('../scripts/getModelInputs.R')
  
  # get appropriate model code
  modelCode <- get(paste0('SIR_', alarmFit, '_', infPeriod))
  
  # for reproducibility so inits are always the same
  set.seed(seed + 3)

  # model-specific constants, data, and inits
  modelInputs <- getModelInput(alarmFit, incData, infPeriod, smoothWindow, 
                               N, I0, R0, Rstar0, lengthI)
  
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
  
  # if gaussian process model, use slice sampling
  if (alarmFit == 'gp') {
    
    paramsForSlice <- c('beta', 'l', 'sigma')
    myConfig$removeSampler(paramsForSlice)
    myConfig$addSampler(target = paramsForSlice[1], type = "slice")
    myConfig$addSampler(target = paramsForSlice[2], type = "slice")
    myConfig$addSampler(target = paramsForSlice[3], type = "slice")
    
  } else if (alarmFit == 'thresh') {
    
    # block sampler for transmission parameters
    paramsForBlock <- c('beta', 'delta')
    myConfig$removeSampler(paramsForBlock)
    myConfig$addSampler(target = paramsForBlock, type = "RW_block",
                        control = list(adaptInterval = 100,
                                       propCov = diag(c(0.1, 0.1))))
    
    myConfig$removeSampler('H')
    myConfig$addSampler(target = 'H', type = "slice")
    
  } else if (alarmFit == 'hill') {
    
    # block sampler for transmission parameters
    paramsForBlock <- c('beta', 'delta', 'nu', 'x0')
    myConfig$removeSampler(paramsForBlock)
    myConfig$addSampler(target = paramsForBlock, type = "RW_block",
                        control = list(adaptInterval = 100,
                                       propCov = diag(c(0.2, 0.2, 1, 800))))
    
  }
  myConfig$addMonitors(c('Rstar'))
  # browser()
  
  myMCMC <- buildMCMC(myConfig)
  compiled <- compileNimble(myModel, myMCMC) 
  

  
  runMCMC(compiled$myMCMC, 
          niter = niter, 
          nburnin = nburn,
          thin = nthin,
          setSeed  = seed)
  
}
