################################################################################
# function to fit models 
# inputs: 
#   incData - observed incidence data
#   alarmFit - type of function to use to describe the alarm
################################################################################

fitAlarmModel <- function(incData, deathData, smoothI, N, E0, I0, R0, 
                          smoothWindow, alarmFit, seed) {
    
    source('./scripts/modelCodes.R')
    source('./scripts/getModelInputs.R')
    
    # get appropriate model code
    modelCode <- get(paste0('SEIR_', alarmFit))
    
    # for reproducibility so inits are always the same
    set.seed(seed + 3)
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(alarmFit = alarmFit, 
                                 incData = incData, deathData = deathData,
                                 smoothI = smoothI,
                                 N = N, E0 = E0, I0 = I0, R0 = R0)
    
    
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
    
    # proposal for Estar
    myConfig$removeSamplers('Estar') # Nodes will be expanded
    myConfig$addSampler(target = 'Estar',
                        type = "EstarUpdate")
    
    # proposal for Istar
    myConfig$removeSamplers('Istar') # Nodes will be expanded
    myConfig$addSampler(target = 'Istar',
                        type = "transUpdate",
                        control = list(fixed = incData))
    
    # proposal for Rstar
    myConfig$removeSamplers('Rstar') # Nodes will be expanded
    myConfig$addSampler(target = 'Rstar',
                        type = "transUpdate",
                        control = list(fixed = deathData))
    
    # customize samplers depending on model being fitted
    if (alarmFit == 'power') {
        
        # block sampler for beta, rateI parameters
        paramsForBlock <- c('beta', 'rateE', 'rateI')
        myConfig$removeSampler(paramsForBlock)
        myConfig$addSampler(target = paramsForBlock, 
                            type = "AF_slice")
        
        # slice sampler for k
        paramsForSlice <- c('k')
        myConfig$removeSampler(paramsForSlice)
        myConfig$addSampler(target = paramsForSlice[1], type = "slice")
        
    } else if (alarmFit == 'thresh') {
        
        # block sampler for beta, rateI parameters
        paramsForBlock <- c('beta', 'rateE', 'rateI')
        myConfig$removeSampler(paramsForBlock)
        myConfig$addSampler(target = paramsForBlock, 
                            type = "AF_slice")
        
        # slice sampler for delta and H
        paramsForSlice <- c('delta', 'H')
        myConfig$removeSampler(paramsForSlice)
        myConfig$addSampler(target = paramsForSlice[1], type = "slice")
        myConfig$addSampler(target = paramsForSlice[2], type = "slice")
        
        
    } else if (alarmFit == 'hill') {
        
        # block sampler for beta, rateI parameters
        paramsForBlock <- c('beta', 'rateE', 'rateI')
        myConfig$removeSampler(paramsForBlock)
        myConfig$addSampler(target = paramsForBlock, 
                            type = "AF_slice")
        
        # slice samplers for delta, nu and x0
        paramsForSlice <- c('delta', 'nu', 'x0')
        myConfig$removeSampler(paramsForSlice)
        myConfig$addSampler(target = paramsForSlice[1], type = "slice")
        myConfig$addSampler(target = paramsForSlice[2], type = "slice")
        myConfig$addSampler(target = paramsForSlice[3], type = "slice")
        
    } else if (alarmFit == 'spline') {
        
        # block sampler for transmission parameters
        paramsForBlock <- c('beta', 'b', 'rateE', 'rateI')
        myConfig$removeSampler(paramsForBlock)
        myConfig$addSampler(target = paramsForBlock, 
                            type = "AF_slice")
        
        # block slice sampler for knots
        paramsForBlock <- c('knots')
        myConfig$removeSampler(paramsForBlock)
        myConfig$addSampler(target = paramsForBlock, 
                            type = "AF_slice")
        
    } else if (alarmFit == 'splineFixKnot') {
        
        # block sampler for transmission parameters
        paramsForBlock <- c('beta', 'b', 'rateE', 'rateI')
        myConfig$removeSampler(paramsForBlock)
        myConfig$addSampler(target = paramsForBlock, 
                            type = "AF_slice")
        
        
    } else if (alarmFit == 'gp') {
        
        # if gaussian process model, use slice sampling for GP parameters
        paramsForSlice <- c('l', 'sigma')
        myConfig$removeSampler(paramsForSlice)
        myConfig$addSampler(target = paramsForSlice[1], type = "slice")
        myConfig$addSampler(target = paramsForSlice[2], type = "slice")
        
        # sample beta and rateI jointly
        paramsForBlock <- c('beta', 'rateE', 'rateI')
        myConfig$removeSampler(paramsForBlock)
        myConfig$addSampler(target = paramsForBlock, 
                            type = "AF_slice")
        
        # multivariate normal sampler for logit_yAlarm
        myConfig$addSampler(target = 'logit_yAlarm[2:10]',
                            type = 'RW_block')
        
        
    } else if (alarmFit == 'betatSpline') {
        
        # block sampler for transmission parameters
        paramsForBlock <- c('b', 'rateE', 'rateI')
        myConfig$removeSampler(paramsForBlock)
        myConfig$addSampler(target = paramsForBlock, 
                            type = "AF_slice")
        
        # block slice sampler for knots
        paramsForBlock <- c('knots')
        myConfig$removeSampler(paramsForBlock)
        myConfig$addSampler(target = paramsForBlock, 
                            type = "AF_slice")
        
    } else if (alarmFit == 'basic') {
        
        # block sampler for transmission parameters
        paramsForBlock <- c('beta', 'rateE', 'rateI')
        myConfig$removeSampler(paramsForBlock)
        myConfig$addSampler(target = paramsForBlock, 
                            type = "AF_slice")
        
        
    }
    
    myConfig$addMonitors(c('Estar', 'Istar', 'Rstar', 'R0_eff'))
    print(myConfig)
    print(myConfig$getUnsampledNodes())
    nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
    myMCMC <- buildMCMC(myConfig)
    compiled <- compileNimble(myModel, myMCMC) 
    
    runMCMC(compiled$myMCMC, 
            niter = niter, 
            nburnin = nburn,
            thin = nthin,
            setSeed  = seed)
    
}
