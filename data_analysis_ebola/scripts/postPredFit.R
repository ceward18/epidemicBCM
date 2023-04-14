################################################################################
# function to do posterior predictive distribution of observed epidemic curve 
################################################################################


postPredFit <- function(incData, deathData, smoothI, smoothWindow,
                        N, E0, I0, R0, 
                        alarmFit, paramsPost, alarmSamples) {
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(alarmFit = alarmFit, 
                                 incData = incData, deathData = deathData,
                                 smoothI = smoothI,
                                 N = N, E0 = E0, I0 = I0, R0 = R0)
    
    # model code
    if (!alarmFit %in% c('basic', 'betatSpline')) {
        modelCode <- get(paste0('SEIR_', alarmFit, '_sim'))
        modelInputs$constantsList$bw <- smoothWindow
    } else {
        modelCode <- get(paste0('SEIR_', alarmFit))
    }
    
    # compile model and simulator
    if (alarmFit %in% c('spline', 'splineFixKnot', 'betatSpline')) {
        # spline model needs initial values to avoid warning from NA knots
        myModelPred <- nimbleModel(modelCode, 
                                   constants = modelInputs$constantsList,
                                   inits = modelInputs$initsList)
    } else {
        myModelPred <- nimbleModel(modelCode, 
                                   constants = modelInputs$constantsList)
    }
    
    compiledPred  <- compileNimble(myModelPred) 
    
    tau <- modelInputs$constantsList$tau
    dataNodes <- paste0('Istar[', 1:tau, ']')
    dataNodes <- c(dataNodes, paste0('Estar[', 1:tau, ']'), paste0('Rstar[', 1:tau, ']'))
    
    sim_R <- simulator(myModelPred, dataNodes)
    sim_C <- compileNimble(sim_R)
    
    # get order of parameters
    parentNodes <- myModelPred$getParents(dataNodes, stochOnly = TRUE)
    parentNodes <- parentNodes[-which(parentNodes %in% dataNodes)]
    parentNodes <- myModelPred$expandNodeNames(parentNodes, returnScalarComponents = TRUE)
    
    nPost <- 10000
    postPredInc <- matrix(NA, nrow = tau, ncol = nPost)
    set.seed(1)
    for (j in 1:nPost) {
        
        postIdx <- sample(1:nrow(paramsPost), 1)
        
        if (alarmFit != 'betatSpline') {
            betaPost <- paramsPost[postIdx,'beta']
        }
        
        # model specific parameters
        if (alarmFit == 'power') {
            
            alarmParamPost <- paramsPost[postIdx, 'k']
            trueVals <- c(betaPost, alarmParamPost)
            
            
        } else if (alarmFit == 'thresh') {
            
            alarmParamPost <- paramsPost[postIdx, c('delta', 'H')]
            trueVals <- c(betaPost, alarmParamPost)
            
            
        } else if (alarmFit == 'hill') {
            
            alarmParamPost <- paramsPost[postIdx, c('delta', 'nu', 'x0')]
            trueVals <- c(betaPost, alarmParamPost)
            
            
        } else if (alarmFit == 'spline') {
            
            bPost <- paramsPost[postIdx, grep('b\\[', colnames(paramsPost))]
            knotsPost <- paramsPost[postIdx, grep('knots\\[', colnames(paramsPost))]
            trueVals <- c(betaPost, bPost, knotsPost)
            
            
        }  else if (alarmFit == 'splineFixKnot') {
            
            bPost <- paramsPost[postIdx, grep('b\\[', colnames(paramsPost))]
            trueVals <- c(betaPost, bPost)
            
            
        } else if (alarmFit == 'gp') {
            
            logitAlarmPost <- logit(alarmSamples[,postIdx])[-1]
            names(logitAlarmPost) <- paste0('logit_', names(logitAlarmPost))
            trueVals <- c(betaPost, logitAlarmPost)
            
        } else if (alarmFit == 'basic') {
            
            trueVals <- c(betaPost)
            
        } else if (alarmFit == 'betatSpline') {
            
            bPost <- paramsPost[postIdx, grep('b\\[', colnames(paramsPost))]
            knotsPost <- paramsPost[postIdx, grep('knots\\[', colnames(paramsPost))]
            trueVals <- c(bPost, knotsPost)
            
        }
        
        # for exponential exposure period
        rateEPost <- paramsPost[postIdx, 'rateE']
        trueVals <- c(trueVals, rateEPost)
        
        # for exponential infectious period
        rateIPost <- paramsPost[postIdx, 'rateI']
        trueVals <- c(trueVals, rateIPost)
        
        
        trueVals <- trueVals[parentNodes]
       
        postPredInc[,j] <- apply(sim_C$run(trueVals, 10), 2, median)[grep('Istar', dataNodes)]
    }
    
    postPredInc
}
