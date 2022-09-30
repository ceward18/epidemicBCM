################################################################################
# function to do posterior predictive distribution of observed epidemic curve 
################################################################################


postPredFit <- function(incData, smoothI, N, I0, R0, Rstar0, lengthI,
                        alarmFit, prior, peak, smoothWindow, 
                        paramsPost, alarmSamples) {
    
    # model code
    if (!alarmFit %in% c('basic', 'betatSpline')) {
        modelCode <- get(paste0('SIR_', alarmFit, '_sim'))
    } else {
        modelCode <- get(paste0('SIR_', alarmFit))
    }
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(alarmFit = alarmFit, 
                                 incData = incData, smoothI = smoothI, 
                                 prior = prior, peak = peak,
                                 smoothWindow = smoothWindow,
                                 N = N, I0 = I0, R0 = R0,
                                 Rstar0 = Rstar0, lengthI = lengthI)
    
    modelInputs$constantsList$bw <- smoothWindow
    
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
    dataNodes <- c(dataNodes, paste0('Rstar[', 1:tau, ']'))
    
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
        
        # for exponential infectious period
        rateIPost <- paramsPost[postIdx, 'rateI']
        trueVals <- c(trueVals, rateIPost)
        
        # initial conditions 
        SIR_initPost <- paramsPost[postIdx, grep('SIR_init', colnames(paramsPost))]
        SIR_initPost[3] <- N - sum(SIR_initPost)
        names(SIR_initPost)[3] <- "SIR_init[3]"
        trueVals <- c(trueVals, SIR_initPost)
        
        trueVals <- trueVals[parentNodes]
       
        postPredInc[,j] <- apply(sim_C$run(trueVals, 10), 2, median)[grep('Istar', dataNodes)]
    }
    
    postPredInc
}
