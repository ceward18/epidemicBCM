################################################################################
# Calculate WAIC from matrix of posterior samples and model object
# need to recompile model outside of the parallel environment that was used 
#   to run multiple chains
# called from summarize post after samples have been combined across chains
################################################################################

getWAIC <- function(samples, incData, smoothI, N, I0, R0, Rstar0, lengthI,
                    prior, peak, alarmFit) {

    # get appropriate model code
    modelCode <- get(paste0('SIR_', alarmFit))
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(alarmFit = alarmFit, 
                                 incData = incData, smoothI = smoothI,
                                 prior = prior, peak = peak,
                                 N = N, I0 = I0, R0 = R0,
                                 Rstar0 = Rstar0, lengthI = lengthI)

    ### create nimble model
    myModel <- nimbleModel(modelCode, 
                           data = modelInputs$dataList, 
                           constants = modelInputs$constantsList,
                           inits = modelInputs$initsList)
    
    compiled <- compileNimble(myModel) 
    
    
    incDataSamples <- matrix(rep(incData), NROW(samples),
                             ncol = length(incData), byrow = T)
    colnames(incDataSamples) <- paste0('Istar[', 1:ncol(incDataSamples), ']')
   
    samplesIstar <- cbind(samples, 
                          incDataSamples)
    
    waicList <- calculateWAIC(samplesIstar, compiled)
 
    data.frame(waic = waicList$WAIC,
               lppd = waicList$lppd,
               pWAIC = waicList$pWAIC)
    
}
