################################################################################
# Calculate WAIC from matrix of posterior samples and model object
# need to recompile model outside of the parallel environment that was used 
#   to run multiple chains
# called from summarize post after samples have been combined across chains
################################################################################

getWAIC <- function(samples, incData, N, I0, R0, lengthI,
                    infPeriod, alarmFit, smoothWindow) {
    
    source('./scripts/modelCodes.R')

    # get appropriate model code
    modelCode <- get(paste0('SIR_', alarmFit, '_', infPeriod))
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(alarmFit, incData, infPeriod, smoothWindow, 
                                 N, I0, R0)

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
