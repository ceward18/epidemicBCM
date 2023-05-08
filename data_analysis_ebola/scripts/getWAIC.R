################################################################################
# Calculate WAIC from matrix of posterior samples and model object
# need to recompile model outside of the parallel environment that was used 
#   to run multiple chains
# called from summarize post after samples have been combined across chains
################################################################################

# use my own function to calculate likelihood

getLL <- nimbleFunction(   
    setup = function(model, dataNodes) {
        
        parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
        
    }, 
    
    run = function(samples = double(2), tau = double(0)) {
        returnType(double(2))
        
        nSamples <- dim(samples)[1]
        nParams <- dim(samples)[2]
        ll <- matrix(0, nrow = nSamples, ncol = tau)
        
        for (i in 1:nSamples) {
            
            values(model, parentNodes) <<- samples[i, 1:nParams] 
            model$calculate()
            
            for (j in 1:tau) {
                
                ll[i, j] <- model$logProb_Estar[j] + model$logProb_Istar[j] + model$logProb_Rstar[j]
                if (ll[i, j] == -Inf) ll[i, j] <- 0
            }
        }
        
        return(ll)
    })



getWAICFun <- function(ll) {
    
    pWAIC2 <- sum(apply(ll, 1, var))
    
    lppd <- sum(log(rowMeans(exp(ll))))
    
    -2 * (lppd - pWAIC2)
    
}

getWAIC <- function(samples, incData, deathData, smoothI,
                    N, E0, I0, R0, intTime, alarmFit) {
    
    # get appropriate model code
    modelCode <- get(paste0('SEIR_', alarmFit))
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(alarmFit = alarmFit, 
                                 incData = incData, deathData = deathData,
                                 smoothI = smoothI,
                                 N = N, E0 = E0, I0 = I0, R0 = R0, intTime = intTime)
    
    ### create nimble model
    myModel <- nimbleModel(modelCode, 
                           data = modelInputs$dataList, 
                           constants = modelInputs$constantsList,
                           inits = modelInputs$initsList)
    
    compiled <- compileNimble(myModel) 
   
    parentNodes <- myModel$getParents(c('Estar', 'Istar', 'Rstar'), stochOnly = TRUE)
    parentNodes <- myModel$expandNodeNames(parentNodes, returnScalarComponents = TRUE)
    
    samples <- samples[,parentNodes]
    
    ll_R <- getLL(myModel, c('Estar', 'Istar', 'Rstar'))
    ll_C <- compileNimble(ll_R)
    
    ll_IterTime <- ll_C$run(samples, modelInputs$constantsList$tau)
    
    
    data.frame(waic = getWAICFun(ll_IterTime))
    
}
