################################################################################
# Calculate WAIC from matrix of posterior samples and model object
# need to recompile model outside of the parallel environment that was used 
#   to run multiple chains
# called from summarize post after samples have been combined across chains
################################################################################

getWAIC <- function(samples, incData, N, I0, R0, lengthI,
                    infPeriod, alarmFit, smoothWindow) {
    
    source('./scripts/modelCodes.R')
    
    # constants that are the same for all models
    S0 <- N - I0 - R0
    tau <- length(incData)
    
    # get appropriate model code
    modelCode <- get(paste0('SIR_', alarmFit, '_', infPeriod))
    
    # model-specific constants, data, and inits
    
    if (alarmFit == 'thresh') {
        
        ### constants
        n <- 50
        maxI <- ceiling(max(movingAverage(incData, smoothWindow)) )
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              probRstar = rep(1/lengthI, lengthI),
                              bw = smoothWindow,
                              lengthI = lengthI,
                              n = n,
                              xAlarm = xAlarm)
        
        ### data
        dataList <- list(Istar = incData)
        
        ### inits
        initsList <- list(beta = runif(1, 0, 1),
                          delta = runif(1, 0, 1),
                          H = runif(1, 0, maxI))
        
    } else if (alarmFit == 'hill') {
        
        ### constants
        n <- 50
        maxI <- ceiling(max(movingAverage(incData, smoothWindow)) )
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              probRstar = rep(1/lengthI, lengthI),
                              bw = smoothWindow,
                              lengthI = lengthI,
                              n = n,
                              xAlarm = xAlarm,
                              maxI = maxI)
        
        ### data
        dataList <- list(Istar = incData)
        
        ### inits
        initsList <- list(beta = runif(1, 0, 1),
                          delta = runif(1, 0, 1),
                          nu = runif(1, 0, 10),
                          x0 = max(rnorm(1, maxI/2, 10), 1))
        
    } else if (alarmFit == 'power') {
        
        ### constants 
        n <- 50
        maxI <- ceiling(max(movingAverage(incData, smoothWindow)) )
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              probRstar = rep(1/lengthI, lengthI),
                              bw = smoothWindow,
                              lengthI = lengthI,
                              n = n,
                              xAlarm = xAlarm)
        
        ### data
        dataList <- list(Istar = incData)
        
        ### inits
        initsList <- list(beta = runif(1, 0, 1),
                          k = runif(1, 0, 1))
        
    } else if (alarmFit == 'spline') {
        
        ### constants
        n <- 50
        maxI <- ceiling(max(movingAverage(incData, smoothWindow)))
        xAlarm <- seq(0, maxI, length.out = n)
        nb <- 3
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              probRstar = rep(1/lengthI, lengthI),
                              bw = smoothWindow,
                              xAlarm = xAlarm,
                              n = n,
                              maxI = maxI,
                              lengthI = lengthI,
                              nb = nb)
        
        ### data
        dataList <- list(Istar = incData,
                         constrain_knots = 1,
                         constrain_min = 1,
                         constrain_max = 1)
        
        ### inits (must satisfy constraint)
        repeat {
            initsList <- list(beta = runif(1, 0, 1),
                              b = rnorm(nb, 0, 4),
                              knots = as.vector(quantile(xAlarm, 
                                                         probs = sort(runif(nb - 1, 0.1, 0.8)))))
            
            cond <- all(splineAlarm(xAlarm, initsList$b, initsList$knots) >= 0) & 
                all(splineAlarm(xAlarm, initsList$b, initsList$knots) <= 1)
            
            if (cond) break
        }
        
    } else if (alarmFit == 'gp') {
        
        ### constants
        n <- 10
        maxI <- ceiling(max(movingAverage(incData, smoothWindow)))
        xAlarm <- seq(0, maxI, length.out = n)
        distMat <- as.matrix(dist(matrix(xAlarm)))
        
        uniqueDists <- distMat[lower.tri(distMat)]
        minDist <- min(uniqueDists)
        maxDist <- max(uniqueDists)
        midDist <- getl(maxDist)
        
        # parameters of inverse gamma distribution for prior on lengthscale
        vals <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
                            min = minDist, mid = midDist, max = maxDist)$par, 2)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              probRstar = rep(1/lengthI, lengthI),
                              bw = smoothWindow,
                              lengthI = lengthI,
                              dists = distMat,
                              mu0 = 1,
                              ones = logit(seq(0.0001, 0.9999, length.out= n)),
                              n = n,
                              xAlarm = xAlarm,
                              c = vals[1],
                              d = vals[2])
        
        ### data
        dataList <- list(Istar = incData)
        
        ### inits 
        initsList <- list(beta = runif(1, 0, 1),
                          l = rinvgamma(1, vals[1], vals[2]),
                          sigma = rgamma(1, 100, 50))
    
        
    } else if (alarmFit == 'betat') {
        
        ### set up grid for gaussian process of beta over epidemic time
        distMat <- as.matrix(dist(matrix(1:tau)))
        
        uniqueDists <- distMat[lower.tri(distMat)]
        minDist <- min(uniqueDists)
        maxDist <- max(uniqueDists)
        midDist <- getl(maxDist)
        
        
        vals <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
                            min = minDist, mid = midDist, max = maxDist)$par, 2)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              probRstar = rep(1/lengthI, lengthI),
                              lengthI = lengthI,
                              dists = distMat,
                              zeros = rep(0, tau),
                              c = vals[1],
                              d = vals[2])
        
        ### data
        dataList <- list(Istar = incData)
        
        ### inits 
        initsList <- list(l = rinvgamma(1, vals[1], vals[2]),
                          sigma = rgamma(1, 100, 50))
        
        
    } else if (alarmFit == 'basic') {
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              probRstar = rep(1/lengthI, lengthI),
                              lengthI = lengthI)
        
        ### data
        dataList <- list(Istar = incData)
        
        ### inits 
        initsList <- list(beta = runif(1, 0, 1))
        
        
    }
    

    ### create nimble model
    myModel <- nimbleModel(modelCode, 
                           data = dataList, 
                           constants = constantsList,
                           inits = initsList)
    
    compiled <- compileNimble(myModel) 
    
    
    incDataSamples <- matrix(rep(incData), NROW(samples),
                             ncol = length(incData), byrow = T)
    colnames(incDataSamples) <- paste0('Istar[', 1:tau, ']')
   
    samplesIstar <- cbind(samples, 
                          incDataSamples)
    
    waicList <- calculateWAIC(samplesIstar, compiled)
 
    data.frame(waic = waicList$WAIC,
               lppd = waicList$lppd,
               pWAIC = waicList$pWAIC)
    
}
