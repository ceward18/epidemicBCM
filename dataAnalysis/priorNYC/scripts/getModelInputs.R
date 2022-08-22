################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(alarmFit, incData, smoothI, infPeriod, prior,
                          N, I0, R0, Rstar0, lengthI) {
    
    # constants that are the same for all models
    S0 <- N - I0 - R0
    tau <- length(incData)
    
    # for exponential infectious period
    if (prior == 1) {
        # strong centered on truth
        bb <- 1350
        aa <- 1/7*bb
    } else if (prior == 2) {
        # strong but misspecified
        bb <- 4725
        aa <- 1/2*bb
    } else if (prior == 3) {
        # vague centered on truth
        bb <- 100
        aa <- 1/7*bb
    } else if (prior == 4) {
        # vague and misspecified
        bb <- 350
        aa <- 1/2*bb
    }
    
    if (alarmFit == 'thresh') {
        
        ### constants
        minI <- floor(min(smoothI))
        maxI <- ceiling(max(smoothI))
        n <- 50
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              Rstar0 = Rstar0,
                              lengthI = lengthI,
                              n = n,
                              xAlarm = xAlarm,
                              maxI = maxI,
                              minI = minI)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI= smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 1),
                          delta = runif(1, 0, 1),
                          H = runif(1, 0, maxI/N/4))
        
        ### MCMC specifications
        niter <- 700000
        nburn <- 500000
        nthin <- 10
        
    } else if (alarmFit == 'hill') {
        
        ### constants
        minI <- floor(min(smoothI))
        maxI <- ceiling(max(smoothI))
        n <- 50
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              Rstar0 = Rstar0,
                              lengthI = lengthI,
                              n = n,
                              xAlarm = xAlarm,
                              minI = minI,
                              maxI = maxI)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI= smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 1),
                          delta = runif(1, 0, 1),
                          nu = runif(1, 0, 10),
                          x0 = max(rnorm(1, maxI/4, 10), 1))
        
        ### MCMC specifications
        niter <- 700000
        nburn <- 500000
        nthin <- 10
        
    } else if (alarmFit == 'power') {
        
        ### constants 
        maxI <- ceiling(max(smoothI))
        n <- 50
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              Rstar0 = Rstar0,
                              lengthI = lengthI,
                              n = n,
                              xAlarm = xAlarm)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI= smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 1),
                          k = runif(1, 0, 1))
        
        ### MCMC specifications
        niter <- 700000
        nburn <- 500000
        nthin <- 10
        
    } else if (alarmFit == 'spline') {
        
        ### constants
        n <- 50
        minI <- floor(min(smoothI))
        maxI <- ceiling(max(smoothI))
        xAlarm <- seq(0, maxI, length.out = n)
        nb <- 3
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              Rstar0 = Rstar0,
                              xAlarm = xAlarm,
                              n = n,
                              minI = minI,
                              maxI = maxI,
                              lengthI = lengthI,
                              nb = nb)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI = smoothI,
                         constrain_knots = 1,
                         constrain_min = 1,
                         constrain_max = 1)
        
        ### inits (must satisfy constraint)
        repeat {
            initsList <- list(beta = runif(1, 1/7, 1),
                              b = rnorm(nb, 0, 4),
                              knots = as.vector(quantile(xAlarm, 
                                                         probs = sort(runif(nb - 1, 0.2, 0.8)))))
            
            cond <- all(splineAlarm(xAlarm, initsList$b, initsList$knots) >= 0) & 
                all(splineAlarm(xAlarm, initsList$b, initsList$knots) <= 1) & 
                all(initsList$knots > minI)
            
            if (cond) break
        }
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 600000
        nthin <- 10
        
    } else if (alarmFit == 'gp') {
        
        ### constants
        n <- 10
        maxI <- ceiling(max(smoothI))
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
                              Rstar0 = Rstar0,
                              lengthI = lengthI,
                              dists = distMat,
                              mu0 = 1,
                              ones = logit(seq(0.0001, 0.9999, length.out= n)),
                              n = n,
                              xAlarm = xAlarm,
                              c = vals[1],
                              d = vals[2])
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI= smoothI)
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 1),
                          l = rinvgamma(1, vals[1], vals[2]),
                          sigma = rgamma(1, 100, 50))
        
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 600000
        nthin <- 10
        
    } else if (alarmFit == 'betatSpline') {
        
        ### constants
        timeVec <- 1:tau
        nb <- 4
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              timeVec = timeVec,
                              lengthI = lengthI,
                              Rstar0 = Rstar0,
                              nb = nb)
        
        ### data
        dataList <- list(Istar = incData,
                         constrain_knots = 1)
        
        ### inits
        initsList <- list(b = rnorm(nb, 0, 4),
                          knots = as.vector(quantile(timeVec, 
                                                     probs = sort(runif(nb - 1, 0.2, 0.5)))))
        
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 600000
        nthin <- 10
        
        xAlarm <- NULL
        
    } else if (alarmFit == 'basic') {
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              Rstar0 = Rstar0,
                              lengthI = lengthI)
        
        ### data
        dataList <- list(Istar = incData)
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 1))
        
        
        ### MCMC specifications
        niter <- 300000
        nburn <- 100000
        nthin <- 10
        
        xAlarm <- NULL
        
    }
    
    # adjust specs if model is exponential infectious period
    if (infPeriod == 'exp') {
        
        # adjust constants
        constantsList$lengthI <- NULL
        constantsList$Rstar0 <- NULL
        constantsList$aa <- aa
        constantsList$bb <- bb
        
        # add initial value for rateI
        initsList$rateI <- rgamma(1, aa, bb)
        
        # add initial value for Rstar (everyone removed lengthI days later)
        initsList$Rstar = c(Rstar0, 
                            dataList$Istar[1:(tau-lengthI)])
        
    }
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin,
         xAlarm = xAlarm)
    
    
}