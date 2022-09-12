################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(alarmFit, incData, smoothI, prior, peak,
                          N, I0, R0, Rstar0, lengthI) {
    
    # constants that are the same for all models
    S0 <- N - I0 - R0
    tau <- length(incData)

    # parameters for prior on infectious period rate
    if (prior == 1) {
        # strong centered on 1/5
        bb <- 5633
        aa <- 1/5 * bb
    } else if (prior == 2) {
        # strong centered on 1/2
        bb <- 419
        aa <- 1/2 * bb
    } else if (prior == 3) {
        # vague centered on 1/5
        bb <- 437
        aa <- 1/5 * bb
    } else if (prior == 4) {
        # vague centered on 1/2
        bb <- 28
        aa <- 1/2 * bb
    }
    
    if (alarmFit == 'power') {
        
        ### constants 
        maxI <- ceiling(max(smoothI))
        n <- 50
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              n = n,
                              xAlarm = xAlarm,
                              aa = aa,
                              bb = bb)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI= smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 1),
                          k = runif(1, 0, 1),
                          rateI = rgamma(1, aa, bb),
                          Rstar = c(Rstar0, 
                                    dataList$Istar[1:(tau-lengthI)]))
        
        ### MCMC specifications
        niter <- 600000
        nburn <- 300000
        nthin <- 10
        
    } else if (alarmFit == 'thresh') {
        
        ### constants
        minI <- floor(min(smoothI))
        maxI <- ceiling(max(smoothI))
        n <- 50
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              n = n,
                              xAlarm = xAlarm,
                              maxI = maxI,
                              minI = minI,
                              aa = aa,
                              bb = bb)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI= smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 1),
                          delta = runif(1, 0, 1),
                          H = runif(1, 0, maxI/N/4),
                          rateI = rgamma(1, aa, bb),
                          Rstar = c(Rstar0, 
                                    dataList$Istar[1:(tau-lengthI)]))
        
        ### MCMC specifications
        niter <- 600000
        nburn <- 300000
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
                              n = n,
                              xAlarm = xAlarm,
                              minI = minI,
                              maxI = maxI,
                              aa = aa,
                              bb = bb)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI= smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 1),
                          delta = runif(1, 0, 1),
                          nu = runif(1, 0, 10),
                          x0 = max(rnorm(1, maxI/4, 10), 1),
                          rateI = rgamma(1, aa, bb),
                          Rstar = c(Rstar0, 
                                    dataList$Istar[1:(tau-lengthI)]))
        
        ### MCMC specifications
        niter <- 600000
        nburn <- 300000
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
                              xAlarm = xAlarm,
                              n = n,
                              minI = minI,
                              maxI = maxI,
                              nb = nb,
                              aa = aa,
                              bb = bb)
        
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
                                                         probs = sort(runif(nb - 1, 0.2, 0.8)))),
                              rateI = rgamma(1, aa, bb),
                              Rstar = c(Rstar0, 
                                        dataList$Istar[1:(tau-lengthI)]))
            
            cond <- all(splineAlarm(xAlarm, initsList$b, initsList$knots) >= 0) & 
                all(splineAlarm(xAlarm, initsList$b, initsList$knots) <= 1) & 
                all(initsList$knots > minI)
            
            if (cond) break
        }
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 500000
        nthin <- 10
        
    }  else if (alarmFit == 'splineFixKnot') {
        
        ### constants
        n <- 50
        minI <- floor(min(smoothI))
        maxI <- ceiling(max(smoothI))
        xAlarm <- seq(0, maxI, length.out = n)
        nb <- 3
        
        # use posterior means from full run as knots
        paramPost <- readRDS('./resultsCombined/paramsPostAll.rds')
        knotsPost <- paramPost[paramPost$alarmFit == 'spline' &
                                   paramPost$param %in% c('knots[1]', 'knots[2]'), 
                               c('param', 'mean', 'peak', 'prior')]
        knots <- knotsPost[knotsPost$peak == peak & knotsPost$prior == prior, 'mean']
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              xAlarm = xAlarm,
                              n = n,
                              nb = nb,
                              knots = knots,
                              aa = aa,
                              bb = bb)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI = smoothI,
                         constrain_min = 1,
                         constrain_max = 1)
        
        ### inits (must satisfy constraint)
        repeat {
            initsList <- list(beta = runif(1, 1/7, 1),
                              b = rnorm(nb, 0, 4),
                              rateI = rgamma(1, aa, bb),
                              Rstar = c(Rstar0, 
                                        dataList$Istar[1:(tau-lengthI)]))
            
            cond <- all(splineAlarm(xAlarm, initsList$b, constantsList$knots) >= 0) & 
                all(splineAlarm(xAlarm, initsList$b, constantsList$knots) <= 1) 
            
            if (cond) break
        }
        
        ### MCMC specifications
        niter <- 600000
        nburn <- 300000
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
                              dists = distMat,
                              mu0 = 1,
                              ones = logit(seq(0.0001, 0.9999, length.out= n)),
                              n = n,
                              xAlarm = xAlarm,
                              c = vals[1],
                              d = vals[2],
                              aa = aa,
                              bb = bb)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI= smoothI)
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 1),
                          l = rinvgamma(1, vals[1], vals[2]),
                          sigma = rgamma(1, 100, 50),
                          rateI = rgamma(1, aa, bb),
                          Rstar = c(Rstar0, 
                                    dataList$Istar[1:(tau-lengthI)]))
        
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 500000
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
                              nb = nb,
                              aa = aa,
                              bb = bb)
        
        ### data
        dataList <- list(Istar = incData,
                         constrain_knots = 1)
        
        ### inits
        initsList <- list(b = rnorm(nb, 0, 4),
                          knots = as.vector(quantile(timeVec, 
                                                     probs = sort(runif(nb - 1, 0.2, 0.5)))),
                          rateI = rgamma(1, aa, bb),
                          Rstar = c(Rstar0, 
                                    dataList$Istar[1:(tau-lengthI)]))
        
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 500000
        nthin <- 10
        
        xAlarm <- NULL
        
    } else if (alarmFit == 'basic') {
        
        ### constants
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              aa = aa,
                              bb = bb)
        
        ### data
        dataList <- list(Istar = incData)
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 1),
                          rateI = rgamma(1, aa, bb),
                          Rstar = c(Rstar0, 
                                    dataList$Istar[1:(tau-lengthI)]))
        
        
        ### MCMC specifications
        niter <- 400000
        nburn <- 100000
        nthin <- 10
        
        xAlarm <- NULL
        
    }
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin,
         xAlarm = xAlarm)
    
    
}