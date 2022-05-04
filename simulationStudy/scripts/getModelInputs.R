################################################################################
#
################################################################################

getModelInput <- function(alarmFit, incData, infPeriod, smoothWindow) {
    
    # constants that are the same for all models
    N <- 1e6
    I0 <- 5
    tau <- length(incData)
    
    # for fixed infectious period
    lengthI <- 7
    
    # for exponential infectious period
    # puts 95% probability of mean infectious period between 6 and 8 days
    bb <- 1350
    aa <- 1/7*bb
    
    if (alarmFit == 'thresh') {
        
        ### constants
        smoothI <- c(0, movingAverage(incData, smoothWindow))
        n <- 100
        maxI <- ceiling(max(smoothI))
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              I0 = I0,
                              lengthI = lengthI,
                              n = n,
                              xAlarm = xAlarm,
                              maxI = maxI)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI = smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 5/7),
                          delta = runif(1, 0, 1),
                          H = runif(1, 0, maxI/N/2))
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 600000
        nthin <- 10
        
    } else if (alarmFit == 'hill') {
        
        ### constants
        smoothI <- c(0, movingAverage(incData, smoothWindow))
        n <- 50
        maxI <- ceiling(max(smoothI))
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              I0 = I0,
                              lengthI = lengthI,
                              n = n,
                              xAlarm = xAlarm,
                              maxI = maxI)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI = smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 5/7),
                          delta = runif(1, 0, 1),
                          nu = runif(1, 0, 10),
                          x0 = max(rnorm(1, maxI/2, 10), 1))
        
        ### MCMC specifications
        niter <- 600000
        nburn <- 400000
        nthin <- 10
        
    } else if (alarmFit == 'power') {
        
        ### constants 
        smoothI <- c(0, movingAverage(incData, smoothWindow))
        n <- 50
        maxI <- ceiling(max(smoothI))
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              I0 = I0,
                              lengthI = lengthI,
                              n = n,
                              xAlarm = xAlarm)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI = smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 5/7),
                          k = runif(1, 0, 1))
        
        ### MCMC specifications
        niter <- 600000
        nburn <- 400000
        nthin <- 10
        
    } else if (alarmFit == 'spline') {
        
        ### constants
        smoothI <- c(0, movingAverage(incData, smoothWindow))
        n <- 50
        maxI <- ceiling(max(smoothI))
        xAlarm <- seq(0, maxI, length.out = n)
        nb <- 3
        
        constantsList <- list(tau = tau,
                              N = N,
                              I0 = I0,
                              xAlarm = xAlarm,
                              n = n,
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
            initsList <- list(beta = runif(1, 1/7, 5/7),
                              b = rnorm(nb, 0, 4),
                              knots = as.vector(quantile(xAlarm, 
                                                         probs = sort(runif(nb - 1, 0, 0.4)))))
            
            cond <- all(splineAlarm(xAlarm, initsList$b, initsList$knots) >= 0) & 
                all(splineAlarm(xAlarm, initsList$b, initsList$knots) <= 1)
            
            if (cond) break
        }
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 600000
        nthin <- 10
        
    } else if (alarmFit == 'gp') {
        
        ### constants
        smoothI <- c(0, movingAverage(incData, smoothWindow))
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
                              I0 = I0,
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
                         smoothI = smoothI)
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 5/7),
                          l = rinvgamma(1, vals[1], vals[2]),
                          sigma = rgamma(1, 100, 50))
        
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 600000
        nthin <- 10
        
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
                              I0 = I0,
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
        
        
        ### MCMC specifications
        niter <- 1200000
        nburn <- 1000000
        nthin <- 10
        
        xAlarm <- NULL
        
    } else if (alarmFit == 'betatSpline') {
        
        ### constants
        timeVec <- 1:tau
        nb <- 4
        
        constantsList <- list(tau = tau,
                              N = N,
                              I0 = I0,
                              timeVec = timeVec,
                              lengthI = lengthI,
                              nb = nb)
        
        ### data
        dataList <- list(Istar = incData,
                         constrain_knots = 1)
        
        ### inits
        initsList <- list(b = rnorm(nb, 0, 4),
                          knots = as.vector(quantile(timeVec, 
                                                     probs = sort(runif(nb - 1, 0, 0.4)))))
        
        
        ### MCMC specifications
        niter <- 600000
        nburn <- 400000
        nthin <- 10
        
        xAlarm <- NULL
        
    } else if (alarmFit == 'basic') {
        
        constantsList <- list(tau = tau,
                              N = N,
                              I0 = I0,
                              lengthI = lengthI)
        
        ### data
        dataList <- list(Istar = incData)
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 5/7))
        
        
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
        constantsList$aa <- aa
        constantsList$bb <- bb
        
        # add initial value for rateI
        initsList$rateI <- rgamma(1, aa, bb)
        
    }
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin,
         xAlarm = xAlarm)
    
}