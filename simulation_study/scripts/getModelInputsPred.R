################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################

getModelInputPred <- function(alarmFit, incData, smoothWindow, obsData) {
    
    # constants that are the same for all models
    N <- 1e6
    I0 <- 5
    tau <- length(incData)
    
    # prior for exponential infectious period
    # puts 95% probability of mean infectious period between 6 and 8 days
    bb <- 1350
    aa <- 1/7*bb
    
    if (alarmFit == 'power') {
        
        ### constants 
        smoothI <- c(0, movingAverage(incData, smoothWindow))
        smoothIObs <- c(0, movingAverage(obsData, smoothWindow))
        n <- 50
        maxI <- ceiling(max(smoothIObs))
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              I0 = I0,
                              n = n,
                              xAlarm = xAlarm,
                              aa = aa,
                              bb = bb)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI = smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 5/7),
                          k = runif(1, 0, 1),
                          rateI = rgamma(1, aa, bb))
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 500000
        nthin <- 10
        
    } else  if (alarmFit == 'thresh') {
        
        ### constants
        smoothI <- c(0, movingAverage(incData, smoothWindow))
        smoothIObs <- c(0, movingAverage(obsData, smoothWindow))
        n <- 100
        minI <- floor(min(smoothIObs))
        maxI <- ceiling(max(smoothIObs))
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              I0 = I0,
                              n = n,
                              xAlarm = xAlarm,
                              minI = minI,
                              maxI = maxI,
                              aa = aa,
                              bb = bb)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI = smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 5/7),
                          delta = runif(1, 0, 1),
                          H = runif(1, 0, maxI/N/2),
                          rateI = rgamma(1, aa, bb))
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 500000
        nthin <- 10
        
    } else if (alarmFit == 'hill') {
        
        ### constants
        smoothI <- c(0, movingAverage(incData, smoothWindow))
        smoothIObs <- c(0, movingAverage(obsData, smoothWindow))
        n <- 50
        minI <- floor(min(smoothIObs))
        maxI <- ceiling(max(smoothIObs))
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              I0 = I0,
                              n = n,
                              xAlarm = xAlarm,
                              minI = minI,
                              maxI = maxI,
                              aa = aa,
                              bb = bb)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothI = smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 5/7),
                          delta = runif(1, 0, 1),
                          nu = runif(1, 0, 10),
                          x0 = max(rnorm(minI, maxI/2, 10), 1),
                          rateI = rgamma(1, aa, bb))
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 500000
        nthin <- 10
        
    } else if (alarmFit == 'spline') {
        
        ### constants
        smoothI <- c(0, movingAverage(incData, smoothWindow))
        smoothIObs <- c(0, movingAverage(obsData, smoothWindow))
        n <- 50
        minI <- floor(min(smoothIObs))
        maxI <- ceiling(max(smoothIObs))
        xAlarm <- seq(0, maxI, length.out = n)
        nb <- 3
        
        constantsList <- list(tau = tau,
                              N = N,
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
            initsList <- list(beta = runif(1, 1/7, 5/7),
                              b = rnorm(nb, 0, 4),
                              knots = as.vector(quantile(xAlarm, 
                                                         probs = sort(runif(nb - 1, 
                                                                            0.4, 
                                                                            0.7)))),
                              rateI = rgamma(1, aa, bb))
            
            cond <- all(splineAlarm(xAlarm, initsList$b, initsList$knots) >= 0) & 
                all(splineAlarm(xAlarm, initsList$b, initsList$knots) <= 1)
            
            if (cond) break
        }
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 500000
        nthin <- 10
        
    } else if (alarmFit == 'gp') {
        
        ### constants
        smoothI <- c(0, movingAverage(incData, smoothWindow))
        smoothIObs <- c(0, movingAverage(obsData, smoothWindow))
        n <- 10
        maxI <- ceiling(max(smoothIObs))
        xAlarm <- seq(0, maxI, length.out = n)
        distMat <- as.matrix(dist(matrix(xAlarm)))
        
        uniqueDists <- distMat[lower.tri(distMat)]
        maxDist <- max(uniqueDists)
        midDist <- getl(maxDist)
        
        # parameters of inverse gamma distribution for prior on lengthscale
        vals <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), 
                            method = 'L-BFGS-B', mid = midDist)$par, 2)
        
        constantsList <- list(tau = tau,
                              N = N,
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
                         smoothI = smoothI)
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 5/7),
                          l = rinvgamma(1, vals[1], vals[2]),
                          sigma = rgamma(1, 150, 50),
                          rateI = rgamma(1, aa, bb))
        
        
        ### MCMC specifications
        niter <- 800000
        nburn <- 500000
        nthin <- 10
        
    } else if (alarmFit == 'basic') {
        
        constantsList <- list(tau = tau,
                              N = N,
                              I0 = I0,
                              aa = aa,
                              bb = bb)
        
        ### data
        dataList <- list(Istar = incData)
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 5/7),
                          rateI = rgamma(1, aa, bb))
        
        
        ### MCMC specifications
        niter <- 350000
        nburn <- 50000
        nthin <- 10
        
        xAlarm <- NULL
    }
    
    # initialize removals 3 days after infections
    initsList$Rstar <- c(rep(0, 3), I0, dataList$Istar[1:(tau-4)]) 
    names(initsList$Rstar) <- paste0('Rstar[', 1:tau, ']')
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin,
         xAlarm = xAlarm)
    
}
