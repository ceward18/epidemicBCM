################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(alarmFit, incData, deathData, smoothI,
                          N, E0, I0, R0, intTime) {
    
    # Mean of prior for initial values
    S0 <- N - E0 - I0 - R0
    
    # constants that are the same for all models
    tau <- length(incData)
    
    # initial value of Istar (unknown when 22 cases became infectious)
    # initialize with unknown infections happening early in the epidemic
    firstIstar <- 10
    Istar <- incData + c(rep(0, E0 + 1), 
                         rmultinom(1, 22, rep(1/firstIstar, firstIstar)),
                         rep(0, tau - firstIstar - E0 - 1))
    
    # initial value of Estar (exposed 2 days before infectious)
    Estar <- c(Istar[3:tau], rep(0, 2))
    
    # initial value of Rstar (unknown when 78 cases were removed)
    # initialize with all removals happening on the last time point
    Rstar <- deathData + c(rep(0, tau-1), 78)
    
    if (alarmFit == 'power') {
        
        ### constants 
        maxI <- ceiling(max(smoothI))
        n <- 50
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              E0 = E0,
                              I0 = I0,
                              n = n,
                              xAlarm = xAlarm)
        
        ### data
        dataList <- list(smoothI = smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 1),
                          k = runif(1, 0, 1),
                          rateE = rgamma(1, 20, 100),
                          rateI = rgamma(1, 20, 100),
                          Estar = Estar,
                          Istar = Istar,
                          Rstar = Rstar)
        
    } else if (alarmFit == 'thresh') {
        
        ### constants
        minI <- floor(min(smoothI))
        maxI <- ceiling(max(smoothI))
        n <- 50
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              E0 = E0,
                              I0 = I0,
                              n = n,
                              xAlarm = xAlarm,
                              maxI = maxI,
                              minI = minI)
        
        ### data
        dataList <- list(smoothI= smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 1),
                          delta = runif(1, 0, 1),
                          H = runif(1, 0, maxI/N/3),
                          rateE = rgamma(1, 20, 100),
                          rateI = rgamma(1, 20, 100),
                          Estar = Estar,
                          Istar = Istar,
                          Rstar = Rstar)
        
    } else if (alarmFit == 'hill') {
        
        ### constants
        minI <- floor(min(smoothI))
        maxI <- ceiling(max(smoothI))
        n <- 50
        xAlarm <- seq(0, maxI, length.out = n)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              E0 = E0,
                              I0 = I0,
                              n = n,
                              xAlarm = xAlarm,
                              minI = minI,
                              maxI = maxI)
        
        ### data
        dataList <- list(smoothI= smoothI)
        
        ### inits
        initsList <- list(beta = runif(1, 1/7, 1),
                          delta = runif(1, 0, 1),
                          nu = runif(1, 0, 10),
                          x0 = max(rnorm(1, maxI/2, 10), 1),
                          rateE = rgamma(1, 20, 100),
                          rateI = rgamma(1, 20, 100),
                          Estar = Estar,
                          Istar = Istar,
                          Rstar = Rstar)
        
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
                              E0 = E0,
                              I0 = I0,
                              xAlarm = xAlarm,
                              n = n,
                              minI = minI,
                              maxI = maxI,
                              nb = nb)
        
        ### data
        dataList <- list(smoothI = smoothI,
                         constrain_knots = 1,
                         constrain_min = 1,
                         constrain_max = 1)
        
        ### inits (must satisfy constraint)
        repeat {
            initsList <- list(beta = runif(1, 1/7, 1),
                              b = rnorm(nb, 0, 4),
                              knots = as.vector(quantile(xAlarm, 
                                                         probs = sort(runif(nb - 1, 0.5, 0.8)))),
                              rateE = rgamma(1, 20, 100),
                              rateI = rgamma(1, 20, 100),
                              Estar = Estar,
                              Istar = Istar,
                              Rstar = Rstar)
            
            cond <- all(splineAlarm(xAlarm, initsList$b, initsList$knots) >= 0) & 
                all(splineAlarm(xAlarm, initsList$b, initsList$knots) <= 1) & 
                all(initsList$knots > minI)
            
            if (cond) break
        }
        
    }  else if (alarmFit == 'splineFixKnot') {
        
        ### constants
        n <- 50
        minI <- floor(min(smoothI))
        maxI <- ceiling(max(smoothI))
        xAlarm <- seq(0, maxI, length.out = n)
        nb <- 3
        
        # use posterior means from full run as knots
        paramPost <- readRDS('./results/paramsPostAll.rds')
        knotsPost <- paramPost[paramPost$alarmFit == 'spline' &
                                   paramPost$param %in% c('knots[1]', 'knots[2]'), 
                               c('param', 'mean', 'peak', 'prior', 'smoothWindow')]
        knots <- knotsPost[knotsPost$peak == peak & 
                               knotsPost$prior == prior & 
                               knotsPost$smoothWindow == smoothWindow, 'mean']
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              E0 = E0,
                              I0 = I0,
                              xAlarm = xAlarm,
                              n = n,
                              nb = nb,
                              knots = knots)
        
        ### data
        dataList <- list(smoothI = smoothI,
                         constrain_min = 1,
                         constrain_max = 1)
        
        ### inits (must satisfy constraint)
        repeat {
            initsList <- list(beta = runif(1, 1/7, 1),
                              b = rnorm(nb, 0, 4),
                              rateE = rgamma(1, 20, 100),
                              rateI = rgamma(1, 20, 100),
                              Estar = Estar,
                              Istar = Istar,
                              Rstar = Rstar)
            
            cond <- all(splineAlarm(xAlarm, initsList$b, constantsList$knots) >= 0) & 
                all(splineAlarm(xAlarm, initsList$b, constantsList$knots) <= 1) 
            
            if (cond) break
        }
        
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
                              E0 = E0,
                              I0 = I0,
                              dists = distMat,
                              mu0 = 1,
                              ones = logit(seq(0.0001, 0.9999, length.out= n)),
                              n = n,
                              xAlarm = xAlarm,
                              c = vals[1],
                              d = vals[2])
        
        ### data
        dataList <- list(smoothI= smoothI)
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 1),
                          l = rinvgamma(1, vals[1], vals[2]),
                          sigma = rgamma(1, 150, 50),
                          rateE = rgamma(1, 20, 100),
                          rateI = rgamma(1, 20, 100),
                          Estar = Estar,
                          Istar = Istar,
                          Rstar = Rstar)
        
    } else if (alarmFit == 'betatSpline') {
        
        ### constants
        timeVec <- 1:tau
        nb <- 4
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              E0 = E0,
                              I0 = I0,
                              timeVec = timeVec,
                              nb = nb)
        
        ### data
        dataList <- list(constrain_knots = 1)
        
        ### inits
        initsList <- list(b = rnorm(nb, 0, 4),
                          knots = as.vector(quantile(timeVec, 
                                                     probs = sort(runif(nb - 1, 0.2, 0.8)))),
                          rateE = rgamma(1, 20, 100),
                          rateI = rgamma(1, 20, 100),
                          Estar = Estar,
                          Istar = Istar,
                          Rstar = Rstar)
        
        xAlarm <- NULL
        
    } else if (alarmFit == 'basic') {
        
        ### constants
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              E0 = E0,
                              I0 = I0)
        
        ### data
        dataList <- NULL
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 1),
                          rateE = rgamma(1, 20, 100),
                          rateI = rgamma(1, 20, 100),
                          Estar = Estar,
                          Istar = Istar,
                          Rstar = Rstar)
        
        xAlarm <- NULL
        
    } else if (alarmFit == 'basicInt') {
        
        # design matrix for intervention
        X <- cbind(1, cumsum(1:tau >= intTime))
        
        
        ### constants
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              E0 = E0,
                              I0 = I0,
                              X = X)
        
        ### data
        dataList <- NULL
        
        ### inits 
        initsList <- list(transParams = c(runif(1, -1, 0), 
                                          runif(1, -0.1, 0)), # a priori assume intervention is successful
                          rateE = rgamma(1, 20, 100),
                          rateI = rgamma(1, 20, 100),
                          Estar = Estar,
                          Istar = Istar,
                          Rstar = Rstar)
        
        xAlarm <- NULL
        
    }
    
    ### MCMC specifications
    niter <- 500000
    nburn <- 250000
    nthin <- 20
    
    ### MCMC specifications
    niter <- 500
    nburn <- 0
    nthin <- 1
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin,
         xAlarm = xAlarm)
    
    
}
