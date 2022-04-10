################################################################################
# function to fit models 
# inputs: 
#   incData - observed incidence data
#   infPeriod - type of infectious period (fixed or exponential)
#   alarmFit - type of function to use to describe the alarm
#   smoothWindow - width of smoothing window 
################################################################################

fitAlarmModel <- function(incData, N, I0, R0,
                          infPeriod, alarmFit, smoothWindow, seed) {
  
  source('./scripts/modelCodes.R')
  
  # constants that are the same for all models
  S0 <- N - I0 - R0
  tau <- length(incData)
  lengthI <- 7
  
  # get appropriate model code
  modelCode <- get(paste0('SIR_', alarmFit, '_', infPeriod))
  
  # for reproducibility so inits are always the same
  set.seed(seed + 3)

  # model-specific constants, data, and inits
  
  if (alarmFit == 'thresh') {
    
    ### constants
    n <- 50
    maxI <- ceiling(max(movingAverage(incData, smoothWindow)) )
    xAlarm <- seq(0, maxI, length.out = n)
    
    constantsList <- list(tau = tau,
                          N = N,
                          I0 = I0,
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
                      H = runif(1, 0, maxI/N))
    
    ### MCMC specifications
    niter <- 800000
    nburn <- 600000
    nthin <- 10
    
  } else if (alarmFit == 'hill') {
    
    ### constants
    n <- 50
    maxI <- ceiling(max(movingAverage(incData, smoothWindow)) )
    xAlarm <- seq(0, maxI, length.out = n)
    
    constantsList <- list(tau = tau,
                          N = N,
                          I0 = I0,
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
    
    ### MCMC specifications
    niter <- 600000
    nburn <- 300000
    nthin <- 15
    
  } else if (alarmFit == 'power') {
    
    ### constants 
    n <- 50
    maxI <- ceiling(max(movingAverage(incData, smoothWindow)) )
    xAlarm <- seq(0, maxI, length.out = n)
    
    constantsList <- list(tau = tau,
                          N = N,
                          I0 = I0,
                          bw = smoothWindow,
                          lengthI = lengthI,
                          n = n,
                          xAlarm = xAlarm)
    
    ### data
    dataList <- list(Istar = incData)
    
    ### inits
    initsList <- list(beta = runif(1, 0, 1),
                      k = runif(1, 0, 1))
    
    ### MCMC specifications
    niter <- 600000
    nburn <- 300000
    nthin <- 15
    
  } else if (alarmFit == 'spline') {
    
    ### constants
    n <- 50
    maxI <- ceiling(max(movingAverage(incData, smoothWindow)))
    xAlarm <- seq(0, maxI, length.out = n)
    nb <- 3
    
    constantsList <- list(tau = tau,
                          N = N,
                          I0 = I0,
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
    
    ### MCMC specifications
    niter <- 800000
    nburn <- 600000
    nthin <- 15
    
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
                          I0 = I0,
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
    
    
    ### MCMC specifications
    niter <- 1000000
    nburn <- 600000
    nthin <- 20
    
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
    niter <- 1400000
    nburn <- 1000000
    nthin <- 10
    
  } else if (alarmFit == 'basic') {
    
    constantsList <- list(tau = tau,
                          N = N,
                          I0 = I0,
                          lengthI = lengthI)
    
    ### data
    dataList <- list(Istar = incData)
    
    ### inits 
    initsList <- list(beta = runif(1, 0, 1))
    
    
    ### MCMC specifications
    niter <- 200000
    nburn <- 50000
    nthin <- 10
    
  }
  
   
  
  ### create nimble model
  myModel <- nimbleModel(modelCode, 
                         data = dataList, 
                         constants = constantsList,
                         inits = initsList)
  myConfig <- configureMCMC(myModel)
  
  # need to ensure all stochastic nodes are monitored for WAIC calculation
  
  if (!alarmFit %in% c('betat', 'basic')) {
    myConfig$addMonitors(c('yAlarm', 'alarm'))
  }
  if (alarmFit == 'betat') {
    myConfig$addMonitors(c('beta'))
  }
  
  # if exponential infectious period, need to use special proposal
  if (infPeriod == 'exp') {
    myConfig$removeSamplers('Rstar') # Nodes will be expanded
    myConfig$addSampler(target = c('Rstar'),
                        type = "RstarUpdate")
  }
  
  # if gaussian process model, use slice sampling
  if (alarmFit == 'gp') {
    
    paramsForSlice <- c('beta', 'l', 'sigma')
    myConfig$removeSampler(paramsForSlice)
    myConfig$addSampler(target = paramsForSlice[1], type = "slice")
    myConfig$addSampler(target = paramsForSlice[2], type = "slice")
    myConfig$addSampler(target = paramsForSlice[3], type = "slice")
    
  } else if (alarmFit == 'betat') {
    
    paramsForSlice <- c('l', 'sigma')
    myConfig$removeSampler(paramsForSlice)
    myConfig$addSampler(target = paramsForSlice[1], type = "slice")
    myConfig$addSampler(target = paramsForSlice[2], type = "slice")
    
  } else if (alarmFit == 'thresh') {
    
    # block sampler for transmission parameters
    paramsForBlock <- c('beta', 'delta', 'H')
    myConfig$removeSampler(paramsForBlock)
    myConfig$addSampler(target = paramsForBlock, type = "RW_block",
                        control = list(adaptInterval = 100,
                                       propCov = diag(c(0.2, 0.2, 0.002))))
    
  } else if (alarmFit == 'hill') {
    
    # block sampler for transmission parameters
    paramsForBlock <- c('beta', 'delta', 'nu', 'x0')
    myConfig$removeSampler(paramsForBlock)
    myConfig$addSampler(target = paramsForBlock, type = "RW_block",
                        control = list(adaptInterval = 100,
                                       propCov = diag(c(0.2, 0.2, 1, 20))))
    
  }
  
  myMCMC <- buildMCMC(myConfig)
  compiled <- compileNimble(myModel, myMCMC) 
  
  runMCMC(compiled$myMCMC, 
          niter = niter, 
          nburnin = nburn,
          thin = nthin,
          setSeed  = seed)
  
}
