################################################################################
# NIMBLE Models to Simulate from
# 14 and 30 day smoothing
# incidence-based alarms
# threshold, hill, power alarms
# fixed infectious periods
################################################################################

library(nimble)

# number of simulated epidemics in each scenario
nSim <- 50


################################################################################
### Helper functions

# threshold alarm function
thresholdAlarm <- nimbleFunction(     
  run = function(x = double(0), N = double(0), delta = double(0), H = double(0)) {
    returnType(double(0))
    
    result <- delta * (x / N > H)
    
    return(result)
  })

# power alarm function
powerAlarm <- nimbleFunction(     
  run = function(x = double(0), N = double(0), k = double(0)) {
    returnType(double(0))
    
    result <- 1 - (1 - x / N)^(1 / k)
    
    return(result)
  })

# Hill-Langmuir alarm function
hillAlarm <- nimbleFunction(     
  run = function(x = double(0), nu = double(0), x0 = double(0), delta = double(0)) {
    returnType(double(0))
    
    result <- delta / (1 + (x0 / x) ^ nu)
    
    return(result)
  })

# calculate moving average for smoothing
movingAverage <- nimbleFunction(     
  run = function(x = double(1), bw = double(0)) {
    returnType(double(1))
    n <- length(x)
    
    out <- rep(0, n)
    for (i in 1:n) {
      
      if (i < bw) {
        t1 = 1
        t2 = i
      } else {
        t1 = i - bw + 1
        t2 = i
      }
      
      out[i] <- mean(x[t1:t2])
    }
    
    return(out)
  })

################################################################################
# simulator function

simulator <- nimbleFunction(
  setup = function(model, dataNodes) {
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
    # exclude data from parent nodes
    parentNodes <- parentNodes[-which(parentNodes %in% dataNodes)]
    cat("Stochastic parents of data are: ", paste(parentNodes, sep = ','), ".\n")
    simNodes <- model$getDependencies(parentNodes, self = FALSE,
                                      downstream = T)
    
    nData <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
  },
  run = function(params = double(1), nSim = double()) {
    simDat <- matrix(nrow = nSim, ncol = nData)   
    for(i in 1:nSim) {
      values(model, parentNodes) <<- params
      model$simulate(simNodes, includeData = TRUE)
      simDat[i, ] <- values(model, dataNodes)
    }
    return(simDat)
    returnType(double(2))
  })

################################################################################
# parameters for all

N <- 1e6
I0 <- 5
S0 <- N - I0
tau <- 100
lengthI <- 7
beta <- 0.36

dataNodes <- c(paste0('Istar[', 1:tau, ']'))


################################################################################
# Threshold alarms
################################################################################

SIR_thresh_fixed <-  nimbleCode({
  
  S[1] <- N - I0 
  I[1] <- I0
  Rstar[1:(lengthI-1)] <- 0
  Rstar[lengthI] <- I0
  
  ### first time point for incidence based alarm
  # compute alarm
  smoothI[1] <- 0
  alarm[1] <- thresholdAlarm(smoothI[1],  N, delta, H)
  
  probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
  
  Istar[1] ~ dbin(probSI[1], S[1])
  Rstar[1 + lengthI] <- Istar[1]
  
  # update S and I
  S[2] <- S[1] - Istar[1]
  I[2] <- I[1] + Istar[1] - Rstar[1]
  
  ### rest of time points
  for(t in 2:tau) {
    
    # compute alarm
    smoothI[t] <- movingAverage(Istar[1:(t-1)], bw)[t-1]
    alarm[t] <- thresholdAlarm(smoothI[t],  N, delta, H)
    
    probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
    
    Istar[t] ~ dbin(probSI[t], S[t])
    
    Rstar[t + lengthI] <- Istar[t]
    
    # update S and I
    S[t + 1] <- S[t] - Istar[t]
    I[t + 1] <- I[t] + Istar[t] - Rstar[t]
    
  }
  
  # priors
  beta ~ dgamma(0.1, 0.1)
  delta ~ dbeta(1, 1)
  H ~ dgamma(1, 10)
  rateI ~ dgamma(20, 100)
  
})

################################################################################
### 14-day Incidence
SIR_thresh_fixed_14 <- nimbleModel(code = SIR_thresh_fixed,
                                     constants = list(N = N, 
                                                      tau = tau,
                                                      I0 = I0,
                                                      bw = 14,
                                                      lengthI = lengthI))
cModel <- compileNimble(SIR_thresh_fixed_14)
sim_R <- simulator(SIR_thresh_fixed_14, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_thresh_fixed_14)

trueVals <- c(beta = beta, 
              delta = 0.66,
              H = 0.0001)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- dataNodes

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]


saveRDS(toSave, './Data/thresh_fixed_14.rds')

################################################################################
### 30-day Incidence
SIR_thresh_fixed_30 <- nimbleModel(code = SIR_thresh_fixed,
                                     constants = list(N = N, 
                                                      tau = tau,
                                                      I0 = I0,
                                                      bw = 30,
                                                      lengthI = lengthI))

cModel <- compileNimble(SIR_thresh_fixed_30)
sim_R <- simulator(SIR_thresh_fixed_30, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_thresh_fixed_30)

trueVals <- c(beta = beta, 
              delta = 0.66,
              H = 0.00008)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- dataNodes

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]


saveRDS(toSave, './Data/thresh_fixed_30.rds')


################################################################################
# Hill Alarm
################################################################################

SIR_hill_fixed <-  nimbleCode({
  
  S[1] <- N - I0 
  I[1] <- I0
  Rstar[1:(lengthI-1)] <- 0
  Rstar[lengthI] <- I0
  
  ### first time point for incidence based alarm
  # compute alarm
  smoothI[1] <- 0
  alarm[1] <- hillAlarm(smoothI[1], nu, x0, delta)
  
  probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
  
  Istar[1] ~ dbin(probSI[1], S[1])
  Rstar[1 + lengthI] <- Istar[1]
  
  # update S and I
  S[2] <- S[1] - Istar[1]
  I[2] <- I[1] + Istar[1] - Rstar[1]
  
  ### rest of time points
  for(t in 2:tau) {
    
    # compute alarm
    smoothI[t] <- movingAverage(Istar[1:(t-1)], bw)[t-1]
    alarm[t] <- hillAlarm(smoothI[t], nu, x0, delta)
    
    probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
    
    Istar[t] ~ dbin(probSI[t], S[t])
    
    Rstar[t + lengthI] <- Istar[t]
    
    # update S and I
    S[t + 1] <- S[t] - Istar[t]
    I[t + 1] <- I[t] + Istar[t] - Rstar[t]
    
  }
  
  # priors
  beta ~ dgamma(0.1, 0.1)
  nu ~ dgamma(0.1, 0.1)
  x0 ~ dnorm(100, sd=10)
  delta ~ dbeta(1, 1)
  
})

################################################################################
### 14-day Incidence
SIR_hill_fixed_14 <- nimbleModel(code = SIR_hill_fixed,
                                       constants = list(N = N, 
                                                        tau = tau,
                                                        I0 = I0,
                                                        bw = 14,
                                                        lengthI = lengthI))
cModel <- compileNimble(SIR_hill_fixed_14)
sim_R <- simulator(SIR_hill_fixed_14, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_hill_fixed_14)

trueVals <- c(beta = beta, 
              nu = 5,
              x0 = 140,
              delta = 0.65)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- dataNodes

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

saveRDS(toSave, './Data/hill_fixed_14.rds')

################################################################################
### 30-day Incidence
SIR_hill_fixed_30 <- nimbleModel(code = SIR_hill_fixed,
                                       constants = list(N = N, 
                                                        tau = tau,
                                                        I0 = I0,
                                                        bw = 30,
                                                        lengthI = lengthI))

cModel <- compileNimble(SIR_hill_fixed_30)
sim_R <- simulator(SIR_hill_fixed_30, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_hill_fixed_30)

trueVals <- c(beta = beta, 
              nu = 5,
              x0 = 80,
              delta = 0.7)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- dataNodes

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

saveRDS(toSave, './Data/hill_fixed_30.rds')


################################################################################
# Power Alarm
################################################################################

SIR_power_fixed <-  nimbleCode({
  
  S[1] <- N - I0 
  I[1] <- I0
  Rstar[1:(lengthI-1)] <- 0
  Rstar[lengthI] <- I0
  
  ### first time point for incidence based alarm
  # compute alarm
  smoothI[1] <- 0
  alarm[1] <- powerAlarm(smoothI[1], N, k)
  
  probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
  
  Istar[1] ~ dbin(probSI[1], S[1])
  Rstar[1 + lengthI] <- Istar[1]
  
  # update S and I
  S[2] <- S[1] - Istar[1]
  I[2] <- I[1] + Istar[1] - Rstar[1]
  
  ### rest of time points
  for(t in 2:tau) {
    
    # compute alarm
    smoothI[t] <- movingAverage(Istar[1:(t-1)], bw)[t-1]
    alarm[t] <- powerAlarm(smoothI[t], N, k)
    
    probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
    
    Istar[t] ~ dbin(probSI[t], S[t])
    
    Rstar[t + lengthI] <- Istar[t]
    
    # update S and I
    S[t + 1] <- S[t] - Istar[t]
    I[t + 1] <- I[t] + Istar[t] - Rstar[t]
    
  }
  
  # priors
  beta ~ dgamma(0.1, 0.1)
  k ~ dgamma(0.1, 0.1)
  
})

################################################################################
### 14-day Incidence
SIR_power_fixed_14 <- nimbleModel(code = SIR_power_fixed,
                                     constants = list(N = N, 
                                                      tau = tau,
                                                      I0 = I0,
                                                      bw = 14,
                                                      lengthI = lengthI))
cModel <- compileNimble(SIR_power_fixed_14)
sim_R <- simulator(SIR_power_fixed_14, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_power_fixed_14)

trueVals <- c(beta = beta, 
              k = 0.0003)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- dataNodes

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

saveRDS(toSave, './Data/power_fixed_14.rds')

################################################################################
### 30-day Incidence
SIR_power_fixed_30 <- nimbleModel(code = SIR_power_fixed,
                                     constants = list(N = N, 
                                                      tau = tau,
                                                      I0 = I0,
                                                      bw = 30,
                                                      lengthI = lengthI))

cModel <- compileNimble(SIR_power_fixed_30)
sim_R <- simulator(SIR_power_fixed_30, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_power_fixed_30)

trueVals <- c(beta = beta, 
              k = 0.00018)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- dataNodes

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

saveRDS(toSave, './Data/power_fixed_30.rds')
