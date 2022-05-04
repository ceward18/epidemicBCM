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

# load model scripts and simulator function
source('./scripts/modelCodesSim.R')

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
