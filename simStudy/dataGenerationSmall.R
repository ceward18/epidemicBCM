################################################################################
# NIMBLE Models to Simulate from
# 14 and 30 day smoothing
# incidence-based alarms
# threshold, hill, power alarms
# exponential infectious periods
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
rateI <- 1/7
beta <- 0.42   #R0 = 3

dataNodes <- c(paste0('Istar[', 1:tau, ']'),
               paste0('Rstar[', 1:tau, ']'))

################################################################################
# Threshold alarms
################################################################################

################################################################################
### 14-day Incidence
SIR_thresh_exp_14 <- nimbleModel(code = SIR_thresh_exp,
                                   constants = list(N = N, 
                                                    tau = tau,
                                                    I0 = I0,
                                                    bw = 14,
                                                    a = 1,
                                                    b = 1,
                                                    maxI = 1000))
cModel <- compileNimble(SIR_thresh_exp_14)
sim_R <- simulator(SIR_thresh_exp_14, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_thresh_exp_14)

trueVals <- c(beta = beta,
              delta = 0.7,
              H = 0.00012,
              rateI = rateI)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- dataNodes

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

saveRDS(toSave, './Data/thresh_exp_14.rds')

################################################################################
### 30-day Incidence
SIR_thresh_exp_30 <- nimbleModel(code = SIR_thresh_exp,
                                   constants = list(N = N, 
                                                    tau = tau,
                                                    I0 = I0,
                                                    bw = 30,
                                                    a = 1,
                                                    b = 1,
                                                    maxI = 1000))

cModel <- compileNimble(SIR_thresh_exp_30)
sim_R <- simulator(SIR_thresh_exp_30, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_thresh_exp_30)

trueVals <- c(beta = beta, 
              delta = 0.78,
              H = 0.00007,
              rateI = rateI)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- dataNodes

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

saveRDS(toSave, './Data/thresh_exp_30.rds')

################################################################################
# Hill Alarm
################################################################################

################################################################################
### 14-day Incidence
SIR_hill_exp_14 <- nimbleModel(code = SIR_hill_exp,
                                 constants = list(N = N, 
                                                  tau = tau,
                                                  I0 = I0,
                                                  bw = 14,
                                                  a = 1,
                                                  b = 1,
                                                  maxI = 1000))
cModel <- compileNimble(SIR_hill_exp_14)
sim_R <- simulator(SIR_hill_exp_14, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_hill_exp_14)

trueVals <- c(beta = beta, 
              delta = 0.75,
              nu = 5,
              x0 = 150,
              rateI = rateI)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- dataNodes

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

saveRDS(toSave, './Data/hill_exp_14.rds')

################################################################################
### 30-day Incidence
SIR_hill_exp_30 <- nimbleModel(code = SIR_hill_exp,
                                 constants = list(N = N, 
                                                  tau = tau,
                                                  I0 = I0,
                                                  bw = 30,
                                                  a = 1,
                                                  b = 1,
                                                  maxI = 1000))

cModel <- compileNimble(SIR_hill_exp_30)
sim_R <- simulator(SIR_hill_exp_30, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_hill_exp_30)

trueVals <- c(beta = beta, 
              delta = 0.85,
              nu = 5,
              x0 = 100,
              rateI = rateI)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- dataNodes

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

saveRDS(toSave, './Data/hill_exp_30.rds')


################################################################################
# Power Alarm
################################################################################

################################################################################
### 14-day Incidence
SIR_power_exp_14 <- nimbleModel(code = SIR_power_exp,
                                  constants = list(N = N, 
                                                   tau = tau,
                                                   I0 = I0,
                                                   bw = 14,
                                                   a = 1,
                                                   b = 1))
cModel <- compileNimble(SIR_power_exp_14)
sim_R <- simulator(SIR_power_exp_14, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_power_exp_14)

trueVals <- c(beta = beta, 
              k = 0.0003,
              rateI = rateI)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- dataNodes

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

saveRDS(toSave, './Data/power_exp_14.rds')

################################################################################
### 30-day Incidence
SIR_power_exp_30 <- nimbleModel(code = SIR_power_exp,
                                  constants = list(N = N, 
                                                   tau = tau,
                                                   I0 = I0,
                                                   bw = 30,
                                                   a = 1,
                                                   b = 1))

cModel <- compileNimble(SIR_power_exp_30)
sim_R <- simulator(SIR_power_exp_30, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_power_exp_30)

trueVals <- c(beta = beta, 
              k = 0.00018,
              rateI = rateI)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- dataNodes

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

saveRDS(toSave, './Data/power_exp_30.rds')
