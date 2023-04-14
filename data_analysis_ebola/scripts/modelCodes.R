################################################################################
# Script containing model code and associated nimbleFunctions
# to be sourced by model fitting script
################################################################################


################################################################################
### Helper functions

# threshold alarm function
thresholdAlarm <- nimbleFunction(     
    run = function(x = double(0), N = double(0), delta = double(0), H = double(0)) {
        returnType(double(0))
        
        result <- delta * (x / N > H)
        
        return(result)
    })
assign('thresholdAlarm', thresholdAlarm, envir = .GlobalEnv)

# power alarm function
powerAlarm <- nimbleFunction(     
    run = function(x = double(0), N = double(0), k = double(0)) {
        returnType(double(0))
        
        result <- 1 - (1 - x / N)^(1 / k)
        
        return(result)
    })
assign('powerAlarm', powerAlarm, envir = .GlobalEnv)

# Hill-Langmuir alarm function
hillAlarm <- nimbleFunction(     
    run = function(x = double(0), nu = double(0), x0 = double(0), delta = double(0)) {
        returnType(double(0))
        
        result <- delta / (1 + (x0 / x) ^ nu)
        
        return(result)
    })
assign('hillAlarm', hillAlarm, envir = .GlobalEnv)

# calculate moving average for smoothing
movingAverage <- nimbleFunction(     
    run = function(x = double(1), bw = double(0)) {
        returnType(double(1))
        n <- length(x)
        bw <- floor(bw)
        
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
assign('movingAverage', movingAverage, envir = .GlobalEnv)

# get effective reproductive number at time t using a forward sum
get_R0 <- nimbleFunction(     
    run = function(betat = double(1), rateI = double(0), N = double(0), S = double(1)) {
        returnType(double(1))
        
        maxInf <- 14
        
        # probability of transition given 1 infectious individual (vector length t)
        pi_SI <- 1 - exp(- betat / N )
        
        # infectious probability of removal
        pi_IR <- 1 - exp(-rateI)
        
        # for exponential periods
        multVec <- c(1, rep(NA, maxInf))
        for (i in 2:length(multVec)) {
            multVec[i] <- multVec[i-1] * (1 - pi_IR)
        }
        
        nTime <- length(pi_SI)
        bw <- length(multVec)
        sumSmooth <- rep(NA, nTime - bw)
        for(k in 1:(nTime - bw)){
            t1 <- k
            t2 <- k + bw - 1
            sumSmooth[k] <- sum(pi_SI[t1:t2] * multVec) * S[k]
        }
        
        return(sumSmooth)
    })
assign('get_R0', get_R0, envir = .GlobalEnv)

# get effective reproductive number for no BC model using a forward sum
get_R0_basic <- nimbleFunction(     
    run = function(beta = double(0), rateI = double(0), N = double(0), S = double(1)) {
        returnType(double(1))
        
        maxInf <- 14
        
        # probability of transition given 1 infectious individual (vector length t)
        pi_SI <- rep(1 - exp(- beta / N ), length(S))
        
        # infectious probability of removal
        pi_IR <- 1 - exp(-rateI)
        
        # for exponential periods
        multVec <- c(1, rep(NA, maxInf))
        for (i in 2:length(multVec)) {
            multVec[i] <- multVec[i-1] * (1 - pi_IR)
        }
        
        nTime <- length(pi_SI)
        bw <- length(multVec)
        sumSmooth <- rep(NA, nTime - bw)
        for(k in 1:(nTime - bw)){
            t1 <- k
            t2 <- k + bw - 1
            sumSmooth[k] <- sum(pi_SI[t1:t2] * multVec) * S[k]
        }
        
        return(sumSmooth)
    })
assign('get_R0_basic', get_R0_basic, envir = .GlobalEnv)

# squared exponential covariance for gaussian process
sqExpCov <- nimbleFunction(     
    run = function(dists = double(2), sigma = double(0), l = double(0)) {
        returnType(double(2))
        n <- dim(dists)[1]
        result <- matrix(nrow = n, ncol = n, init = FALSE)
        sigma2 <- sigma*sigma
        l2 <- 2 * l^2
        for(i in 1:n) {
            for(j in 1:n) {
                
                result[i, j] <- sigma2*exp(- dists[i,j]^2 / l2)
                
                if (i == j) {
                    result[i, j] <- result[i, j] + 1e-6
                }
            }
        }
        
        return(result)
    })
assign('sqExpCov', sqExpCov, envir = .GlobalEnv)

# spline alarm define first as an R function
splineAlarmR <- function(x, b, knots) {
    xBasis <- splines::ns(x, knots = knots)
    c(xBasis %*% b)
    
}
# then convert to something that can be compiled in nimble
splineAlarm <- nimbleRcall(function(x = double(1), b = double(1), knots= double(1)){}, 
                           Rfun = 'splineAlarmR',
                           returnType = double(1)
)
assign('splineAlarm', splineAlarm, envir = .GlobalEnv)

# spline for beta[t] define first as an R function
splineBetaR <- function(x, b, knots) {
    xBasis <- splines::ns(x, knots = knots, Boundary.knots = c(min(x) - 50, max(x)))
    c(xBasis %*% b)
    
}
# then convert to something that can be compiled in nimble
splineBeta <- nimbleRcall(function(x = double(1), b = double(1), knots= double(1)){}, 
                          Rfun = 'splineBetaR',
                          returnType = double(1)
)
assign('splineBeta', splineBeta, envir = .GlobalEnv)

# linear interpolation function to get alarm values for each observed incidence value
nim_approx <- nimbleFunction(     
    run = function(x = double(1), y = double(1), xout = double(0)) {
        returnType(double(0))
        
        # if xout is > max(x), return the closest value
        if (xout >= max(x)) {
            return(y[which(x == max(x))[1]])
        }
        
        # x values on either side of xout
        xPlacement <- 1 * (xout < x)
        # last 0 and first 1
        leftIdx <- max(which(xPlacement == 0))
        rightIdx <- min(which(xPlacement == 1))
        
        x0 <- x[leftIdx]
        y0 <- y[leftIdx]
        
        x1 <- x[rightIdx]
        y1 <- y[rightIdx]
        
        # linear interpolation for x
        out <- y0  + (xout - x0) * ((y1 - y0)/(x1 - x0))
        
        return(out)
    })
assign('nim_approx', nim_approx, envir = .GlobalEnv)


# function to simulate from myModel
simulator <- nimbleFunction(
    setup = function(model, dataNodes) {
        parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
        # exclude data from parent nodes
        parentNodes <- parentNodes[-which(parentNodes %in% dataNodes)]
        parentNodes <- model$expandNodeNames(parentNodes, returnScalarComponents = TRUE)
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
assign('simulator', simulator, envir = .GlobalEnv)

# determine mean for prior on lengthscale parameter
getl <- function(maxDist) {
    sqrt(- (maxDist/2) ^2 / (2 * log(0.025)))
}

# determine shape and scale for prior on rho
myF <- function(x, min, mid,  max) {
    a <- x[1]; b <- x[2]
    
    # mean is b/(a - 1)
    distMean <-  b/(a - 1) 
    
    # variance is b^2 / ((a-1)^2 * (a-2))
    distSD <- sqrt(b^2 / ((a-1)^2 * (a-2)))
    sdEst <- 4
    
    summat <- c(distMean - mid,
                distSD - sdEst)
    sum(summat^2)
}

################################################################################
### Special proposal function for exposure times in exponential model

EstarUpdate <- nimbleFunction(
    name = 'Estar',                              
    contains = sampler_BASE,                     
    setup = function(model, mvSaved, target, control) {                 # REQUIRED setup arguments
        calcNodes <- model$getDependencies(target) 
        
        # number of update attempts 
        nUpdates <- 30
    },  # setup can't return anything
    run = function() {
        currentValue <- model[[target]]                                   
        currentLogProb <- model$getLogProb(calcNodes)                    
        
        # repeat proposal many times 
        for (it in 1:nUpdates) {
            
            proposalValue <- currentValue
            
            nTimePoints <- length(currentValue)
            
            # move a removal time
            possibleSubtract <- which(currentValue > 0)
            subtractIdx <- possibleSubtract[runif(1, 1, length(possibleSubtract) + 1)]
            addIdx <- runif(1, 1, nTimePoints + 1)
            
            proposalValue[subtractIdx] <- proposalValue[subtractIdx] - 1
            proposalValue[addIdx] <- proposalValue[addIdx] + 1
            
            # g(old|new) - g(new|old)
            # possibly have different number of values to subtract from 
            newPossibleSubtract <- which(proposalValue > 0)
            g <- -log(length(newPossibleSubtract)) +log(length(possibleSubtract))
            
            # put proposal value in model
            model[[target]] <<- proposalValue                                
            proposalLogProb <- model$calculate(calcNodes)                     
            logAcceptanceRatio <- proposalLogProb - currentLogProb + g            
            
            accept <- decide(logAcceptanceRatio)                              
            
            if (accept) {
                # no changes to model object needed
                currentLogProb <- proposalLogProb
                currentValue <- proposalValue
                
            } else {
                # reject proposal and revert model to current state
                model[[target]] <<- currentValue
                
                # current full conditional (calculate overwrites the stored value)
                currentLogProb <- model$calculate(calcNodes) 
            }
            
        } # end loop
        
        # synchronize model -> mvSaved after nUpdates
        copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        
    },
    methods = list(                              # required method for sampler_BASE base class
        reset = function() {}
    )
)
assign('EstarUpdate', EstarUpdate, envir = .GlobalEnv)

# update for transition vectors Istar and Rstar which are partially observed
transUpdate <- nimbleFunction(
    name = 'transUpdate',                              
    contains = sampler_BASE,                     
    setup = function(model, mvSaved, target, control) {                 # REQUIRED setup arguments
        calcNodes <- model$getDependencies(target) 
        
        # observed part of vector is passed as control argument
        fixedValues <- control$fixed
        
        # number of update attempts 
        nUpdates <- 30
    },  # setup can't return anything
    run = function() {
        
        currentValue <- model[[target]]                                   
        currentLogProb <- model$getLogProb(calcNodes)       
        
        updatePart <- currentValue - fixedValues
        
        # repeat proposal many times 
        for (it in 1:nUpdates) {
            
            updatePartProposal <- updatePart
            
            nTimePoints <- length(updatePart)
            
            # move a transition time
            possibleSubtract <- which(updatePart > 0)
            subtractIdx <- possibleSubtract[runif(1, 1, length(possibleSubtract) + 1)]
            addIdx <- runif(1, 1, nTimePoints + 1)
            
            updatePartProposal[subtractIdx] <- updatePartProposal[subtractIdx] - 1
            updatePartProposal[addIdx] <- updatePartProposal[addIdx] + 1
            
            # g(old|new) - g(new|old)
            # possibly have different number of values to subtract from 
            newPossibleSubtract <- which(updatePartProposal > 0)
            g <- -log(length(newPossibleSubtract)) +log(length(possibleSubtract))
            
            # put proposal value in model
            proposalValue <- updatePartProposal + fixedValues
            
            model[[target]] <<- proposalValue                               
            proposalLogProb <- model$calculate(calcNodes)                     
            logAcceptanceRatio <- proposalLogProb - currentLogProb + g      
            
            
            accept <- decide(logAcceptanceRatio)                              
            
            if (accept) {
                # no changes to model object needed
                currentLogProb <- proposalLogProb
                currentValue <- proposalValue
                
            } else {
                # reject proposal and revert model to current state
                model[[target]] <<- currentValue
                
                # current full conditional (calculate overwrites the stored value)
                currentLogProb <- model$calculate(calcNodes) 
            }
            
        } # end loop
        
        # synchronize model -> mvSaved after nUpdates
        copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        
    },
    methods = list(                              # required method for sampler_BASE base class
        reset = function() {}
    )
)
assign('transUpdate', transUpdate, envir = .GlobalEnv)

################################################################################
### Power alarm models

# for model fitting (smoothI known)

SEIR_power <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # compute alarm
        alarm[t] <- powerAlarm(smoothI[t], N, k)
        
        probSE[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
        
    }
    
    # estimated effective R0
    R0_eff[1:(tau-15)] <- get_R0(beta * (1 - alarm[1:tau]), rateI, N, S[1:tau])
    
    for (i in 1:n) {
        yAlarm[i] <- powerAlarm(xAlarm[i], N, k)
    }
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    k ~ dgamma(0.1, 0.1)
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
})

# for forward simulation (smoothI unknown)


SEIR_power_sim <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### first time point
    smoothI[1] <- 0
    alarm[1] <- powerAlarm(smoothI[1], N, k)
    
    probSE[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    Estar[1] ~ dbin(probSE[1], S[1])
    Istar[1] ~ dbin(probEI, E[1])
    Rstar[1] ~ dbin(probIR, I[1])
    
    # update S E I
    S[2] <- S[1] - Estar[1]
    E[2] <- E[1] + Estar[1] - Istar[1]
    I[2] <- I[1] + Istar[1] - Rstar[1]
    
    ### rest of time points
    for(t in 2:tau) {
        
        # compute alarm
        smoothI[t] <-  movingAverage(Istar[1:(t-1)], bw)[t-1]
        alarm[t] <- powerAlarm(smoothI[t], N, k)
        
        probSE[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    k ~ dgamma(0.1, 0.1)
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
})



################################################################################
### Threshold alarm models

# for model fitting (smoothI known)

SEIR_thresh <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # compute alarm
        alarm[t] <- thresholdAlarm(smoothI[t],  N, delta, H)
        
        probSE[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    # estimated effective R0
    R0_eff[1:(tau-15)] <- get_R0(beta * (1 - alarm[1:tau]), rateI, N, S[1:tau])
    
    for (i in 1:n) {
        yAlarm[i] <- thresholdAlarm(xAlarm[i],  N, delta, H)
    }
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    delta ~ dbeta(1, 1)
    H ~ dunif(minI/N, maxI/N)
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
})

SEIR_thresh_sim <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### first time point
    smoothI[1] <- 0
    alarm[1] <- powerAlarm(smoothI[1], N, k)
    
    probSE[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    Estar[1] ~ dbin(probSE[1], S[1])
    Istar[1] ~ dbin(probEI, E[1])
    Rstar[1] ~ dbin(probIR, I[1])
    
    # update S E I
    S[2] <- S[1] - Estar[1]
    E[2] <- E[1] + Estar[1] - Istar[1]
    I[2] <- I[1] + Istar[1] - Rstar[1]
    
    ### rest of time points
    for(t in 2:tau) {
        
        # compute alarm
        smoothI[t] <-  movingAverage(Istar[1:(t-1)], bw)[t-1]
        alarm[t] <- thresholdAlarm(smoothI[t],  N, delta, H)
        
        probSE[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    delta ~ dbeta(1, 1)
    H ~ dunif(minI/N, maxI/N)
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
})

################################################################################
### Hill alarm models

# for model fitting (smoothI known)

SEIR_hill <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # compute alarm
        alarm[t] <- hillAlarm(smoothI[t], nu, x0, delta)
        
        probSE[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
    }
    
    # estimated effective R0
    R0_eff[1:(tau-15)] <- get_R0(beta * (1 - alarm[1:tau]), rateI, N, S[1:tau])
    
    
    for (i in 1:n) {
        yAlarm[i] <- hillAlarm(xAlarm[i],  nu, x0, delta)
    }
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    delta ~ dbeta(1, 1)
    nu ~ dunif(0, 50)
    x0 ~ dunif(minI, maxI)
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
})

SEIR_hill_sim <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### first time point
    smoothI[1] <- 0
    alarm[1] <- powerAlarm(smoothI[1], N, k)
    
    probSE[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    Estar[1] ~ dbin(probSE[1], S[1])
    Istar[1] ~ dbin(probEI, E[1])
    Rstar[1] ~ dbin(probIR, I[1])
    
    # update S E I
    S[2] <- S[1] - Estar[1]
    E[2] <- E[1] + Estar[1] - Istar[1]
    I[2] <- I[1] + Istar[1] - Rstar[1]
    
    ### rest of time points
    for(t in 2:tau) {
        
        # compute alarm
        smoothI[t] <-  movingAverage(Istar[1:(t-1)], bw)[t-1]
        alarm[t] <- hillAlarm(smoothI[t], nu, x0, delta)
        
        probSE[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    delta ~ dbeta(1, 1)
    nu ~ dunif(0, 50)
    x0 ~ dunif(minI, maxI)
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
})

################################################################################
### Spline alarm models (knots estimated)

# for model fitting (smoothI known)

SEIR_spline <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # compute alarm
        alarm[t] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[t])
        
        probSE[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    # estimated effective R0
    R0_eff[1:(tau-15)] <- get_R0(beta * (1 - alarm[1:tau]), rateI, N, S[1:tau])
    
    yAlarm[1:n] <- splineAlarm(xAlarm[1:n], b[1:nb], knots[1:(nb - 1)])
    
    # constrain yAlarm to be between 0 and 1
    minYAlarm <- min(yAlarm[1:n])
    maxYAlarm <- max(yAlarm[1:n])
    constrain_min ~ dconstraint(minYAlarm >= 0)
    constrain_max ~ dconstraint(maxYAlarm <= 1)
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    for (i in 1:nb) {
        b[i] ~ dnorm(0, sd = 100)
    }
    for (i in 1:(nb - 1)) {
        knots[i] ~ dunif(min = minI, max = maxI - 1)
    }
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
    # constrain knots to be ordered
    constrain_knots ~ dconstraint(knots[1] < knots[2])
    
})


SEIR_spline_sim <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### first time point
    smoothI[1] <- 0
    alarm[1] <- powerAlarm(smoothI[1], N, k)
    
    probSE[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    Estar[1] ~ dbin(probSE[1], S[1])
    Istar[1] ~ dbin(probEI, E[1])
    Rstar[1] ~ dbin(probIR, I[1])
    
    # update S E I
    S[2] <- S[1] - Estar[1]
    E[2] <- E[1] + Estar[1] - Istar[1]
    I[2] <- I[1] + Istar[1] - Rstar[1]
    
    ### rest of time points
    for(t in 2:tau) {
        
        # compute alarm
        smoothI[t] <-  movingAverage(Istar[1:(t-1)], bw)[t-1]
        alarm[t] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[t])
        
        probSE[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    yAlarm[1:n] <- splineAlarm(xAlarm[1:n], b[1:nb], knots[1:(nb - 1)])
    
    # constrain yAlarm to be between 0 and 1
    minYAlarm <- min(yAlarm[1:n])
    maxYAlarm <- max(yAlarm[1:n])
    constrain_min ~ dconstraint(minYAlarm >= 0)
    constrain_max ~ dconstraint(maxYAlarm <= 1)
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    for (i in 1:nb) {
        b[i] ~ dnorm(0, sd = 100)
    }
    for (i in 1:(nb - 1)) {
        knots[i] ~ dunif(min = minI, max = maxI - 1)
    }
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
    # constrain knots to be ordered
    constrain_knots ~ dconstraint(knots[1] < knots[2])
    
})


################################################################################
### Spline alarm models (knots fixed)

# for model fitting (smoothI known)

SEIR_splineFixKnot <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # compute alarm
        alarm[t] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[t])
        
        probSE[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    # estimated effective R0
    R0_eff[1:(tau-15)] <- get_R0(beta * (1 - alarm[1:tau]), rateI, N, S[1:tau])
    
    yAlarm[1:n] <- splineAlarm(xAlarm[1:n], b[1:nb], knots[1:(nb - 1)])
    
    # constrain yAlarm to be between 0 and 1
    minYAlarm <- min(yAlarm[1:n])
    maxYAlarm <- max(yAlarm[1:n])
    constrain_min ~ dconstraint(minYAlarm >= 0)
    constrain_max ~ dconstraint(maxYAlarm <= 1)
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    for (i in 1:nb) {
        b[i] ~ dnorm(0, sd = 100)
    }
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
    
})

SEIR_splineFixKnot_sim <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### first time point
    smoothI[1] <- 0
    alarm[1] <- powerAlarm(smoothI[1], N, k)
    
    probSE[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    Estar[1] ~ dbin(probSE[1], S[1])
    Istar[1] ~ dbin(probEI, E[1])
    Rstar[1] ~ dbin(probIR, I[1])
    
    # update S E I
    S[2] <- S[1] - Estar[1]
    E[2] <- E[1] + Estar[1] - Istar[1]
    I[2] <- I[1] + Istar[1] - Rstar[1]
    
    ### rest of time points
    for(t in 2:tau) {
        
        # compute alarm
        smoothI[t] <-  movingAverage(Istar[1:(t-1)], bw)[t-1]
        alarm[t] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[t])
        
        probSE[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    yAlarm[1:n] <- splineAlarm(xAlarm[1:n], b[1:nb], knots[1:(nb - 1)])
    
    # constrain yAlarm to be between 0 and 1
    minYAlarm <- min(yAlarm[1:n])
    maxYAlarm <- max(yAlarm[1:n])
    constrain_min ~ dconstraint(minYAlarm >= 0)
    constrain_max ~ dconstraint(maxYAlarm <= 1)
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    for (i in 1:nb) {
        b[i] ~ dnorm(0, sd = 100)
    }
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
})

################################################################################
### Gaussian process alarm models

# for model fitting (smoothI known)

SEIR_gp <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # compute alarm
        alarm[t] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[t])
        
        probSE[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
    }
    
    # estimated effective R0
    R0_eff[1:(tau-15)] <- get_R0(beta * (1 - alarm[1:tau]), rateI, N, S[1:tau])
    
    yAlarm[1] <- 0
    mu[2:n] <- mu0 * ones[2:n]
    cov[2:n, 2:n] <- sqExpCov(dists[2:n, 2:n], sigma, l)
    logit(yAlarm[2:n]) ~ dmnorm(mu[2:n], cov = cov[2:n, 2:n])
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    sigma ~ dgamma(150, 50)
    l ~ dinvgamma(c, d)
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
})

SEIR_gp_sim <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### first time point
    smoothI[1] <- 0
    alarm[1] <- powerAlarm(smoothI[1], N, k)
    
    probSE[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    Estar[1] ~ dbin(probSE[1], S[1])
    Istar[1] ~ dbin(probEI, E[1])
    Rstar[1] ~ dbin(probIR, I[1])
    
    # update S E I
    S[2] <- S[1] - Estar[1]
    E[2] <- E[1] + Estar[1] - Istar[1]
    I[2] <- I[1] + Istar[1] - Rstar[1]
    
    ### rest of time points
    for(t in 2:tau) {
        
        # compute alarm
        smoothI[t] <-  movingAverage(Istar[1:(t-1)], bw)[t-1]
        alarm[t] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[t])
        
        probSE[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    yAlarm[1] <- 0
    mu[2:n] <- mu0 * ones[2:n]
    cov[2:n, 2:n] <- sqExpCov(dists[2:n, 2:n], sigma, l)
    logit(yAlarm[2:n]) ~ dmnorm(mu[2:n], cov = cov[2:n, 2:n])
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    sigma ~ dgamma(150, 50)
    l ~ dinvgamma(c, d)
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
})

################################################################################
### Standard SIR model with no alarm

SEIR_basic <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### rest of time points
    for(t in 1:tau) {
        
        probSE[t] <- 1 - exp(- beta * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
    }
    
    # estimated effective R0
    R0_eff[1:(tau-15)] <- get_R0_basic(beta, rateI, N, S[1:tau])
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
})


################################################################################
### SIR model with time-varying beta (no alarm)
# Spline

SEIR_betatSpline <-  nimbleCode({
    
    S[1] <- S0
    E[1] <- E0
    I[1] <- I0
    
    probEI <- 1 - exp(-rateE)
    probIR <- 1 - exp(-rateI)
    
    ### rest of time points
    for(t in 1:tau) {
        
        probSE[t] <- 1 - exp(- beta[t] * I[t] / N)
        
        Estar[t] ~ dbin(probSE[t], S[t])
        Istar[t] ~ dbin(probEI, E[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S E I
        S[t + 1] <- S[t] - Estar[t]
        E[t + 1] <- E[t] + Estar[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
    }
    
    # estimated effective R0
    R0_eff[1:(tau-15)] <- get_R0(beta[1:tau], rateI, N, S[1:tau])
    
    log(beta[1:tau]) <- splineBeta(timeVec[1:tau], b[1:nb], knots[1:(nb - 1)])
    
    # priors
    for (i in 1:nb) {
        b[i] ~ dnorm(0, sd = 100)
    }
    for (i in 1:(nb - 1)) {
        knots[i] ~ dunif(min = 1, max = tau)
    }
    rateI ~ dgamma(20, 100)
    rateE ~ dgamma(20, 100)
    
    # constrain knots to be ordered
    constrain_knots ~ dconstraint(knots[1] < knots[2] & knots[2] < knots[3] )
    
})


