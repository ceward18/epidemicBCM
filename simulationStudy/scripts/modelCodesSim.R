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

################################################################################
### Threshold alarm models

### Fixed infectious period

SIR_thresh_fixed <-  nimbleCode({
    
    S[1] <- N - I0 
    I[1] <- I0
    Rstar[1:(lengthI-1)] <- 0
    Rstar[lengthI] <- I0
    
    ### rest of time points
    for(t in 1:tau) {
        
        # compute alarm
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
    H ~ dunif(0, maxI/N)
    
})

### Exponential infectious period

SIR_thresh_exp <-  nimbleCode({
    
    S[1] <- N - I0 
    I[1] <- I0
    
    probIR <- 1 - exp(-rateI)
    
    ### first time point for incidence based alarm
    # compute alarm
    smoothI[1] <- 0
    alarm[1] <- thresholdAlarm(smoothI[1],  N, delta, H)
    
    probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    Istar[1] ~ dbin(probSI[1], S[1])
    Rstar[1] ~ dbin(probIR, I[1])
    
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
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    delta ~ dbeta(1, 1)
    H ~ dunif(0, maxI/N)
    rateI ~ dgamma(aa, bb)
    
})

################################################################################
### Hill alarm models

### Fixed infectious period

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
    delta ~ dbeta(1, 1)
    nu ~ dgamma(0.1, 0.1)
    x0 ~ dunif(1, maxI)
    
})

### Exponential infectious period

SIR_hill_exp <-  nimbleCode({
    
    S[1] <- N - I0 
    I[1] <- I0
    
    probIR <- 1 - exp(-rateI)
    
    ### first time point for incidence based alarm
    # compute alarm
    smoothI[1] <- 0
    alarm[1] <- hillAlarm(smoothI[1], nu, x0, delta)
    
    probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    Istar[1] ~ dbin(probSI[1], S[1])
    Rstar[1] ~ dbin(probIR, I[1])
    
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
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    delta ~ dbeta(1, 1)
    nu ~ dgamma(0.1, 0.1)
    x0 ~ dunif(1, maxI)
    rateI ~ dgamma(aa, bb)
    
})

################################################################################
### Power alarm models

### Fixed infectious period

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

### Exponential infectious period

SIR_power_exp <-  nimbleCode({
    
    S[1] <- N - I0 
    I[1] <- I0
    
    probIR <- 1 - exp(-rateI)
    
    ### first time point for incidence based alarm
    # compute alarm
    smoothI[1] <- 0
    alarm[1] <- powerAlarm(smoothI[1], N, k)
    
    probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    Istar[1] ~ dbin(probSI[1], S[1])
    Rstar[1] ~ dbin(probIR, I[1])
    
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
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    k ~ dgamma(0.1, 0.1)
    rateI ~ dgamma(aa, bb)
    
})

################################################################################
### Spline alarm models

### Fixed infectious period

SIR_spline_fixed <-  nimbleCode({
    
    S[1] <- N - I0 
    I[1] <- I0
    Rstar[1:(lengthI-1)] <- 0
    Rstar[lengthI] <- I0
    
    ### first time point 
    # compute alarm
    smoothI[1] <- 0
    alarm[1] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[1])
    
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
        alarm[t] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[t])
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Istar[t] ~ dbin(probSI[t], S[t])
        
        Rstar[t + lengthI] <- Istar[t]
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
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
        knots[i] ~ dunif(min = 1, max = maxI)
    }
    
    # constrain knots to be ordered
    constrain_knots ~ dconstraint(knots[1] < knots[2])
    
})

### Exponential infectious period

SIR_spline_exp <-  nimbleCode({
    
    S[1] <- N - I0 
    I[1] <- I0
    
    probIR <- 1 - exp(-rateI)
    
    ### first time point 
    # compute alarm
    smoothI[1] <- 0
    alarm[1] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[1])
    
    probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    Istar[1] ~ dbin(probSI[1], S[1])
    Rstar[1] ~ dbin(probIR, I[1])
    
    # update S and I
    S[2] <- S[1] - Istar[1]
    I[2] <- I[1] + Istar[1] - Rstar[1]
    
    ### rest of time points
    for(t in 2:tau) {
        
        # compute alarm
        smoothI[t] <- movingAverage(Istar[1:(t-1)], bw)[t-1]
        alarm[t] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[t])
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Istar[t] ~ dbin(probSI[t], S[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
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
        knots[i] ~ dunif(min = 1, max = maxI)
    }
    rateI ~ dgamma(aa, bb)
    
    # constrain knots to be ordered
    constrain_knots ~ dconstraint(knots[1] < knots[2])
    
})


################################################################################
### Gaussian process alarm models

### Fixed infectious period

SIR_gp_fixed <-  nimbleCode({
    
    S[1] <- N - I0 
    I[1] <- I0
    Rstar[1:(lengthI-1)] <- 0
    Rstar[lengthI] <- I0
    
    ### first time point 
    # compute alarm
    smoothI[1] <- 0
    alarm[1] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[1])
    
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
        alarm[t] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[t])
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Istar[t] ~ dbin(probSI[t], S[t])
        
        Rstar[t + lengthI] <- Istar[t]
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    yAlarm[1] <- 0
    mu[2:n] <- mu0 * ones[2:n]
    cov[2:n, 2:n] <- sqExpCov(dists[2:n, 2:n], sigma, l)
    logit(yAlarm[2:n]) ~ dmnorm(mu[2:n], cov = cov[2:n, 2:n])
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    sigma ~ dgamma(100, 50)
    l ~ dinvgamma(c, d)
    
})

### Exponential infectious period

SIR_gp_exp <-  nimbleCode({
    
    S[1] <- N - I0 
    I[1] <- I0
    
    probIR <- 1 - exp(-rateI)
    
    ### first time point 
    # compute alarm
    smoothI[1] <- 0
    alarm[1] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[1])
    
    probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    Istar[1] ~ dbin(probSI[1], S[1])
    Rstar[1] ~ dbin(probIR, I[1])
    
    # update S and I
    S[2] <- S[1] - Istar[1]
    I[2] <- I[1] + Istar[1] - Rstar[1]
    
    ### rest of time points
    for(t in 2:tau) {
        
        # compute alarm
        smoothI[t] <- movingAverage(Istar[1:(t-1)], bw)[t-1]
        alarm[t] <- nim_approx(xAlarm[1:n], yAlarm[1:n], smoothI[t])
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        Istar[t] ~ dbin(probSI[t], S[t])
        Rstar[t] ~ dbin(probIR, I[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Rstar[t]
        
    }
    
    yAlarm[1] <- 0
    mu[2:n] <- mu0 * ones[2:n]
    cov[2:n, 2:n] <- sqExpCov(dists[2:n, 2:n], sigma, l)
    logit(yAlarm[2:n]) ~ dmnorm(mu[2:n], cov = cov[2:n, 2:n])
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    sigma ~ dgamma(100, 50)
    l ~ dinvgamma(c, d)
    rateI ~ dgamma(aa, bb)
    
})
