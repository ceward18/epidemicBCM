################################################################################
# summarize posterior samples from three chains
# gelman-rubin to assess convergence
# summarize posterior means and credible intervals:
#   model parameters   
#   alarm function
# posterior prediction mean and credible intervals
################################################################################

library(coda)
library(nimble)

# source relevant scripts
source('./scripts/modelCodes.R')
source('./scripts/getModelInputs.R')
source('./scripts/postPredFit.R')
source('./scripts/getWAIC.R')

summarizePost <- function(resThree, incData, deathData, smoothI, smoothWindow,
                          N, E0, I0, R0, alarmFit) {
    
    if (!alarmFit %in% c('betatSpline', 'basic')) {
        paramSamples1 <- resThree[[1]][,-grep('alarm|yAlarm|Estar|Istar|Rstar|R0|SIR_init\\[3\\]',
                                              colnames(resThree[[1]]))]
        paramSamples2 <- resThree[[2]][,-grep('alarm|yAlarm|Estar|Istar|Rstar|R0|SIR_init\\[3\\]', 
                                              colnames(resThree[[2]]))]
        paramSamples3 <- resThree[[3]][,-grep('alarm|yAlarm|Estar|Istar|Rstar|R0|SIR_init\\[3\\]', 
                                              colnames(resThree[[3]]))]
    } else if (alarmFit == 'betatSpline') {
        paramSamples1 <- resThree[[1]][,-grep('beta|Estar|Istar|Rstar|R0|SIR_init\\[3\\]', 
                                              colnames(resThree[[1]]))]
        paramSamples2 <- resThree[[2]][,-grep('beta|Estar|Istar|Rstar|R0|SIR_init\\[3\\]', 
                                              colnames(resThree[[2]]))]
        paramSamples3 <- resThree[[3]][,-grep('beta|Estar|Istar|Rstar|R0|SIR_init\\[3\\]', 
                                              colnames(resThree[[3]]))]
    } else if (alarmFit == 'basic') {
        paramSamples1 <- resThree[[1]][,-grep('Estar|Istar|Rstar|Rstar|R0|SIR_init\\[3\\]', 
                                              colnames(resThree[[1]])), drop = F]
        paramSamples2 <- resThree[[2]][,-grep('Estar|Istar|Rstar|Rstar|R0|SIR_init\\[3\\]', 
                                              colnames(resThree[[2]])), drop = F]
        paramSamples3 <- resThree[[3]][,-grep('Estar|Istar|Rstar|Rstar|R0|SIR_init\\[3\\]', 
                                              colnames(resThree[[3]])), drop = F]
    }
    
    # Estar Posterior
    EstarSamples1 <-  resThree[[1]][,grep('Estar', colnames(resThree[[1]]))]
    EstarSamples2 <-  resThree[[2]][,grep('Estar', colnames(resThree[[2]]))]
    EstarSamples3 <-  resThree[[3]][,grep('Estar', colnames(resThree[[3]]))]
    EstarPost <- rbind(EstarSamples1, EstarSamples2, EstarSamples3)
    
    # Istar Posterior
    IstarSamples1 <-  resThree[[1]][,grep('Istar', colnames(resThree[[1]]))]
    IstarSamples2 <-  resThree[[2]][,grep('Istar', colnames(resThree[[2]]))]
    IstarSamples3 <-  resThree[[3]][,grep('Istar', colnames(resThree[[3]]))]
    IstarPost <- rbind(IstarSamples1, IstarSamples2, IstarSamples3)
    
    # Rstar Posterior
    RstarSamples1 <-  resThree[[1]][,grep('Rstar', colnames(resThree[[1]]))]
    RstarSamples2 <-  resThree[[2]][,grep('Rstar', colnames(resThree[[2]]))]
    RstarSamples3 <-  resThree[[3]][,grep('Rstar', colnames(resThree[[3]]))]
    RstarPost <- rbind(RstarSamples1, RstarSamples2, RstarSamples3)
    
    ##############################################################################
    ### gelman-rubin
    res_mcmc <- mcmc.list(mcmc(paramSamples1),
                          mcmc(paramSamples2),
                          mcmc(paramSamples3))
    gdiag <- data.frame(gelman.diag(res_mcmc, multivariate = F)$psrf)
    colnames(gdiag) <- c('gr', 'grUpper')
    gdiag$param <- rownames(gdiag)
    rownames(gdiag) <- NULL
    
    ##############################################################################
    ### posterior mean and 95% CI for parameters
    paramsPost <- rbind(paramSamples1, paramSamples2, paramSamples3)
    postMeans <- colMeans(paramsPost)
    postCI <- apply(paramsPost, 2, quantile, probs = c(0.025, 0.975))
    postParams <- data.frame(param = names(postMeans),
                             mean = postMeans,
                             lower = postCI[1,],
                             upper = postCI[2,])
    rownames(postParams) <- NULL
    
    ##############################################################################
    ### posterior distribution of alarm function 
    if (!alarmFit %in% c('betatSpline', 'basic')) {
        alarmSamples1 <- t(resThree[[1]][,grep('yAlarm', colnames(resThree[[1]]))])
        alarmSamples2 <- t(resThree[[2]][,grep('yAlarm', colnames(resThree[[2]]))])
        alarmSamples3 <- t(resThree[[3]][,grep('yAlarm', colnames(resThree[[3]]))])
        alarmSamples <- cbind(alarmSamples1, alarmSamples2, alarmSamples3)
        
        postMeans <- rowMeans(alarmSamples)
        postCI <- apply(alarmSamples, 1, quantile, probs = c(0.025, 0.975))
        
        # get xAlarm
        if (alarmFit == 'gp') {
            n <- 10
        } else {
            n <- 50
        }
        maxI <- ceiling(max(smoothI))
        xAlarm <- seq(0, maxI, length.out = n)
        
        postAlarm <- data.frame(xAlarm = xAlarm, 
                                mean = postMeans,
                                lower = postCI[1,],
                                upper = postCI[2,])
        
        rownames(postAlarm) <- NULL
        
    } else {
        postAlarm <- data.frame(xAlarm = NA, 
                                mean = NA,
                                lower = NA,
                                upper = NA)
    }
    
    ##############################################################################
    ### posterior distribution of beta[t] when estimated directly
    if (alarmFit == 'betatSpline') {
        betaSamples1 <- t(resThree[[1]][,grep('beta', colnames(resThree[[1]]))])
        betaSamples2 <- t(resThree[[2]][,grep('beta', colnames(resThree[[2]]))])
        betaSamples3 <- t(resThree[[3]][,grep('beta', colnames(resThree[[3]]))])
        betaSamples <- cbind(betaSamples1, betaSamples2, betaSamples3)
        
        postMeans <- rowMeans(betaSamples)
        postCI <- apply(betaSamples, 1, quantile, probs = c(0.025, 0.975))
        postBeta <- data.frame(time = 1:length(postMeans),
                               mean = postMeans,
                               lower = postCI[1,],
                               upper = postCI[2,])
    } else {
        postBeta <- data.frame(time = NA, 
                               mean = NA,
                               lower = NA,
                               upper = NA)
    }
    
    ##############################################################################
    ### Posterior distribution of effective reproductive number
    
    R0Samples1 <-  resThree[[1]][,grep('R0_eff', colnames(resThree[[1]]))]
    R0Samples2 <-  resThree[[2]][,grep('R0_eff', colnames(resThree[[2]]))]
    R0Samples3 <-  resThree[[3]][,grep('R0_eff', colnames(resThree[[3]]))]
    
    # combine posterior parameters with posterior R0
    R0Samples <- rbind(R0Samples1, R0Samples2, R0Samples3)
    
    postMeans <- colMeans(R0Samples)
    postCI <- apply(R0Samples, 2, quantile, probs = c(0.025, 0.975))
    postR0 <- data.frame(time = 1:length(postMeans),
                         mean = postMeans,
                         lower = postCI[1,],
                         upper = postCI[2,])
    rownames(postR0) <- NULL

    ##############################################################################
    ### posterior predictive model fit
    
    postPredObs <- postPredFit(incData = incData, deathData = deathData,
                               smoothI = smoothI, smoothWindow = smoothWindow, 
                               N = N, E0 = E0, I0 = I0, R0 = R0,  
                               alarmFit = alarmFit, 
                               paramsPost = paramsPost, alarmSamples = alarmSamples)
    
    postMean <- rowMeans(postPredObs)
    postCI <- apply(postPredObs, 1, quantile, probs = c(0.025, 0.975))
    postPredFit <- data.frame(time = 1:length(incData),
                              mean = postMean,
                              lower = postCI[1,],
                              upper = postCI[2,])
   
    
    ##############################################################################
    ### WAIC values
    
    # samples to use for WAIC calculation differ by model
    if (alarmFit %in% c('thresh', 'hill', 'power', 'basic', 'betatSpline')) {
        
        samples <- rbind(resThree[[1]], resThree[[2]], resThree[[3]])
        
    } else if (alarmFit %in% c('gp', 'spline', 'splineFixKnot')) {
        
        # includes Rstar, so don't need to add later
        samps1 <- resThree[[1]][,-grep('alarm', colnames(resThree[[1]]))]
        samps2 <- resThree[[2]][,-grep('alarm', colnames(resThree[[2]]))]
        samps3 <- resThree[[3]][,-grep('alarm', colnames(resThree[[3]]))]
        
        samples <- rbind(samps1, samps2, samps3)
        
        if (alarmFit == 'gp') {
            # need yAlarm on the logit scale
            yAlarmCols <- grep('yAlarm', colnames(samples))
            samples[,yAlarmCols] <- logit(samples[,yAlarmCols])
            colnames(samples)[yAlarmCols] <- paste0('logit_yAlarm[', 1:length(yAlarmCols), ']')
        }
        
    } 
    
    waic <- getWAIC(samples = samples, incData = incData, deathData = deathData,
                    smoothI = smoothI, smoothWindow = smoothWindow,
                    N = N, E0 = E0, I0 = I0, R0 = R0, alarmFit = alarmFit)
    
    ##############################################################################
    
    ### output
    list(gdiag = gdiag,
         postParams = postParams,
         postAlarm = postAlarm,
         postPredFit = postPredFit,
         postBeta = postBeta,
         postR0 = postR0,
         waic = waic)
    
    
    
}





