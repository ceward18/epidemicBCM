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
source('./scripts/postPred.R')
source('./scripts/getWAIC.R')

summarizePost <- function(resThree, incData, alarmBase, alarmFit, smoothWindow, prior) {
    
    
    if (!alarmFit %in% c('betat', 'basic', 'betatSpline')) {
        paramsamples1 <- resThree[[1]][,-grep('alarm|yAlarm|Rstar', colnames(resThree[[1]]))]
        paramsamples2 <- resThree[[2]][,-grep('alarm|yAlarm|Rstar', colnames(resThree[[2]]))]
        paramsamples3 <- resThree[[3]][,-grep('alarm|yAlarm|Rstar', colnames(resThree[[3]]))]
    } else if (alarmFit %in% c('betat', 'betatSpline')) {
        paramsamples1 <- resThree[[1]][,-grep('beta\\[|Rstar', colnames(resThree[[1]]))]
        paramsamples2 <- resThree[[2]][,-grep('beta\\[|Rstar', colnames(resThree[[2]]))]
        paramsamples3 <- resThree[[3]][,-grep('beta\\[|Rstar', colnames(resThree[[3]]))]
    } else if (alarmFit == 'basic') {
        paramsamples1 <- resThree[[1]][,-grep('Rstar', colnames(resThree[[1]]))]
        paramsamples2 <- resThree[[2]][,-grep('Rstar', colnames(resThree[[2]]))]
        paramsamples3 <- resThree[[3]][,-grep('Rstar', colnames(resThree[[3]]))]
    }
    
    RstarSamples1 <-  resThree[[1]][,grep('Rstar', colnames(resThree[[1]]))]
    RstarSamples2 <-  resThree[[2]][,grep('Rstar', colnames(resThree[[2]]))]
    RstarSamples3 <-  resThree[[3]][,grep('Rstar', colnames(resThree[[3]]))]
    
    # combine posterior parameters with posterior Rstar
    RstarPost <- rbind(RstarSamples1, RstarSamples2, RstarSamples3)
    
    ##############################################################################
    ### gelman-rubin
    res_mcmc <- mcmc.list(mcmc(paramsamples1), 
                          mcmc(paramsamples2), 
                          mcmc(paramsamples3))
    gdiag <- data.frame(gelman.diag(res_mcmc, multivariate = F)$psrf)
    colnames(gdiag) <- c('gr', 'grUpper')
    gdiag$param <- rownames(gdiag)
    rownames(gdiag) <- NULL
    
    ##############################################################################
    ### posterior mean and 95% CI for parameters
    paramsPost <- rbind(paramsamples1, paramsamples2, paramsamples3)
    postMeans <- colMeans(paramsPost)
    postCI <- apply(paramsPost, 2, quantile, probs = c(0.025, 0.975))
    postParams <- data.frame(param = names(postMeans),
                             mean = postMeans,
                             lower = postCI[1,],
                             upper = postCI[2,])
    rownames(postParams) <- NULL
    
    ##############################################################################
    ### posterior distribution of alarm function
    if (!alarmFit %in% c('betat', 'basic', 'betatSpline')) {
        alarmSamples1 <- t(resThree[[1]][,grep('yAlarm', colnames(resThree[[1]]))])
        alarmSamples2 <- t(resThree[[2]][,grep('yAlarm', colnames(resThree[[2]]))])
        alarmSamples3 <- t(resThree[[3]][,grep('yAlarm', colnames(resThree[[3]]))])
        alarmSamples <- cbind(alarmSamples1, alarmSamples2, alarmSamples3)
        
        postMeans <- rowMeans(alarmSamples)
        postCI <- apply(alarmSamples, 1, quantile, probs = c(0.025, 0.975))
        
        # get xAlarm
        modelInputs <- getModelInput(alarmFit, incData, smoothWindow, prior)
        
        xAlarm <- modelInputs$xAlarm
        
        postAlarm <- data.frame(xAlarm = xAlarm, 
                                mean = postMeans,
                                lower = postCI[1,],
                                upper = postCI[2,])
    } else {
        postAlarm <- data.frame(xAlarm = NA, 
                                mean = NA,
                                lower = NA,
                                upper = NA)
    }
    
    ##############################################################################
    ### posterior distribution of beta[t] when estimated directly
    if (alarmFit %in% c('betat', 'betatSpline')) {
        betaSamples1 <- t(resThree[[1]][,grep('beta\\[', colnames(resThree[[1]]))])
        betaSamples2 <- t(resThree[[2]][,grep('beta\\[', colnames(resThree[[2]]))])
        betaSamples3 <- t(resThree[[3]][,grep('beta\\[', colnames(resThree[[3]]))])
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
    ### posterior predictive forecasting assuming first 50 days have been observed
    if (!alarmFit %in% c('betat', 'betatSpline')) {
        
        postPredInc <- postPred(incData = incData, alarmFit = alarmFit, 
                                smoothWindow = smoothWindow, prior = prior,
                                paramsPost = paramsPost, alarmSamples = alarmSamples, 
                                RstarPost = RstarPost)
        
        postMean <- rowMeans(postPredInc)
        postCI <- apply(postPredInc, 1, quantile, probs = c(0.025, 0.975))
        postEpiPred <- data.frame(time = 51:100,
                                  mean = postMean,
                                  lower = postCI[1,],
                                  upper = postCI[2,])
    } else {
        postEpiPred <- data.frame(time = NA, 
                                  mean = NA,
                                  lower = NA,
                                  upper = NA)
    }
    
    ##############################################################################
    ### WAIC values
    
    # samples to use for WAIC calculation differ by model
    if (alarmFit %in% c('thresh', 'hill', 'power', 'basic')) {
        
        samples <- paramsPost
        
    } else if (alarmFit %in% c('gp', 'spline')) {
        
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
        
        
    } else if (alarmFit %in% c('betat', 'betatSpline')) {
        
        samples <- rbind(resThree[[1]], resThree[[2]], resThree[[3]])
        
        # samples needs to be log_beta to get WAIC
        if (alarmFit == 'betat') {
            betaCols <- grep('beta\\[', colnames(samples))
            samples[,betaCols] <- log(samples[,betaCols])
            colnames(samples)[betaCols] <- paste0('log_beta[', 1:length(betaCols), ']')
        }
        
    }
    
    if (!alarmFit %in% c('gp', 'spline', 'betat', 'betatSpline')) {
        
        samples <- cbind(samples, RstarPost)
        
    }
    
    waic <- getWAIC(samples = samples, incData = incData, alarmFit = alarmFit,
                    smoothWindow = smoothWindow, prior = prior)
    
    ### output
    list(gdiag = gdiag,
         postParams = postParams,
         postAlarm = postAlarm,
         postEpiPred = postEpiPred,
         postBeta = postBeta,
         waic= waic)
    
    
    
}





