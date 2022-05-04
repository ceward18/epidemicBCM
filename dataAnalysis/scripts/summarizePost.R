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

summarizePost <- function(resThree, incData, N, I0, R0, lengthI, 
                          alarmBase, alarmFit, infPeriod, smoothWindow) {
  
  if (!alarmFit %in% c('betatSpline', 'basic')) {
    paramSamples1 <- resThree[[1]][,-grep('alarm|yAlarm|Rstar', colnames(resThree[[1]]))]
    paramSamples2 <- resThree[[2]][,-grep('alarm|yAlarm|Rstar', colnames(resThree[[2]]))]
    paramSamples3 <- resThree[[3]][,-grep('alarm|yAlarm|Rstar', colnames(resThree[[3]]))]
  } else if (alarmFit == 'betatSpline') {
    paramSamples1 <- resThree[[1]][,-grep('beta|Rstar', colnames(resThree[[1]]))]
    paramSamples2 <- resThree[[2]][,-grep('beta|Rstar', colnames(resThree[[2]]))]
    paramSamples3 <- resThree[[3]][,-grep('beta|Rstar', colnames(resThree[[3]]))]
  } else if (alarmFit == 'basic') {
    paramSamples1 <- resThree[[1]][,-grep('Rstar', colnames(resThree[[1]])), drop = F]
    paramSamples2 <- resThree[[2]][,-grep('Rstar', colnames(resThree[[2]])), drop = F]
    paramSamples3 <- resThree[[3]][,-grep('Rstar', colnames(resThree[[3]])), drop = F]
  }
  
  RstarSamples1 <-  resThree[[1]][,grep('Rstar', colnames(resThree[[1]]))]
  RstarSamples2 <-  resThree[[2]][,grep('Rstar', colnames(resThree[[2]]))]
  RstarSamples3 <-  resThree[[3]][,grep('Rstar', colnames(resThree[[3]]))]
  
  # combine posterior parameters with posterior Rstar
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
    maxI <- ceiling(max(movingAverage(incData, smoothWindow)) + 1)
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
  ### Posterior distribution of reproductive number
  ### R0 = beta(1-a) * lengthI, beta * lengthI
  
  if (!alarmFit %in% c('betatSpline', 'basic')) {
    
    # want R0 by xAlarm and epidemic time
    betaSamples <- paramsPost[,'beta']
    
    # by xAlarm
    R0SamplesAlarm <- sapply(1:ncol(alarmSamples), function(k)
      betaSamples[k] * (1 - alarmSamples[,k]) * lengthI)
    
    postMeans <- rowMeans(R0SamplesAlarm)
    postCI <- apply(R0SamplesAlarm, 1, quantile, probs = c(0.025, 0.975))
    postR0Alarm <- data.frame(xAlarm = xAlarm,
                              mean = postMeans,
                              lower = postCI[1,],
                              upper = postCI[2,])
    rownames(postR0Alarm) <- NULL
    
    # by time
    alarmTimeSamples1 <- t(resThree[[1]][,grep('alarm', colnames(resThree[[1]]))])
    alarmTimeSamples2 <- t(resThree[[2]][,grep('alarm', colnames(resThree[[2]]))])
    alarmTimeSamples3 <- t(resThree[[3]][,grep('alarm', colnames(resThree[[3]]))])
    alarmTimeSamples <- cbind(alarmTimeSamples1, alarmTimeSamples2, alarmTimeSamples3)
    
    R0Samples <- sapply(1:ncol(alarmTimeSamples), function(k)
      betaSamples[k] * (1 - alarmTimeSamples[,k]) * lengthI)
    
  } else if (alarmFit == 'betatSpline') {
    
    postR0Alarm <- data.frame(xAlarm = NA, 
                              mean = NA,
                              lower = NA,
                              upper = NA)
    
    # beta samples already calculated (matrix)
    R0Samples <- betaSamples * lengthI
    
  } else if (alarmFit == 'basic') {
    
    postR0Alarm <- data.frame(xAlarm = NA, 
                              mean = NA,
                              lower = NA,
                              upper = NA)
    
    betaSamples <- paramsPost[,'beta']
    R0Samples <- betaSamples * lengthI
    
    R0Samples <- matrix(rep(R0Samples, length(incData)),
                        nrow = length(incData),
                        ncol = length(R0Samples),
                        byrow = T)
    
  }
  
  postMeans <- rowMeans(R0Samples)
  postCI <- apply(R0Samples, 1, quantile, probs = c(0.025, 0.975))
  postR0 <- data.frame(time = 1:length(incData),
                       mean = postMeans,
                       lower = postCI[1,],
                       upper = postCI[2,])
  rownames(postR0) <- NULL
  
  ##############################################################################
  ### posterior predictive forecasting
  
  # not including for now as it is time consuming
  # if (!alarmFit %in% c('betatSpline')) {
  #   
  #   postPredInc <- postPred(incData, N, I0, R0, lengthI,
  #                           alarmFit, infPeriod, smoothWindow, 
  #                           cbind(paramsPost, RstarPost), alarmSamples)
  # 
  #   postMean <- rowMeans(postPredInc)
  #   postCI <- apply(postPredInc, 1, quantile, probs = c(0.025, 0.975))
  #   postEpiPred <- data.frame(time = 1:length(postMean) + length(incData),
  #                             mean = postMean,
  #                             lower = postCI[1,],
  #                             upper = postCI[2,])
  # } else {
  #   postEpiPred <- data.frame(time = NA, 
  #                             mean = NA,
  #                             lower = NA,
  #                             upper = NA)
  # }
  
  postEpiPred <- data.frame(time = NA, 
                            mean = NA,
                            lower = NA,
                            upper = NA)
  
  ##############################################################################
  ### WAIC values
  
  # samples to use for WAIC calculation differ by model
  if (alarmFit %in% c('thresh', 'hill', 'power', 'basic', 'betatSpline')) {
    
    samples <- rbind(resThree[[1]], resThree[[2]], resThree[[3]])
    
  } else if (alarmFit %in% c('gp', 'spline')) {
    
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
  
  waic <- getWAIC(samples = samples, incData = incData, 
                  N = N, I0 = I0, R0 = R0, lengthI = lengthI,
                  infPeriod = infPeriod, 
                  alarmFit = alarmFit, smoothWindow = smoothWindow)
  
  ### output
  list(gdiag = gdiag,
       postParams = postParams,
       postAlarm = postAlarm,
       postEpiPred = postEpiPred,
       postBeta = postBeta,
       postR0Alarm = postR0Alarm,
       postR0 = postR0,
       waic= waic)
  
  
  
}





