################################################################################
# Figures in paper
################################################################################

# load libraries
library(nimble)
library(ggplot2)

source('./scripts/modelCodes.R')

################################################################################
### Figure 1 - Example alarm functions
################################################################################


x <- seq(0, 1000, 1)
N <- 10000


### power alarm 

kVals <- c(0.001, 0.01, 0.1)

yPower1 <- powerAlarm(x, N, k = kVals[1])
yPower2 <- powerAlarm(x, N, k = kVals[2])
yPower3 <- powerAlarm(x, N, k = kVals[3])

### threshold alarm 

greek1 <- "delta"
greek2 <- 'nu'
greek3 <- 'x[0]'

deltaVals <- c(0.9, 0.6, 0.4)
hVals <- c(100, 600, 250)

yThresh1 <- thresholdAlarm(x, N, delta = deltaVals[1], H = hVals[1]/N)
yThresh2 <- thresholdAlarm(x, N, delta = deltaVals[2], H = hVals[2]/N)
yThresh3 <- thresholdAlarm(x, N, delta = deltaVals[3], H = hVals[3]/N)

.expressions <- mapply(sprintf, greek1, "=", deltaVals, 
                       "H =", hVals ,
                       MoreArgs = list(fmt = '%s~"%s %s,"~"%s %s"'))
legend_expressions_thresh <-parse(text = .expressions)

### hill alarm
deltaVals <- c(1, 0.7, 0.3)
nuVals <- c(1, 8, 4)
x0Vals <- c(100, 600, 250)

.expressions <- mapply(sprintf, greek1, "=", deltaVals, 
                       greek2, "=", nuVals, 
                       greek3, "=", x0Vals,
                       MoreArgs = list(fmt = '%s~"%s %.1f,"~%s~"%s %s,"~%s~"%s %s"'))
legend_expressions_hill <-parse(text = .expressions)

yHill1 <- hillAlarm(x, nu = nuVals[1], x0 = x0Vals[1], delta = deltaVals[1])
yHill2 <- hillAlarm(x, nu = nuVals[2], x0 = x0Vals[2], delta = deltaVals[2])
yHill3 <- hillAlarm(x, nu = nuVals[3], x0 = x0Vals[3], delta = deltaVals[3])


pal <- c('#1E88E5', '#FFC107', '#D81B60' )
lineTypes <- c(1,2,4)
cexMain <- 2.4
cexAxis <- 1.8
cexLab <- 2
cexLeg <- 1.8
lineWidths <- 4
lineLength <- 4

layoutMat <- matrix(1:6, nrow = 2, ncol = 3)

pdf('./figures/fig1_ex_alarms.pdf', width = 12, height = 5)
layout(mat = layoutMat,
       heights = c(3, 1), 
       widths = c(1, 1, 1)) 

# power
par(mar=c(5.1, 5.1, 4.1, 1.1))
plot(x, yPower1, type = 'l', col = pal[1], lwd = lineWidths, ylim = c(0, 1), lty = lineTypes[1],
     main = 'Power Alarm', xlab = 'Incidence', ylab= 'Alarm Value',
     cex.main = cexMain, cex.axis = cexAxis, cex.lab = cexLab)
lines(x, yPower2, col = pal[2], lwd = lineWidths, lty = lineTypes[2])
lines(x, yPower3, col = pal[3], lwd = lineWidths, lty = lineTypes[3])

par(mar=c(0, 0, 0, 0))
plot.new()
legend('top', bty = 'n',
       paste0('k = ', kVals),
       col = pal, lty = lineTypes,
       lwd = lineWidths,
       cex = cexLeg,
       seg.len= lineLength)

# threshold
par(mar=c(5.1, 5.1, 4.1, 1.1))
plot(x, yThresh1, type = 'l', col = pal[1], lwd = lineWidths, ylim = c(0, 1), lty = lineTypes[1],
     main = 'Threshold Alarm', xlab = 'Incidence', ylab= 'Alarm Value',
     cex.main = cexMain, cex.axis = cexAxis, cex.lab = cexLab)
lines(x, yThresh2, col = pal[2], lwd = lineWidths, lty = lineTypes[2])
lines(x, yThresh3, col = pal[3], lwd = lineWidths, lty = lineTypes[3])

par(mar=c(0, 0, 0, 0))
plot.new()
legend('top', legend=legend_expressions_thresh,
       bty = 'n',
       col = pal, lty = lineTypes,
       lwd = lineWidths,
       cex = cexLeg,
       seg.len= lineLength)

# hill
par(mar=c(5.1, 5.1, 4.1, 1.1))
plot(x, yHill1, type = 'l', col = pal[1], lwd = lineWidths, ylim = c(0, 1), lty = lineTypes[1],
     main = 'Hill Alarm', xlab = 'Incidence', ylab= 'Alarm Value',
     cex.main = cexMain, cex.axis = cexAxis, cex.lab = cexLab)
lines(x, yHill2, col = pal[2], lwd = lineWidths, lty = lineTypes[2])
lines(x, yHill3, col = pal[3], lwd = lineWidths, lty = lineTypes[3])

par(mar=c(0, 0, 0, 0))
plot.new()
legend('top', legend=legend_expressions_hill, 
       bty = 'n',
       col = pal, lty = lineTypes,
       lwd = lineWidths,
       cex = cexLeg,
       seg.len= lineLength)


dev.off()


################################################################################
### Figure 2 - Example epidemics from alarm functions
################################################################################

### Want to show:

# Threshold/Hill flatten curve
# Power lowers peak location
# smoothing affects smoothness/timing of peaks


### constants for all
N <- 1e6
I0 <- 5
S0 <- N - I0
tau <- 150
rateI <- 1/5
beta <- 3/5   #R0 = 3

dataNodes <- c(paste0('Istar[', 1:tau, ']'),
               paste0('Rstar[', 1:tau, ']'))


### Define simulator functions for each alarm/smoothing combination

constants1 <- list(N = N, tau = tau, I0 = I0, bw = 1, aa = 1, bb = 1)
constants14 <- list(N = N, tau = tau, I0 = I0, bw = 14, aa = 1, bb = 1)
constants30 <- list(N = N, tau = tau, I0 = I0, bw = 30, aa = 1, bb = 1)

## Power alarm

# 1 day average
powerModel <- nimbleModel(code = SIR_power_sim,
                           constants = constants1)
simPower1 <- simulator(powerModel, dataNodes)

# 14 day average
powerModel <- nimbleModel(code = SIR_power_sim,
                           constants = constants14)
simPower14 <- simulator(powerModel, dataNodes)

# 30 day average
powerModel <- nimbleModel(code = SIR_power_sim,
                           constants = constants30)
simPower30 <- simulator(powerModel, dataNodes)


## Threshold alarm

# 1 day average
threshModel <- nimbleModel(code = SIR_thresh_sim,
                           constants = constants1)
simThresh1 <- simulator(threshModel, dataNodes)

# 14 day average
threshModel <- nimbleModel(code = SIR_thresh_sim,
                           constants = constants14)
simThresh14 <- simulator(threshModel, dataNodes)

# 30 day average
threshModel <- nimbleModel(code = SIR_thresh_sim,
                           constants = constants30)
simThresh30 <- simulator(threshModel, dataNodes)

## Hill alarm

# 1 day average
hillModel <- nimbleModel(code = SIR_hill_sim,
                           constants = constants1)
simHill1 <- simulator(threshModel, dataNodes)

# 14 day average
hillModel <- nimbleModel(code = SIR_hill_sim,
                           constants = constants14)
simHill14 <- simulator(hillModel, dataNodes)

# 30 day average
hillModel <- nimbleModel(code = SIR_hill_sim,
                           constants = constants30)
simHill30 <- simulator(hillModel, dataNodes)


### Threshold flattening curve

# using smoothing of 1 day
# compared to no BC model


# no BC Model
trueVals <- c(beta = beta, 
              delta = 0,
              rateI = rateI,
              H = 1)

set.seed(123)
noBC <- simThresh1$run(trueVals, 1)[1:tau]

# delta = 0.2
trueVals <- c(beta = beta, 
              delta = 0.2,
              rateI = rateI,
              H = 100/N)

set.seed(123)
thresh1 <- simThresh1$run(trueVals, 1)[1:tau]

# delta = 0.4
trueVals <- c(beta = beta, 
              delta = 0.4,
              rateI = rateI,
              H = 100/N)

set.seed(123)
thresh2 <- simThresh1$run(trueVals, 1)[1:tau]

# delta = 0.6
trueVals <- c(beta = beta, 
              delta = 0.6,
              rateI = rateI,
              H = 100/N)

set.seed(123)
thresh3 <- simThresh1$run(trueVals, 1)[1:tau]

toPlot <- data.frame(alarm = c(rep('No behavior change', tau * 3),
                               rep('Threshold', tau * 3)),
                     time = rep(1:tau, 6),
                     delta = c(rep(rep(c(0.2, 0.4, 0.6), each = tau), 1)),
                     Istar = c(noBC, noBC, noBC,
                               thresh1, thresh2, thresh3))


ggplot(toPlot, aes(x = time, y = Istar, group = alarm, color = alarm)) + 
    geom_line(size = 1) +
    facet_wrap(~delta) + 
    theme_bw()


# threshold vs hill smoothness


# threshold 1, 14, 30
trueVals <- c(beta = beta, 
              delta = 0.9,
              rateI = rateI,
              H = 100/N)

set.seed(123)
thresh1 <- simThresh1$run(trueVals, 1)[1:tau]
set.seed(123)
thresh2 <- simThresh14$run(trueVals, 1)[1:tau]
set.seed(123)
thresh3 <- simThresh30$run(trueVals, 1)[1:tau]

# hill 1, 14, 30
trueVals <- c(beta = beta, 
              delta = 0.9,
              nu = 3,
              rateI = rateI,
              x0 = 250)

set.seed(123)
hill1 <- simHill1$run(trueVals, 1)[1:tau]
set.seed(123)
hill2 <- simHill14$run(trueVals, 1)[1:tau]
set.seed(123)
hill3 <- simHill30$run(trueVals, 1)[1:tau]


toPlot <- data.frame(alarm = c(rep('Threshold', tau * 3),
                               rep('Hill', tau * 3)),
                     time = rep(1:tau, 6),
                     smoothWindow = rep(rep(c(1, 14, 30), each = tau), 2),
                     Istar = c(thresh1, thresh2, thresh3,
                               hill1, hill2, hill3))

ggplot(subset(toPlot, time <=100), 
       aes(x = time, y = Istar, group = alarm, color = alarm)) + 
    geom_line(size = 1) +
    facet_wrap(~smoothWindow) + 
    theme_bw()



### Power
# generate power alarms for 1 day, 14 days, 30 days

powerModel1 <- nimbleModel(code = SIR_power_sim,
                           constants = list(N = N, 
                                            tau = tau,
                                            I0 = I0,
                                            bw = 1,
                                            aa = 1,
                                            bb = 1))
powerModel14 <- nimbleModel(code = SIR_power_sim,
                            constants = list(N = N, 
                                             tau = tau,
                                             I0 = I0,
                                             bw = 14,
                                             aa = 1,
                                             bb = 1))
powerModel30 <- nimbleModel(code = SIR_power_sim,
                            constants = list(N = N, 
                                             tau = tau,
                                             I0 = I0,
                                             bw = 30,
                                             aa = 1,
                                             bb = 1))

simPower1 <- simulator(powerModel1, dataNodes)
simPower14 <- simulator(powerModel14, dataNodes)
simPower30 <- simulator(powerModel30, dataNodes)

kVals <- c(0.001, 0.01, 0.1)

powerSims <- data.frame(alarm = 'power',
                        time = rep(1:tau, 3 * length(kVals)),
                        smoothWindow = rep(rep(c(1, 14, 30), each = tau), length(kVals)),
                        k = rep(kVals, each = tau * 3),
                        Istar = NA)

for (i in 1:length(kVals)) {
    trueVals <- c(beta = beta, 
                  k = kVals[i],
                  rateI = rateI)
    
    set.seed(123)
    powerSims$Istar[which(powerSims$k == kVals[i] & 
                              powerSims$smoothWindow == 1)] <- simPower1$run(trueVals, 1)[1:tau]
    set.seed(123)
    powerSims$Istar[which(powerSims$k == kVals[i] & 
                              powerSims$smoothWindow == 14)] <- simPower14$run(trueVals, 1)[1:tau]
    set.seed(123)
    powerSims$Istar[which(powerSims$k == kVals[i] & 
                              powerSims$smoothWindow == 30)] <- simPower30$run(trueVals, 1)[1:tau]
    
    
}


ggplot(powerSims, aes(x = time,  y = Istar, 
                      group = factor(smoothWindow), col = factor(smoothWindow))) + 
    geom_line() + 
    facet_wrap(~k, scales = 'free')

ggplot(powerSims, aes(x = time,  y = Istar, 
                      group = factor(k), col = factor(k))) + 
    geom_line() + 
    facet_grid(k~smoothWindow, scales = 'free')


### Threshold
# generate threshold alarms for 1 day, 14 days, 30 days

threshModel1 <- nimbleModel(code = SIR_thresh_sim,
                            constants = list(N = N, 
                                             tau = tau,
                                             I0 = I0,
                                             bw = 1,
                                             aa = 1,
                                             bb = 1))
threshModel14 <- nimbleModel(code = SIR_thresh_sim,
                             constants = list(N = N, 
                                              tau = tau,
                                              I0 = I0,
                                              bw = 14,
                                              aa = 1,
                                              bb = 1))
threshModel30 <- nimbleModel(code = SIR_thresh_sim,
                             constants = list(N = N, 
                                              tau = tau,
                                              I0 = I0,
                                              bw = 30,
                                              aa = 1,
                                              bb = 1))

simThresh1 <- simulator(threshModel1, dataNodes)
simThresh14 <- simulator(threshModel14, dataNodes)
simThresh30 <- simulator(threshModel30, dataNodes)

deltaVals <- c(0.9, 0.6, 0.4)
hVals <- c(100, 600, 250)/N

threshSims <- data.frame(alarm = 'thresh',
                         time = rep(1:tau, 3 * length(deltaVals)),
                         smoothWindow = rep(rep(c(1, 14, 30), each = tau), length(deltaVals)),
                         delta = rep(deltaVals, each = tau * 3),
                         H = rep(hVals, each = tau * 3),
                         Istar = NA)

for (i in 1:length(deltaVals)) {
    trueVals <- c(beta = beta, 
                  delta = deltaVals[i],
                  rateI = rateI,
                  H = 100/N)
    
    
    set.seed(123)
    threshSims$Istar[which(threshSims$delta == deltaVals[i] & 
                               threshSims$smoothWindow == 1)] <- simThresh1$run(trueVals, 1)[1:tau]
    set.seed(123)
    threshSims$Istar[which(threshSims$delta == deltaVals[i] & 
                               threshSims$smoothWindow == 14)] <- simThresh14$run(trueVals, 1)[1:tau]
    set.seed(123)
    threshSims$Istar[which(threshSims$delta == deltaVals[i] & 
                               threshSims$smoothWindow == 30)] <- simThresh30$run(trueVals, 1)[1:tau]
    
    
}


ggplot(threshSims, aes(x = time,  y = Istar, 
                       group = factor(smoothWindow), col = factor(smoothWindow))) + 
    geom_line() + 
    facet_wrap(~delta, scales = 'free')

ggplot(threshSims, aes(x = time,  y = Istar, 
                       group = factor(delta), col = factor(delta))) + 
    geom_line() + 
    facet_grid(delta~smoothWindow, scales = 'free')







