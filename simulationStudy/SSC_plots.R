################################################################################
# Output plots as pdf images for SSC PRESENTATION
# ALSO ONE PLOT FOR BAYESM
################################################################################

library(openxlsx)
library(ggplot2)
library(grid)
library(gridExtra)
library(nimble)
library(plyr)
library(knitr)
library(scales)
library(kableExtra)

# functions to calculate the alarms
source('./scripts/modelCodes.R')

theme_set(theme_bw() + 
              theme(strip.background = element_rect(fill = 'white'),
                    strip.text = element_text(size = 18),
                    axis.title = element_text(size = 17),
                    axis.text = element_text(size = 15),
                    plot.title = element_text(size = 19, h = 0.5)))

infPeriodSpec <- 'exp'



################################################################################

library(nimble)

source('./simulationStudy/scripts/modelCodes.R')

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

layoutMat <- matrix(1:6, nrow = 3, ncol = 2, byrow = T)

pdf('./Figures/exampleAlarms_SSC.pdf', width = 7, height = 9)
layout(mat = layoutMat,
       heights = c(1, 1, 1), 
       widths = c(6, 5)) 

# power
par(mar=c(5.1, 5.1, 4.1, 1.1))
plot(x, yPower1, type = 'l', col = pal[1], lwd = lineWidths, ylim = c(0, 1), lty = lineTypes[1],
     main = 'Power Alarm', xlab = 'Incidence', ylab= 'Alarm Value',
     cex.main = cexMain, cex.axis = cexAxis, cex.lab = cexLab)
lines(x, yPower2, col = pal[2], lwd = lineWidths, lty = lineTypes[2])
lines(x, yPower3, col = pal[3], lwd = lineWidths, lty = lineTypes[3])

par(mar=c(0, 0, 0, 0))
plot.new()
legend('left', bty = 'n',
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
legend('left', legend=legend_expressions_thresh,
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
legend('left', legend=legend_expressions_hill, 
       bty = 'n',
       col = pal, lty = lineTypes,
       lwd = lineWidths,
       cex = cexLeg,
       seg.len= lineLength)


dev.off()

################################################################################

### get parameters for true alarms
paramsTruth <- read.xlsx('simParamsSummary.xlsx')
paramsTruth <- paramsTruth[paramsTruth$infPeriod == infPeriodSpec,]

N <- 1e6

xAlarm <- 0:400
trueAlarmThresh14 <- thresholdAlarm(xAlarm, N = N, 
                                    delta = paramsTruth$delta[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 14)],
                                    H = paramsTruth$H[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 14)])
trueAlarmThresh30 <- thresholdAlarm(xAlarm, N = N, 
                                    delta = paramsTruth$delta[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 30)],
                                    H = paramsTruth$H[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 30)])
trueAlarmHill14 <- hillAlarm(xAlarm, 
                             nu = paramsTruth$nu[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 14)],
                             x0 = paramsTruth$x0[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 14)],
                             delta = paramsTruth$delta[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 14)])
trueAlarmHill30 <- hillAlarm(xAlarm, 
                             nu = paramsTruth$nu[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 30)],
                             x0 = paramsTruth$x0[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 30)],
                             delta = paramsTruth$delta[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 30)])
trueAlarmPower14 <- powerAlarm(xAlarm, N = N, 
                               k = paramsTruth$k[
                                   (paramsTruth$alarmGen == 'power' & 
                                        paramsTruth$smoothWindow == 14)])
trueAlarmPower30 <- powerAlarm(xAlarm, N = N, 
                               k = paramsTruth$k[
                                   (paramsTruth$alarmGen == 'power' & 
                                        paramsTruth$smoothWindow == 30)])

trueAlarms <- data.frame(xAlarm = rep(xAlarm, 6),
                         trueAlarm = c(trueAlarmThresh14,
                                       trueAlarmThresh30,
                                       trueAlarmHill14,
                                       trueAlarmHill30,
                                       trueAlarmPower14,
                                       trueAlarmPower30),
                         alarmGen = c(rep('thresh', length(xAlarm)*2),
                                      rep('hill', length(xAlarm)*2),
                                      rep('power', length(xAlarm)*2)),
                         smoothWindow = rep(rep(c(14, 30), each = length(xAlarm)), 3))


# get example data
set.seed(1)
simNumber <- round(runif(1, 0.5, 50.5))

nDays <- 100

# get true epidemic curves for each scenario
trueCurveThresh14 <- readRDS(paste0('./Data/thresh_', infPeriodSpec, '_14.rds'))[simNumber,1:nDays]
trueCurveThresh30 <- readRDS(paste0('./Data/thresh_', infPeriodSpec, '_30.rds'))[simNumber,1:nDays]
trueCurveHill14 <- readRDS(paste0('./Data/hill_', infPeriodSpec, '_14.rds'))[simNumber,1:nDays]
trueCurveHill30 <- readRDS(paste0('./Data/hill_', infPeriodSpec, '_30.rds'))[simNumber,1:nDays]
trueCurvePower14 <- readRDS(paste0('./Data/power_', infPeriodSpec, '_14.rds'))[simNumber,1:nDays]
trueCurvePower30 <- readRDS(paste0('./Data/power_', infPeriodSpec, '_30.rds'))[simNumber,1:nDays]

trueCurves <- data.frame(time = rep(1:length(trueCurveThresh14), 6),
                         truth = c(trueCurveThresh14, trueCurveThresh30,
                                   trueCurveHill14, trueCurveHill30,
                                   trueCurvePower14, trueCurvePower30),
                         alarmGen = c(rep('thresh', length(trueCurveThresh14)*2),
                                      rep('hill', length(trueCurveThresh14)*2),
                                      rep('power', length(trueCurveThresh14)*2)),
                         smoothWindow = rep(rep(c(14, 30), each = length(trueCurveThresh14)), 3))

trueAlarms$alarmGen <- factor(trueAlarms$alarmGen,
                              levels = c('power', 'thresh', 'hill'),
                              labels = c('Power', 'Threshold', 'Hill'))

trueCurves$alarmGen <- factor(trueCurves$alarmGen,
                              levels = c('power', 'thresh', 'hill'),
                              labels = c('Power', 'Threshold', 'Hill'))


p1 <- ggplot(subset(trueAlarms, smoothWindow == 30), 
             aes(x = xAlarm, y = trueAlarm)) +
    geom_line(size = 1) +
    facet_wrap(~alarmGen) +
    labs(x = '30-day average incidence', y = 'Alarm',
         title = 'True alarm functions') + 
    ylim(0, 1)

p2 <- ggplot(subset(trueCurves, smoothWindow == 30), 
             aes(x = time, y = truth)) +
    geom_line(size = 1) +
    facet_wrap(~alarmGen) +
    labs(x = 'Epidemic Time', y = 'Incidence',
         title = 'Example simulated epidemic curves') +
    geom_vline(xintercept = 50, linetype = 2) +
    annotate('text', x = 25, y = 530, label='Train', hjust = 0.5, size = 6) +
    annotate('text', x = 75, y = 530, label='Test', hjust = 0.5, size = 6) +
    ylim(0, 550)

pdf('./Figures/sim_setup_SSC.pdf', height = 7, width = 11)
grid.arrange(p1, p2, ncol = 1)
dev.off()

################################################################################



### get parameters for true alarms
paramsTruth <- read.xlsx('simParamsSummary.xlsx')
paramsTruth <- paramsTruth[paramsTruth$infPeriod == infPeriodSpec,]

N <- 1e6


xAlarm <- 0:400
trueAlarmThresh14 <- thresholdAlarm(xAlarm, N = N, 
                                    delta = paramsTruth$delta[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 14)],
                                    H = paramsTruth$H[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 14)])
trueAlarmThresh30 <- thresholdAlarm(xAlarm, N = N, 
                                    delta = paramsTruth$delta[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 30)],
                                    H = paramsTruth$H[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 30)])
trueAlarmHill14 <- hillAlarm(xAlarm, 
                             nu = paramsTruth$nu[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 14)],
                             x0 = paramsTruth$x0[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 14)],
                             delta = paramsTruth$delta[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 14)])
trueAlarmHill30 <- hillAlarm(xAlarm, 
                             nu = paramsTruth$nu[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 30)],
                             x0 = paramsTruth$x0[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 30)],
                             delta = paramsTruth$delta[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 30)])
trueAlarmPower14 <- powerAlarm(xAlarm, N = N, 
                               k = paramsTruth$k[
                                   (paramsTruth$alarmGen == 'power' & 
                                        paramsTruth$smoothWindow == 14)])
trueAlarmPower30 <- powerAlarm(xAlarm, N = N, 
                               k = paramsTruth$k[
                                   (paramsTruth$alarmGen == 'power' & 
                                        paramsTruth$smoothWindow == 30)])

trueAlarms <- data.frame(xAlarm = rep(xAlarm, 6),
                         trueAlarm = c(trueAlarmThresh14,
                                       trueAlarmThresh30,
                                       trueAlarmHill14,
                                       trueAlarmHill30,
                                       trueAlarmPower14,
                                       trueAlarmPower30),
                         alarmGen = c(rep('thresh', length(xAlarm)*2),
                                      rep('hill', length(xAlarm)*2),
                                      rep('power', length(xAlarm)*2)),
                         smoothWindow = rep(rep(c(14, 30), each = length(xAlarm)), 3))



### load posterior estimates
alarmAll <- readRDS('./resultsFinal/alarmPostAll.rds')
alarmAll <- alarmAll[alarmAll$infPeriod == infPeriodSpec,]

# format for better plotting
alarmAll$alarmFit <- factor(alarmAll$alarmFit,
                            levels = c('thresh', 'hill', 'power', 'spline', 'gp'),
                            labels = c('Threshold', 'Hill', 'Power',
                                       'Spline', 'Gaussian Process'))



theme_set(theme_bw() + 
              theme(strip.background = element_rect(fill = 'white'),
                    strip.text = element_text(size = 14),
                    axis.title = element_text(size = 14),
                    axis.text = element_text(size = 12),
                    plot.title = element_text(size = 14, h = 0.5)))




alarmAll$alarmGen <- factor(alarmAll$alarmGen,
                            levels = c('power', 'thresh', 'hill'),
                            labels = c('Power', 'Threshold', 'Hill'))

trueAlarms$alarmGen <- factor(trueAlarms$alarmGen,
                            levels = c('power', 'thresh', 'hill'),
                            labels = c('Power', 'Threshold', 'Hill'))

trueAlarms <- trueAlarms[-which(trueAlarms$alarmGen == 'Threshold' & trueAlarms$xAlarm>170),]
trueAlarms <- trueAlarms[-which(trueAlarms$alarmGen == 'Hill' & trueAlarms$xAlarm>220),]
trueAlarms <- trueAlarms[-which(trueAlarms$alarmGen == 'Power' & trueAlarms$xAlarm>300),]

pdf('./Figures/sim_alarms30_SSC.pdf', height = 6, width = 8)
ggplot() +  
    geom_line(data = subset(alarmAll, smoothWindow == 30 & 
                                alarmFit %in% c('Spline', 'Gaussian Process')), 
              aes(x = xAlarm, y = mean, group = simNumber), 
              col = adjustcolor('grey50', alpha = 0.4)) +
    geom_line(data = subset(trueAlarms, smoothWindow == 30), 
              aes(x = xAlarm, y = trueAlarm), col = 'red', size = 0.7) +
    facet_wrap(~alarmFit + alarmGen, scales = 'free_x') + 
    labs(x = '30-day average incidence', y = 'Alarm')+
    ylim(0, 1)  + 
    ggtitle('Posterior mean estimates of the alarm function') 
dev.off()


################################################################################
# BAYESM PLOT


# for a randomly selected simulation
set.seed(1)
simNumber <- round(runif(1, 0.5, 50.5))

nDays <- 100

# get true epidemic curves for each scenario
trueCurveThresh14 <- readRDS(paste0('./Data/thresh_', infPeriodSpec, '_14.rds'))[simNumber,1:nDays]
trueCurveThresh30 <- readRDS(paste0('./Data/thresh_', infPeriodSpec, '_30.rds'))[simNumber,1:nDays]
trueCurveHill14 <- readRDS(paste0('./Data/hill_', infPeriodSpec, '_14.rds'))[simNumber,1:nDays]
trueCurveHill30 <- readRDS(paste0('./Data/hill_', infPeriodSpec, '_30.rds'))[simNumber,1:nDays]
trueCurvePower14 <- readRDS(paste0('./Data/power_', infPeriodSpec, '_14.rds'))[simNumber,1:nDays]
trueCurvePower30 <- readRDS(paste0('./Data/power_', infPeriodSpec, '_30.rds'))[simNumber,1:nDays]

trueCurves <- data.frame(time = rep(1:length(trueCurveThresh14), 6),
                         truth = c(trueCurveThresh14, trueCurveThresh30,
                                   trueCurveHill14, trueCurveHill30,
                                   trueCurvePower14, trueCurvePower30),
                         alarmGen = c(rep('thresh', length(trueCurveThresh14)*2),
                                      rep('hill', length(trueCurveThresh14)*2),
                                      rep('power', length(trueCurveThresh14)*2)),
                         smoothWindow = rep(rep(c(14, 30), each = length(trueCurveThresh14)), 3))

# merge with posterior predictions
postPredAll <- readRDS('./resultsFinal/postPredAll.rds')
postPredAll <- postPredAll[postPredAll$simNumber == simNumber & 
                               postPredAll$infPeriod == infPeriodSpec,]

# format for better plotting
postPredAll$alarmFit <- factor(postPredAll$alarmFit,
                               levels = c('basic', 'power', 'thresh', 'hill', 'spline', 'gp'),
                               labels = c('No Behavioral Change', 'Power', 'Threshold', 'Hill',
                                          'Spline', 'Behavioral Change'))

myCol <- 'blue'

theme_set(theme_bw() + 
              theme(strip.background = element_rect(fill = 'white'),
                    strip.text = element_text(size = 28),
                    axis.title = element_text(size = 26),
                    axis.text = element_text(size = 24),
                    plot.title = element_text(size = 36, h = 0.5),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),))

pdf('../../nycBCMAnalysis/figures/simPostPred.pdf', height = 4.5, width = 12)
ggplot() +
    geom_line(data = subset(trueCurves, alarmGen == 'hill' & smoothWindow == 30),
              aes(x = time, y = truth)) + 
    geom_line(data = subset(postPredAll, alarmGen == 'hill' & smoothWindow == 30 
                            & alarmFit %in% c('No Behavioral Change', 'Behavioral Change')),
              aes(x = time, y = mean), col = myCol, size = 1) + 
    geom_ribbon(data=subset(postPredAll, alarmGen == 'hill' & smoothWindow == 30 
                            & alarmFit %in% c('No Behavioral Change', 'Behavioral Change')),
                aes(x = time, ymin=lower, ymax=upper), alpha=0.3, fill = myCol) +
    facet_wrap(~alarmFit, nrow = 1) +
    labs(x = 'Epidemic time', y = 'Incidence') 
dev.off()

