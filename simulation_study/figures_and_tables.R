################################################################################
# Figures in paper
################################################################################

# load libraries
library(nimble)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(RColorBrewer)
library(ggpubr)
library(openxlsx)
library(plyr)
library(kableExtra)

source('./scripts/modelCodes.R')


# check convergence

grAll <- readRDS('./results/grAll.rds')

# which didn't converge
notConverge <- grAll[which(grAll$gr > 1.1),  ]
notConvergeModels <-  notConverge[!duplicated(notConverge[,-which(colnames(notConverge) %in% 
                                                                      c('gr', 'grUpper', 'param'))]),
                                  c('alarmGen', 'alarmFit', 'smoothWindow', 'simNumber')]
notConvergeModels$noConverge <- 1


table(notConvergeModels$alarmGen, 
      notConvergeModels$alarmFit,
      notConvergeModels$smoothWindow)


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

incRange <- 1:300

dataNodes <- c(paste0('Istar[', 1:tau, ']'),
               paste0('Rstar[', 1:tau, ']'))


### Define simulator functions for each alarm/smoothing combination

constants1 <- list(N = N, tau = tau, I0 = I0, bw = 1, aa = 1, bb = 1)
constants14 <- list(N = N, tau = tau, I0 = I0, bw = 14, aa = 1, bb = 1)
constants30 <- list(N = N, tau = tau, I0 = I0, bw = 30, aa = 1, bb = 1)

## Power alarm

# 1 day average
powerModel1 <- nimbleModel(code = SIR_power_sim,
                           constants = constants1)
simPower1 <- simulator(powerModel1, dataNodes)

# 14 day average
powerModel14 <- nimbleModel(code = SIR_power_sim,
                            constants = constants14)
simPower14 <- simulator(powerModel14, dataNodes)

# 30 day average
powerModel30 <- nimbleModel(code = SIR_power_sim,
                            constants = constants30)
simPower30 <- simulator(powerModel30, dataNodes)


## Threshold alarm

# 1 day average
threshModel1 <- nimbleModel(code = SIR_thresh_sim,
                            constants = constants1)
simThresh1 <- simulator(threshModel1, dataNodes)

# 14 day average
threshModel14 <- nimbleModel(code = SIR_thresh_sim,
                             constants = constants14)
simThresh14 <- simulator(threshModel14, dataNodes)

# 30 day average
threshModel30 <- nimbleModel(code = SIR_thresh_sim,
                             constants = constants30)
simThresh30 <- simulator(threshModel30, dataNodes)

## Hill alarm

# 1 day average
hillModel1 <- nimbleModel(code = SIR_hill_sim,
                          constants = constants1)
simHill1 <- simulator(hillModel1, dataNodes)

# 14 day average
hillModel14 <- nimbleModel(code = SIR_hill_sim,
                           constants = constants14)
simHill14 <- simulator(hillModel14, dataNodes)

# 30 day average
hillModel30 <- nimbleModel(code = SIR_hill_sim,
                           constants = constants30)
simHill30 <- simulator(hillModel30, dataNodes)


### Threshold/hill flattening curve, power not

# using smoothing of 30 days
# compared to no BC model

### no BC Model
trueVals <- c(beta = beta, 
              delta = 0,
              rateI = rateI,
              H = 1)

set.seed(123)
noBC <- simThresh1$run(trueVals, 1)[1:tau]

### power
trueVals <- c(beta = beta, 
              k = 0.05,
              rateI = rateI)

set.seed(123)
power1 <- simPower30$run(trueVals, 1)[1:tau]
powerAlarm1 <- powerAlarm(incRange, N = N, k = trueVals['k'])

trueVals['k'] <- 0.01

set.seed(123)
power2 <- simPower30$run(trueVals, 1)[1:tau]
powerAlarm2 <- powerAlarm(incRange, N = N, k = trueVals['k'])

trueVals['k'] <- 0.003

set.seed(123)
power3 <- simPower30$run(trueVals, 1)[1:tau]
powerAlarm3 <- powerAlarm(incRange, N = N, k = trueVals['k'])


### threshold
# delta = 0.2
trueVals <- c(beta = beta, 
              delta = 0.2,
              rateI = rateI,
              H = 100/N)

set.seed(123)
thresh1 <- simThresh30$run(trueVals, 1)[1:tau]
threshAlarm1 <- thresholdAlarm(incRange, N = N, 
                               delta = trueVals['delta'], H = trueVals['H'])

# delta = 0.4
trueVals['delta'] <- 0.4

set.seed(123)
thresh2 <- simThresh30$run(trueVals, 1)[1:tau]
threshAlarm2 <- thresholdAlarm(incRange, N = N, 
                               delta = trueVals['delta'], H = trueVals['H'])

# delta = 0.6
trueVals['delta'] <- 0.6

set.seed(123)
thresh3 <- simThresh30$run(trueVals, 1)[1:tau]
threshAlarm3 <- thresholdAlarm(incRange, N = N, 
                               delta = trueVals['delta'], H = trueVals['H'])

### Hill
# delta = 0.2
trueVals <- c(beta = beta, 
              delta = 0.2,
              nu = 5,
              rateI = rateI,
              x0 = 100)


set.seed(123)
hill1 <- simHill30$run(trueVals, 1)[1:tau]
hillAlarm1 <- hillAlarm(incRange, 
                        nu = trueVals['nu'], x0 = trueVals['x0'],
                        delta = trueVals['delta'])

# delta = 0.4
trueVals['delta'] <- 0.4

set.seed(123)
hill2 <- simHill30$run(trueVals, 1)[1:tau]
hillAlarm2 <- hillAlarm(incRange, 
                        nu = trueVals['nu'], x0 = trueVals['x0'],
                        delta = trueVals['delta'])

# delta = 0.6
trueVals['delta'] <- 0.6

set.seed(123)
hill3 <- simHill30$run(trueVals, 1)[1:tau]
hillAlarm3 <- hillAlarm(incRange, 
                        nu = trueVals['nu'], x0 = trueVals['x0'],
                        delta = trueVals['delta'])

toPlot_epi_A <- data.frame(alarm = c(rep('Power', tau * 4),
                                     rep('Threshold', tau * 4),
                                     rep('Hill', tau * 4)),
                           time = rep(1:tau, 12),
                           param = c(rep(c(0, 0.2, 0.4, 0.6), each = tau),
                                     rep(c(0, 0.05, 0.01, 0.003), each = tau),
                                     rep(c(0, 0.2, 0.4, 0.6), each = tau)),
                           Istar = c(noBC, power1, power2, power3,
                                     noBC, thresh1, thresh2, thresh3,
                                     noBC, hill1, hill2, hill3))

toPlot_alarm_A <- data.frame(alarm = c(rep('Power', length(incRange) * 3),
                                       rep('Threshold', length(incRange) * 3),
                                       rep('Hill', length(incRange) * 3)),
                             incRange = rep(incRange, 9),
                             param = c(rep(c(0.2, 0.4, 0.6), each = length(incRange)),
                                       rep(c(0.05, 0.01, 0.003), each = length(incRange)),
                                       rep(c(0.2, 0.4, 0.6), each = length(incRange))),
                             alarmVal = c(powerAlarm1, powerAlarm2, powerAlarm3,
                                       threshAlarm1, threshAlarm2, threshAlarm3,
                                       hillAlarm1, hillAlarm2, hillAlarm3))

toPlot_epi_A$paramLab <- factor(toPlot_epi_A$param, 
                            levels = c(0, 0.2, 0.4, 0.6, 0.05, 0.01, 0.003))
toPlot_alarm_A$paramLab <- factor(toPlot_alarm_A$param, 
                            levels = c(0.2, 0.4, 0.6, 0.05, 0.01, 0.003))

plotLabsDelta <- c('No behavioral change',
                   ~delta~'= 0.2',
                   ~delta~'= 0.4',
                   ~delta~'= 0.6')

plotLabsK <- c('No behavioral change',
               'k = 0.05', 'k = 0.01', 'k = 0.003')

pal <- c('black', 'royalblue3', 'cyan3', 'darkgoldenrod2')

legendLocation1 <- c(0.2, 0.80)
legendLocation2 <- c(0.65, 0.77)

myTheme1 <- theme(legend.position = legendLocation1,
                 legend.title= element_blank(),
                 plot.title = element_text(h = 0.5, size = 15),
                 axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 legend.text = element_text(size = 12))

myTheme2 <- theme(legend.position = legendLocation2,
                 legend.title= element_blank(),
                 plot.title = element_text(h = 0.5, size = 14),
                 axis.title = element_text(size = 12),
                 axis.text = element_text(size = 10),
                 legend.text = element_text(size = 12))

# alarms
p1 <- ggplot(subset(toPlot_alarm_A, alarm == 'Power'), 
             aes(x = incRange, y = alarmVal, color = paramLab)) + 
    geom_line(size = 1, alpha = 0.8) +
    theme_bw() + 
    labs(x = '30-day Incidence', y = 'Alarm',
         title = 'Power Alarm') + 
    myTheme1 +
    scale_color_manual(values = pal[-1],
                       labels = plotLabsK[-1]) +
    ylim(0, 1)

p2 <- ggplot(subset(toPlot_alarm_A, alarm == 'Threshold'), 
             aes(x = incRange, y = alarmVal, color = paramLab)) + 
    geom_line(size = 1, alpha = 0.8) +
    theme_bw() + 
    labs(x = '30-day Incidence', y = 'Alarm',
         title = 'Threshold Alarm') + 
    myTheme1 + 
    scale_color_manual(values = pal[-1],
                       labels = plotLabsDelta[-1]) +
    ylim(0, 1)

p3 <- ggplot(subset(toPlot_alarm_A, alarm == 'Hill'), 
             aes(x = incRange, y = alarmVal, color = paramLab)) + 
    geom_line(size = 1, alpha = 0.8) +
    theme_bw() + 
    labs(x = '30-day Incidence', y = 'Alarm',
         title = 'Hill Alarm') + 
    myTheme1 + 
    scale_color_manual(values = pal[-1],
                       labels = plotLabsDelta[-1]) +
    ylim(0, 1)

# epidemics
p4 <- ggplot(subset(toPlot_epi_A, alarm == 'Power'), 
             aes(x = time, y = Istar, color = paramLab)) + 
    geom_line(size = 1, alpha = 0.8) +
    theme_bw() + 
    labs(x = 'Epidemic Time', y = 'Incidence',
         title = 'Power Alarm') + 
    myTheme2 +
    scale_color_manual(values = pal,
                       labels = plotLabsK)

p5 <- ggplot(subset(toPlot_epi_A, alarm == 'Threshold'), 
             aes(x = time, y = Istar, color = paramLab)) + 
    geom_line(size = 1, alpha = 0.8) +
    theme_bw() + 
    labs(x = 'Epidemic Time', y = 'Incidence',
         title = 'Threshold Alarm') + 
    myTheme2 + 
    scale_color_manual(values = pal,
                       labels = plotLabsDelta)

p6 <- ggplot(subset(toPlot_epi_A, alarm == 'Hill'), 
             aes(x = time, y = Istar, color = paramLab)) + 
    geom_line(size = 1, alpha = 0.8) +
    theme_bw() + 
    labs(x = 'Epidemic Time', y = 'Incidence',
         title = 'Hill Alarm') + 
    myTheme2 + 
    scale_color_manual(values = pal,
                       labels = plotLabsDelta)

# pdf('./figures/fig2_ex_epidemics.pdf', height = 8, width = 13)
# grid.arrange(p1, p2, p3,
#              p4, p5, p6, nrow = 2)
# dev.off()




# threshold vs hill smoothness

# power 1, 14, 30
trueVals <- c(beta = beta, 
              k = 0.0005,
              rateI = rateI)

set.seed(123)
power1 <- simPower1$run(trueVals, 1)[1:tau]
set.seed(123)
power2 <- simPower14$run(trueVals, 1)[1:tau]
set.seed(123)
power3 <- simPower30$run(trueVals, 1)[1:tau]

# threshold 1, 14, 30
trueVals <- c(beta = beta, 
              delta = 0.8,
              rateI = rateI,
              H = 350/N)

set.seed(123)
thresh1 <- simThresh1$run(trueVals, 1)[1:tau]
set.seed(123)
thresh2 <- simThresh14$run(trueVals, 1)[1:tau]
set.seed(123)
thresh3 <- simThresh30$run(trueVals, 1)[1:tau]

# hill 1, 14, 30
trueVals <- c(beta = beta, 
              delta = 0.85,
              nu = 2,
              rateI = rateI,
              x0 = 450)

set.seed(123)
hill1 <- simHill1$run(trueVals, 1)[1:tau]
set.seed(123)
hill2 <- simHill14$run(trueVals, 1)[1:tau]
set.seed(123)
hill3 <- simHill30$run(trueVals, 1)[1:tau]


toPlot_B <- data.frame(alarm = c(rep('Power', tau * 3),
                                 rep('Threshold', tau * 3),
                                 rep('Hill', tau * 3)),
                       time = rep(1:tau, 9),
                       smoothWindow = rep(rep(c(1, 14, 30), each = tau), 3),
                       Istar = c(power1, power2, power3,
                                 thresh1, thresh2, thresh3,
                                 hill1, hill2, hill3))

toPlot_B$alarm <- factor(toPlot_B$alarm,
                         levels = c('Power', 'Threshold', 'Hill'))
toPlot_B$smoothWindow <- paste0(toPlot_B$smoothWindow,
                                '-day average')

### create plots

p7 <- ggplot(subset(toPlot_B, time <=80), 
       aes(x = time, y = Istar, color = alarm)) + 
    geom_line(size = 1, alpha = 0.8) +
    facet_wrap(~smoothWindow) + 
    theme_bw() +
    scale_color_manual(values = brewer.pal(3, 'Dark2')) +
    theme(strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) +
    labs(color = 'Alarm',
         x = 'Epidemic Time', y = 'Incidence') 

pdf('./figures/fig2_ex_epidemics.pdf', height = 7, width = 12)
ggarrange(ggarrange(p4, p5, p6, NULL, nrow = 1,
                    widths = c(1, 1, 1, 0.1)), 
          p7, nrow = 2,
          labels = c('(A)', '(B)'))
dev.off()





################################################################################
# Figure 3 - Simulation set up 30 day smoothing alarms and example epidemics
################################################################################


theme_set(theme_bw() + 
              theme(strip.background = element_rect(fill = 'white'),
                    strip.text = element_text(size = 16),
                    axis.title = element_text(size = 16),
                    axis.text = element_text(size = 14),
                    plot.title = element_text(size = 17, h = 0.5)))


### get parameters for true alarms
paramsTruth <- read.xlsx('simParamsSummary.xlsx')

N <- 1e6

xAlarm <- 0:400

# True alarm functions used for data generation
trueAlarmPower30 <- powerAlarm(xAlarm, N = N, 
                               k = paramsTruth$k[
                                   (paramsTruth$alarmGen == 'power' & 
                                        paramsTruth$smoothWindow == 30)])

trueAlarmThresh30 <- thresholdAlarm(xAlarm, N = N, 
                                    delta = paramsTruth$delta[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 30)],
                                    H = paramsTruth$H[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 30)])

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

trueAlarms <- data.frame(xAlarm = rep(xAlarm, 3),
                         trueAlarm = c(trueAlarmPower30,
                                       trueAlarmThresh30,
                                       trueAlarmHill30),
                         alarmGen = c(rep('power', length(xAlarm)),
                                      rep('thresh', length(xAlarm)),
                                      rep('hill', length(xAlarm))))


# get example data
set.seed(1)
simNumber <- round(runif(1, 0.5, 50.5))

nDays <- 100

# get true epidemic curves for each scenario
trueCurvePower30 <- readRDS(paste0('./Data/power_30.rds'))[simNumber,1:nDays]
trueCurveThresh30 <- readRDS(paste0('./Data/thresh_30.rds'))[simNumber,1:nDays]
trueCurveHill30 <- readRDS(paste0('./Data/hill_30.rds'))[simNumber,1:nDays]

trueCurves <- data.frame(time = rep(1:nDays, 3),
                         truth = c(trueCurvePower30, 
                                   trueCurveThresh30,
                                   trueCurveHill30),
                         alarmGen = c(rep('power', nDays),
                                      rep('thresh', nDays),
                                      rep('hill', nDays)))

trueAlarms$alarmGen <- factor(trueAlarms$alarmGen,
                              levels = c('power', 'thresh', 'hill'),
                              labels = c('Power', 'Threshold', 'Hill'))

trueCurves$alarmGen <- factor(trueCurves$alarmGen,
                              levels = c('power', 'thresh', 'hill'),
                              labels = c('Power', 'Threshold', 'Hill'))


p1 <- ggplot(trueAlarms, aes(x = xAlarm, y = trueAlarm)) +
    geom_line(size = 1) +
    facet_wrap(~alarmGen) +
    labs(x = '30-day average incidence', y = 'Alarm',
         title = 'True alarm functions') + 
    ylim(0, 1)

p2 <- ggplot(trueCurves, aes(x = time, y = truth)) +
    geom_line(size = 1) +
    facet_wrap(~alarmGen) +
    labs(x = 'Epidemic Time', y = 'Incidence',
         title = 'Example simulated epidemic curves') +
    geom_vline(xintercept = 50, linetype = 2) +
    annotate('text', x = 25, y = 530, label='Train', hjust = 0.5, size = 6) +
    annotate('text', x = 75, y = 530, label='Test', hjust = 0.5, size = 6) +
    ylim(0, 550)

pdf('./figures/fig3_sim_setup.pdf', height = 7, width = 9)
grid.arrange(p1, p2, ncol = 1)
dev.off()




################################################################################
# Figure 4 - Posterior estimation of the alarm functions for 30-day smoothing
# Supplemental Figure 5 - for 14-day smoothing
################################################################################

theme_set(theme_bw() + 
              theme(strip.background = element_rect(fill = 'white'),
                    strip.text = element_text(size = 14),
                    axis.title = element_text(size = 14),
                    axis.text = element_text(size = 12),
                    plot.title = element_text(size = 14, h = 0.5)))

### get parameters for true alarms
paramsTruth <- read.xlsx('simParamsSummary.xlsx')
paramsTruth <- paramsTruth[paramsTruth$infPeriod == infPeriodSpec,]

N <- 1e6

xAlarm <- 0:500
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
alarmAll <- readRDS('./results/alarmPostAll.rds')

# format for better plotting
alarmAll$alarmFit <- factor(alarmAll$alarmFit,
                            levels = c('thresh', 'hill', 'power', 'spline', 'gp'),
                            labels = c('Threshold', 'Hill', 'Power',
                                       'Spline', 'Gaussian Process'))


### smoothing window 14 days
p1 <- ggplot() +  
    geom_line(data = subset(alarmAll, smoothWindow == 14 & alarmGen == 'power'), 
              aes(x = xAlarm, y = mean, group = simNumber), 
              col = adjustcolor('grey50', alpha = 0.4)) +
    geom_line(data = subset(trueAlarms, smoothWindow == 14 & alarmGen == 'power'), 
              aes(x = xAlarm, y = trueAlarm), col = 'red', size = 0.7) +
    facet_wrap(~alarmFit) + 
    labs(x = '14-day average incidence', y = 'Alarm')+
    ylim(0, 1) + 
    ggtitle('Data generation: power alarm')+
    xlim(0, 370)

p2 <- ggplot() +  
    geom_line(data = subset(alarmAll, smoothWindow == 14 &  alarmGen == 'thresh'), 
              aes(x = xAlarm, y = mean, group = simNumber), 
              col = adjustcolor('grey50', alpha = 0.4)) +
    geom_line(data = subset(trueAlarms, smoothWindow == 14 & alarmGen == 'thresh'), 
              aes(x = xAlarm, y = trueAlarm), col = 'red', size = 0.7) +
    facet_wrap(~alarmFit) + 
    labs(x = '
         14-day average incidence', y = 'Alarm') +
    ylim(0, 1)  + 
    ggtitle('Data generation: threshold alarm')+
    xlim(0, 250)

p3 <- ggplot() +  
    geom_line(data = subset(alarmAll, smoothWindow == 14 & alarmGen == 'hill'), 
              aes(x = xAlarm, y = mean, group = simNumber), 
              col = adjustcolor('grey50', alpha = 0.4)) +
    geom_line(data = subset(trueAlarms, smoothWindow == 14 & alarmGen == 'hill'), 
              aes(x = xAlarm, y = trueAlarm), col = 'red', size = 0.7) +
    facet_wrap(~alarmFit) + 
    labs(x = '14-day average incidence', y = 'Alarm')+
    ylim(0, 1)  + 
    ggtitle('Data generation: Hill alarm') +
    xlim(0, 300)



pdf('./figures/supp_fig5_sim_alarms14.pdf', height = 8, width = 7)
grid.arrange(p1, p2, p3, nrow = 3,
             top=textGrob("         Posterior mean of alarm functions", gp = gpar(fontsize = 16, font = 1)))

dev.off()


### smoothing window 30 days
p1 <- ggplot() +  
    geom_line(data = subset(alarmAll, smoothWindow == 30 & alarmGen == 'power'), 
              aes(x = xAlarm, y = mean, group = simNumber), 
              col = adjustcolor('grey50', alpha = 0.4)) +
    geom_line(data = subset(trueAlarms, smoothWindow == 30 & alarmGen == 'power'), 
              aes(x = xAlarm, y = trueAlarm), col = 'red', size = 0.7) +
    facet_wrap(~alarmFit) + 
    labs(x = '30-day average incidence', y = 'Alarm')+
    ylim(0, 1) + 
    ggtitle('Data generation: power alarm')+
    xlim(0, 300)

p2 <- ggplot() +  
    geom_line(data = subset(alarmAll, smoothWindow == 30 &  alarmGen == 'thresh'), 
              aes(x = xAlarm, y = mean, group = simNumber), 
              col = adjustcolor('grey50', alpha = 0.4)) +
    geom_line(data = subset(trueAlarms, smoothWindow == 30 & alarmGen == 'thresh'), 
              aes(x = xAlarm, y = trueAlarm), col = 'red', size = 0.7) +
    facet_wrap(~alarmFit) + 
    labs(x = '
         30-day average incidence', y = 'Alarm') +
    ylim(0, 1)  + 
    ggtitle('Data generation: threshold alarm')+
    xlim(0, 170)

p3 <- ggplot() +  
    geom_line(data = subset(alarmAll, smoothWindow == 30 & alarmGen == 'hill'), 
              aes(x = xAlarm, y = mean, group = simNumber), 
              col = adjustcolor('grey50', alpha = 0.4)) +
    geom_line(data = subset(trueAlarms, smoothWindow == 30 & alarmGen == 'hill'), 
              aes(x = xAlarm, y = trueAlarm), col = 'red', size = 0.7) +
    facet_wrap(~alarmFit) + 
    labs(x = '30-day average incidence', y = 'Alarm')+
    ylim(0, 1)  + 
    ggtitle('Data generation: Hill alarm') +
    xlim(0, 220)



pdf('./figures/fig4_sim_alarms30.pdf', height = 8, width = 7)
grid.arrange(p1, p2, p3, nrow = 3,
             top=textGrob("         Posterior mean of alarm functions", gp = gpar(fontsize = 16, font = 1)))

dev.off()


################################################################################
# Figure 5 - Posterior prediction for 30-day smoothing
# Supplemental Figure 6 - for 14-day smoothing
################################################################################


theme_set(theme_bw() + 
              theme(strip.background = element_rect(fill = 'white'),
                    strip.text = element_text(size = 14),
                    axis.title = element_text(size = 14),
                    axis.text = element_text(size = 12),
                    plot.title = element_text(size = 14, h = 0.5)))

# for a randomly selected simulation
set.seed(1)
simNumber <- round(runif(1, 0.5, 50.5))

nDays <- 100

# get true epidemic curves for each scenario
trueCurveThresh14 <- readRDS(paste0('./Data/thresh_14.rds'))[simNumber,1:nDays]
trueCurveThresh30 <- readRDS(paste0('./Data/thresh_30.rds'))[simNumber,1:nDays]
trueCurveHill14 <- readRDS(paste0('./Data/hill_14.rds'))[simNumber,1:nDays]
trueCurveHill30 <- readRDS(paste0('./Data/hill_30.rds'))[simNumber,1:nDays]
trueCurvePower14 <- readRDS(paste0('./Data/power_14.rds'))[simNumber,1:nDays]
trueCurvePower30 <- readRDS(paste0('./Data/power_30.rds'))[simNumber,1:nDays]

trueCurves <- data.frame(time = rep(1:length(trueCurveThresh14), 6),
                         truth = c(trueCurveThresh14, trueCurveThresh30,
                                   trueCurveHill14, trueCurveHill30,
                                   trueCurvePower14, trueCurvePower30),
                         alarmGen = c(rep('thresh', length(trueCurveThresh14)*2),
                                      rep('hill', length(trueCurveThresh14)*2),
                                      rep('power', length(trueCurveThresh14)*2)),
                         smoothWindow = rep(rep(c(14, 30), each = length(trueCurveThresh14)), 3))

# merge with posterior predictions
postPredAll <- readRDS('./results/postPredAll.rds')
postPredAll <- postPredAll[postPredAll$simNumber == simNumber,]

# format for better plotting
postPredAll$alarmFit <- factor(postPredAll$alarmFit,
                               levels = c('basic', 'power', 'thresh', 'hill', 'spline', 'gp'),
                               labels = c('No Behavioral Change', 'Power', 'Threshold', 'Hill',
                                          'Spline', 'Gaussian Process'))

myCol <- 'blue'

### 14-day Smoothing Window
p1 <- ggplot() +
    geom_line(data = subset(trueCurves, alarmGen == 'power' & smoothWindow == 14),
              aes(x = time, y = truth)) + 
    geom_line(data = subset(postPredAll, alarmGen == 'power' & smoothWindow == 14),
              aes(x = time, y = mean), col = myCol, size = 0.5) + 
    geom_ribbon(data=subset(postPredAll, alarmGen == 'power' & smoothWindow == 14),
                aes(x = time, ymin = lower, ymax = upper), alpha = 0.3, fill = myCol) +
    facet_wrap(~alarmFit, nrow = 1) +
    labs(x = 'Epidemic time', y = 'Incidence') + 
    ggtitle('Data generation: power alarm')

p2 <- ggplot() +
    geom_line(data = subset(trueCurves, alarmGen == 'thresh' & smoothWindow == 14),
              aes(x = time, y = truth)) + 
    geom_line(data = subset(postPredAll, alarmGen == 'thresh' & smoothWindow == 14),
              aes(x = time, y = mean), col = myCol, size = 0.5) + 
    geom_ribbon(data=subset(postPredAll, alarmGen == 'thresh' & smoothWindow == 14),
                aes(x = time, ymin=lower, ymax=upper), alpha=0.3, fill = myCol) +
    facet_wrap(~alarmFit, nrow = 1) +
    labs(x = 'Epidemic time', y = 'Incidence') + 
    ggtitle('Data generation: threshold alarm')

p3 <- ggplot() +
    geom_line(data = subset(trueCurves, alarmGen == 'hill' & smoothWindow == 14),
              aes(x = time, y = truth)) + 
    geom_line(data = subset(postPredAll, alarmGen == 'hill' & smoothWindow == 14),
              aes(x = time, y = mean), col = myCol, size = 0.5) + 
    geom_ribbon(data=subset(postPredAll, alarmGen == 'hill' & smoothWindow == 14),
                aes(x = time, ymin=lower, ymax=upper), alpha=0.3, fill = myCol) +
    facet_wrap(~alarmFit, nrow = 1) +
    labs(x = 'Epidemic time', y = 'Incidence')  + 
    ggtitle('Data generation: Hill alarm')


pdf('./figures/supp_fig6_sim_postPred14.pdf', height = 8, width = 9)
grid.arrange(p1, p2, p3, nrow = 3,
             top=textGrob("        Posterior predictive forecasting", gp = gpar(fontsize = 18, font = 1)))

dev.off()



### 30-day Smoothing Window
p1 <- ggplot() +
    geom_line(data = subset(trueCurves, alarmGen == 'power' & smoothWindow == 30),
              aes(x = time, y = truth)) + 
    geom_line(data = subset(postPredAll, alarmGen == 'power' & smoothWindow == 30),
              aes(x = time, y = mean), col = myCol, size = 0.5) + 
    geom_ribbon(data=subset(postPredAll, alarmGen == 'power' & smoothWindow == 30),
                aes(x = time, ymin = lower, ymax = upper), alpha = 0.3, fill = myCol) +
    facet_wrap(~alarmFit, nrow = 1) +
    labs(x = 'Epidemic time', y = 'Incidence') + 
    ggtitle('Data generation: power alarm')

p2 <- ggplot() +
    geom_line(data = subset(trueCurves, alarmGen == 'thresh' & smoothWindow == 30),
              aes(x = time, y = truth)) + 
    geom_line(data = subset(postPredAll, alarmGen == 'thresh' & smoothWindow == 30),
              aes(x = time, y = mean), col = myCol, size = 0.5) + 
    geom_ribbon(data=subset(postPredAll, alarmGen == 'thresh' & smoothWindow == 30),
                aes(x = time, ymin=lower, ymax=upper), alpha=0.3, fill = myCol) +
    facet_wrap(~alarmFit, nrow = 1) +
    labs(x = 'Epidemic time', y = 'Incidence') + 
    ggtitle('Data generation: threshold alarm') +
    ylim(0, 600)

p3 <- ggplot() +
    geom_line(data = subset(trueCurves, alarmGen == 'hill' & smoothWindow == 30),
              aes(x = time, y = truth)) + 
    geom_line(data = subset(postPredAll, alarmGen == 'hill' & smoothWindow == 30),
              aes(x = time, y = mean), col = myCol, size = 0.5) + 
    geom_ribbon(data=subset(postPredAll, alarmGen == 'hill' & smoothWindow == 30),
                aes(x = time, ymin=lower, ymax=upper), alpha=0.3, fill = myCol) +
    facet_wrap(~alarmFit, nrow = 1) +
    labs(x = 'Epidemic time', y = 'Incidence')  + 
    ggtitle('Data generation: Hill alarm')


pdf('./figures/fig5_sim_postPred30.pdf', height = 8, width = 9)
grid.arrange(p1, p2, p3, nrow = 3,
             top=textGrob("        Posterior predictive forecasting", gp = gpar(fontsize = 18, font = 1)))

dev.off()







################################################################################
# Table 1 - Summaries of WAIC
# Supplemental Table 3 - for 14-day smoothing
################################################################################

### Summary table with mean (SD) and % selected WAIC

waicAll <- readRDS('./results/waicAll.rds')

waicMin <- ddply(waicAll, .(alarmGen, smoothWindow, simNumber), summarize,
                 minWAIC = min(waic))

waicAll <- merge(waicAll, waicMin, 
                 by = c('alarmGen', 'smoothWindow', 'simNumber'),
                 all.x = T)

# get mean and sd across simulations
waicSummary <- ddply(waicAll, .(alarmGen, alarmFit, smoothWindow), summarize,
                     meanWAIC = mean(waic),
                     sdWAIC = sd(waic))


# get % selected
waicSelected <- ddply(waicAll, .(alarmGen, smoothWindow, simNumber), summarize,
                      selected = alarmFit[which(waic == minWAIC)])
waicSelected$selected <- factor(waicSelected$selected,
                                levels = c('power', 'thresh', 'hill',
                                           'spline', 'gp', 'betatSpline', 'basic'))

waicTab14 <- with(subset(waicSelected, smoothWindow == 14), 
                  table(selected, alarmGen, exclude = NULL))
waicTab14 <- as.data.frame(prop.table(waicTab14, 2))
waicTab14$smoothWindow <- 14

waicTab30 <- with(subset(waicSelected, smoothWindow == 30), 
                  table(selected, alarmGen, exclude = NULL))
waicTab30 <- as.data.frame(prop.table(waicTab30, 2))
waicTab30$smoothWindow <- 30

waicTab <- rbind.data.frame(waicTab14, waicTab30)
colnames(waicTab)[1] <- 'alarmFit'

# merge with other summaries
waicFinal <- merge(waicSummary, waicTab, by = c('alarmGen', 'alarmFit', 'smoothWindow'))

### 14-day smoothing
waicFinal14 <- subset(waicFinal, smoothWindow == 14)
waicFinal14$WAIC <- paste0(sprintf("%.2f", round(waicFinal14$meanWAIC, 2)), 
                           ' (',
                           sprintf("%.2f", round(waicFinal14$sdWAIC, 2)),
                           ')')

waicFinal14$Perc <- paste0(round(waicFinal14$Freq * 100, 2), '\\%')

waicFinal14 <- waicFinal14[-which(colnames(waicFinal14) %in% 
                                      c('smoothWindow', 'sdWAIC',
                                        'Freq'))]
waicFinal14$alarmGen <- factor(waicFinal14$alarmGen,
                               levels = c('power', 'thresh', 'hill'),
                               labels = c('Power',
                                          'Threshold', 
                                          'Hill'))


waicFinal14$alarmFit <- factor(waicFinal14$alarmFit,
                               levels = c('power', 'thresh', 'hill',
                                          'spline', 'gp', 'betatSpline', 'basic'),
                               labels = c('Power', 'Threshold', 'Hill',
                                          'Spline', 'Gaussian Process', 
                                          '$\\beta_t$', 'No Behavioral Change'))

waicFinal14 <- waicFinal14[order(waicFinal14$alarmGen, 
                                 waicFinal14$meanWAIC),]
waicFinal14 <- waicFinal14[,-which(colnames(waicFinal14) == 'meanWAIC')]


kable(waicFinal14, row.names = F, format = 'latex', align = 'llcc', 
      booktabs = T, escape = F, 
      col.names = linebreak(c('\\textbf{Data generation}',
                              '\\textbf{Model fitted}', 
                              '\\textbf{WAIC}\n\\textbf{Mean (SD)}',
                              '\\textbf{\\% selected}'), align = 'c')) %>% 
    collapse_rows(columns = 1, latex_hline = 'major') 

### 30-day smoothing
waicFinal30 <- subset(waicFinal, smoothWindow == 30)
waicFinal30$WAIC <- paste0(sprintf("%.2f", round(waicFinal30$meanWAIC, 2)), 
                           ' (',
                           sprintf("%.2f", round(waicFinal30$sdWAIC, 2)),
                           ')')

waicFinal30$Perc <- paste0(round(waicFinal30$Freq * 100, 2), '\\%')

waicFinal30 <- waicFinal30[-which(colnames(waicFinal30) %in% 
                                      c('smoothWindow', 'sdWAIC',
                                        'Freq'))]
waicFinal30$alarmGen <- factor(waicFinal30$alarmGen,
                               levels = c('power', 'thresh', 'hill'),
                               labels = c('Power',
                                          'Threshold', 
                                          'Hill'))


waicFinal30$alarmFit <- factor(waicFinal30$alarmFit,
                               levels = c('power', 'thresh', 'hill',
                                          'spline', 'gp', 'betatSpline', 'basic'),
                               labels = c('Power', 'Threshold', 'Hill',
                                          'Spline', 'Gaussian Process', 
                                          '$\\beta_t$', 'No Behavioral Change'))

waicFinal30 <- waicFinal30[order(waicFinal30$alarmGen, 
                                 waicFinal30$meanWAIC),]
waicFinal30 <- waicFinal30[,-which(colnames(waicFinal30) == 'meanWAIC')]

kable(waicFinal30, row.names = F, format = 'latex', align = 'llcc', 
      booktabs = T, escape = F, 
      col.names = linebreak(c('\\textbf{Data generation}',
                              '\\textbf{Model fitted}', 
                              '\\textbf{WAIC}\n\\textbf{Mean (SD)}',
                              '\\textbf{\\% selected}'), align = 'c')) %>% 
    collapse_rows(columns = 1, latex_hline = 'major') 


################################################################################
# Supplemental Figures 1-3 - Posterior dists for parameters in true models
################################################################################

theme_set(theme_bw() + 
              theme(strip.background = element_rect(fill = 'white'),
                    strip.text = element_text(size = 16),
                    axis.title = element_text(size = 16),
                    axis.text = element_text(size = 14),
                    plot.title = element_text(size = 17, h = 0.5)))


### load truth
paramsTruth <- read.xlsx('simParamsSummary.xlsx')

### wide to long
paramsTruth <- reshape(paramsTruth, 
                       varying = c("beta", "delta", "H", "nu", "x0", "k", "rateI"), 
                       v.names = "truth",
                       timevar = "param", 
                       times = c("beta", "delta", "H", "nu", "x0", "k", "rateI"), 
                       new.row.names = 1:1000,
                       direction = "long")

paramsTruth <- paramsTruth[-which(is.na(paramsTruth$truth)),]
paramsTruth <- paramsTruth[order(paramsTruth$alarmGen, paramsTruth$smoothWindow, paramsTruth$param),]
paramsTruth <- paramsTruth[-which(colnames(paramsTruth) %in% c('infPeriod', 'id'))]

paramsPostAll <- readRDS('./results/paramsPostAll.rds')

# only keep rows where alarmGen == alarmFit
paramsPostAll <- paramsPostAll[paramsPostAll$alarmGen == paramsPostAll$alarmFit,]

# merge with truth
paramsPostAll <- merge(paramsPostAll, paramsTruth, by = c('param', 'alarmGen', 'smoothWindow'),
                       all.x = T)

paramsPostAll <- paramsPostAll[order(paramsPostAll$alarmGen, 
                                     paramsPostAll$smoothWindow,
                                     paramsPostAll$simNumber, 
                                     paramsPostAll$param),]

paramsPostAll$param <- factor(paramsPostAll$param, 
                              levels=c('beta', 'delta', 'H', 'k', 'nu', 'x0', 'rateI'),
                              labels=c('beta', 'delta', 'H', 'k', 'nu', 'x[0]', 'gamma'))

paramsPostAll$smoothWindow <- factor(paramsPostAll$smoothWindow , 
                                     labels=c('14-day smoothing', '30-day smoothing'))


### power alarm
pdf('./figures/supp_fig1_sim_powerParams.pdf', height = 6, width = 9)
ggplot(data = subset(paramsPostAll, alarmFit == 'power'), 
       aes(x = simNumber, y = mean)) +
    geom_point(size = 1) + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.9,
                  position = position_dodge(width = 0.9)) +
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, size = 1) +
    facet_wrap(smoothWindow ~ param, nrow = 2, scales = 'free_y',
               labeller =  labeller(smoothWindow = label_value, param = label_parsed)) +
    labs(x = 'Simulation Number', y = '')
dev.off()

### threshold alarm
pdf('./figures/supp_fig2_sim_threshParams.pdf', height = 6, width = 11)
ggplot(data = subset(paramsPostAll, alarmFit == 'thresh'), 
       aes(x = simNumber, y = mean)) +
    geom_point(size = 1) + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.9,
                  position = position_dodge(width = 0.9)) +
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, size = 1) +
    facet_wrap(smoothWindow ~ param, nrow = 2, scales = 'free_y',
               labeller =  labeller(smoothWindow = label_value, param = label_parsed)) +
    labs(x = 'Simulation Number', y = '') 
dev.off()

### hill alarm
pdf('./figures/supp_fig3_sim_hillParams.pdf', height = 6, width = 13)
ggplot(data = subset(paramsPostAll, alarmFit == 'hill'), 
       aes(x = simNumber, y = mean)) +
    geom_point(size = 1) + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.9,
                  position = position_dodge(width = 0.9)) +
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, size = 1) +
    facet_wrap(smoothWindow ~ param, nrow = 2, scales = 'free_y',
               labeller =  labeller(smoothWindow = label_value, param = label_parsed)) +
    labs(x = 'Simulation Number', y = '')
dev.off()


################################################################################
# Supplemental Figure 4 - Posterior for beta_t model
################################################################################

theme_set(theme_bw() + 
              theme(strip.background = element_rect(fill = 'white'),
                    strip.text = element_text(size = 14),
                    axis.title = element_text(size = 14),
                    axis.text = element_text(size = 12),
                    plot.title = element_text(size = 14, h = 0.5)))

set.seed(1)
simNumber <- round(runif(1, 0.5, 50.5))

times <- 1:50

# get true epidemic curves for each scenario
trueCurveThresh14 <- readRDS(paste0('./Data/thresh_14.rds'))[simNumber,times]
trueCurveThresh30 <- readRDS(paste0('./Data/thresh_30.rds'))[simNumber,times]
trueCurveHill14 <- readRDS(paste0('./Data/hill_14.rds'))[simNumber,times]
trueCurveHill30 <- readRDS(paste0('./Data/hill_30.rds'))[simNumber,times]
trueCurvePower14 <- readRDS(paste0('./Data/power_14.rds'))[simNumber,times]
trueCurvePower30 <- readRDS(paste0('./Data/power_30.rds'))[simNumber,times]


# using epidemic trajectory and true parameters, find value of alarm[t]/beta[t]
paramsTruth <- read.xlsx('simParamsSummary.xlsx')

N <- 1e6
trueAlarmThresh14 <- thresholdAlarm(movingAverage(trueCurveThresh14, 14),
                                    N = N, 
                                    delta = paramsTruth$delta[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 14)],
                                    H = paramsTruth$H[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 14)])
trueAlarmThresh30 <- thresholdAlarm(movingAverage(trueCurveThresh30, 30), 
                                    N = N, 
                                    delta = paramsTruth$delta[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 30)],
                                    H = paramsTruth$H[
                                        (paramsTruth$alarmGen == 'thresh' & 
                                             paramsTruth$smoothWindow == 30)])
trueAlarmHill14 <- hillAlarm(movingAverage(trueCurveHill14, 14), 
                             nu = paramsTruth$nu[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 14)],
                             x0 = paramsTruth$x0[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 14)],
                             delta = paramsTruth$delta[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 14)])
trueAlarmHill30 <- hillAlarm(movingAverage(trueCurveHill30, 30), 
                             nu = paramsTruth$nu[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 30)],
                             x0 = paramsTruth$x0[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 30)],
                             delta = paramsTruth$delta[
                                 (paramsTruth$alarmGen == 'hill' & 
                                      paramsTruth$smoothWindow == 30)])
trueAlarmPower14 <- powerAlarm(movingAverage(trueCurvePower14, 14), 
                               N = N, 
                               k = paramsTruth$k[
                                   (paramsTruth$alarmGen == 'power' & 
                                        paramsTruth$smoothWindow == 14)])
trueAlarmPower30 <- powerAlarm(movingAverage(trueCurvePower30, 30),
                               N = N, 
                               k = paramsTruth$k[
                                   (paramsTruth$alarmGen == 'power' & 
                                        paramsTruth$smoothWindow == 30)])

trueAlarmsSim <- data.frame(time = rep(times, 6),
                            trueAlarm = c(trueAlarmThresh14,
                                          trueAlarmThresh30,
                                          trueAlarmHill14,
                                          trueAlarmHill30,
                                          trueAlarmPower14,
                                          trueAlarmPower30),
                            alarmGen = c(rep('thresh', length(times)*2),
                                         rep('hill', length(times)*2),
                                         rep('power', length(times)*2)),
                            smoothWindow = rep(rep(c(14, 30), each = length(times)), 3))

# from true alarms, get true beta[t]
trueBeta <- trueAlarmsSim
trueBeta$trueBeta <- paramsTruth$beta[1] * (1 - trueBeta$trueAlarm)

### load posterior estimates
betaPostAll <- readRDS('./results/betaPostAll.rds')
betaPostAll <- betaPostAll[betaPostAll$simNumber == simNumber,]


betaPostAll$smoothWindow <- factor(betaPostAll$smoothWindow , 
                                   labels=c('14-day smoothing', '30-day smoothing'))
betaPostAll$alarmGen <- factor(betaPostAll$alarmGen, 
                               levels = c('power', 'thresh', 'hill'),
                               labels=c('Power', 'Threshold', 'Hill'))
trueBeta$smoothWindow <- factor(trueBeta$smoothWindow , 
                                labels=c('14-day smoothing', '30-day smoothing'))
trueBeta$alarmGen <- factor(trueBeta$alarmGen, 
                            levels = c('power', 'thresh', 'hill'),
                            labels=c('Power', 'Threshold', 'Hill'))

pdf('./figures/supp_fig4_sim_betatCurves.pdf', height = 5, width = 7)
ggplot() +  
    geom_line(data = subset(betaPostAll,), 
              aes(x = time, y = mean, group = simNumber), 
              col = 'orangered', size = 1) +
    geom_ribbon(data=subset(betaPostAll),
                aes(x = time, ymin = lower, ymax = upper), 
                alpha = 0.3, fill = 'orangered') +
    geom_line(data = subset(trueBeta), 
              aes(x = time, y = trueBeta), col = 'black') +
    facet_wrap(~smoothWindow + alarmGen, nrow = 2, 
               labeller =  labeller(smoothWindow = label_value, alarmGen = label_parsed)) +
    labs(x = 'Epidemic time', y = expression(beta[t])) +
    ggtitle(expression(paste('Posterior mean and 95% credible intervals for ', beta[t], ' model')))
dev.off()

