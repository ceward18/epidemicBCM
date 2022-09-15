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


