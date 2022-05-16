################################################################################
# Output plots as pdf images for publication
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
                  strip.text = element_text(size = 16),
                  axis.title = element_text(size = 16),
                  axis.text = element_text(size = 14),
                  plot.title = element_text(size = 17, h = 0.5)))

infPeriodSpec <- 'exp'


################################################################################
# describing simulation study true alarm function and data plots

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

pdf('./Figures/sim_setup.pdf', height = 7, width = 9)
grid.arrange(p1, p2, ncol = 1)
dev.off()

################################################################################
# For supplemental material - posterior distributions of parameters from each 
#   epidemic


### load truth
paramsTruth <- read.xlsx('simParamsSummary.xlsx')
paramsTruth <- paramsTruth[paramsTruth$infPeriod == infPeriodSpec,]

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

paramsPostAll <- readRDS('./resultsFinal/paramsPostAll.rds')
paramsPostAll <- paramsPostAll[paramsPostAll$infPeriod == infPeriodSpec,]

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
pdf('./Figures/sim_powerParams.pdf', height = 6, width = 9)
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
pdf('./Figures/sim_threshParams.pdf', height = 6, width = 11)
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
pdf('./Figures/sim_hillParams.pdf', height = 6, width = 13)
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
# posterior means of alarm functions
# 14-day in supplemental, 30-day in paper

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



pdf('./Figures/sim_alarms14.pdf', height = 8, width = 7)
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



pdf('./Figures/sim_alarms30.pdf', height = 8, width = 7)
grid.arrange(p1, p2, p3, nrow = 3,
             top=textGrob("         Posterior mean of alarm functions", gp = gpar(fontsize = 16, font = 1)))

dev.off()

################################################################################
# posterior prediction 
# 14-day in supplemental, 30-day in paper


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


pdf('./Figures/sim_postPred14.pdf', height = 8, width = 9)
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


pdf('./Figures/sim_postPred30.pdf', height = 8, width = 9)
grid.arrange(p1, p2, p3, nrow = 3,
             top=textGrob("        Posterior predictive forecasting", gp = gpar(fontsize = 18, font = 1)))

dev.off()


################################################################################
# beta[t] figure
# supplemental only


set.seed(1)
simNumber <- round(runif(1, 0.5, 50.5))

times <- 1:50

# get true epidemic curves for each scenario
trueCurveThresh14 <- readRDS(paste0('./Data/thresh_', infPeriodSpec, '_14.rds'))[simNumber,times]
trueCurveThresh30 <- readRDS(paste0('./Data/thresh_', infPeriodSpec, '_30.rds'))[simNumber,times]
trueCurveHill14 <- readRDS(paste0('./Data/hill_', infPeriodSpec, '_14.rds'))[simNumber,times]
trueCurveHill30 <- readRDS(paste0('./Data/hill_', infPeriodSpec, '_30.rds'))[simNumber,times]
trueCurvePower14 <- readRDS(paste0('./Data/power_', infPeriodSpec, '_14.rds'))[simNumber,times]
trueCurvePower30 <- readRDS(paste0('./Data/power_', infPeriodSpec, '_30.rds'))[simNumber,times]


# using epidemic trajectory and true parameters, find value of alarm[t]/beta[t]
paramsTruth <- read.xlsx('simParamsSummary.xlsx')
paramsTruth <- paramsTruth[paramsTruth$infPeriod == infPeriodSpec,]

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
betaPostAll <- readRDS('./resultsFinal/betaPostAll.rds')
betaPostAll <- betaPostAll[betaPostAll$simNumber == simNumber &
                             betaPostAll$infPeriod == infPeriodSpec,]


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

pdf('./Figures/sim_betatCurves.pdf', height = 5, width = 7)
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

################################################################################
# WAIC Table or Figure
# 14-day in supplemental, 30-day in paper

waicAll <- readRDS('./resultsFinal/waicAll.rds')
waicAll<- waicAll[waicAll$infPeriod == infPeriodSpec,]

waicMin <- ddply(waicAll, .(alarmGen, smoothWindow, simNumber), summarize,
                 minWAIC = min(waic))

waicAll <- merge(waicAll, waicMin, 
                 by = c('alarmGen', 'smoothWindow', 'simNumber'),
                 all.x = T)

waicSelected <- ddply(waicAll, .(alarmGen, smoothWindow, simNumber), summarize,
                      selected = alarmFit[which(waic == minWAIC)])
waicSelected$selected <- factor(waicSelected$selected,
                                levels = c('power', 'thresh', 'hill',
                                           'spline', 'gp', 'betatSpline', 'basic'))


# format for better plotting
alarmFitLabs <- c('Power', 'Threshold', 'Hill',
                  'Spline', 'Gaussian~Process', 
                  'beta[t]', 'No~Behavioral~Change')

waicAll$alarmFit <- factor(waicAll$alarmFit,
                           levels = rev(c('power', 'thresh', 'hill',
                                          'spline', 'gp', 'betatSpline', 'basic')),
                           labels = rev(alarmFitLabs))

waicAll$alarmGen <- factor(waicAll$alarmGen,
                           levels = c('power', 'thresh', 'hill'),
                           labels = c('Power', 'Threshold', 'Hill'))

# 14 day smoothing selections
waicTab <- with(subset(waicSelected, smoothWindow == 14), 
                table(selected, alarmGen, exclude = NULL))
waicTab <- as.data.frame(prop.table(waicTab, 2))
waicTab$smoothWindow <- 14
colnames(waicTab)[1] <- 'alarmFit'

waicTab <- waicTab[-which(waicTab$alarmGen == 'power' & 
                            waicTab$alarmFit %in% c('thresh', 'hill')),]
waicTab <- waicTab[-which(waicTab$alarmGen == 'thresh' & 
                            waicTab$alarmFit %in% c('power', 'hill')),]
waicTab <- waicTab[-which(waicTab$alarmGen == 'hill' & 
                            waicTab$alarmFit %in% c('power', 'thresh')),]


waicTab$alarmFit <- factor(waicTab$alarmFit,
                           levels = rev(c('power', 'thresh', 'hill',
                                          'spline', 'gp', 'betatSpline', 'basic')),
                           labels = rev(c(alarmFitLabs)))

waicTab$alarmGen <- factor(waicTab$alarmGen,
                           levels = c('power', 'thresh', 'hill'),
                           labels = c('Power', 'Threshold', 'Hill'))

myLabFun <- function(x) {
  parse(text=x)
}


# boxplots?
pdf('./Figures/sim_waic14.pdf', height = 7, width = 6)
ggplot(subset(waicAll, smoothWindow == 14), aes(x = alarmFit, y = waic, fill = alarmFit)) +
  geom_boxplot(outlier.size  = 1) +
  geom_text(data = waicTab, 
            aes(x = alarmFit, y = 1500, 
                label = scales::percent(Freq)), size = 5) + 
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1, 5, 1, 1), "lines"),
        strip.background = element_blank()) +
  scale_x_discrete(labels = myLabFun) +
  coord_flip(clip = "off", ylim = c(min(waicAll$waic), 1250)) + 
  facet_wrap(~alarmGen, scales = 'free', ncol = 1) + 
  labs(x = '', y = 'WAIC') +
  geom_text(aes(x = Inf, y = 1480), label = c('% chosen', rep(NA, 749)),
            vjust = -31.2, size = 5, fontface = 2) 
dev.off()



# 30 day smoothing selections
waicTab <- with(subset(waicSelected, smoothWindow == 30), 
                table(selected, alarmGen, exclude = NULL))
waicTab <- as.data.frame(prop.table(waicTab, 2))
waicTab$smoothWindow <- 30
colnames(waicTab)[1] <- 'alarmFit'

waicTab <- waicTab[-which(waicTab$alarmGen == 'power' & 
                            waicTab$alarmFit %in% c('thresh', 'hill')),]
waicTab <- waicTab[-which(waicTab$alarmGen == 'thresh' & 
                            waicTab$alarmFit %in% c('power', 'hill')),]
waicTab <- waicTab[-which(waicTab$alarmGen == 'hill' & 
                            waicTab$alarmFit %in% c('power', 'thresh')),]


waicTab$alarmFit <- factor(waicTab$alarmFit,
                           levels = rev(c('power', 'thresh', 'hill',
                                          'spline', 'gp', 'betatSpline', 'basic')),
                           labels = rev(c(alarmFitLabs)))

waicTab$alarmGen <- factor(waicTab$alarmGen,
                           levels = c('power', 'thresh', 'hill'),
                           labels = c('Power', 'Threshold', 'Hill'))

endLim <- 1550

pdf('./Figures/sim_waic30.pdf', height = 7, width = 6)
ggplot(subset(waicAll, smoothWindow == 30), aes(x = alarmFit, y = waic, fill = alarmFit)) +
  geom_boxplot(outlier.size  = 1) +
  geom_text(data = waicTab, 
            aes(x = alarmFit, y = endLim + 280, 
                label = scales::percent(Freq)), size = 5) + 
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(1, 5, 1, 1), "lines"),
        strip.background = element_blank()) +
  scale_x_discrete(labels = myLabFun) +
  coord_flip(clip = "off", ylim = c(min(waicAll$waic), endLim)) + 
  facet_wrap(~alarmGen, scales = 'free', ncol = 1) + 
  labs(x = '', y = 'WAIC') +
  geom_text(aes(x = Inf, y = endLim + 290), label = c('% chosen', rep(NA, 749)),
            vjust = -31.2, size = 5, fontface = 2) 
dev.off()


### Summary table with mean (SD) and % selected WAIC

waicAll <- readRDS('./resultsFinal/waicAll.rds')
waicAll<- waicAll[waicAll$infPeriod == infPeriodSpec,]

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











