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

# functions to calculate the alarms
source('./scripts/modelCodes.R')

theme_set(theme_bw() + 
              theme(strip.background = element_rect(fill = 'white'),
                    strip.text = element_text(size = 16),
                    axis.title = element_text(size = 16),
                    axis.text = element_text(size = 14),
                    plot.title = element_text(size = 17, h = 0.5)))

################################################################################
# For supplemental material - posterior distributions of parameters from each 
#   epidemic


### load truth
paramsTruth <- read.xlsx('simParamsSummary.xlsx')
paramsTruth <- paramsTruth[paramsTruth$infPeriod == 'fixed',]

### wide to long
paramsTruth <- reshape(paramsTruth, 
                       varying = c("beta", "delta", "H", "nu", "x0", "k"), 
                       v.names = "truth",
                       timevar = "param", 
                       times = c("beta", "delta", "H", "nu", "x0", "k"), 
                       new.row.names = 1:1000,
                       direction = "long")

paramsTruth <- paramsTruth[-which(is.na(paramsTruth$truth)),]
paramsTruth <- paramsTruth[order(paramsTruth$alarmGen, paramsTruth$smoothWindow, paramsTruth$param),]
paramsTruth <- paramsTruth[-which(colnames(paramsTruth) %in% c('infPeriod', 'id'))]

paramsPostAll <- readRDS('./resultsFinal/paramsPostAll.rds')

# only keep rows where alarmGen == alarmFit
paramsPostAll <- paramsPostAll[paramsPostAll$alarmGen == paramsPostAll$alarmFit,]

# merge with truth
paramsPostAll <- merge(paramsPostAll, paramsTruth, by = c('param', 'alarmGen', 'smoothWindow'),
                       all.x = T)

paramsPostAll <- paramsPostAll[order(paramsPostAll$alarmGen, 
                                     paramsPostAll$smoothWindow,
                                     paramsPostAll$simNumber, 
                                     paramsPostAll$param),]

paramsPostAll$param <- factor(paramsPostAll$param , 
                              labels=c('beta', 'delta', 'H', 'k', 'nu', 'x[0]'))

paramsPostAll$smoothWindow <- factor(paramsPostAll$smoothWindow , 
                                     labels=c('14-day smoothing', '30-day smoothing'))

### threshold alarm
pdf('./Figures/sim_threshParams.pdf', height = 6, width = 9)
ggplot(data = subset(paramsPostAll, alarmFit == 'thresh'), 
       aes(x = simNumber, y = mean)) +
    geom_point() + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.9,
                  position = position_dodge(width = 0.9)) +
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2) +
    facet_wrap(smoothWindow ~ param, nrow = 2, scales = 'free_y',
               labeller =  labeller(smoothWindow = label_value, param = label_parsed)) +
    labs(x = 'Simulation Number', y = '') 
dev.off()

### hill alarm
pdf('./Figures/sim_hillParams.pdf', height = 6, width = 12)
ggplot(data = subset(paramsPostAll, alarmFit == 'hill'), 
       aes(x = simNumber, y = mean)) +
    geom_point() + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.9,
                  position = position_dodge(width = 0.9)) +
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2) +
    facet_wrap(smoothWindow ~ param, nrow = 2, scales = 'free_y',
               labeller =  labeller(smoothWindow = label_value, param = label_parsed)) +
    labs(x = 'Simulation Number', y = '')
dev.off()

### power alarm
pdf('./Figures/sim_powerParams.pdf', height = 6, width = 6)
ggplot(data = subset(paramsPostAll, alarmFit == 'power'), 
       aes(x = simNumber, y = mean)) +
    geom_point() + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.9,
                  position = position_dodge(width = 0.9)) +
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2) +
    facet_wrap(smoothWindow ~ param, nrow = 2, scales = 'free_y',
               labeller =  labeller(smoothWindow = label_value, param = label_parsed)) +
    labs(x = 'Simulation Number', y = '')
dev.off()



################################################################################
# posterior means of alarm functions
# 14-day in supplemental, 30-day in paper

### get parameters for true alarms
paramsTruth <- read.xlsx('simParamsSummary.xlsx')
paramsTruth <- paramsTruth[paramsTruth$infPeriod == 'fixed',]

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
              col = adjustcolor('grey40', alpha = 0.5)) +
    geom_line(data = subset(trueAlarms, smoothWindow == 14 & alarmGen == 'power'), 
              aes(x = xAlarm, y = trueAlarm), col = 'red', size = 1) +
    facet_wrap(~alarmFit) + 
    labs(x = '14-day average incidence', y = 'Alarm')+
    ylim(0, 1) + 
    ggtitle('Power Alarm')

p2 <- ggplot() +  
    geom_line(data = subset(alarmAll, smoothWindow == 14 &  alarmGen == 'thresh'), 
              aes(x = xAlarm, y = mean, group = simNumber), 
              col = adjustcolor('grey40', alpha = 0.5)) +
    geom_line(data = subset(trueAlarms, smoothWindow == 14 & alarmGen == 'thresh'), 
              aes(x = xAlarm, y = trueAlarm), col = 'red', size = 1) +
    facet_wrap(~alarmFit) + 
    labs(x = '14-day average incidence', y = 'Alarm') +
    ylim(0, 1)  + 
    ggtitle('Threshold Alarm')

p3 <- ggplot() +  
    geom_line(data = subset(alarmAll, smoothWindow == 14 & alarmGen == 'hill'), 
              aes(x = xAlarm, y = mean, group = simNumber), 
              col = adjustcolor('grey40', alpha = 0.5)) +
    geom_line(data = subset(trueAlarms, smoothWindow == 14 & alarmGen == 'hill'), 
              aes(x = xAlarm, y = trueAlarm), col = 'red', size = 1) +
    facet_wrap(~alarmFit) + 
    labs(x = '14-day average incidence', y = 'Alarm')+
    ylim(0, 1)  + 
    ggtitle('Hill Alarm') 



pdf('./Figures/sim_alarms14.pdf', height = 7, width = 7)
grid.arrange(p1, p2, p3, nrow = 3,
             top=textGrob("Alarm Functions Posterior Mean", gp = gpar(fontsize = 16, font = 1)))

dev.off()




