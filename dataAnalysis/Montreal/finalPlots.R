################################################################################
# Plots for Montreal Analysis
################################################################################

library(ggplot2)
library(grid)
library(gridExtra)
library(nimble)
library(plyr)
library(knitr)
library(scales)
library(kableExtra)

# functions to calculate the alarms
source('../scripts/modelCodes.R')

theme_set(theme_bw() + 
              theme(strip.background = element_rect(fill = 'white'),
                    strip.text = element_text(size = 18),
                    axis.title = element_text(size = 18),
                    axis.text = element_text(size = 16),
                    legend.title = element_text(size = 18),
                    legend.text = element_text(size = 16)))

infPeriodSpec <- 'fixed'

### read data
montreal <- read.csv('./Data/montrealClean.csv')
montreal$smoothedCases <- round(movingAverage(montreal$dailyCases, 3))
montreal$cumulativeCases <- cumsum(montreal$smoothedCases)
montreal$date <- as.Date(montreal$date)

### set up ggplot theme

theme_set(theme_bw() + 
              theme(strip.background = element_rect(fill = 'white'),
                    strip.text = element_text(size = 16),
                    axis.title = element_text(size = 16),
                    axis.text = element_text(size = 14),
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14),
                    plot.title = element_text(size = 18, h = 0.5)))

infPeriodSpec <- 'fixed'


################################################################################
# Plot 1 - Montreal data indicating peak locations

montreal$Peak <- factor(montreal$peak)

pal <- c('mediumpurple2', 'darkolivegreen3', 'darkorange2', 'steelblue1')

pdf('./Figures/montreal_data.pdf', width = 7, height = 4)
ggplot(montreal, aes(x = date, y = smoothedCases)) + 
    geom_line(linetype = 2) + 
    geom_line(data = subset(montreal, peak == 1), col = pal[1], lwd = 1.2) + 
    geom_line(data = subset(montreal, peak == 2), col = pal[2], lwd = 1.2) + 
    geom_line(data = subset(montreal, peak == 4), col = pal[3], lwd = 1.2) +
    geom_line(data = subset(montreal, peak == 5), col = pal[4], lwd = 1.2) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = 'Date', y = 'Incidene', title = 'Montreal COVID-19 Case Counts') +
    annotate('text', x = as.Date('2020-04-25'), y = 1000, label = 'Peak 1') +
    annotate('text', x = as.Date('2021-01-15'), y = 1500, label = 'Peak 2') +
    annotate('text', x = as.Date('2021-12-25'), y = 5000, label = 'Peak 3') +
    annotate('text', x = as.Date('2022-04-10'), y = 1100, label = 'Peak 4')
dev.off()

################################################################################
# Plot 2 - all alarm Functions estimated by peak

### load posterior estimates of alarm functions
alarmAll <- readRDS('./resultsFinal/alarmPostAll.rds')
alarmAll <- alarmAll[alarmAll$infPeriod == infPeriodSpec,]

# format for better plotting
alarmAll$alarmFit <- factor(alarmAll$alarmFit,
                            levels = c('spline', 'gp', 'thresh', 'hill', 'power'),
                            labels = c('Spline', 'Gaussian Process', 
                                       'Threshold', 'Hill', 'Power'))

alarmAll$xAlarmScale <- alarmAll$xAlarm / montreal$Population[1] * 100000

# peak 1


pal <- c('mediumpurple3', 'chartreuse3', 'darkorange2', 'steelblue3')

p1 <- ggplot(subset(alarmAll, peak == 1), aes(x = xAlarmScale, y = mean)) +  
    geom_line(col = pal[1]) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill = pal[1]) +
    facet_wrap(~alarmFit, nrow = 1) +
    labs(x = '30-day incidence per 100,000', y = 'Alarm') + 
    ggtitle('Peak 1') 

p2 <- ggplot(subset(alarmAll, peak == 2), aes(x = xAlarmScale, y = mean)) +  
    geom_line(col = pal[2]) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill = pal[2]) +
    facet_wrap(~alarmFit, nrow = 1) +
    labs(x = '30-day incidence per 100,000', y = 'Alarm') + 
    ggtitle('Peak 2') 


p3 <- ggplot(subset(alarmAll, peak == 4), aes(x = xAlarmScale, y = mean)) +  
    geom_line(col = pal[3]) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill = pal[3]) +
    facet_wrap(~alarmFit, nrow = 1) +
    labs(x = '30-day incidence per 100,000', y = 'Alarm') + 
    ggtitle('Peak 3') 


p4 <- ggplot(subset(alarmAll, peak == 5), aes(x = xAlarmScale, y = mean)) +  
    geom_line(col = pal[4]) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill = pal[4]) +
    facet_wrap(~alarmFit, nrow = 1) +
    labs(x = '30-day incidence per 100,000', y = 'Alarm') + 
    ggtitle('Peak 4') 

pdf('./Figures/montreal_allAlarms.pdf', width = 10, height = 12)
grid.arrange(p1, p2, p3, p4, ncol = 1)
dev.off()


################################################################################
# Plot 3 - spline alarm functions by peak

pal <- c('mediumpurple3', 'limegreen', 'darkorange2', 'royalblue1')

alarmAll <- alarmAll[alarmAll$peak %in% c(1,2,4,5),]
alarmAll$peak <- factor(alarmAll$peak,
                        labels = paste0('Peak ', 1:4))


pdf('./Figures/montreal_splineAlarms.pdf', height = 3.5, width = 12)
ggplot(subset(alarmAll, alarmFit == 'Spline' & smoothWindow == 14), 
       aes(x = xAlarmScale, y = mean)) +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill = peak), alpha=0.3) +
    geom_line(aes(color = peak), size = 1) +
    facet_wrap(~peak, nrow = 1, scales = 'free_x') +
    labs(x = '30-day incidence per 100,000', y = 'Alarm') + 
    ggtitle('Spline Alarm Functions') +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) + 
    guides(fill = 'none', color = 'none')
dev.off()

# nyc$peakFac <- factor(nyc$peak, 
#                       labels = paste0('Peak ', 1:4))
# 
# p1 <- ggplot(subset(nyc, !is.na(peak)), aes(x = date, y = smoothedCases, color = peakFac)) + 
#     geom_line(size = 1) + 
#     facet_wrap(~peakFac, nrow = 1, scales = 'free_x') +
#     scale_y_continuous(labels = scales::comma) +
#     labs(x = 'Date', y = 'Incidene', title = 'NYC COVID-19 Case Counts') +
#     scale_color_manual(values = pal)+ 
#     guides(color = 'none')
# 
# grid.arrange(p1, p2, ncol = 1)


################################################################################
# R0 over time

### load posterior estimates of R0 over epidemic time
r0All <- readRDS('./resultsFinal/r0PostAll.rds')
r0All <- r0All[r0All$infPeriod == infPeriodSpec,]

# format for better plotting
r0All$alarmFit <- factor(r0All$alarmFit,
                         levels = c('spline', 'gp', 'thresh', 'hill', 'power',
                                    'betatSpline', 'basic'),
                         labels = c('Spline', 'Gaussian Process', 
                                    'Threshold', 'Hill', 'Power',
                                    'Beta[t]', 'No Behavioral Change'))
r0All <- r0All[order(r0All$alarmFit, 
                     r0All$peak, 
                     r0All$time),]


datNew <- montreal

idxStart <- 15

datNew$timeFull <- 1:nrow(datNew)-idxStart
datNew$timeFull[1:idxStart] <- NA
datNew$timePeak5 <- datNew$timePeak4 <- datNew$timePeak3 <- datNew$timePeak2 <- datNew$timePeak1 <- NA

# peak 1 time 1 is row 6
datNew$timePeak1[(idxStart + 1):(max(which(datNew$peak == 1)))] <-
    1:(max(which(datNew$peak == 1)) - idxStart)

datNew$timePeak2[(min(which(datNew$peak == 2)) + 1):(max(which(datNew$peak == 2)))] <-
    (min(which(datNew$peak == 2)) + 1):(max(which(datNew$peak == 2))) -
    min(which(datNew$peak == 2))

datNew$timePeak3[(min(which(datNew$peak == 3)) + 1):(max(which(datNew$peak == 3)))] <-
    (min(which(datNew$peak == 3)) + 1):(max(which(datNew$peak == 3))) -
    min(which(datNew$peak == 3))

datNew$timePeak4[(min(which(datNew$peak == 4)) + 1):(max(which(datNew$peak == 4)))] <-
    (min(which(datNew$peak == 4)) + 1):(max(which(datNew$peak == 4))) -
    min(which(datNew$peak == 4))

datNew$timePeak5[(min(which(datNew$peak == 5)) + 1):(max(which(datNew$peak == 5)))] <-
    (min(which(datNew$peak == 5)) + 1):(max(which(datNew$peak == 5))) -
    min(which(datNew$peak == 5))

r0Peak1 <- subset(r0All, peak == '1')
r0Peak1 <- merge(r0Peak1, datNew, 
                 by.x = 'time', by.y = 'timePeak1', 
                 all.x = T)

r0Peak2 <- subset(r0All, peak == '2')
r0Peak2 <- merge(r0Peak2, datNew, 
                 by.x = 'time', by.y = 'timePeak2', 
                 all.x = T)

r0Peak3 <- subset(r0All, peak == '3')
r0Peak3 <- merge(r0Peak3, datNew, 
                 by.x = 'time', by.y = 'timePeak3', 
                 all.x = T)

r0Peak4 <- subset(r0All, peak == '4')
r0Peak4 <- merge(r0Peak4, datNew, 
                 by.x = 'time', by.y = 'timePeak4', 
                 all.x = T)

r0Peak5 <- subset(r0All, peak == '5')
r0Peak5 <- merge(r0Peak5, datNew, 
                 by.x = 'time', by.y = 'timePeak5', 
                 all.x = T)

r0AllPeaks <- rbind.data.frame(r0Peak1[,-c(17:20)], 
                               r0Peak2[,-c(17:20)], 
                               r0Peak3[,-c(17:20)], 
                               r0Peak4[,-c(17:20)], 
                               r0Peak5[,-c(17:20)])


ggplot(r0AllPeaks, aes(x = date, y = mean)) +  
    geom_line() +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill = 'blue') +
    facet_wrap(~peak.x + alarmFit, nrow = 4,
               scales = 'free') +
    labs(x = 'Date', y = 'R0') + 
    geom_hline(yintercept = 1, linetype = 2)

# R0 spline only

pal <- c('mediumpurple3', 'limegreen', 'darkorange2', 'royalblue1')

r0AllPeaks <- r0AllPeaks[r0AllPeaks$peak.x != 3,]
r0AllPeaks$peak.x <- factor(r0AllPeaks$peak.x,
                        labels = paste0('Peak ', 1:4))


pdf('./Figures/montreal_splineR0.pdf', height = 3.5, width = 12)
ggplot(subset(r0AllPeaks, alarmFit == 'Spline' & smoothWindow == 30), 
       aes(x = date, y = mean)) +  
  geom_line(aes(col = peak.x), size = 1) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = peak.x), alpha=0.3) +
  facet_wrap(~peak.x , nrow = 1,
             scales = 'free') +
  labs(x = 'Date', y = 'R0') + 
  ggtitle('Reproductive number over time') +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) + 
  guides(fill = 'none', color = 'none') +
    scale_x_date(date_labels = '%b')
dev.off()

################################################################################
# WAIC table by peak

options(scipen=20)
waicAll <- readRDS('./resultsFinal/waicAll.rds')
waicAll <- waicAll[waicAll$infPeriod == infPeriodSpec,]

waicAll <- waicAll[order(waicAll$peak, waicAll$waic),]
# summaries of WAIC values by model
waicAll$waic <- round(waicAll$waic, 2)
waicAll$lppd <- round(waicAll$lppd, 2)
waicAll$pWAIC <- round(waicAll$pWAIC, 2)

waicAll$alarmFit <- factor(waicAll$alarmFit,
                            levels = c('spline', 'gp', 'thresh', 'hill', 'power',
                                       'basic', 'betatSpline'),
                            labels = c('Spline', 'Gaussian Process', 
                                       'Threshold', 'Hill', 'Power',
                                       'No Behavioral Change', '$\\beta_t$'))

kable(waicAll[,c('peak', 'alarmFit', 'smoothWindow', 'waic')], 
      format = 'latex', align = 'llc', 
      row.names = F, booktabs = T, escape = F, 
      col.names = linebreak(c('\\textbf{Peak}', 
                              '\\textbf{Model fitted}', 
                              '\\textbf{WAIC}'), align = 'c')) %>% 
    collapse_rows(columns = 1, latex_hline = 'major') 



