################################################################################
# Figures in paper
# NYC Analysis
################################################################################

# load libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(nimble)
library(plyr)
library(dplyr)
library(knitr)
library(kableExtra)
library(ggh4x)

# for movingAverage function
source('./scripts/modelCodes.R')

################################################################################
### Pull in and set up all results needed for figures/tables
################################################################################

################################################################################
### set up NYC  data
dat <- read.csv('./Data/nycClean.csv')
dat <- dat[-c((max(which(dat$peak == 2)) + 1):nrow(dat)),]

dat$date <- as.Date(dat$date)
dat$smoothedCases <- round(movingAverage(dat$dailyCases, 7))
dat <- dat[dat$date < as.Date('2022-03-15'), ]
dat$Peak <- factor(dat$peak)

################################################################################
### Gelman-rubin

grAll <- readRDS('./results/grAll.rds')
grAll <- subset(grAll, peak == c(1,2))

# which didn't converge
notConverge <- grAll[which(grAll$gr > 1.1),  ]
notConvergeModels <- notConverge[
    !duplicated(notConverge
                [,-which(colnames(notConverge) %in% c('gr', 'grUpper', 'param'))]),
    c('alarmFit', 'smoothWindow', 'prior', 'peak')]
notConvergeModels$noConverge <- 1

################################################################################
### WAIC

### load WAIC values
waicAll <- readRDS('./results/waicAll.rds')
waicAll <- subset(waicAll, peak == c(1,2))

### flag those that did not converge
waicAll <- merge(waicAll, notConvergeModels,
                 by = c('alarmFit', 'smoothWindow', 'prior', 'peak'),
                 all.x = T)
waicAll$noConverge[is.na(waicAll$noConverge)] <- 0
waicAll <- waicAll[order(waicAll$peak, waicAll$waic),]
waicAll$waic[waicAll$noConverge == 1] <- NA

# minimum for each peak
# minimum WAIC within each prior/peak combination
minWAIC <- ddply(subset(waicAll,  noConverge == 0 ), 
                 .(peak, prior), summarize,
                 alarmFit = alarmFit[which.min(waic)],
                 smoothWindow = smoothWindow[which.min(waic)])
minWAIC$isMin <- 1

################################################################################
### Posterior alarms

### load posterior estimates of alarm functions
alarmAll <- readRDS('./results/alarmPostAll.rds')
alarmSub <- subset(alarmAll, prior == 5)
alarmSub <- subset(alarmSub, peak == c(1,2))

### identify those that did not converge 
alarmSub <- merge(alarmSub, notConvergeModels,
                  by = c('alarmFit', 'smoothWindow', 'prior', 'peak'),
                  all.x = T)
alarmSub$noConverge[is.na(alarmSub$noConverge)] <- 0

# NA for those that did not converge
alarmSub$mean[alarmSub$noConverge == 1] <- NA
alarmSub$lower[alarmSub$noConverge == 1] <- NA
alarmSub$upper[alarmSub$noConverge == 1] <- NA

# wave 1 and 3 60-day incidence
# wave 2 30-day incidence
alarmSub <- alarmSub[-which(alarmSub$peak == 1 & alarmSub$smoothWindow == 30),]
alarmSub <- alarmSub[-which(alarmSub$peak == 2 & alarmSub$smoothWindow == 60),]

# format for better plotting
alarmSub$alarmFit <- factor(alarmSub$alarmFit,
                            levels = c('spline', 'gp',
                                       'thresh', 'hill', 'power'),
                            labels = c('Spline', 'Gaussian Process', 
                                       'Threshold', 'Hill', 'Power'))

alarmSub$Peak <- factor(alarmSub$peak, labels = paste0('Wave ', 1:2))

################################################################################
### Posterior R0

### load posterior estimates of R0 over epidemic time
r0All <- readRDS('./results/r0PostAll.rds')
r0All <- subset(r0All, peak == c(1,2))
r0All <- subset(r0All, prior == 5)

### identify those that did not converge 
r0All <- merge(r0All, notConvergeModels,
               by = c('alarmFit', 'smoothWindow', 'prior', 'peak'),
               all.x = T)
r0All$noConverge[is.na(r0All$noConverge)] <- 0

# NA for those that did not converge
r0All$mean[r0All$noConverge == 1] <- NA
r0All$lower[r0All$noConverge == 1] <- NA
r0All$upper[r0All$noConverge == 1] <- NA

# wave 1 and 3 60-day incidence
# wave 2 30-day incidence
r0All <- r0All[-which(r0All$peak == 1 & r0All$smoothWindow %in% c(30)),]
r0All <- r0All[-which(r0All$peak == 2 & r0All$smoothWindow %in% c(60)),]

### create columns numbering timing of each wave from 0
dat$timePeak2 <- dat$timePeak1 <- NA

# peak 1 time 1 is row 6
dat$timePeak1[6:(max(which(dat$peak == 1)))] <- 1:(max(which(dat$peak == 1)) - 5)

dat$timePeak2[(min(which(dat$peak == 2)) + 1):(max(which(dat$peak == 2)))] <-
    (min(which(dat$peak == 2)) + 1):(max(which(dat$peak == 2))) -
    min(which(dat$peak == 2)) 


r0Peak1 <- subset(r0All, peak == '1')
r0Peak1 <- merge(r0Peak1, dat, 
                 by.x = 'time', by.y = 'timePeak1', 
                 all.x = T)
r0Peak1 <- r0Peak1[,-which(colnames(r0Peak1) %in% paste0('timePeak', 1:2))]

r0Peak2 <- subset(r0All, peak == '2')
r0Peak2 <- merge(r0Peak2, dat, 
                 by.x = 'time', by.y = 'timePeak2', 
                 all.x = T)
r0Peak2 <- r0Peak2[,-which(colnames(r0Peak2) %in% paste0('timePeak', 1:2))]

r0PeakAll <- rbind.data.frame(r0Peak1, r0Peak2)

# format for better plotting
r0PeakAll$alarmFit <- factor(r0PeakAll$alarmFit,
                             levels = c('spline', 'gp',
                                        'thresh', 'hill', 'power',
                                        'betatSpline', 'basic'),
                             labels = c('Spline', 'Gaussian Process', 
                                        'Threshold', 'Hill', 'Power',
                                        'beta[t]', 'No Behavioral Change'))

################################################################################
### Posterior predictive fit

### load posterior predictive fit
postPredFitAll <- readRDS('./results/postPredFitAll.rds')
postPredFitAll <- subset(postPredFitAll, peak == c(1,2))
postPredFitAll <- subset(postPredFitAll, prior == 5)


### remove those that did not converge 
postPredFitAll <- merge(postPredFitAll, notConvergeModels,
                        by = c('alarmFit', 'smoothWindow', 'prior', 'peak'),
                        all.x = T)
postPredFitAll$noConverge[is.na(postPredFitAll$noConverge)] <- 0

# NA for those that did not converge
postPredFitAll$mean[postPredFitAll$noConverge == 1] <- NA
postPredFitAll$lower[postPredFitAll$noConverge == 1] <- NA
postPredFitAll$upper[postPredFitAll$noConverge == 1] <- NA

# wave 1 and 3 60-day incidence
# wave 2 30-day incidence
postPredFitAll <- postPredFitAll[
    -which(postPredFitAll$peak == 1 & 
               postPredFitAll$smoothWindow %in% c(30)),]
postPredFitAll <- postPredFitAll[
    -which(postPredFitAll$peak == 2 &
               postPredFitAll$smoothWindow %in% c(60)),]


# merge with dates from dataset
postPredFitPeak1 <- subset(postPredFitAll, peak == '1')
postPredFitPeak1 <- merge(postPredFitPeak1, dat, 
                          by.x = 'time', by.y = 'timePeak1', 
                          all.x = T)
postPredFitPeak1 <- postPredFitPeak1[,-which(colnames(postPredFitPeak1) %in%
                                                 paste0('timePeak', 1:2))]

postPredFitPeak2 <- subset(postPredFitAll, peak == '2')
postPredFitPeak2 <- merge(postPredFitPeak2, dat, 
                          by.x = 'time', by.y = 'timePeak2', 
                          all.x = T)
postPredFitPeak2 <- postPredFitPeak2[,-which(colnames(postPredFitPeak2) %in%
                                                 paste0('timePeak', 1:2))]

postPredFitPeakAll <- rbind.data.frame(postPredFitPeak1, 
                                       postPredFitPeak2)


################################################################################
# Figure 5: NYC Data
################################################################################

pal <- c('darkorchid2', 'darkorange2')

jpeg('./figures/fig5_nyc_data.jpg', units = 'in', res = 500, width = 6, height = 3)
ggplot(dat, aes(x = date, y = smoothedCases)) + 
    geom_line(linetype = 2) + 
    geom_line(data = subset(dat, peak == 1), col = pal[1], lwd = 1.2) + 
    geom_line(data = subset(dat, peak == 2), col = pal[2], lwd = 1.2) + 
    scale_y_continuous(labels = scales::comma, limits = c(0, 6000)) +
    scale_x_date(date_breaks = "3 month", date_minor_breaks = "1 month",
                 date_labels = "%b '%y") +
    labs(x = 'Date', y = 'Incidence', title = 'NYC COVID-19 Case Counts') +
    annotate('text', x = as.Date('2020-04-15'), y = 5800, label = 'Wave 1') +
    annotate('text', x = as.Date('2021-01-15'), y = 5800, label = 'Wave 2') +
    theme_bw() + 
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          plot.title = element_text(size = 11, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
dev.off()

################################################################################
# Table 2: All WAIC values for converged models
################################################################################

### WAIC table
waicTab <- merge(waicAll, minWAIC, by = c('peak', 'prior',
                                          'alarmFit', 'smoothWindow'),
                 all.x = T)
waicTab <- waicTab[waicTab$prior == 5,
                   c('peak', 'alarmFit', 'smoothWindow', 'waic')]

waicTab <- waicTab[order(waicTab$peak, waicTab$waic),]

waicTab$alarmFit <- factor(waicTab$alarmFit,
                           levels = c('thresh', 'hill', 'power',
                                      'spline', 'gp',
                                      'betatSpline',  'basic'),
                           labels = c('Threshold', 'Hill', 'Power',
                                      'Spline', 
                                      'Gaussian Process',
                                      '$\\beta_t$', 'No Behavior Change'))

waicTab$smoothWindow <- factor(waicTab$smoothWindow,
                               levels = c(1, 60, 30),
                               labels = c('None', '60-day', '30-day'))

waicTab$peak <- factor(waicTab$peak,
                       labels = paste0('Wave ', 1:2))
waicTab$waic <- sprintf("%.2f", round(waicTab$waic, 2))

kable(waicTab, row.names = F, format = 'latex', align = 'lccc', 
      booktabs = T, escape = F, 
      col.names = linebreak(c('\\textbf{Wave}',
                              '\\textbf{Model fitted}', 
                              '\\textbf{Smoothing}', 
                              '\\textbf{WAIC}'), align = 'c')) %>% 
    collapse_rows(columns = 1, latex_hline = 'major') 


################################################################################
# Figure 6: Alarms 
################################################################################

myTheme <- theme_bw() + 
    theme(plot.title = element_text(h = 0.5),
          strip.background = element_rect(color = 'white',fill = 'white'),
          strip.text = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 8),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())


pal <- c('darkorchid2', 'darkorange2')

p1 <- ggplot(subset(alarmSub, peak == 1),
             aes(x = xAlarm, y = mean, ymin=lower, ymax=upper)) +  
    geom_line(col = pal[1]) +
    geom_ribbon(alpha=0.3, fill = pal[1]) +
    facet_grid(~alarmFit) +
    labs(x = '60-day incidence', y = 'Alarm') + 
    ylim(0, 0.7) +
    labs(title = 'Wave 1') +
    myTheme 

p2 <- ggplot(subset(alarmSub, peak == 2),
             aes(x = xAlarm, y = mean, ymin=lower, ymax=upper)) +  
    geom_line(col = pal[2]) +
    geom_ribbon(alpha=0.3, fill = pal[2]) +
    facet_grid(~alarmFit) +
    labs(x = '30-day incidence', y = 'Alarm') + 
    ylim(0, 0.7) +
    labs(title = 'Wave 2') +
    myTheme 

jpeg('./figures/fig6_nyc_alarms.jpg', units = 'in', res = 500, width = 8, height = 5)
grid.arrange(p1, p2, nrow = 2)
dev.off()

################################################################################
# Figure 7 - Posterior prediction and R0
################################################################################

### Merge posterior predictive and R0

r0Peak1$output <- 'R_0'
r0Peak2$output <- 'R_0'

postPredFitPeak1$output <- 'Posterior Predictive Fit'
postPredFitPeak2$output <- 'Posterior Predictive Fit'

allOut1 <- rbind.data.frame(r0Peak1, postPredFitPeak1)
allOut2 <- rbind.data.frame(r0Peak2, postPredFitPeak2)

allOut <- rbind.data.frame(allOut1, allOut2)

allOut$smoothedCases[allOut$output == 'R_0'] <- NA
allOut$r0Thresh <- 1
allOut$r0Thresh[allOut$output == 'Posterior Predictive Fit'] <- NA

# format for better plotting
allOut$output <- factor(allOut$output, 
                        levels = c('R_0', 'Posterior Predictive Fit'), 
                        labels = c('R[0](t)', "atop('Posterior', 'Predictive Fit')"))

allOut$Peak <- factor(allOut$Peak, 
                      labels = paste0('Wave ', 1:2))



allOut$alarmFitLab<- factor(allOut$alarmFit,
                          levels = c('spline', 'gp',
                                     'thresh', 'hill', 'power',
                                     'betatSpline', 'basic'),
                          labels = c('Spline', "atop('BC Model', 'Gaussian Process')", 
                                     'Threshold', 'Hill', 'Power',
                                     'Fleixble~beta[t]', "atop('No Behavioral', 'Change')"))


pal <- c('darkorchid2',  'darkorange2')

scale_r01 <- scale_y_continuous(limits = c(0.8, 3.8))
scale_r02 <- scale_y_continuous(limits = c(0.45, 1.2))
scale_pp1 <- scale_y_continuous(limits = c(0, 6000), labels = scales::comma)
scale_pp1_b <- scale_y_continuous(limits = c(0, 50000), labels = scales::comma)
scale_pp2 <- scale_y_continuous(limits = c(0, 5500), labels = scales::comma)

scales <- list(
    scale_r01, scale_pp1, scale_r02, scale_pp2, 
    scale_r01, scale_pp1_b, scale_r02, scale_pp2, 
    scale_r01, scale_pp1, scale_r02, scale_pp2
)

my_strips <- strip_nested(
    # Horizontal strips
    text_x = elem_list_text(size = c(28, 22)),
    by_layer_x = TRUE,
    # Vertical strips
    text_y = elem_list_text(size = 22),
    by_layer_y = FALSE
)

jpeg('./figures/fig7_nyc_r0_postPred.jpg', units = 'in', res = 500, height = 11, width = 17)
ggplot(subset(allOut, 
                  alarmFit %in% c('gp', 'basic', 'betatSpline')), 
       aes(x = date, y = mean, ymin=lower, ymax=upper)) +  
    geom_line(aes(y = smoothedCases), color = 'black', linewidth = 0.5) +
    geom_line(aes(col = Peak), size = 0.8) +
    geom_ribbon(aes(fill = Peak), alpha=0.3) +
    geom_line(aes(y = r0Thresh), color = 'black', linewidth = 0.5, linetype = 2) +
    facet_nested(alarmFitLab~Peak + output, scales = 'free', independent = 'y',
                        labeller =  labeller(Peak = label_value, 
                                             alarmFitLab = label_parsed, 
                                             output = label_parsed),
                 nest_line = element_line(linetype = 1),
                 switch = 'y',
                 strip = my_strips) +
    facetted_pos_scales(y = scales) +
    labs(x = 'Date', y = '') + 
    scale_color_manual(values = pal) + 
    scale_fill_manual(values = pal) + 
    scale_x_date(date_labels = "%b '%y") +
    guides(col = 'none', fill = 'none') +
    theme_bw() + 
    theme(strip.background = element_rect(color = 'white',fill = 'white'),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),
          axis.title = element_text(size = 22),
          strip.placement = 'outside',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
dev.off()

