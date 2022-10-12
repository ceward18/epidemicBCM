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

# functions to calculate the alarms
source('./scripts/modelCodes.R')

# load data
### read data
dat <- read.csv('./Data/nycClean.csv')
dat$date <- as.Date(dat$date)
dat$smoothedCases <- round(movingAverage(dat$dailyCases, 7))
dat <- dat[dat$date < as.Date('2022-03-15'), ]

################################################################################
# Figure XX: NYC Data
################################################################################

dat$Peak <- factor(dat$peak)

pal <- c('mediumpurple2', 'darkolivegreen3', 'darkorange2', 'steelblue1')

pdf('./figures/fig6_nyc_data.pdf', width = 7, height = 4)
ggplot(dat, aes(x = date, y = smoothedCases)) + 
    geom_line(linetype = 2) + 
    geom_line(data = subset(dat, peak == 1), col = pal[1], lwd = 1.2) + 
    geom_line(data = subset(dat, peak == 2), col = pal[2], lwd = 1.2) + 
    geom_line(data = subset(dat, peak == 3), col = pal[3], lwd = 1.2) + 
    geom_line(data = subset(dat, peak == 4), col = pal[4], lwd = 1.2) +
    scale_y_continuous(labels = scales::comma) +
    scale_x_date(date_breaks = "2 month", date_minor_breaks = "1 month",
                 date_labels = "%b %y") +
    labs(x = 'Date', y = 'Incidene', title = 'NYC COVID-19 Case Counts') +
    annotate('text', x = as.Date('2020-04-15'), y = 8000, label = 'Wave 1') +
    annotate('text', x = as.Date('2021-02-10'), y = 8000, label = 'Wave 2') +
    annotate('text', x = as.Date('2021-09-01'), y = 4000, label = 'Wave 3') +
    annotate('text', x = as.Date('2021-11-15'), y = 30000, label = 'Wave 4') + 
    theme_bw() + 
    theme(strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
          plot.title = element_text(size = 14, h = 0.5))
dev.off()

################################################################################
# Figure XX: All WAIC values for converged models
################################################################################

grAll <- readRDS('./results/grAll.rds')

# which didn't converge
notConverge <- grAll[which(grAll$gr > 1.2),  ]
notConvergeModels <- notConverge[
    !duplicated(notConverge
                [,-which(colnames(notConverge) %in% c('gr', 'grUpper', 'param'))]),
    c('alarmFit', 'smoothWindow', 'prior', 'peak')]
notConvergeModels$noConverge <- 1


### load WAIC values
waicAll <- readRDS('./results/waicAll.rds')

### flag those that did not converge
waicAll <- merge(waicAll, notConvergeModels,
                 by = c('alarmFit', 'smoothWindow', 'prior', 'peak'),
                 all.x = T)
waicAll$noConverge[is.na(waicAll$noConverge)] <- 0
waicAll <- waicAll[order(waicAll$peak, waicAll$waic),]
waicAll <- waicAll[waicAll$alarmFit != 'spline',]


waicAll$alarmFit <- factor(waicAll$alarmFit,
                           levels = c('thresh', 'hill', 'power',
                                      'splineFixKnot', 'gp',
                                      'betatSpline',  'basic'),
                           labels = c('Threshold', 'Hill', 'Power',
                                      'Spline', 
                                      'Gaussian Process',
                                      'Beta[t]', 'No Behavior Change'))


waicAll$smoothWindow <- factor(waicAll$smoothWindow,
                               levels = c(1, 60, 30),
                               labels = c('None', '60-day', '30-day'))

waicAll$waic[waicAll$noConverge == 1] <- NA

# minimum for each peak
# minimum WAIC within each prior/peak combination
minWAIC <- ddply(subset(waicAll,  noConverge == 0 ), 
                 .(peak, prior), summarize,
                 alarmFit = alarmFit[which.min(waic)],
                 smoothWindow = smoothWindow[which.min(waic)])
minWAIC$isMin <- 1

waicPlot <- merge(waicAll, minWAIC, by = c('peak', 'prior',
                                           'alarmFit', 'smoothWindow'),
                  all.x = T)

waicPlot$isMin <- factor(waicPlot$isMin, 
                         labels = c('Minimum WAIC'))

waicPlot$prior <- factor(waicPlot$prior, 
                         levels = c(1,3,2,4, 5),
                         labels=c('Mean 5, strong', 'Mean 5, weaker', 
                                  'Mean 2, strong','Mean 2, weaker',
                                  'Mean 3'))

waicPlot$Peak <- factor(waicPlot$peak,
                        labels = paste0('Wave ', 1:4))


pdf('./figures/supp_fig7_nyc_waic.pdf', width = 12, height = 10.5)
ggplot(waicPlot, 
       aes(x = alarmFit,  y = waic, fill = smoothWindow, col = isMin)) +
    geom_bar(stat = 'identity', position=position_dodge2(padding = 0.1),
             size = 0.7) + 
    facet_grid(prior ~ Peak, scales = 'free' ) + 
    coord_flip() + 
    scale_fill_manual(values = c('grey30', 'grey70', 'slateblue2'),
                      limits = c('None', '30-day', '60-day'))+
    scale_color_manual(values = c('red'), na.translate = FALSE) +
    theme_bw() + 
    theme(strip.background = element_rect(color = 'white',
                                          fill = 'white'),
          strip.text.x = element_text(size = 17),
          strip.text.y = element_text(size = 14),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 8),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          plot.title = element_text(size = 18, h = 0.5)) + 
    guides(color = guide_legend(override.aes = list(fill = "white"),
                                title = element_blank())) + 
    labs(fill = 'Smoothing',
         x = '', y = 'WAIC')
dev.off()


# just prior 5 for main results

pdf('./figures/fig7_nyc_waic.pdf', width = 12, height = 4)
ggplot(subset(waicPlot, prior == 'Mean 3'), 
       aes(x = alarmFit,  y = waic, fill = smoothWindow, col = isMin)) +
    geom_bar(stat = 'identity', position=position_dodge2(padding = 0.1),
             size = 0.7) + 
    facet_grid( ~ Peak, scales = 'free' ) + 
    coord_flip() + 
    scale_fill_manual(values = c('grey30', 'grey70', 'slateblue2'),
                      limits = c('None', '30-day', '60-day'))+
    scale_color_manual(values = c('red'), na.translate = FALSE) +
    theme_bw() + 
    theme(strip.background = element_rect(color = 'white',
                                          fill = 'white'),
          strip.text.x = element_text(size = 17),
          strip.text.y = element_text(size = 14),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 8),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          plot.title = element_text(size = 18, h = 0.5)) + 
    guides(color = guide_legend(override.aes = list(fill = "white"),
                                title = element_blank())) + 
    labs(fill = 'Smoothing',
         title = "WAIC values by wave",
         x = '', y = 'WAIC')
dev.off()

################################################################################
# Tables XX: All posterior mean and 95% CI's
################################################################################

## load posterior estimates of parameters
paramsPostAll <- readRDS('./results/paramsPostAll.rds')

# remove those that did not converge
paramsPostAll <- merge(paramsPostAll, notConvergeModels,
                       by = c('alarmFit', 'smoothWindow', 'prior', 'peak'),
                       all.x = T)
paramsPostAll$noConverge[is.na(paramsPostAll$noConverge)] <- 0


# remove spline initial run (knots estimated)
paramsPostAll <- paramsPostAll[paramsPostAll$alarmFit != 'spline',]

rateIParams <- paramsPostAll[paramsPostAll$param == 'rateI',]
meanIParams <- rateIParams
meanIParams$mean <- 1/meanIParams$mean
meanIParams$lower <- 1/meanIParams$lower
meanIParams$upper <- 1/meanIParams$upper
meanIParams$param <- 'meanI'

paramsPostAll <- rbind.data.frame(paramsPostAll, meanIParams)
paramsPostAll <- paramsPostAll[-which(paramsPostAll$param == 'rateI'),]


paramsPostAll$Peak <- factor(paramsPostAll$peak,
                             labels = paste0('Wave ', 1:4))

paramsPostAll$param <- factor(paramsPostAll$param, 
                              levels = c("beta", "meanI", 
                                         "b[1]", "b[2]", "b[3]", "b[4]", 
                                         "delta", "H", "k",  
                                         "knots[1]", "knots[2]", "knots[3]",
                                         "l", "nu", "sigma", "x0"))

paramsPostAll$mean <- sprintf('%.2f',paramsPostAll$mean)
paramsPostAll$lower <- sprintf('%.2f',paramsPostAll$lower)
paramsPostAll$upper <- sprintf('%.2f',paramsPostAll$upper)

paramsPostAll$mean[paramsPostAll$noConverge == 1] <- '-'
paramsPostAll$lower[paramsPostAll$noConverge == 1] <- '-'
paramsPostAll$upper[paramsPostAll$noConverge == 1] <- '-'


#### Tables instead of figures

postParamTableAlarm <- function(alarmFits, smoothWindows) {
    
    tab <- subset(paramsPostAll,
                  alarmFit == alarmFits & smoothWindow == smoothWindows)
    
    tab$Prior <- factor(tab$prior, 
                        levels = c(1,3,2,4),
                        labels=c('Mean 5 \n Strong', 
                                 'Mean 5 \n Weaker', 
                                 'Mean 2 \n Strong',
                                 'Mean 2 \n Weaker'))
    
    tab$val <- paste0(tab$mean, ' \n (',
                      tab$lower, ', ',
                      tab$upper, ')')
    tab <- tab[,c('Peak', 'Prior',  'param', 'val')]
    tab <- tab[order(tab$Peak, tab$Prior, tab$param),]
    
    tabWide <- reshape(tab, 
                       timevar = "param",
                       idvar = c('Peak', 'Prior'),
                       direction = "wide")
    
    tabWide <- tabWide[order(tabWide$Peak, tabWide$Prior),]
    
    if (alarmFits == 'power') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$\\beta$', '$1/\\gamma$', '$k$') 
    } else if (alarmFits == 'thresh') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$\\beta$', '$1/\\gamma$', '$\\delta$', '$H$') 
    } else if (alarmFits == 'hill') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$\\beta$', '$1/\\gamma$', '$\\delta$', '$\\nu$', '$x_0$') 
    } else if (alarmFits == 'splineFixKnot') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$\\beta$', '$1/\\gamma$', '$b_1$', '$b_2$', '$b_3$') 
    } else if (alarmFits == 'gp') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$\\beta$', '$1/\\gamma$', '$\\ell$', '$\\sigma$') 
    } else if (alarmFits == 'betatSpline') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$1/\\gamma$', 
                                                '$b_1$', '$b_2$', '$b_3$', '$b_4$', 
                                                '$\\text{knot}_1$', '$\\text{knot}_2$', '$\\text{knot}_3$') 
    } else if (alarmFits == 'basic') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$\\beta$', '$1/\\gamma$') 
    } 
    
    tabWide %>%
        mutate_all(linebreak, align = 'c') %>%
        kable(format = 'latex', booktabs = T, 
              row.names = F, escape = F, align = c(rep('l', 2), rep('c', ncol(tabWide) - 2))) %>%
        kable_styling()  %>%
        collapse_rows(columns = 1, valign = "top")
}


postParamTableAlarm(alarmFits = 'power', smoothWindows = 30)
postParamTableAlarm(alarmFits = 'power', smoothWindows = 60)

postParamTableAlarm(alarmFits = 'thresh', smoothWindows = 30)
postParamTableAlarm(alarmFits = 'thresh', smoothWindows = 60)

postParamTableAlarm(alarmFits = 'hill', smoothWindows = 30)
postParamTableAlarm(alarmFits = 'hill', smoothWindows = 60)

postParamTableAlarm(alarmFits = 'splineFixKnot', smoothWindows = 30)
postParamTableAlarm(alarmFits = 'splineFixKnot', smoothWindows = 60)

postParamTableAlarm(alarmFits = 'gp', smoothWindows = 30)
postParamTableAlarm(alarmFits = 'gp', smoothWindows = 60)

# tables for 
postParamTableAlarm(alarmFits = 'basic', smoothWindows = 1)
postParamTableAlarm(alarmFits = 'betatSpline', smoothWindows = 1)


################################################################################
# Supplemental Figures 8-12: Alarm functions by model
################################################################################

### load posterior estimates of alarm functions
alarmAll <- readRDS('./results/alarmPostAll.rds')

### remove those that did not converge (TEMPORARY)
alarmAll <- merge(alarmAll, notConvergeModels,
                  by = c('alarmFit', 'smoothWindow', 'prior', 'peak'),
                  all.x = T)
alarmAll$noConverge[is.na(alarmAll$noConverge)] <- 0


# NA for those that did not converge
alarmAll$mean[alarmAll$noConverge == 1] <- NA
alarmAll$lower[alarmAll$noConverge == 1] <- NA
alarmAll$upper[alarmAll$noConverge == 1] <- NA

# remove spline without fixed knots
alarmAll <- alarmAll[alarmAll$alarmFit != 'spline',]

alarmAll$prior <- factor(alarmAll$prior, 
                         levels = c(1,3,2,4),
                         labels=c('Mean 5, strong', 'Mean 5, weaker', 
                                  'Mean 2, strong', 'Mean 2, weaker'))

alarmAll$Smoothing <- factor(alarmAll$smoothWindow,
                             labels = c( '30-day', '60-day'))

alarmAll$Peak <- paste0('Wave ', alarmAll$peak)

myTheme <- theme_bw() +
    theme(strip.background = element_rect(color = 'white',
                                          fill = 'white'),
          strip.text = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 8),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12))

pal <- c('orangered', 'blue')

pdf('./figures/supp_fig8_nyc_alarm_power.pdf', width = 8, height = 7)
ggplot(subset(alarmAll, alarmFit == 'power'),
       aes(x = xAlarm, y = mean, col = Smoothing, fill = Smoothing)) +  
    geom_line(size = 0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
    facet_grid(prior~Peak, scales = 'free_x') +
    labs(x = 'Smoothed incidence', y = 'Alarm') +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    myTheme
dev.off()


pdf('./figures/supp_fig9_nyc_alarm_thresh.pdf', width = 8, height = 7)
ggplot(subset(alarmAll, alarmFit == 'thresh'),
       aes(x = xAlarm, y = mean, col = Smoothing, fill = Smoothing)) +  
    geom_line(size = 0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
    facet_grid(prior~Peak, scales = 'free_x') +
    labs(x = 'Smoothed incidence', y = 'Alarm') +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    myTheme
dev.off()

pdf('./figures/supp_fig10_nyc_alarm_hill.pdf', width = 8, height = 7)
ggplot(subset(alarmAll, alarmFit == 'hill'),
       aes(x = xAlarm, y = mean, col = Smoothing, fill = Smoothing)) +  
    geom_line(size = 0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
    facet_grid(prior~Peak, scales = 'free_x') +
    labs(x = 'Smoothed incidence', y = 'Alarm') +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    myTheme
dev.off()


pdf('./figures/supp_fig11_nyc_alarm_spline.pdf', width = 8, height = 7)
ggplot(subset(alarmAll, alarmFit == 'splineFixKnot'),
       aes(x = xAlarm, y = mean, col = Smoothing, fill = Smoothing)) +  
    geom_line(size = 0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
    facet_grid(prior~Peak, scales = 'free_x') +
    labs(x = 'Smoothed incidence', y = 'Alarm') +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    myTheme
dev.off()


pdf('./figures/supp_fig12_nyc_alarm_gp.pdf', width = 8, height = 7)
ggplot(subset(alarmAll, alarmFit == 'gp'),
       aes(x = xAlarm, y = mean, col = Smoothing, fill = Smoothing)) +  
    geom_line(size = 0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
    facet_grid(prior~Peak, scales = 'free_x') +
    labs(x = 'Smoothed incidence', y = 'Alarm') +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    myTheme
dev.off()


################################################################################
# Supplemental Figures 8-12: R0(t) by model
################################################################################

################################################################################
### load posterior estimates of R0 over epidemic time
r0All <- readRDS('./results/r0PostAll.rds')

### remove those that did not converge (TEMPORARY)
r0All <- merge(r0All, notConvergeModels,
               by = c('alarmFit', 'smoothWindow', 'prior', 'peak'),
               all.x = T)
r0All$noConverge[is.na(r0All$noConverge)] <- 0


# NA for those that did not converge
r0All$mean[r0All$noConverge == 1] <- NA
r0All$lower[r0All$noConverge == 1] <- NA
r0All$upper[r0All$noConverge == 1] <- NA

# remove spline without fixed knots
r0All <- r0All[r0All$alarmFit != 'spline',]


r0All <- r0All[order(r0All$alarmFit, 
                     r0All$smoothWindow,
                     r0All$peak, 
                     r0All$time),]

r0All$prior <- factor(r0All$prior, 
                         levels = c(1,3,2,4),
                         labels=c('Mean 5, strong', 'Mean 5, weaker', 
                                  'Mean 2, strong', 'Mean 2, weaker'))

r0All$Smoothing <- factor(r0All$smoothWindow,
                             labels = c('None', '30-day', '60-day'))

r0All$Peak <- paste0('Wave ', r0All$peak)

pdf('./figures/supp_figX_nyc_r0_power.pdf', width = 8, height = 7)
ggplot(subset(r0All, alarmFit == 'power'),
       aes(x = time, y = mean, col = Smoothing, fill = Smoothing)) +  
    geom_line(size = 0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
    facet_grid(prior~Peak, scales = 'free_x') +
    labs(x = 'Epidemic Time', y = expression(R[0](t))) +
    geom_hline(yintercept = 1, linetype = 2) + 
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    myTheme
dev.off()

pdf('./figures/supp_figX_nyc_r0_thresh.pdf', width = 8, height = 7)
ggplot(subset(r0All, alarmFit == 'thresh'),
       aes(x = time, y = mean, col = Smoothing, fill = Smoothing)) +  
    geom_line(size = 0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
    facet_grid(prior~Peak, scales = 'free_x') +
    labs(x = 'Epidemic Time', y = expression(R[0](t))) +
    geom_hline(yintercept = 1, linetype = 2) + 
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    myTheme
dev.off()

pdf('./figures/supp_figX_nyc_r0_hill.pdf', width = 8, height = 7)
ggplot(subset(r0All, alarmFit == 'hill'),
       aes(x = time, y = mean, col = Smoothing, fill = Smoothing)) +  
    geom_line(size = 0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
    facet_grid(prior~Peak, scales = 'free_x') +
    labs(x = 'Epidemic Time', y = expression(R[0](t))) +
    geom_hline(yintercept = 1, linetype = 2) + 
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    myTheme
dev.off()

pdf('./figures/supp_figX_nyc_r0_spline.pdf', width = 8, height = 7)
ggplot(subset(r0All, alarmFit == 'splineFixKnot'),
       aes(x = time, y = mean, col = Smoothing, fill = Smoothing)) +  
    geom_line(size = 0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
    facet_grid(prior~Peak, scales = 'free_x') +
    labs(x = 'Epidemic Time', y = expression(R[0](t))) +
    geom_hline(yintercept = 1, linetype = 2) + 
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    myTheme
dev.off()

pdf('./figures/supp_figX_nyc_r0_gp.pdf', width = 8, height = 7)
ggplot(subset(r0All, alarmFit == 'gp'),
       aes(x = time, y = mean, col = Smoothing, fill = Smoothing)) +  
    geom_line(size = 0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
    facet_grid(prior~Peak, scales = 'free_x') +
    labs(x = 'Epidemic Time', y = expression(R[0](t))) +
    geom_hline(yintercept = 1, linetype = 2) + 
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    myTheme
dev.off()
