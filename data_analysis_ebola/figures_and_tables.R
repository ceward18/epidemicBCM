################################################################################
# Figures and tables in paper
# Ebola Analysis
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


################################################################################
### Pull in and set up all results needed for figures/tables
################################################################################

################################################################################
### set up Ebola data

dat <- read.csv('./data/ebolaClean.csv')
dat$date <- as.Date(dat$date)
dat$time <- 1:nrow(dat)

################################################################################
### Gelman-rubin - to assess convergence
grAll <- readRDS('./results/grAll.rds')
any(grAll$gr > 1.1) # all have converged

################################################################################
### WAIC

waicAll <- readRDS('./results/waicAll.rds')
waicAll$modelLab <- c(rep('BC Model', 5), rep('Standard Approach', 3))
waicAll <- waicAll[order(waicAll$waic),]
waicAll <- waicAll[,3:1]

# format for better labels
waicAll$alarmFit <- factor(waicAll$alarmFit,
                           levels = c('thresh', 'hill', 'gp',
                                      'power', 'spline', 
                                      'betatSpline', 'basicInt', 'basic'),
                           labels = c('Threshold', 'Hill', 'Gaussian Process',
                                      'Power', 'Spline', 
                                      '$\\beta_t$', 'Intervention',
                                      'No Behavior Change'))

################################################################################
### Posterior alarms

alarmPostAll <- readRDS('./results/alarmPostAll.rds')

# format labels for plotting
alarmPostAll$alarmFit <- factor(alarmPostAll$alarmFit, 
                                levels = c('power', 'thresh', 'hill',
                                           'spline', 'gp'),
                                labels = c('Power', 'Threshold', 'Hill', 
                                           'Spline', 'Gaussian Process'))




################################################################################
### Posterior R0

R0PostAll <- readRDS('./results/R0PostAll.rds')

# get the actual date
R0PostAll <- merge(R0PostAll, dat, by = 'time', all.x = T)

# compare only to best BC model (threshold)
R0PostAll <- subset(R0PostAll, alarmFit %in% c('basic', 'basicInt',
                                               'thresh', 'betatSpline'))

# format labels for plotting
R0PostAll$alarmFit <- factor(R0PostAll$alarmFit, 
                             levels = c('thresh', 'basicInt', 
                                        'betatSpline', 'basic'),
                             labels = c("atop('BC Model', 'Threshold')", 'Intervention', 
                                        'Flexible~beta[t]',
                                        "atop('No Behavioral', 'Change')"))

################################################################################
### Posterior predictive fit

postPredFitAll <- readRDS('./results/postPredFitAll.rds')

# get the actual date
postPredFitAll <- merge(postPredFitAll, dat, by = 'time', all.x = T)

# compare only to best BC model (threshold)
postPredFitAll <- subset(postPredFitAll, alarmFit %in% c('basic', 'basicInt',
                                                         'thresh', 'betatSpline'))

# convert to cumulative scale
postPredFitAllCum <- postPredFitAll %>%
    group_by(alarmFit) %>%
    mutate(mean = cumsum(mean),
           lower = cumsum(lower),
           upper = cumsum(upper),
           onset = cumsum(onset)) %>%
    as.data.frame()

# format labels for plotting
postPredFitAllCum$alarmFit <- factor(postPredFitAllCum$alarmFit, 
                                     levels = c('thresh', 'basicInt', 
                                                'betatSpline', 'basic'),
                                     labels = c("atop('BC Model', 'Threshold')", 
                                                'Intervention', 
                                                'Flexible~beta[t]',
                                                "atop('No Behavioral', 'Change')"))


################################################################################
# Figure 5: Ebola Data
################################################################################

# convert to long fomat for plotting
datWide <- reshape(dat, 
                   varying = c("onset", "death"), 
                   v.names = "count",
                   timevar = "type", 
                   times = c("onset", "death"), 
                   new.row.names = 1:300,
                   direction = "long")

datWide$type <- factor(datWide$type, 
                       levels = c('onset', 'death'),
                       labels = c('Cases', 'Deaths'))

jpeg('./figures/fig5_ebola_data.jpg', units = 'in', res = 500, width = 8, height = 5)
ggplot(datWide, aes(x=date, y=count)) +
    geom_bar(stat="identity") + 
    facet_wrap(~type) +
    labs(x = 'Date', y = 'Count', fill = '',
         title = 'DRC EVD Case and Death Counts') + 
    theme_bw() + 
    theme(strip.background = element_rect(color = 'white',fill = 'white'),
          strip.text = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust=1),
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 18, h = 0.5),
          strip.placement = 'outside',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
dev.off()


################################################################################
# Table 2: All WAIC values for converged models
################################################################################

### WAIC table

kable(waicAll, row.names = F, format = 'latex', align = 'lcc', 
      booktabs = T, escape = F, format.args = list(big.mark = ","),
      col.names = linebreak(c('\\textbf{Type}',
                              '\\textbf{Model fitted}', 
                              '\\textbf{WAIC}'), align = 'c')) %>% 
    collapse_rows(columns = 1, latex_hline = 'major') 



################################################################################
# Figure 6: Alarms
################################################################################

jpeg('./figures/fig6_ebola_alarms.jpg', units = 'in', res = 500, width = 8, height = 2.5)
ggplot(alarmPostAll,
       aes(x = xAlarm, y = mean, ymin=lower, ymax=upper)) +  
    geom_line() +
    geom_ribbon(alpha=0.3) +
    facet_grid(~alarmFit) +
    labs(x = 'Observed cumulative incidence', y = 'Alarm') + 
    theme_bw() + 
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
dev.off()

################################################################################
# Figure 7: Posterior prediction and R0
################################################################################

R0PostAll$intDate <- as.Date('1995-05-09')


# convert to long format for plotting
postPredFitLong <- reshape(postPredFitAllCum, 
                           varying = c("mean", "onset"), 
                           v.names = "y",
                           timevar = "type", 
                           times = c("mean", "onset"), 
                           new.row.names = 1:1056,
                           direction = "long")
postPredFitLong$type <- factor(postPredFitLong$type,
                               levels = c('onset', 'mean'))

p1 <- ggplot(R0PostAll,
             aes(x = date, y = mean, ymin=lower, ymax=upper)) +  
    geom_vline(aes(xintercept = intDate,  col = factor(intDate)), linewidth = 1) + 
    geom_line(linewidth = 0.8) +
    geom_ribbon(alpha=0.3) +
    geom_hline(yintercept = 1, linetype = 2) +
    facet_grid(~alarmFit, labeller =  labeller(alarmFit = label_parsed)) +
    labs(x = '', y =  bquote(paste('\n', R[0](t))), col = '') + 
    theme_bw() + 
    theme(plot.title = element_text(h = 0.5),
          strip.background = element_rect(color = 'white',fill = 'white'),
          strip.text = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          # axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(unit(c(0, 5.5, 0, 20.5), "points"))) +
    scale_color_manual(values = 'tomato', labels = 'Intervention Date,\nMay 9, 1995')


p2 <- ggplot(postPredFitLong,
             aes(x = date, y = y, ymin=lower, ymax=upper)) +  
    geom_line( aes(linetype = type, linewidth = type)) +
    geom_ribbon(alpha=0.3) +
    facet_grid(~alarmFit, labeller =  labeller(alarmFit = label_parsed)) +
    labs(x = '', y =  'Cumulative Incidence', col = '', 
         linetype = '', linewidth = '') + 
    theme_bw() + 
    theme(plot.title = element_text(h = 0.5),
          strip.background = element_rect(color = 'white',fill = 'white'),
          strip.text = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.key.width = unit(2, "line"),
          plot.margin = margin(unit(c(0, 29, 5.5, 5.5), "points"))) + 
    scale_linetype_manual(values = 1:2, 
                          labels = c( '\nObserved\n', 'Posterior\npredictive\nmean')) +
    scale_linewidth_manual(values = c(0.8, 0.8), 
                           labels = c( '\nObserved\n', 'Posterior\npredictive\nmean'))

jpeg('./figures/fig7_ebola_r0_postPred.jpg', units = 'in', res = 500, width = 9, height = 5)
grid.arrange(p1, p2, nrow = 2)
dev.off()



