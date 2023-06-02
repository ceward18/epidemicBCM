################################################################################
# Helper functions to generate the supplemental material document
################################################################################

################################################################################
# functions to plot chains/alarms/beta[t]

plotParamChains <- function(chains, param) {
    ylims <- c(min(chains[[1]][,param],
                   chains[[2]][,param],
                   chains[[3]][,param]),
               max(chains[[1]][,param],
                   chains[[2]][,param],
                   chains[[3]][,param]))
    
    if (param == 'rateI') {
        paramName <- 'gamma'
    } else {
        paramName <- param
    }
    
    plot(chains[[1]][,param], type = 'l', 
         ylim = ylims, main = paramName, ylab = '')
    lines(chains[[2]][,param], col = 'red')
    lines(chains[[3]][,param], col = 'blue')
}

plotAlarmChains <- function(chains) {
    yAlarm1 <- chains[[1]][,grep('yAlarm', colnames(chains[[1]]))]
    yAlarm2 <- chains[[2]][,grep('yAlarm', colnames(chains[[2]]))]
    yAlarm3 <- chains[[3]][,grep('yAlarm', colnames(chains[[3]]))]
    
    yAlarmPost <- cbind.data.frame(x = rep(1:ncol(yAlarm1), 3),
                                   Chain = rep(1:3, each = ncol(yAlarm1)),
                                   mean = c(colMeans(yAlarm1),
                                            colMeans(yAlarm2),
                                            colMeans(yAlarm3)),
                                   lower = c(apply(yAlarm1, 2, quantile, probs = 0.025),
                                             apply(yAlarm2, 2, quantile, probs = 0.025),
                                             apply(yAlarm3, 2, quantile, probs = 0.025)),
                                   upper = c(apply(yAlarm1, 2, quantile, probs = 0.975),
                                             apply(yAlarm2, 2, quantile, probs = 0.975),
                                             apply(yAlarm3, 2, quantile, probs = 0.975)))
    
    yAlarmPost$Chain <- factor(yAlarmPost$Chain)
    
    ggplot(yAlarmPost, aes(x = x, y = mean,
                           ymin = lower, ymax = upper,
                           group = Chain, col = Chain, fill = Chain)) +
        geom_line(size = 1) + 
        geom_ribbon(alpha = 0.2) +
        theme_bw() +
        ggtitle('Alarm function') + 
        scale_color_manual(values = c('black', 'red', 'blue')) + 
        scale_fill_manual(values = c('black', 'red', 'blue')) +
        labs(x = '',  y = 'Alarm') +
        theme(plot.title = element_text(h = 0.5, size = 14),
              legend.title = element_text(h = 0.5, size = 13),
              legend.text = element_text(h = 0.5, size = 12),
              axis.text.y = element_text(h = 0.5, size = 12),
              axis.title.y = element_text(h = 0.5, size = 13),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank())
    
}

plotPowerChains <- function(chains) {
    
    par(mfrow = c(1,3))
    
    param <- 'beta'
    plotParamChains(chains, param)
    
    param <- 'rateI'
    plotParamChains(chains, param)
    
    param <- 'k'
    plotParamChains(chains, param)
    
    
}

plotThreshChains <- function(chains) {
    
    par(mfrow = c(1,4))
    
    param <- 'beta'
    plotParamChains(chains, param)
    
    param <- 'rateI'
    plotParamChains(chains, param)
    
    param <- 'delta'
    plotParamChains(chains, param)
    
    param <- 'H'
    plotParamChains(chains, param)
}

plotHillChains <- function(chains) {
    
    par(mfrow = c(2,3))
    
    param <- 'beta'
    plotParamChains(chains, param)
    
    param <- 'rateI'
    plotParamChains(chains, param)
    
    param <- 'delta'
    plotParamChains(chains, param)
    
    param <- 'nu'
    plotParamChains(chains, param)
    
    param <- 'x0'
    plotParamChains(chains, param)
    
}

plotSplineChains <- function(chains) {
    
    par(mfrow = c(2,4))
    
    param <- 'beta'
    plotParamChains(chains, param)
    
    param <- 'rateI'
    plotParamChains(chains, param)
    
    param <- 'b[1]'
    plotParamChains(chains, param)
    
    param <- 'b[2]'
    plotParamChains(chains, param)
    
    param <- 'b[3]'
    plotParamChains(chains, param)
    
    param <- 'knots[1]'
    plotParamChains(chains, param)
    
    param <- 'knots[2]'
    plotParamChains(chains, param)
    
}

plotSplineFixedChains <- function(chains) {
    
    par(mfrow = c(2,4))
    
    param <- 'beta'
    plotParamChains(chains, param)
    
    param <- 'rateI'
    plotParamChains(chains, param)
    
    param <- 'b[1]'
    plotParamChains(chains, param)
    
    param <- 'b[2]'
    plotParamChains(chains, param)
    
    param <- 'b[3]'
    plotParamChains(chains, param)
    
}

plotGPChains <- function(chains) {
    
    par(mfrow = c(1,4))
    
    param <- 'beta'
    plotParamChains(chains, param)
    
    param <- 'rateI'
    plotParamChains(chains, param)
    
    param <- 'l'
    plotParamChains(chains, param)
    
    param <- 'sigma'
    plotParamChains(chains, param)
    
}

plotBetatChains <- function(chains) {
    
    par(mfrow = c(2,4))
    
    param <- 'rateI'
    plotParamChains(chains, param)
    
    param <- 'b[1]'
    plotParamChains(chains, param)
    
    param <- 'b[2]'
    plotParamChains(chains, param)
    
    param <- 'b[3]'
    plotParamChains(chains, param)
    
    param <- 'b[4]'
    plotParamChains(chains, param)
    
    param <- 'knots[1]'
    plotParamChains(chains, param)
    
    param <- 'knots[2]'
    plotParamChains(chains, param)
    
    param <- 'knots[3]'
    plotParamChains(chains, param)
    
}

plotBetaChains <- function(chains) {
    yAlarm1 <- chains[[1]][,grep('beta', colnames(chains[[1]]))]
    yAlarm2 <- chains[[2]][,grep('beta', colnames(chains[[2]]))]
    yAlarm3 <- chains[[3]][,grep('beta', colnames(chains[[3]]))]
    
    yAlarmPost <- cbind.data.frame(x = rep(1:ncol(yAlarm1), 3),
                                   chain = rep(1:3, each = ncol(yAlarm1)),
                                   mean = c(colMeans(yAlarm1),
                                            colMeans(yAlarm2),
                                            colMeans(yAlarm3)),
                                   lower = c(apply(yAlarm1, 2, quantile, probs = 0.025),
                                             apply(yAlarm2, 2, quantile, probs = 0.025),
                                             apply(yAlarm3, 2, quantile, probs = 0.025)),
                                   upper = c(apply(yAlarm1, 2, quantile, probs = 0.975),
                                             apply(yAlarm2, 2, quantile, probs = 0.975),
                                             apply(yAlarm3, 2, quantile, probs = 0.975)))
    
    yAlarmPost$Chain <- factor(yAlarmPost$chain)
    
    ggplot(yAlarmPost, aes(x = x, y = mean,
                           ymin = lower, ymax = upper,
                           group = Chain, col = Chain, fill = Chain)) +
        geom_line(size = 1) + 
        geom_ribbon(alpha = 0.2) +
        theme_bw() +
        ggtitle(expression(beta[t]~'function')) +
        scale_color_manual(values = c('black', 'red', 'blue')) + 
        scale_fill_manual(values = c('black', 'red', 'blue')) +
        labs(x = 'Epidemic Time',  y = expression(beta[t])) +
        theme(plot.title = element_text(h = 0.5, size = 14),
              legend.title = element_text(h = 0.5, size = 13),
              legend.text = element_text(h = 0.5, size = 12),
              axis.text.y = element_text(h = 0.5, size = 12),
              axis.title.y = element_text(h = 0.5, size = 13),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank())
    
}

################################################################################
# create table of posterior parameters - Ebola

postParamTableEbola <- function(paramsPostAll, alarmFits, scale = FALSE) {
    
    tab <- paramsPostAll
    
    tab$val <- paste0(tab$mean, ' \n (',
                      tab$lower, ', ',
                      tab$upper, ')')
    tab <- tab[,c('alarmFit', 'param', 'val')]
    tab <- tab[order(tab$alarmFit, tab$param),]
    
    tabWide <- reshape(tab, 
                       timevar = "param",
                       idvar = 'alarmFit',
                       direction = "wide")
    
    colnames(tabWide) <- c('Model', 
                           '$\\beta$', '$1/\\lambda$', '$1/\\gamma$', 
                           '$\\beta_1$', '$\\beta_2$',
                           '$b_1$', '$b_2$', '$b_3$', '$b_4$', 
                           '$\\text{knot}_1$', '$\\text{knot}_2$', '$\\text{knot}_3$',
                           '$\\ell$', '$\\sigma$',
                           '$\\delta$', '$\\nu$', '$x_0$',
                           '$k$', '$H$')
    
    tabWide <- tabWide[tabWide$Model == alarmFits,]
    
    # remove NA columns
    tabWide <- tabWide[, colSums(is.na(tabWide)) < nrow(tabWide)]
    
    if (scale) {
        tabWide[,-1] %>%
            mutate_all(linebreak, align = 'c') %>%
            kable(format = 'latex', booktabs = T, 
                  row.names = F, escape = F, 
                  align = rep('c', ncol(tabWide))) %>%
            kable_styling(latex_options = c("HOLD_position", "scale_down"))   
    } else {
        tabWide[,-1] %>%
            mutate_all(linebreak, align = 'c') %>%
            kable(format = 'latex', booktabs = T, 
                  row.names = F, escape = F, 
                  align = rep('c', ncol(tabWide))) %>%
            kable_styling(latex_options = c("HOLD_position"))  
    }
   
}


################################################################################
# create table of posterior parameters - COVID-19

postParamTableAlarm <- function(paramsPostAll, alarmFits, smoothWindows) {
    
    tab <- subset(paramsPostAll,
                  alarmFit == alarmFits & smoothWindow == smoothWindows)
    
    tab$Prior <- factor(tab$prior, 
                        levels = c(5,2,4,1,3),
                        labels=c('Mean 3', 
                                 'Mean 2, strong','Mean 2, weaker', 
                                 'Mean 4, strong', 'Mean 4, weaker'))
    
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
        colnames(tabWide)[3:ncol(tabWide)] <- c('$S_0$', '$I_0$', '$\\beta$', '$1/\\gamma$', '$k$') 
    } else if (alarmFits == 'thresh') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$S_0$', '$I_0$', '$\\beta$', '$1/\\gamma$', '$\\delta$', '$H$') 
    } else if (alarmFits == 'hill') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$S_0$', '$I_0$', '$\\beta$', '$1/\\gamma$', '$\\delta$', '$\\nu$', '$x_0$') 
    } else if (alarmFits == 'splineFixKnot') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$S_0$', '$I_0$', '$\\beta$', '$1/\\gamma$', '$b_1$', '$b_2$', '$b_3$') 
    } else if (alarmFits == 'spline') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$S_0$', '$I_0$', '$\\beta$', '$1/\\gamma$', 
                                                '$b_1$', '$b_2$', '$b_3$', '$\\text{knot}_1$', '$\\text{knot}_2$') 
    } else if (alarmFits == 'gp') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$S_0$', '$I_0$', '$\\beta$', '$1/\\gamma$', '$\\ell$', '$\\sigma$') 
    } else if (alarmFits == 'betatSpline') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$S_0$', '$I_0$', '$1/\\gamma$', 
                                                '$b_1$', '$b_2$', '$b_3$', '$b_4$', 
                                                '$\\text{knot}_1$', '$\\text{knot}_2$', '$\\text{knot}_3$') 
    } else if (alarmFits == 'basic') {
        colnames(tabWide)[3:ncol(tabWide)] <- c('$S_0$', '$I_0$', '$\\beta$', '$1/\\gamma$') 
    } 
    
    tabWide %>%
        mutate_all(linebreak, align = 'c') %>%
        kable(format = 'latex', booktabs = T, 
              row.names = F, escape = F, 
              align = c(rep('l', 2), rep('c', ncol(tabWide) - 2))) %>%
        kable_styling(latex_options = c("HOLD_position", "scale_down"))  %>%
        collapse_rows(columns = 1, valign = "top")
}


