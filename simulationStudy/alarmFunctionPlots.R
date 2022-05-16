
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
lineWidths <- 4
lineLength <- 4

layoutMat <- matrix(1:6, nrow = 2, ncol = 3, byrow = F)

pdf('./simulationStudy/Figures/exampleAlarms.pdf', height = 4.5, width = 11.5)
layout(mat = layoutMat,
       widths = c(1, 1, 1), 
       heights = c(7, 2)) 

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
       cex = cexLab,
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
       cex = cexLab,
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
       cex = cexLab,
       seg.len= lineLength)


dev.off()