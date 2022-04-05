################################################################################
# combine batch runs
################################################################################

library(nimble)
source('./scripts/modelCodes.R')

outputFiles <- list.files('./Output')

################################################################################

################################################################################
# Gelman rubin

grFiles <- outputFiles[grep('grBatch', outputFiles)]

grAll <- readRDS(paste0('./Output/', grFiles[1]))

for (i in 2:length(grFiles)) {
    gr_i <- readRDS(paste0('./Output/', grFiles[i]))
    grAll <-rbind.data.frame(grAll, gr_i)
}

# which didn't converge
notConverge <- grAll[which(grAll$gr > 1.1),  ]
notConvergeModels <-  notConverge[!duplicated(notConverge),
                                  c('alarmBase', 'alarmGen', 'alarmFit',
                                    'infPeriod', 'smoothWindow', 'simNumber')]
notConvergeModels$noConverge <- 1


saveRDS(grAll, './resultsFinal/grAll.rds')

################################################################################
# posterior alarms - for models that estimate the alarm function

alarmFiles <- outputFiles[grep('alarmPostBatch', outputFiles)]

alarmAll <- readRDS(paste0('./Output/', alarmFiles[1]))

for (i in 2:length(alarmFiles)) {
    alarm_i <- readRDS(paste0('./Output/', alarmFiles[i]))
    alarmAll <-rbind.data.frame(alarmAll, alarm_i)
}


################################################################################
# posterior alarms - for models that estimate the alarm function

paramsPostFiles <- outputFiles[grep('paramsPostBatch', outputFiles)]

paramsPostAll <- readRDS(paste0('./Output/', paramsPostFiles[1]))

for (i in 2:length(paramsPostFiles)) {
    paramsPost_i <- readRDS(paste0('./Output/', paramsPostFiles[i]))
    paramsPostAll <-rbind.data.frame(paramsPostAll, paramsPost_i)
}

paramsPostAll$param <- rownames(paramsPostAll)