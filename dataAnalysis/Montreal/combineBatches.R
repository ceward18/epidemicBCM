################################################################################
# combine batch runs
################################################################################

library(nimble)

outputFiles <- list.files('./Output')
outputFiles <- outputFiles[grep('Batch', outputFiles)]

################################################################################

################################################################################
# Gelman rubin

grFiles <- outputFiles[grep('gr', outputFiles)]

grAll <- readRDS(paste0('./Output/', grFiles[1]))

for (i in 2:length(grFiles)) {
    gr_i <- readRDS(paste0('./Output/', grFiles[i]))
    grAll <-rbind.data.frame(grAll, gr_i)
}


saveRDS(grAll, './resultsFinal/grAll.rds')

################################################################################
# posterior alarms - for models that estimate the alarm function

alarmFiles <- outputFiles[grep('alarmPost', outputFiles)]

alarmAll <- readRDS(paste0('./Output/', alarmFiles[1]))

for (i in 2:length(alarmFiles)) {
    alarm_i <- readRDS(paste0('./Output/', alarmFiles[i]))
    alarmAll <-rbind.data.frame(alarmAll, alarm_i)
}

# remove models that didn't estimate an alarm (betat and basic)
alarmAll <- alarmAll[!is.na(alarmAll$mean),]
rownames(alarmAll) <- NULL

saveRDS(alarmAll, './resultsFinal/alarmPostAll.rds')

################################################################################
# posterior parameters 

paramsPostFiles <- outputFiles[grep('paramsPost', outputFiles)]

paramsPostAll <- readRDS(paste0('./Output/', paramsPostFiles[1]))

for (i in 2:length(paramsPostFiles)) {
    paramsPost_i <- readRDS(paste0('./Output/', paramsPostFiles[i]))
    paramsPostAll <-rbind.data.frame(paramsPostAll, paramsPost_i)
}

saveRDS(paramsPostAll, './resultsFinal/paramsPostAll.rds')

################################################################################
# WAIC

waicFiles <- outputFiles[grep('waicPost', outputFiles)]

waicAll <- readRDS(paste0('./Output/', waicFiles[1]))

for (i in 2:length(waicFiles)) {
    waic_i <- readRDS(paste0('./Output/', waicFiles[i]))
    waicAll <-rbind.data.frame(waicAll, waic_i)
}

saveRDS(waicAll, './resultsFinal/waicAll.rds')


################################################################################
# beta posterior (only for betat model)

betaPostFiles <- outputFiles[grep('betaPost', outputFiles)]

betaPostAll <- readRDS(paste0('./Output/', betaPostFiles[1]))

for (i in 2:length(betaPostFiles)) {
    betaPost_i <- readRDS(paste0('./Output/', betaPostFiles[i]))
    betaPostAll <-rbind.data.frame(betaPostAll, betaPost_i)
}

# remove models that didn't estimate posterior predictions (all but betat)
betaPostAll <- betaPostAll[!is.na(betaPostAll$mean),]
rownames(betaPostAll) <- NULL

saveRDS(betaPostAll, './resultsFinal/betaPostAll.rds')


################################################################################
# R0 posterior over time


r0PostFiles <- outputFiles[grep('R0Post', outputFiles)]

r0PostAll <- readRDS(paste0('./Output/', r0PostFiles[1]))

for (i in 2:length(r0PostFiles)) {
    r0Post_i <- readRDS(paste0('./Output/', r0PostFiles[i]))
    r0PostAll <-rbind.data.frame(r0PostAll, r0Post_i)
}

rownames(r0PostAll) <- NULL

saveRDS(r0PostAll, './resultsFinal/r0PostAll.rds')

################################################################################
# R0 posterior over xAlarm (not basic or betat)



r0AlarmPostFiles <- outputFiles[grep('R0AlarmPost', outputFiles)]

r0AlarmPostAll <- readRDS(paste0('./Output/', r0AlarmPostFiles[1]))

for (i in 2:length(r0AlarmPostFiles)) {
    r0AlarmPost_i <- readRDS(paste0('./Output/', r0AlarmPostFiles[i]))
    r0AlarmPostAll <-rbind.data.frame(r0AlarmPostAll, r0AlarmPost_i)
}

rownames(r0AlarmPostAll) <- NULL

saveRDS(r0AlarmPostAll, './resultsFinal/r0AlarmPostAll.rds')

