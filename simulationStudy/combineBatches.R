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


saveRDS(grAll, './resultsFinal/grAll.rds')

################################################################################
# posterior alarms - for models that estimate the alarm function

alarmFiles <- outputFiles[grep('alarmPostBatch', outputFiles)]

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
# posterior predictions

postPredFiles <- outputFiles[grep('epiPredPostBatch', outputFiles)]

postPredAll <- readRDS(paste0('./Output/', postPredFiles[1]))

for (i in 2:length(postPredFiles)) {
    postPred_i <- readRDS(paste0('./Output/', postPredFiles[i]))
    postPredAll <-rbind.data.frame(postPredAll, postPred_i)
}

# remove models that didn't estimate posterior predictions (betat)
postPredAll <- postPredAll[!is.na(postPredAll$mean),]


saveRDS(postPredAll, './resultsFinal/postPredAll.rds')

################################################################################
# posterior parameters 

paramsPostFiles <- outputFiles[grep('paramsPostBatch', outputFiles)]

paramsPostAll <- readRDS(paste0('./Output/', paramsPostFiles[1]))

for (i in 2:length(paramsPostFiles)) {
    paramsPost_i <- readRDS(paste0('./Output/', paramsPostFiles[i]))
    paramsPostAll <-rbind.data.frame(paramsPostAll, paramsPost_i)
}

saveRDS(paramsPostAll, './resultsFinal/paramsPostAll.rds')

################################################################################
# WAIC

waicFiles <- outputFiles[grep('waicPostBatch', outputFiles)]

waicAll <- readRDS(paste0('./Output/', waicFiles[1]))

for (i in 2:length(waicFiles)) {
    waic_i <- readRDS(paste0('./Output/', waicFiles[i]))
    waicAll <-rbind.data.frame(waicAll, waic_i)
}

saveRDS(waicAll, './resultsFinal/waicAll.rds')


################################################################################
# beta posterior (only for betat model)

betaPostFiles <- outputFiles[grep('betaPostBatch', outputFiles)]

betaPostAll <- readRDS(paste0('./Output/', betaPostFiles[1]))

for (i in 2:length(betaPostFiles)) {
    betaPost_i <- readRDS(paste0('./Output/', betaPostFiles[i]))
    betaPostAll <-rbind.data.frame(betaPostAll, betaPost_i)
}

# remove models that didn't estimate posterior predictions (all but betat)
betaPostAll <- betaPostAll[!is.na(betaPostAll$mean),]
rownames(betaPostAll) <- NULL

saveRDS(betaPostAll, './resultsFinal/betaPostAll.rds')

