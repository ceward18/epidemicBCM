################################################################################
# combine batch runs
################################################################################

library(nimble)

outputFiles <- list.files('./Output')

################################################################################

################################################################################
# Gelman rubin

grFiles <- outputFiles[grep('grBatch', outputFiles)]

grAll <- readRDS(paste0('./Output/', grFiles[1]))

for (i in 2:length(grFiles)) {
    gr_i <- readRDS(paste0('./Output/', grFiles[i]))
    if (!'epiSize' %in% colnames(gr_i)) {
        gr_i$epiSize <- 'small'
    }
    
    grAll <-rbind.data.frame(grAll, gr_i)
}


saveRDS(grAll, './resultsFinal/grAll.rds')

################################################################################
# posterior alarms - for models that estimate the alarm function

alarmFiles <- outputFiles[grep('alarmPostBatch', outputFiles)]

alarmAll <- readRDS(paste0('./Output/', alarmFiles[1]))

for (i in 2:length(alarmFiles)) {
    alarm_i <- readRDS(paste0('./Output/', alarmFiles[i]))
    if (!'epiSize' %in% colnames(alarm_i)) {
        alarm_i$epiSize <- 'small'
    }
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
    if (!'epiSize' %in% colnames(postPred_i)) {
        postPred_i$epiSize <- 'small'
    }
    postPredAll <-rbind.data.frame(postPredAll, postPred_i)
}

# remove models that didn't estimate posterior predictions (betat)
postPredAll <- postPredAll[!is.na(postPredAll$mean),]


saveRDS(postPredAll, './resultsFinal/postPredAll.rds')

################################################################################
# R0 posterior over time

# r0PostFiles <- outputFiles[grep('R0Post', outputFiles)]
# 
# r0PostAll <- readRDS(paste0('./', outputFolder, '/', r0PostFiles[1]))
# 
# for (i in 2:length(r0PostFiles)) {
#     r0Post_i <- readRDS(paste0('./', outputFolder, '/', r0PostFiles[i]))
#     r0PostAll <-rbind.data.frame(r0PostAll, r0Post_i)
# }
# 
# rownames(r0PostAll) <- NULL
# 
# saveRDS(r0PostAll, './resultsFinal/r0PostAll.rds')

################################################################################
# posterior parameters 

paramsPostFiles <- outputFiles[grep('paramsPostBatch', outputFiles)]

paramsPostAll <- readRDS(paste0('./Output/', paramsPostFiles[1]))

for (i in 2:length(paramsPostFiles)) {
    paramsPost_i <- readRDS(paste0('./Output/', paramsPostFiles[i]))
    if (!'epiSize' %in% colnames(paramsPost_i)) {
        paramsPost_i$epiSize <- 'small'
    }
    paramsPostAll <-rbind.data.frame(paramsPostAll, paramsPost_i)
}

saveRDS(paramsPostAll, './resultsFinal/paramsPostAll.rds')

################################################################################
# WAIC

waicFiles <- outputFiles[grep('waicPostBatch', outputFiles)]

waicAll <- readRDS(paste0('./Output/', waicFiles[1]))

for (i in 2:length(waicFiles)) {
    waic_i <- readRDS(paste0('./Output/', waicFiles[i]))
    if (!'epiSize' %in% colnames(waic_i)) {
        waic_i$epiSize <- 'small'
    }
    waicAll <-rbind.data.frame(waicAll, waic_i)
}

saveRDS(waicAll, './resultsFinal/waicAll.rds')


################################################################################
# beta posterior (only for betat model)

betaPostFiles <- outputFiles[grep('betaPostBatch', outputFiles)]

betaPostAll <- readRDS(paste0('./Output/', betaPostFiles[1]))

for (i in 2:length(betaPostFiles)) {
    betaPost_i <- readRDS(paste0('./Output/', betaPostFiles[i]))
    if (!'epiSize' %in% colnames(betaPost_i)) {
        betaPost_i$epiSize <- 'small'
    }
    betaPostAll <-rbind.data.frame(betaPostAll, betaPost_i)
}

# remove models that didn't estimate posterior predictions (all but betat)
betaPostAll <- betaPostAll[!is.na(betaPostAll$mean),]
rownames(betaPostAll) <- NULL

saveRDS(betaPostAll, './resultsFinal/betaPostAll.rds')

