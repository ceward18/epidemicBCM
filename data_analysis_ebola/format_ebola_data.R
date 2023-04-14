################################################################################
# Ebola data formatting
################################################################################


library(outbreaks)

ebola <- ebola_kikwit_1995

# data collection starts on March 1
ebola <- ebola[ebola$date >= as.Date('1995-03-01'), ]
ebola <- ebola[c('date', 'onset', 'death')]

write.csv(ebola, './data/ebolaClean.csv', quote = F, row.names = F)

