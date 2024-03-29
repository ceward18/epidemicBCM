################################################################################
# Ebola data formatting
################################################################################


library(outbreaks)

ebola <- ebola_kikwit_1995

# first recorded infection on March 6
ebola <- ebola[ebola$date >= as.Date('1995-03-06'), ]

ebola <- ebola[c('date', 'onset', 'death')]

write.csv(ebola, './data/ebolaClean.csv', quote = F, row.names = F)

