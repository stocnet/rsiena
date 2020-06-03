allEffects <- read.csv('allEffects.csv', as.is=TRUE, comment.char='%')
# add default columns needed internally
allEffects$setting <- rep('', nrow(allEffects))
