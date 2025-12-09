##@getActorStatistics. Use as RSiena:::getTheActorStatistics
getTheActorStatistics <- function(algorithm, data, effects) {
    pData <- sienaSetupDataForCpp(algorithm,
                                    data,
                                    includeBehavior = TRUE,
                                    includeBipartite = FALSE)
    effects <- effects[effects$include, ]
    setup <- sienaSetupEffectsForCpp(pData,
                                    data, 
                                    effects)
    ans <- .Call(C_getTargetActorStatistics, PACKAGE=pkgname,
                 pData, setup$pModel, setup$myeffects,
                 parallelrun = TRUE)
    ans
}

