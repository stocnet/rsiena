##@getActorStatistics. Use as RSiena:::getTheActorStatistics
getTheActorStatistics <- function(algorithm, data, effects) {
    pData <- sienaSetupDataForCpp(data,
                                    includeBehavior = TRUE,
                                    includeBipartite = FALSE)
    effects <- effects[effects$include, ]
    setup <- sienaSetupEffectsForCpp(pData,
                                    data, 
                                    effects)
    sienaSetupModelOptionsForCpp(ans = NULL, algorithm = algorithm, data = data, effects = effects,
                                 pData = pData, pModel = setup$pModel,
                                 profileData = FALSE, parallelrun = TRUE)
    ans <- .Call(C_getTargetActorStatistics, PACKAGE=pkgname,
                 pData, setup$pModel, setup$myeffects,
                 parallelrun = TRUE)
    ans
}

