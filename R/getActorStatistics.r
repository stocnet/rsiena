##@getActorStatistics. Use as RSiena:::getTheActorStatistics
getTheActorStatistics <- function(data, effects) {
    setup <- sienaSetupForCpp(data, effects,
                              includeBehavior = TRUE,
                              includeBipartite = TRUE,
                              returnActorStatistics = TRUE,
                              returnStaticChangeContributions = FALSE,
                              parallelrun = TRUE)
    ans <- .Call(C_getTargets, PACKAGE=pkgname,
                 setup$pData, setup$pModel, setup$myeffects,
                 parallelrun = TRUE,
                 returnActorStatistics = TRUE,
                 returnStaticChangeContributions = FALSE)
    ans
}

