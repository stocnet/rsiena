##@getTargets Miscellaneous Written for Krista. Use as RSiena:::getTargets
getTargets <- function(data, effects) {
        pData <- sienaSetupDataForCpp(NULL,
                                      data,
                                      includeBehavior = TRUE,
                                      includeBipartite = FALSE)
        effects <- effects[effects$include, ]
        setup <- sienaSetupEffectsForCpp(pData,
                                       data, 
                                       effects)
        ans <- .Call(C_getTargets, PACKAGE=pkgname,
                     pData, setup$pModel, setup$myeffects,
                     parallelrun = FALSE,
                     returnActorStatistics = FALSE,
                     returnStaticChangeContributions = FALSE)

    ans
}

##@actorTargets Calculates the actor statistics at a wave (which cannot be
## the last wave), for effects associated with a specified behavior.  Optionally
## replaces the actual data for the wave for this behaviour with imputedData.
## For use with multiple imputation of behavior data.
actorTargets <- function(data, effects, algorithm, behaviorName, wave,
                         imputedData = NULL, rate = FALSE) {
    depvarsNo <- length(data$depvars)
    behNo <- c(1:depvarsNo)[names(data$depvars) == behaviorName]
    beh <- data$depvars[[behNo]][,1,]
    n <- dim(data$depvars[[behNo]])[1]
    waves <- dim(data$depvars[[behNo]])[3]
    effects_filtered <- effects[effects$include,]
    effects_filtered <- effects_filtered[effects_filtered$name == behaviorName,]
    effects_filtered <- effects_filtered[!effects_filtered$basicRate,]
    noEff <- dim(effects_filtered)[1]
    ans2 <- array(NA, dim = c(n, noEff))
    w <- wave
    vars <- data$depvars
    if (!is.null(imputedData)) {
        vars[[behNo]][,,w] <- imputedData
        beh[,w] <- imputedData
    }

    for (j in 1:n) {
        if (rate) {
            data$depvars[[behNo]][,1,] <- beh
            data$depvars[[behNo]][,1, setdiff(1:waves, w)] <- NA
            data$depvars[[behNo]][j,1,(w+1):waves] <- rep(1, waves-w)
            data$depvars[[behNo]][j,1,1:w] <- rep(0, w)
        } else {
            for (k in 1:depvarsNo) {
                data$depvars[[k]][,,1] <- vars[[k]][,,waves]
                data$depvars[[k]][,,2:waves] <- vars[[k]][,,1:(waves-1)]
            }
            data$depvars[[behNo]][,1, setdiff(1:waves, w+1)] <- NA
            data$depvars[[behNo]][j,1, setdiff(1:waves, w+1)] <- 0
        }

        pData <- sienaSetupDataForCpp(algorithm,
                                      data,
                                      includeBehavior = TRUE,
                                      includeBipartite = FALSE)
        setup <- sienaSetupEffectsForCpp(pData,
                                       data, effects_filtered)
        ans <- .Call(C_getTargets, PACKAGE=pkgname,
                     pData, setup$pModel, setup$myeffects,
                     parallelrun = FALSE,
                     returnActorStatistics = FALSE,
                     returnStaticChangeContributions = FALSE)
        ans2[j,] <- ans[,w]
    }

    list(effects = effects_filtered$effectName, 
        effectType = effects_filtered$type, 
        targets = ans2)
    # why not returned a named list?
}
