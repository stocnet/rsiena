##@getTargets Miscellaneous Written for Krista. Use as RSiena:::getTargets
getTargets <- function(data, effects) {
    setup <- sienaSetupForCpp(data, effects,
                              includeBehavior = TRUE,
                              includeBipartite = TRUE,
                              returnActorStatistics = FALSE,
                              returnStaticChangeContributions = FALSE,
                              parallelrun = FALSE)
    ans <- .Call(C_getTargets, PACKAGE=pkgname,
                 setup$pData, setup$pModel, setup$myeffects,
                 parallelrun = FALSE,
                 returnActorStatistics = FALSE,
                 returnStaticChangeContributions = FALSE)
    ans
}

##@actorTargets Calculates the actor statistics at a wave (which cannot be
## the last wave), for effects associated with a specified behavior.  Optionally
## replaces the actual data for the wave for this behaviour with imputedData.
## For use with multiple imputation of behavior data.
actorTargets <- function(data, effects, behaviorName, wave,
                         imputedData = NULL, rate = FALSE)
{
    depvarsNo <- length(data$depvars)
    behNo <- c(1:depvarsNo)[names(data$depvars) == behaviorName]
    beh <- data$depvars[[behNo]][,1,]
    n <- dim(data$depvars[[behNo]])[1]
    waves <- dim(data$depvars[[behNo]])[3]

    effects <- effects[effects$include,]
    effects$setting <- rep("", nrow(effects))
    effects <- effects[effects$name == behaviorName,]
    effects <- effects[!effects$basicRate,]
    noEff <- dim(effects)[1]
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

        # Use sienaSetupForCpp for setup
        setup <- sienaSetupForCpp(data, effects,
                                  includeBehavior = TRUE,
                                  includeBipartite = FALSE,
                                  returnActorStatistics = FALSE,
                                  returnStaticChangeContributions = FALSE,
                                  parallelrun = FALSE)
        ans <- .Call(C_getTargets, PACKAGE=pkgname,
                     setup$pData, setup$pModel, setup$myeffects,
                     parallelrun = FALSE,
                     returnActorStatistics = FALSE,
                     returnStaticChangeContributions = FALSE)
        ans2[j,] <- ans[,w]
    }

    list(effects = effects$effectName, effectType = effects$type, targets = ans2)
}
