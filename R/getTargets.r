##@getTargets Miscellaneous Written for Krista. Use as RSiena:::getTargets
getTargets <- function(data, effects)
{
	f <- unpackData(data)
	effects <- effects[effects$include,]
	##
	pData <- .Call(C_setupData, PACKAGE=pkgname,
		list(as.integer(f$observations)),
		list(f$nodeSets))
	## register a finalizer
	ans <- reg.finalizer(pData, clearData, onexit = FALSE)
	ans<- .Call(C_OneMode, PACKAGE=pkgname,
		pData, list(f$nets))
	ans<- .Call(C_Behavior, PACKAGE=pkgname, pData,
		list(f$behavs))
	ans<-.Call(C_ConstantCovariates, PACKAGE=pkgname,
		pData, list(f$cCovars))
	ans<-.Call(C_ChangingCovariates,PACKAGE=pkgname,
		pData,list(f$vCovars))
	ans<-.Call(C_DyadicCovariates,PACKAGE=pkgname,
		pData,list(f$dycCovars))
	ans<-.Call(C_ChangingDyadicCovariates,PACKAGE=pkgname,
		pData, list(f$dyvCovars))
	storage.mode(effects$parm) <- 'integer'
	storage.mode(effects$group) <- 'integer'
	storage.mode(effects$period) <- 'integer'
	effects$effectPtr <- NA
	myeffects <- split(effects, effects$name)
	ans<- .Call(C_effects, PACKAGE=pkgname,
		pData, myeffects)
	pModel <- ans[[1]][[1]]
	for (i in 1:length(ans[[2]])) ## ans[[2]] is a list of lists of
		## pointers to effects. Each list corresponds to one
		## dependent variable
	{
		effectPtr <- ans[[2]][[i]]
		myeffects[[i]]$effectPtr <- effectPtr
	}
	ans <- .Call(C_getTargets, PACKAGE=pkgname,
		pData, pModel, myeffects,
		NULL, returnActorStatistics=FALSE,
		returnStaticChangeContributions=FALSE)
	ans
}

##@actorTargets Calculates the actor statistics at a wave (which cannot be
## the last wave), for effects associated with a specified behavior.  Optionally
## replaces the actual data for the wave for this behaviour with imputedData.
## For use with multiple imputation of behavior data.
actorTargets <- function(data, effects, x, behaviorName, wave,
			 imputedData = NULL,
			 rate = FALSE)
{
    depvarsNo <- length(data$depvars)
    behNo <- c(1:depvarsNo)[names(data$depvars) == behaviorName]
    beh <- data$depvars[[behNo]][,1,]
    n <- dim(data$depvars[[behNo]])[1]
    waves <- dim(data$depvars[[behNo]])[3]

    effects <- effects[effects$include,]
    effects$setting <- rep("", nrow(effects))
    effects <- effects[effects$name==behaviorName,]
    effects <- effects[!effects$basicRate,]
    noEff <- dim(effects)[1]
    ans2 <- array(rep(NA, n*noEff), dim=c(n,noEff))

    w <- wave

    vars <- data$depvars
    if (!is.null(imputedData))
    {
        vars[[behNo]][,,w] <- imputedData
        beh[,w] <- imputedData
    }

    for (j in 1:n)
    {
        if (rate)
        {
            data$depvars[[behNo]][,1,] <- beh
            data$depvars[[behNo]][,1,c(1:waves)[-w]] <- NA
            data$depvars[[behNo]][j,1,(w+1):waves] <- rep(1,waves-w)
            data$depvars[[behNo]][j,1,1:w] <- rep(0,w)

        } else {

            for (k in 1:depvarsNo)
            {
                data$depvars[[k]][,,1] <- vars[[k]][,,waves]
                data$depvars[[k]][,,2:waves] <- vars[[k]][,,1:(waves-1)]
            }
            data$depvars[[behNo]][,1,c(1:waves)[-(w+1)]] <- NA
            data$depvars[[behNo]][j,1,c(1:waves)[-(w+1)]] <-
            rep(0,waves-1)
        }

        f <- unpackData(data,x)

        pData <- .Call(C_setupData, PACKAGE=pkgname,
        list(as.integer(f$observations)),
        list(f$nodeSets))
        ## register a finalizer
        ans <- reg.finalizer(pData, clearData, onexit = FALSE)
        ans<- .Call(C_OneMode, PACKAGE=pkgname,
        pData, list(f$nets))
        ans<- .Call(C_Behavior, PACKAGE=pkgname, pData,
        list(f$behavs))
        ans<-.Call(C_ConstantCovariates, PACKAGE=pkgname,
        pData, list(f$cCovars))
        ans<-.Call(C_ChangingCovariates,PACKAGE=pkgname,
        pData,list(f$vCovars))
        ans<-.Call(C_DyadicCovariates,PACKAGE=pkgname,
        pData,list(f$dycCovars))
        ans<-.Call(C_ChangingDyadicCovariates,PACKAGE=pkgname,
        pData, list(f$dyvCovars))
        storage.mode(effects$parm) <- 'integer'
        storage.mode(effects$group) <- 'integer'
        storage.mode(effects$period) <- 'integer'
        effects$effectPtr <- NA
        myeffects <- split(effects, effects$name)
        ans<- .Call(C_effects, PACKAGE=pkgname,
        pData, myeffects)
        pModel <- ans[[1]][[1]]
        for (i in 1:length(ans[[2]])) ## ans[[2]] is a list of lists of
        ## pointers to effects. Each list corresponds to one
        ## dependent variable
        {
            effectPtr <- ans[[2]][[i]]
            myeffects[[i]]$effectPtr <- effectPtr
        }
        ans <- .Call(C_getTargets, PACKAGE=pkgname,
        pData, pModel, myeffects, NULL, returnActorStatistics=FALSE,
		returnStaticChangeContributions=FALSE)
        ans2[j,] <- ans[,w]

    }

    list(effects=effects$effectName,effectType=effects$type, targets=ans2)
}
