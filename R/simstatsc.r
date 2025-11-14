#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: simstatsc.r
# *
# * Description: This module contains the code for simulating the process,
# * communicating with C++. Only subsidiary routines used for maximum likelihood
# *****************************************************************************/
##@simstats0c siena07 Simulation Module
simstats0c <- function(z, x, data=NULL, effects=NULL, fromFiniteDiff=FALSE,
                       returnDeps=FALSE, returnChains=FALSE, returnChangeContributions = FALSE,
                       byWave=FALSE, returnDataFrame=FALSE, returnLoglik=FALSE)
{
    ## retrieve stored information
    f <- FRANstore()
    ## browser()
    ## fix up the interface so can call from outside robmon framework
    if (is.null(z$Deriv))
    {
        z$Deriv <- FALSE
    }
    if (is.null(z$Phase))
    {
        z$Phase <- 1 ### nb be aware
    }
    if (fromFiniteDiff)
    {
        returnDeps <- FALSE
    }
    else
    {
        returnDeps <- z$returnDeps
    }
    if (is.null(z$returnActorStatistics))
	{
		z$returnActorStatistics <- FALSE
	}
	if (is.null(z$returnChangeContributions))
	{
		z$returnChangeContributions <- FALSE
	}
    if (is.null(f$seeds))
    {
        seeds <- NULL
    }
    else
    {
        seeds <- f$seeds
    }
    if (is.null(f$randomseed2))
    {
        randomseed2 <- NULL
    }
    else
    {
        if (fromFiniteDiff)
        {
            randomseed2 <- as.integer(f$storedseed)
        }
        else
        {
            randomseed2 <- as.integer(f$randomseed2)
            f$storedseed <- randomseed2
        }
        ## cat(randomseed2, '\n')
    }
	## The possibility to use snow now has been dropped
	## because RSiena requires R >= 2.15.0
	## and snow is superseded.
	## Therefore useStreams was dropped.
	#if (getRversion() < "2.14.0") ##fake this to repeat old results
	#	##	if (TRUE)
	#{
	#	useStreams <- TRUE
	#}
	#else
	#{
	useStreams <- FALSE
	#}
    ## z$int2 is the number of processors if iterating by period, so 1 means
    ## we are not. Now have removed option to parallelize by period
	ans <- .Call(C_forwardModel, PACKAGE=pkgname, z$Deriv, f$pData, seeds,
				 fromFiniteDiff, f$pModel, f$myeffects, z$theta,
				 randomseed2, returnDeps, z$FinDiff.method,
				 !is.null(z$cl) && useStreams, z$addChainToStore,
				  z$returnChains, returnLoglik, 
				  z$returnActorStatistics, z$returnChangeContributions)
    if (!fromFiniteDiff)
    {
        if (z$FinDiff.method)
            f$seeds <- ans[[3]]
    }
    if (z$Deriv )
    {
        sc <- t(ans[[2]])
    }
    else
    {
        sc <-  NULL
    }
    ntim <- ans[[4]]
    fra <- t(ans[[1]])
    f$randomseed2 <- ans[[5]]#[c(1,4,3,2)]
    FRANstore(f)
    if (returnDeps)
    {
        sims <- ans[[6]]
    }
    else
    {
        sims <- NULL
    }
    if (z$returnChains)
    {
        chain <- ans[[7]]
    }
    else
    {
        chain <- NULL
    }
	if (returnLoglik)
	{
		loglik <- ans[[8]]
	}
	else
	{
		loglik <- NULL
	}
	if(z$returnChangeContributions)
	{	
		changeContributions <- ans[[9]]
	}
	else
	{
		changeContributions <- NULL
	}
	if(z$returnActorStatistics)
	{
		actorStatistics <- ans[[10]]
	}
	else
	{
		actorStatistics <- NULL
	}
    if (returnDeps)
    {
        ## attach the names
        names(sims) <- f$groupNames
        periodNo <- 1
        for (i in 1:length(sims))
        {
            names(sims[[i]]) <- f$depNames
            for (j in 1:length(sims[[i]]))
            {
                periodNos <- periodNo:(periodNo  + length(sims[[i]][[j]]) - 1)
                names(sims[[i]][[j]]) <- periodNos
            }
            periodNo <- periodNos[length(periodNos)] + 2
        }
    }
		 ## browser()
    list(sc = sc, fra = fra, ntim0 = ntim, feasible = TRUE, OK = TRUE,
         sims=sims, f$seeds, chain=chain, loglik=loglik,
		 actorStatistics = actorStatistics, changeContributions = changeContributions)
}

##@clearData siena07 Finalizer to clear Data object in C++
clearData <- function(pData)
{
    .Call(C_deleteData, PACKAGE=pkgname, pData)
}
##@clearModel siena07 Finalizer to clear Model object in C++
clearModel <- function(pModel)
{
    .Call(C_deleteModel, PACKAGE=pkgname, pModel)
}
