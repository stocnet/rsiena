#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: maxlikec.r
# *
# * Description: This module contains the code for simulating the process,
# * communicating with C++. For use with maximum likelihood method, so
# * never conditional or from finite differences, or parallel testing.
# *****************************************************************************/
##@maxlikec siena07 ML Simulation Module
maxlikec <- function(z, x, data=NULL, effects=NULL,
                     returnChains=FALSE, byGroup=FALSE, byWave=FALSE,
					 returnDataFrame=FALSE,
					 returnLoglik=FALSE, onlyLoglik=FALSE)
{
# note: parameter x is not used. Just for consistency with other possibilities for FRAN.
    ## retrieve stored information
    f <- FRANstore()
    callGrid <- z$callGrid
    if (nrow(callGrid) == 1)
    {
		if (byGroup)
		{
			theta <- z$thetaMat[1,]
		}
		else
		{
			theta <- z$theta
		}
        ans <- .Call(C_mlPeriod, PACKAGE=pkgname, z$Deriv, f$pData,
                     f$pModel, f$myeffects, theta,
					 1, 1, z$nrunMH, z$addChainToStore,
                     z$returnDataFrame, z$returnDeps,
                     z$returnChains, returnLoglik, onlyLoglik)
		if (!onlyLoglik)
		{
			ans[[6]] <- list(ans[[6]])
			ans[[7]] <- list(ans[[7]])
			if (byGroup)
			{
				ans[[8]] <- list(ans[[8]])
				ans[[9]] <- list(ans[[9]])
				ans[[10]] <- list(ans[[10]])
				ans[[11]] <- list(ans[[11]])
			}
		}
		else
		{
			ans[[2]] <- list(ans[[2]])
			ans[[3]] <- list(ans[[3]])
			ans[[4]] <- list(ans[[4]])
		}
		if (z$returnDeps)
		{
			sims <- ans[[12]]
		}
		else
		{
			sims <- 'there are no simulated dependent variables'
		}
	}
    else
    {
		## z$int2 is the number of processors if iterating by period, so 1 means
		## we are not. Can only parallelize by period withmaxlike.
        if (z$int2 == 1)
        {
            anss <- apply(cbind(callGrid, 1:nrow(callGrid)),
						  1, doMLModel, z$Deriv, z$thetaMat,
                          z$nrunMH, z$addChainToStore,
                          z$returnDataFrame, z$returnDeps,
                          z$returnChains, byGroup, z$theta, returnLoglik,
						  onlyLoglik)
        }
        else
        {
            use <- 1:(min(nrow(callGrid), z$int2))
            anss <- parRapply(z$cl[use], cbind(callGrid, 1:nrow(callGrid)),
							  doMLModel, z$Deriv, z$thetaMat,
                              z$nrunMH, z$addChainToStore,
                              z$returnDataFrame,  z$returnDeps, z$returnChains, byGroup,
							  z$theta, returnLoglik, onlyLoglik)
        }
        ## reorganize the anss so it looks like the normal one
        ans <- list()
 		if (!onlyLoglik)
		{
			ans[[1]] <- sapply(anss, "[[", 1) ## statistics
			ans[[2]] <- NULL ## scores
			ans[[3]] <- NULL ## seeds
			ans[[4]] <- NULL ## ntim
			ans[[5]] <- NULL # randomseed
			if (z$returnChains)
			{
				fff <- lapply(anss, function(x) x[[6]][[1]])
				fff <- split(fff, callGrid[, 1 ]) ## split by group
				ans[[6]] <- fff
			}
			ans[[7]] <- lapply(anss, "[[", 7) ## derivative
			ans[[8]] <- lapply(anss, "[[", 8)
			ans[[9]] <- lapply(anss, "[[", 9)
			ans[[10]] <- lapply(anss, "[[", 10)
			if (!byGroup)
			{
				ans[[8]] <- Reduce("+",  ans[[8]]) ## accepts
				ans[[9]] <- Reduce("+",  ans[[9]]) ## rejects
				ans[[10]] <- Reduce("+",  ans[[10]]) ## aborts
			}
			ans[[11]] <- sapply(anss, "[[", 11)
			if (z$returnDeps)
			{
				fff <- lapply(anss, function(x) x[[12]])
				sims <- split(fff, callGrid[, 1 ]) ## split by group
			}
			else
			{
				sims <- 'no simulated dependent variables'
			}
		}
		else ##onlyLoglik is always byGroup (sienaBayes)
		{
			ans[[1]] <- sum(sapply(anss, '[[', 1))## loglik
			ans[[2]] <- lapply(anss, "[[", 2)
			ans[[3]] <- lapply(anss, "[[", 3)
			ans[[4]] <- lapply(anss, "[[", 4)
		}
    }

    FRANstore(f)

    if (z$Deriv && !onlyLoglik) ## need to reformat the derivatives
    {
        resp <- reformatDerivs(z, f, ans[[7]])
        dff <- resp[[1]]
        dff2 <- resp[[2]]
    }
    else
    {
        dff <- NULL
        dff2 <- NULL
    }
	if (!onlyLoglik)
	{
		fra <- -t(ans[[1]]) ##note sign change

		list(fra = fra, ntim0 = NULL, feasible = TRUE, OK = TRUE,
			 sims=sims, dff=dff, dff2=dff2,
			 chain = ans[[6]], accepts=ans[[8]],
			 rejects= ans[[9]], aborts=ans[[10]], loglik=ans[[11]])
	}
	else
	{
		list(loglik=ans[[1]], accepts=ans[[2]],
			 rejects= ans[[3]], aborts=ans[[4]])
	}
}

##@doMLModel Maximum likelihood
doMLModel <- function(x, Deriv, thetaMat, nrunMH, addChainToStore,
                      returnDataFrame, returnDeps, returnChains,
					  byGroup, theta, returnLoglik, onlyLoglik)
{
    f <- FRANstore()
	if (byGroup)
	{
		theta <- thetaMat[x[1], ]
	}
	else
	{
		theta <- theta
	}
    .Call(C_mlPeriod, PACKAGE=pkgname, Deriv, f$pData,
          f$pModel, f$myeffects, theta,
          as.integer(x[1]), as.integer(x[2]), nrunMH[x[3]], addChainToStore,
		  returnDataFrame, returnDeps, returnChains,
		  returnLoglik, onlyLoglik)
}

##@reformatDerivs Maximum likelihood move this back to inline or internal function
reformatDerivs <- function(z, f, derivList)
{
    ## current format is a list of vectors of the lower? triangle
    ## of the matrix. Need to be put into a symmetric matrix.
    ## Tricky part is getting the rates in the right place!
    ## Note that we do not yet deal with rate effects other than the basic
    dff <- as(matrix(0, z$pp, z$pp), "dsyMatrix")
    nPeriods <- length(derivList)
	dff2 <- vector("list", nPeriods)
	if (z$byWave)
	{
		tmp <- as(matrix(0, z$pp, z$pp), "dsyMatrix")
		dff2[] <- tmp
	}
    for (period in 1:nPeriods)
    {
        dffraw <- derivList[[period]]
        ## start indexes row/col for first effect for the variable
        start <- 1
        ## rawsub is subscript in the vector
        rawsub <- 1
        ## f$myeffects is a list of data frames per dependent variable
        ## of selected effects.
        for (i in 1:length(f$myeffects))
        {
            dffPeriod <- matrix(0, nrow=z$pp, ncol=z$pp)
            ## rows is the number of effects for this variable
            rows <- nrow(f$myeffects[[i]])
            ## nRates is the number of rate effects for this variable
            ## at present nRates will be the number of periods.
            nRates <- sum(f$myeffects[[i]]$type == 'rate')
            ## nonRates is the number of rows/cols in the objective function
            ## part of the derivative matrix dff.
            nonRates <- rows - nRates
            ## first put the basic rate for this variable in the right place
            dffPeriod[cbind(start:(start + nRates -1),
                            start:(start + nRates - 1))] <-
                                dffraw[rawsub:(rawsub + nRates - 1)]
            ##
            ## now the matrix of objective function effects
            ##
            start <- start + nRates
            ##
            rawsub <- rawsub + nRates
            ##
            dffPeriod[start : (start + nonRates - 1 ),
                      start : (start + nonRates - 1)] <-
                          dffraw[rawsub:(rawsub + nonRates * nonRates - 1)]
            start <- start + nonRates
            rawsub <- rawsub + nonRates * nonRates
            dffPeriod <- dffPeriod + t(dffPeriod)
            diag(dffPeriod) <- diag(dffPeriod) / 2
			if (z$byWave)
			{
				dff2[[period]] <- dff2[[period]] - dffPeriod
			}
            dff <- dff - dffPeriod
        }
    }
    list(dff, dff2)
}

