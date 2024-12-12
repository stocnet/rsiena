#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: sienaDataCreate.r
# *
# * Description: This module contains the code to create
# * Siena data object and group data objects.
# *****************************************************************************/
##@addAttributes DataCreate method for attaching attributes to objects
# This is used when creating the data object.
addAttributes <- function(x, name, ...) UseMethod("addAttributes")
# Note: this is used when creating the data object;
# the attributes are not part of the variables.

##@lowIntegers utility function DataCreate
lowIntegers <- function(vals, centr){
	is.wholenumber <- function(x, tol = 1e-6){abs(x - round(x)) < tol}
	all(is.wholenumber(vals), na.rm=TRUE) && (min(vals, na.rm=TRUE) >= 0) &&
					(max(vals, na.rm=TRUE) <= 20)  && (!centr)
}


##@addAttributes.coCovar DataCreate
addAttributes.coCovar <- function(x, name, ...)
{
	storage.mode(x) <- 'double'
	varmean <- mean(x, na.rm=TRUE)
	range2 <- range(x, na.rm=TRUE)
	attr(x, 'moreThan2') <- length(table(x)) > 2
	vartotal <- sum(x, na.rm=TRUE)
	nonMissingCount <- sum(!is.na(x))
	if (attr(x, "centered"))
	{
		x <- x - varmean
	}
	else
	{
		x <- x - 0.0
	}
	attr(x, 'mean') <- varmean
	rr <- rangeAndSimilarity(x, range2)
	if (rr$range[2] == rr$range[1] && !any(is.na(x)))
	{
		attr(x, 'poszvar') <- FALSE
	}
	else
	{
		attr(x, 'poszvar') <- TRUE
	}
	attr(x, 'range') <- ifelse((rr$range[2] - rr$range[1]==0), 1, rr$range[2] - rr$range[1])
	storage.mode(attr(x, 'range')) <- 'double'
	attr(x, 'range2') <- range2
	## attr(x, 'simTotal') <- rr$simTotal
	attr(x, 'simMean') <- rr$simMean
	## attr(x, 'simCnt') <- rr$simCnt
	attr(x, "name") <- name
	attr(x, "vartotal") <- vartotal
	attr(x, "nonMissingCount") <- nonMissingCount
	attr(x, "lowIntegers") <- lowIntegers(x, attr(x, "centered"))
	if ((!is.null(attr(x, "imputationValues"))) && (attr(x, "centered")))
	{
		attr(x, "imputationValues") <- attr(x, "imputationValues") - varmean
	}
	x
}

##@addAttributes.varCovar DataCreate
addAttributes.varCovar <- function(x, name, ...)
{
	tmpmat <- x
	varmean <- mean(x, na.rm=TRUE)
	vartotal <- sum(x, na.rm=TRUE)
	nonMissingCount <- sum(!is.na(x))
	attr(x, "rangep") <- apply(x, 2, range, na.rm=TRUE)
	attr(x, "meanp") <- colMeans(x, na.rm=TRUE)
	cr <- range(x, na.rm=TRUE)
	attr(x, 'range') <- cr[2] - cr[1]
	storage.mode(attr(x, 'range')) <- 'double'
	attr(x, 'mean') <- varmean
	if (attr(x, "centered"))
	{
		x <- x - varmean
	}
	else
	{
		x <- x - 0.0
	}
	rr <- rangeAndSimilarity(tmpmat, cr)
	if (rr$range[2] == rr$range[1] && !any(is.na(tmpmat)))
	{
		attr(x, 'poszvar') <- FALSE
	}
	else
	{
		attr(x, 'poszvar') <- TRUE
	}
	attr(x, 'simMean') <- rr$simMean
	attr(x, 'moreThan2') <- length(unique(x)) > 2
	attr(x, 'name') <- name
	attr(x, "vartotal") <- vartotal
	attr(x, "nonMissingCount") <- nonMissingCount
	attr(x, "lowIntegers") <- lowIntegers(x, attr(x, "centered"))
    if ((!is.null(attr(x, "imputationValues"))) && (attr(x, "centered")))
	{
		attr(x, "imputationValues") <- attr(x, "imputationValues") - varmean
	}
	x
}
##@addAttributes.coDyadCovar DataCreate
addAttributes.coDyadCovar <- function(x, name, bipartite, ...)
{
	sparse <- attr(x, "sparse")
	if (!bipartite) ## remove diagonal for calculation of mean
	{
		if (!sparse)
		{
			diag(x) <- NA
		}
		else
		{
			diag(x[[1]])  <-  NA
		}
	}
	if (sparse)
	{
		nonMissingCount <- sum(!is.na(x[[1]]))
		varmean <- sum(x[[1]], na.rm=TRUE) / nonMissingCount
		## sparse mean is incorrect
		rr <-  range(x[[1]], na.rm=TRUE)
	}
	else
	{
		varmean <- mean(x, na.rm=TRUE)
		rr <-  range(x, na.rm=TRUE)
		nonMissingCount <- sum(!is.na(x))
	}
	attr(x,'mean') <- ifelse(attr(x,'centered'), varmean, 0)
	attr(x,'range') <- rr[2] - rr[1]
	storage.mode(attr(x, 'range')) <- 'double'
	attr(x,'range2') <- rr
	attr(x, 'name') <- name
	attr(x, "nonMissingCount") <- nonMissingCount
	if (!bipartite) #zero the diagonal
	{
		if (sparse)
		{
			diag(x[[1]]) <- 0
		}
		else
		{
			diag(x) <- 0
		}
	}
	x
}
##@addAttributes.varDyadCovar DataCreate
addAttributes.varDyadCovar <- function(x, name, bipartite, ...)
{
	sparse <- attr(x, "sparse")
	vardims <- attr(x, "vardims")
	if (!bipartite) ## remove the diagonal before calculating the mean
	{
		for (obs in 1:vardims[3])
		{
			if (sparse)
			{
				diag(x[[obs]]) <- NA
			}
			else
			{
				diag(x[, , obs]) <- NA
			}
		}
	}
	if (sparse)
	{
		totalValue <- 0
		totalCount <- 0
		meanp <- rep(NA, vardims[3])
		nonMissingCounts <- rep(NA, vardims[3])
		for (obs in 1:vardims[3])
		{
			totalValue <- totalValue + sum(x[[obs]], na.rm=TRUE)
			nonMissingCounts[obs] <- sum(!is.na(x[[obs]]))
			totalCount <- totalCount +	nonMissingCounts[obs]
			meanp[obs] <- sum(x[[obs]], na.rm=TRUE) /
				nonMissingCounts[obs]
	  }
		varmean <- totalValue / totalCount
		rr <- range(sapply(x, range, na.rm=TRUE), na.rm=TRUE)
		attr(x, "meanp") <- meanp
	}
	else
	{
		varmean <- mean(x, na.rm=TRUE)
		rr <-  range(x, na.rm=TRUE)
		attr(x, "meanp") <- colMeans(x, dims=2, na.rm=TRUE)
		nonMissingCounts <- colSums(!is.na(x), dims=2)
	}
	attr(x, "mean") <- ifelse(attr(x, "centered"), varmean, 0)
	attr(x, "range") <- ifelse((rr[2] - rr[1]==0), 1, rr[2] - rr[1])
	storage.mode(attr(x, "range")) <- "double"
	attr(x, "name") <- name
	attr(x, "nonMissingCount") <- nonMissingCounts
	if (!bipartite) ## put diagonal to zero
	{
		for (obs in 1:vardims[3])
		{
			if (!sparse)
			{
				diag(x[, , obs]) <- 0
			}
			else
			{
				diag(x[[obs]]) <- 0
			}
		}
	}
	x
}
##@sienaDataCreate DataCreate
sienaDataCreate<- function(..., nodeSets=NULL, getDocumentation=FALSE)
{
	##@validNodeSet internal sienaDataCreate
	## suppressPackageStartupMessages(require(network))
	validNodeSet <- function(nodeSetName, n)
	{
		sub <- match(nodeSetName, nodeSetNames)
		if (is.na(sub))
		{
			warning(paste("node set",nodeSetName,"not found in",nodeSetNames,'\n'),
					immediate. = TRUE)
		}
		n == length(nodeSets[[sub]])
	}
	if (getDocumentation)
	{
		return(getInternals())
	}
	narg <- nargs()
	## find a set of names for the objects: either the names given in the
	## argument list or the names of the objects in the argument list
	dots <- as.list(substitute(list(...)))[-1] ##first entry is the word 'list'
	if (length(dots) == 0)
	{
		stop('need some objects')
	}
	if (length(dots) == 1)
	{
		ldots <- list(...)
		dotsIsList <- ((is.list(ldots[[1]])) & 
						(! inherits((ldots[[1]]), "sienaDependent")))
# If dotsIsList, it needs to be a list of variables
# The second condition is to rule out the case of a single dependent network
# given as a list of sparse matrices.
		if (dotsIsList) 
		{
			dots <- as.list(substitute(...))[-1]
			narg <- length(ldots)
		}
	}
	else
	{
		dotsIsList <- FALSE
	}
	nm <- names(dots)
	if (is.null(nm))
	{
		fixup <- seq(along=dots)
	}
	else
	{
		fixup <- nm == ''
	}
	dep <- sapply(dots[fixup], function(x) deparse(x)[1])
	if (is.null(nm))
	{
		nm <- dep
	}
	else if (length(dep) > 0)
	{
		nm[fixup] <- dep
	}
	if (!dotsIsList)
	{
		dots <- list(...)
	}
	else
	{
		dots <- (...)
	}
	names(dots) <- nm
	if (any(duplicated(nm)))
	{
		stop('names must be unique')
	}
	## process the inputs: check dimensions,
	## sort out missings and structural zeros and symmetric etc
	## check sizes match the corresponding nodeSets
	observations <- 0
	depvars <- vector('list',narg)
	cCovars <- vector('list',narg)
	vCovars <- vector('list',narg)
	dycCovars <- vector('list',narg)
	dyvCovars <- vector('list',narg)
	compositionChange <- vector('list',narg)
	v1 <- 0; v2 <- 0; v3 <- 0; v4 <- 0; v5 <- 0; v6 <- 0
	for (i in seq(along = dots))
		switch(class(dots[[i]])[1],
			   sienaDependent = {
				   if (attr(dots[[i]],'sparse'))
				   {
					   ##  require(Matrix)
					   netdims <- c(dim(dots[[i]][[1]]), length(dots[[i]]))
				   }
				   else
				   {
					   netdims <- dim(dots[[i]])
				   }
				   if (observations == 0)
				   {
					   observations <- netdims[3]
				   }
				   else if (observations != netdims[3])
				   {
					   stop('differing number of observations')
				   }
				   v1 <- v1 + 1
				   depvars[[v1]] <- dots[[i]]
				   names(depvars)[v1] <- nm[i]
			   },
			   coCovar = {
				   v2 <- v2 + 1
				   cCovars[[v2]] <- dots[[i]]
				   names(cCovars)[v2] <- nm[i]
			   },
			   varCovar = {
				   v3 <- v3 + 1
				   vCovars[[v3]] <- dots[[i]]
				   names(vCovars)[v3] <- nm[i]
			   },
			   coDyadCovar = {
				   v4 <- v4 + 1
				   dycCovars[[v4]] <- dots[[i]]
				   names(dycCovars)[v4] <- nm[i]
			   },
			   varDyadCovar = {
				   v5 <- v5 + 1
				   dyvCovars[[v5]] <- dots[[i]]
				   names(dyvCovars)[v5] <- nm[i]
			   },
			   compositionChange = {
				   v6 <- v6 + 1
				   compositionChange[[v6]] <- dots[[i]]
				   names(compositionChange)[v6] <- nm[i]
			   },
			   stop(paste("invalid object in sienaDataCreate: argument number",
					i, "is of class ", class(dots[[i]]),
					", which is not a valid ... argument."), call.=FALSE)
			   )
	if (v1 == 0)
	{
		stop("need a dependent variable")
	}
	depvars <- depvars[1:v1]
	if (is.null(nodeSets))
	{
		nodeSets <- list(sienaNodeSet(attr(depvars[[1]], "netdims")[1]))
	}
	else
	{
		if (!(inherits(nodeSets, "sienaNodeSet") || inherits(nodeSets[[1]], "sienaNodeSet")))
		{
			stop("nodeSets should be a sienaNodeSet object or a list of such objects")
		}
	}
	nodeSetNames <- sapply(nodeSets,function(x) attr(x,"nodeSetName"))
	names(nodeSets) <- nodeSetNames
	if (v2 == 0)
	{
		cCovars <- list()
	}
	else
	{
		cCovars <- cCovars[1:v2]
	}
	if (v3 == 0)
	{
		vCovars <- list()
	}
	else
	{
		vCovars <- vCovars[1:v3]
	}
	if (v4 == 0)
	{
		dycCovars <- list()
	}
	else
	{
		dycCovars <- dycCovars[1:v4]
	}
	if (v5 == 0)
	{
		dyvCovars <- list()
	}
	else
	{
		dyvCovars <- dyvCovars[1:v5]
	}
	if (v6 == 0)
	{
		compositionChange <- list()
	}
	else
	{
		compositionChange <- compositionChange[1:v6]
	}
	##now can check dimensions and find ranges
	for (i in seq(along = cCovars))
	{
		if (!validNodeSet(attr(cCovars[[i]], 'nodeSet'), length(cCovars[[i]])))
		{
			stop('constant covariate incorrect node set: ', names(cCovars)[i])
		}
		cCovars[[i]] <- addAttributes(cCovars[[i]], names(cCovars)[i])
	}
	for (i in seq(along=vCovars)) ## note that behaviour variables are not here!
	{
		if (observations < 3)
		{
			stop("Changing covariates are not possible with only two waves")
		}
		if (!validNodeSet(attr(vCovars[[i]], 'nodeSet'), nrow(vCovars[[i]])))
			stop('changing covariate incorrect size: ', names(vCovars)[i])
		if (ncol(vCovars[[i]]) < (observations - 1))
			stop('changing covariate not enough columns')
		if (ncol(vCovars[[i]]) != (observations - 1))
		{
			tmpatt <- attributes(vCovars[[i]])
			vCovars[[i]] <- vCovars[[i]][, 1:(observations - 1), drop=FALSE]
			attnames <- names(tmpatt)
			for (att in seq(along=attnames))
			{
				if (!attnames[att] %in% c('dim', 'dimnames'))
				{
					attr(vCovars[[i]], attnames[att]) <- tmpatt[[att]]
				}
			}
		}
		vCovars[[i]] <- addAttributes(vCovars[[i]], names(vCovars)[i])
	}
	for (i in seq(along=dycCovars))
	{
		nattr <- attr(dycCovars[[i]], 'nodeSet')
		bipartite <- attr(dycCovars[[i]], "type") == "bipartite"
		if (attr(dycCovars[[i]], "sparse"))
		{
			thisdycCovar <- dycCovars[[i]][[1]]
		}
		else
		{
			thisdycCovar <- dycCovars[[i]]
		}
		if (!validNodeSet(nattr[1], nrow(thisdycCovar)))
		{
			stop("dyadic covariate incorrect nbr rows ", names(dycCovars)[i])
		}
		if (!validNodeSet(nattr[2], ncol(thisdycCovar)))
		{
			 stop("dyadic covariate incorrect nbr columns ",
				  names(dycCovars)[i])
		 }
		dycCovars[[i]] <- addAttributes(dycCovars[[i]], names(dycCovars)[i],
										bipartite)
	}
	for (i in seq(along=dyvCovars))
	{
		if (observations < 3)
		{
			stop("Changing covariates are not possible with only two waves")
		}
		nattr <- attr(dyvCovars[[i]],'nodeSet')
		sparse <- attr(dyvCovars[[i]], "sparse")
		bipartite <- attr(dyvCovars[[i]], "type") == "bipartite"
		vardims <- attr(dyvCovars[[i]], "vardims")
		if (!validNodeSet(nattr[1], vardims[1]))
		{
			stop('dyadic changing covariate incorrect number of rows ',
				 names(dyvCovars)[i])
		}
		if (!validNodeSet(nattr[2], vardims[2]))
		{
			stop('dyadic changing covariate incorrect number of columns ',
				 names(dyvCovars)[i])
		}
		if (vardims[3] < (observations - 1))
		{
			stop('Dyadic changing covariate not enough observations')
		}
		if (vardims[3] != (observations - 1))
		{
			tmpatt <- attributes(dyvCovars[[i]])
			if (sparse)
			{
				dyvCovars[[i]] <- dyvCovars[[i]][1:(observations - 1)]
			}
			else
			{
				dyvCovars[[i]] <- dyvCovars[[i]][, 1:(observations - 1)]
			}
			attnames <- names(tmpatt)
			for (att in seq(along=attnames))
			{
				if (attnames[att] != "dim")
				{
					attr(dyvCovars[[i]], attnames[att]) <- tmpatt[[att]]
				}
			}
		}
		dyvCovars[[i]] <- addAttributes(dyvCovars[[i]], names(dyvCovars)[i],
										bipartite)
	}
	compnodesets <- sapply(compositionChange, function(x) attr(x, 'nodeSet'))
	if (any(duplicated(compnodesets)))
		stop('Only one composition change allowed for each nodeSet')
	for (i in seq(along = compositionChange))
	{
		thisNodeSet <- attr(compositionChange[[i]], 'nodeSet')
		nodeSetSize <- length(compositionChange[[i]])
		if (!validNodeSet(thisNodeSet, nodeSetSize))
			stop('composition change incorrect size: ',
				 names(compositionChange)[i])
		if (any(sapply(compositionChange[[i]], function(x)
					   any(x < 1.0 | x > observations))))
			stop("invalid times of composition change")
		if (!all(sapply(compositionChange[[i]], length) %% 2 == 0))
			stop(" Each composition change entry must have an ",
				 "even number of digits")
		## generate events and active flags
		activeStart <- matrix(FALSE, nrow=nodeSetSize, ncol=observations)
		action <- matrix(0, nrow=nodeSetSize, ncol=observations)
		events <- vector("list", nodeSetSize * 2 * observations)
		evSubs <- 1
		for (j in 1:nodeSetSize)
		{
			xsubs <- 1
			x <- compositionChange[[i]][[j]]
			repeat
			{
				##process one interval
				##start <- x[xsubs]
				##end <- x[xsubs+1]
				startIndex <- ceiling(x[xsubs])
				endIndex <- trunc(x[xsubs + 1])
			  #	 if (startIndex < observations && startIndex <= activeEndIndex)
			  #	 {
					activeStart[j, startIndex:endIndex] <- TRUE
			  #	 }
				if (x[xsubs] > 1.0)
				{
					period <- trunc(x[xsubs])
					evTime <- x[xsubs] - period
					events[[evSubs]] <- data.frame(event="join",
												   period=period,
												   actor = j, time=evTime)
					evSubs <- evSubs + 1
				}
				if (x[xsubs+1] < observations)
				{
					period <- trunc(x[xsubs+1])
					evTime <- x[xsubs+1] - period
					events[[evSubs]] <- data.frame(event="leave",
												   period=period,
												   actor = j, time=evTime)
					evSubs <- evSubs + 1
				}
				xsubs <- xsubs + 2
				if (xsubs > length(x))
				{
					break
				}
			}
		 #	 cat(j, 'active',activeStart[j,],  x, '\n')
			if (any(!activeStart[j, ]))
			{
				notActive <- which(!activeStart[j, ])
				for (jj in notActive)
				{
					precActive <- jj > 1 && sum(activeStart[j, 1:(jj - 1)]) > 0
					len <- length(activeStart[i,])
					succActive <- (jj < len) &&
					(sum(activeStart[j, (jj + 1):len]) > 0)

					if (!precActive)
					{
						action[j, jj] <- 1
					}
					else if (!succActive)
					{
						action[j, jj] <- 2
					}
					else
					{
						action[j, jj] <- 3
					}
				}
			}
		}
		if (evSubs > 1)
		{
			events <- do.call(rbind, events[1:(evSubs-1)])
		}
		else
		{
			events <- data.frame(event="join", period=1, actor = 1, time=0)[-1,]
		}
		events$event <- factor(events$event, levels=c('join','leave'))
		storage.mode(events$period) <- "integer"
		storage.mode(events$actor) <- "integer"
		attr(compositionChange[[i]], "events") <- events
		attr(compositionChange[[i]], "activeStart") <- activeStart
		attr(compositionChange[[i]], "action") <- action
	}
	## dependent variables. First we sort the list so discrete and then continuous
	## behavior are at the end
	types <- sapply(depvars, function(x)attr(x, "type"))
	depvars <- depvars[c(which(!(types %in% c('behavior', 'continuous'))),
						which(types == 'behavior'), which(types == "continuous"))]
	onemodes <- which(types == "oneMode")
	bipartites <- which(types == "bipartite")
	onemodes.mx <- max(c(onemodes, 0))
	bipartites.mn <- min(c(bipartites, v1+1))
	if (bipartites.mn < onemodes.mx)
	{
		stop("One-mode networks (if any) should be given before bipartite networks (if any).")
	}

	for (i in 1:v1) ## dependent variables
	{
		nattr <- attr(depvars[[i]], 'nodeSet')
		netdims <- attr(depvars[[i]], 'netdims')
		type <- attr(depvars[[i]], 'type')
		sparse <- attr(depvars[[i]], 'sparse')
		myarray <- depvars[[i]]
		if (!validNodeSet(nattr[1], netdims[1]))
			stop('1st net dimension wrong')
		if (type =='bipartite')
			if (!validNodeSet(nattr[2], netdims[2]))
				stop('2nd net dimension wrong')
		attr(depvars[[i]], 'uponly') <- rep(FALSE, observations - 1)
		attr(depvars[[i]], 'downonly') <- rep(FALSE, observations - 1)
		attr(depvars[[i]], 'distance') <- rep(FALSE, observations - 1)
		attr(depvars[[i]], 'vals') <- vector("list", observations)
		attr(depvars[[i]], 'nval') <- rep(NA, observations)
		attr(depvars[[i]], 'noMissing') <- rep(0, observations)
		attr(depvars[[i]], 'noMissingEither') <- rep(0, observations - 1)
		attr(depvars[[i]], 'nonMissingEither') <- rep(0, observations - 1)
		someOnly <- FALSE
		if (type == 'behavior' || type == 'continuous')
		{
			attr(depvars[[i]], 'noMissing') <- FALSE
			attr(depvars[[i]], 'symmetric') <- NA
			for (j in 1:(observations - 1))
			{
				myvector1 <- myarray[, , j]
				myvector2 <- myarray[, , j + 1]
				mydiff <- myvector1 - myvector2
				attr(depvars[[i]], "distance")[j] <- sum(abs(mydiff),
					na.rm=TRUE)
				attr(depvars[[i]], "vals")[[j]] <- table(myvector1,
					useNA="always")
				attr(depvars[[i]], "vals")[[j+1]] <- table(myvector2,
					useNA="always")

				attr(depvars[[i]], "nval")[j] <-  sum(!is.na(myvector1))
				attr(depvars[[i]], "nval")[j + 1] <-  sum(!is.na(myvector2))
				attr(depvars[[i]], 'noMissing')[j] <- sum(is.na(myvector1))
				attr(depvars[[i]], 'noMissing')[j+1] <- sum(is.na(myvector2))
				attr(depvars[[i]], 'noMissingEither')[j] <-
					sum(is.na(myvector2) | is.na(myvector1))
				attr(depvars[[i]], 'nonMissingEither')[j] <-
					sum(!(is.na(myvector2) | is.na(myvector1)))
				if (attr(depvars[[i]], 'allowOnly'))
				{
					if (all(mydiff >= 0, na.rm=TRUE))
					{
						attr(depvars[[i]], 'downonly')[j] <- TRUE
						someOnly <- TRUE
					}
					if (all(mydiff <= 0, na.rm=TRUE))
					{
						attr(depvars[[i]], 'uponly')[j] <- TRUE
						someOnly <- TRUE
					}
				}
			}
			rr <- range(depvars[[i]], na.rm=TRUE)
			if (rr[2] == rr[1] && !any(is.na(depvars[[i]])))
				attr(depvars[[i]], 'poszvar') <- FALSE
			else
				attr(depvars[[i]], 'poszvar') <- TRUE
			crange <- rr[2] - rr[1]
			attr(depvars[[i]],'range') <- crange
			attr(depvars[[i]], "range2") <- rr
			attr(depvars[[i]],'moreThan2') <- length(unique(depvars[[i]])) > 2
			if (type == 'behavior')
			{
			modes <- apply(depvars[[i]][, 1, ], 2, function(z)
					   {
						   taba <- table(round(z))
						   as.numeric(names(which.max(taba)))
					   }
						   )
			}
			else # type == 'continuous', input the medians here; these will be
				 #						 used for missing data imputation
			{
				modes <- apply(depvars[[i]][, 1, ], 2, median, na.rm=TRUE)
			}
			attr(depvars[[i]],'modes') <- modes
			attr(depvars[[i]], 'missing') <- any(is.na(depvars[[i]]))
			tmpmat <- depvars[[i]][, 1, ]
			rr <- rangeAndSimilarity(tmpmat[, - ncol(tmpmat)], rr)
			##	attr(depvars[[i]], 'simTotal') <- rr$simTotal
			##	attr(depvars[[i]], 'simCnt') <- rr$simCnt
			attr(depvars[[i]], 'simMean') <- rr$simMean
            attr(depvars[[i]], 'variance') <- rr$variance
			attr(depvars[[i]], 'structural') <- FALSE
			attr(depvars[[i]], 'balmean') <- NA
			attr(depvars[[i]], 'structmean') <- NA
	   }
		else
		{
			for (j in 1:(observations - 1))
			{
				if (sparse)
				{
					mymat1 <- myarray[[j]]
					mymat2 <- myarray[[j + 1]]
					## remove diagonals if not bipartite
					if (attr(depvars[[i]], "type") != "bipartite")
					{
						diag(mymat1) <- NA
						diag(mymat2) <- NA
					}
					attr(depvars[[i]], 'noMissingEither')[j] <-
						sum(is.na(mymat1) | is.na(mymat2))
					attr(depvars[[i]], 'nonMissingEither')[j] <-
						sum(!(is.na(mymat1) | is.na(mymat2)))
					## remove diagonals if not bipartite
					if (attr(depvars[[i]], "type") != "bipartite")
					{
						attr(depvars[[i]], 'noMissingEither')[j] <-
							attr(depvars[[i]], 'noMissingEither')[j] -
								nrow(mymat1)
					}
					##remove structural values
					x1 <- mymat1@x
					x2 <- mymat2@x
					x1[x1 %in% c(10, 11)] <- NA
					x2[x2 %in% c(10, 11)] <- NA
					mymat1@x <- x1
					mymat2@x <- x2
					mydiff <- mymat2 - mymat1
					attr(depvars[[i]], 'distance')[j] <- sum(mydiff != 0,
															 na.rm = TRUE)
					if (attr(depvars[[i]], 'allowOnly'))
						{
							if (all(mydiff@x >= 0, na.rm=TRUE))
							{
								attr(depvars[[i]], 'uponly')[j] <- TRUE
								someOnly <- TRUE
							}
							if (all(mydiff@x <= 0, na.rm=TRUE))
							{
								attr(depvars[[i]], 'downonly')[j] <- TRUE
								someOnly <- TRUE
							}
						}
				}
				else
				{
					mymat1 <- myarray[, , j]
					mymat2 <- myarray[, , j + 1]
					## remove diagonals if not bipartite
					if (attr(depvars[[i]], "type") != "bipartite")
					{
						diag(mymat1) <- NA
						diag(mymat2) <- NA
					}
					attr(depvars[[i]], 'noMissingEither')[j] <-
						sum(is.na(mymat1) | is.na(mymat2))
					attr(depvars[[i]], 'nonMissingEither')[j] <-
						sum(!(is.na(mymat1) | is.na(mymat2)))
					## remove diagonals if not bipartite as not really missing
					if (attr(depvars[[i]], "type") != "bipartite")
					{
						attr(depvars[[i]], 'noMissingEither')[j] <-
							attr(depvars[[i]], 'noMissingEither')[j] -
								nrow(mymat1)
					}
					##remove structural values
					mymat1[mymat1 %in% c(10,11)] <- NA
					mymat2[mymat2 %in% c(10,11)] <- NA
					mydiff <- mymat2 - mymat1
					attr(depvars[[i]], 'distance')[j] <- sum(mydiff != 0,
															 na.rm = TRUE)
					if (attr(depvars[[i]], 'allowOnly'))
						{
							if (all(mydiff >= 0, na.rm=TRUE))
							{
								attr(depvars[[i]], 'uponly')[j] <- TRUE
								someOnly <- TRUE
							}
							if (all(mydiff <= 0, na.rm=TRUE))
							{
								attr(depvars[[i]], 'downonly')[j] <- TRUE
								someOnly <- TRUE
							}
						}
				}
			}
			if (type == 'oneMode')
			{
				attr(depvars[[i]], 'balmean') <- calcBalmean(depvars[[i]])
				attr(depvars[[i]], 'structmean') <- calcStructmean(depvars[[i]])
				attr(depvars[[i]], 'simMean') <- NA
                attr(depvars[[i]], 'variance') <- NA
				attr(depvars[[i]], 'symmetric') <- TRUE
				attr(depvars[[i]], 'missing') <- FALSE
				attr(depvars[[i]], 'structural') <- FALSE
				maxObsOutDegree <- rep(NA, observations)
				for (j in 1:observations)
				{
					if (sparse)
					{
						mymat <- myarray[[j]]
						diag(mymat) <- NA
					}
					else
					{
						mymat <- myarray[, , j]
						diag(mymat) <- NA
					}
					if (suppressMessages(!isSymmetric(mymat)))
					{
						attr(depvars[[i]], 'symmetric') <- FALSE
					}
					if (sparse)
					{
						if (any(is.na(mymat@x[mymat@i != mymat@j])))
						{
							attr(depvars[[i]], 'missing') <- TRUE
						}
					}
					else if (any(is.na(mymat[row(mymat) != col(mymat)])))
					{
						attr(depvars[[i]], 'missing') <- TRUE
					}
					if (any(!is.na(mymat) & (mymat == 10 | mymat == 11)))
					{
						attr(depvars[[i]], "structural") <- TRUE
					}
					if (sparse)
					{
						nonZeros <- table(mymat@x, useNA="always")
						Zeros <- nrow(mymat) * nrow(mymat - 1) - sum(nonZeros)
						attr(depvars[[i]], "vals")[[j]] <-
							c(table(rep(0, Zeros)), nonZeros)
					}
					else
					{
						attr(depvars[[i]], "vals")[[j]] <-
							table(mymat, useNA="always")
					}
					attr(depvars[[i]], "nval")[j] <-
						sum(!is.na(mymat[row(mymat) != col(mymat)]))
					# Use %% (modulo) to handle structural values:
					maxObsOutDegree[j] <- max(rowSums(mymat %% 10, na.rm=TRUE))
				}
				### need to exclude the structurals here
				if (sparse)
				{
				   vals <- lapply(depvars[[i]], function(x)
								   c(x@x[!(is.na(x@x) |
										   x@x %in% c(10, 11))] , 0))
					attr(depvars[[i]], "range2") <-
						do.call(range, vals)
			   }
				else
				{
					tmp <- depvars[[i]]
					attr(depvars[[i]], "range2") <-
						range(tmp[!(is.na(tmp) | tmp %in% c(10, 11))])
				}
				## average degree etc.
				atts <- attributes(depvars[[i]])
				ones <- sapply(atts$vals, function(x)
						   {
							   if (is.na(x["11"]))
							   {
									if (is.na(x["1"]))
									{
										0
									}
									else
									x["1"]
							   }
							   else
							   {
								   x["1"] + x["11"]
							   }
						   }
							   )
				density <- ones / atts$nval
				degree <- (atts$netdims[1] - 1) * ones / atts$nval
				missings <- 1 - atts$nval/ atts$netdims[1] /
					(atts$netdims[1] - 1)
				noMissing <- atts$netdims[1] * (atts$netdims[1] - 1) -
					atts$nval
				attr(depvars[[i]], "ones") <- ones
				attr(depvars[[i]], "density") <- density
				attr(depvars[[i]], "degree") <- degree
				attr(depvars[[i]], "averageOutDegree") <- mean(degree)
				attr(depvars[[i]], "averageInDegree") <- mean(degree)
				attr(depvars[[i]], "maxObsOutDegree") <- maxObsOutDegree
				attr(depvars[[i]], "missings") <- missings
				attr(depvars[[i]], "noMissing") <- noMissing
		   }
			else #type=='bipartite' not sure what we need here,
				## but include diagonal
			{
				attr(depvars[[i]], 'balmean') <- NA
				attr(depvars[[i]], 'structmean') <- NA
				attr(depvars[[i]], 'simMean') <- NA
                attr(depvars[[i]], 'variance') <- NA
				attr(depvars[[i]], 'symmetric') <- FALSE
				attr(depvars[[i]], 'missing') <- FALSE
				attr(depvars[[i]], 'structural') <- FALSE
				maxObsOutDegree <- rep(NA, observations)
				for (j in 1:observations)
				{
					if (sparse)
					{
						mymat <- myarray[[j]]
					}
					else
					{
						mymat <- myarray[, , j]
					}
					if (any(is.na(mymat)))
					{
						attr(depvars[[i]], 'missing') <- TRUE
					}
					if (any(!is.na(mymat) & (mymat == 10 | mymat == 11)))
					{
						attr(depvars[[i]], "structural") <- TRUE
					}
					if (sparse)
					{
						nonZeros <- table(mymat@x, useNA="always")
						Zeros <- nrow(mymat) * ncol(mymat) - sum(nonZeros)
						attr(depvars[[i]], "vals")[[j]] <-
							c(table(rep(0, Zeros)), nonZeros)
					}
					else
					{
						attr(depvars[[i]], "vals")[[j]] <- table(mymat,
																 useNA="always")
					}
					attr(depvars[[i]], "nval")[j] <- sum(!is.na(mymat))
					# Use %% (modulo) to handle structural values:
					maxObsOutDegree[j] <- max(rowSums(mymat %% 10, na.rm=TRUE))
				}
				### need to exclude the structurals here
				if (sparse)
				{
				   vals <- lapply(depvars[[i]], function(x)
								   c(x@x[!(is.na(x@x) |
										   x@x %in% c(10, 11))] , 0))
					attr(depvars[[i]], "range") <-
						do.call(range, vals)
			   }
				else
				{
					tmp <- depvars[[i]]
					attr(depvars[[i]], "range") <-
						range(tmp[!(is.na(tmp) | tmp %in% c(10, 11))])
				}
				 ## average degree
				atts <- attributes(depvars[[i]])
				ones <- sapply(atts$vals, function(x)
						   {
							   if (is.na(x["11"]))
							   {
								   x["1"]
							   }
							   else
							   {
								   x["1"] + x["11"]
							   }
						   }
							   )
				density <- ones / atts$nval
				degree <- (atts$netdims[2]) * ones / atts$nval
				missings <- 1 - atts$nval/ atts$netdims[1] /
					atts$netdims[2]
				noMissing <- atts$netdims[1] * atts$netdims[2] - atts$nval
				attr(depvars[[i]], "ones") <- ones
				attr(depvars[[i]], "density") <- density
				attr(depvars[[i]], "degree") <- degree
				attr(depvars[[i]], "averageOutDegree") <- mean(degree)
				attr(depvars[[i]], "maxObsOutDegree") <- maxObsOutDegree
				attr(depvars[[i]], "missings") <- missings
				attr(depvars[[i]], "noMissing") <- noMissing
		   }
		}
			
		if (someOnly)
		{
		message("For dependent variable ", names(depvars)[i], ", in some periods,")
		message("there are only increases, or only decreases.")
		message("This will be respected in the simulations. ")
		message("If this is not desired, use allowOnly=FALSE when creating the dependent variable.")
		}
		attr(depvars[[i]], 'name') <- names(depvars)[i]
	}
	## create the object
	z <- NULL
	z$nodeSets <- nodeSets
	z$observations <- observations
	z$depvars <- depvars
	z$cCovars <- cCovars
	z$vCovars <- vCovars
	z$dycCovars <- dycCovars
	z$dyvCovars <- dyvCovars
	z$compositionChange <- compositionChange
	z <- checkConstraints(z)
	z <- covarDist2(z)
	attr(z, "version") <- packageDescription(pkgname, fields = "Version")
	class(z) <- "siena"
	z
}
##@checkConstraints DataCreate
checkConstraints <- function(z)
{
	types <- sapply(z$depvars, function(x)attr(x, "type"))
	symmetrics <- sapply(z$depvars, function(x)attr(x, "symmetric"))
	sparse <- sapply(z$depvars, function(x)attr(x, "sparse"))
	nodeSets <- lapply(z$depvars, function(x)attr(x, "nodeSet"))
	nNets <- length(z$depvars)

	pairsOfNets <- as.matrix(expand.grid(names(z$depvars), names(z$depvars)))
	pairsNames <- paste(pairsOfNets[, 1], pairsOfNets[, 2], sep=",")

	higher <- namedVector(FALSE, pairsNames )
	atLeastOne <- namedVector(FALSE, pairsNames )
	disjoint <- namedVector(FALSE, pairsNames )

	## identify any nets which may relate. These are those that
	## share both node sets and type.
	relates <- data.frame(name=names(z$depvars), type=types,
						  nodeSets=sapply(nodeSets, paste, collapse=","),
						  tn=paste(types, sapply(nodeSets, paste,
						  collapse=",")) , stringsAsFactors=FALSE)

	## just check we have only one row per network
	if (nrow(relates) != nNets)
	{
		stop("Error in checkConstraints")
	}
	## find the ones that do occur more than once
	use <- relates$tn %in% relates$tn[duplicated(relates$tn)]

	## create a working list to be compared, with names for ease of access
	nets <- namedVector(NA, names(z$depvars), listType=TRUE)

	for (net in names(z$depvars)[use])
	{
		if (types[[net]] != "behavior")
		{
			nets[[net]] <- z$depvars[[net]]
		}
	}

	for (i in 1:nrow(pairsOfNets))
	{
		if (pairsOfNets[i, 1] != pairsOfNets[i, 2])
		{
			net1 <- pairsOfNets[i, 1]
			net2 <- pairsOfNets[i, 2]

			type1 <- types[net1]
			type2 <- types[net2]
			nodes1 <- relates[net1, "nodeSets"]
			nodes2 <- relates[net2, "nodeSets"]

			symmetric1 <- symmetrics[net1]
			symmetric2 <- symmetrics[net2]

			if (type1 == type2 && type1 != "behavior" && nodes1 == nodes2
				&& symmetric1 == symmetric2 && type1 != "continuous"
				&& type2 != "continuous")
			{
				higher[i] <- TRUE
				disjoint[i] <- TRUE
				atLeastOne[i] <- TRUE
				depvar1 <- nets[[pairsOfNets[i, 1]]]
				depvar2 <- nets[[pairsOfNets[i, 2]]]
				for (obs in 1:z$observations)
				{
					if (sparse[net1])
					{
						var1 <- depvar1[[obs]]
					}
					else
					{
						var1 <- depvar1[,, obs]
					}
					if (sparse[net2])
					{
						var2 <- depvar2[[obs]]
					}
					else
					{
						var2 <- depvar2[,, obs]
					}
#					var1[var1 %in% c(10, 11)] <- var1[var1 %in% c(10, 11)] - 10
#					var2[var2 %in% c(10, 11)] <- var2[var2 %in% c(10, 11)] - 10
					var1[var1==10] <- 0
					var1[var1==11] <- 1
					var2[var2==10] <- 0
					var2[var2==11] <- 1
					## higher
					if (any(var1 - var2 < 0, na.rm=TRUE))
					{
						higher[i] <- FALSE
					}
					## disjoint
					if (sum(var1 * var2, na.rm=TRUE) > 0)
					{
						disjoint[i] <- FALSE
					}
					##atleastone
					if (any(var1 + var2 == 0, na.rm=TRUE))
					{
						atLeastOne[i] <- FALSE
					}
				}

			}
		}
	}
	attr(z, "higher") <- higher
	attr(z, "disjoint") <- disjoint
	attr(z, "atLeastOne") <- atLeastOne

	## report on constraints; adapted from print01Report
	if (any(higher))
	{
		higherSplit <- strsplit(names(higher)[higher], ",")
		report <- sapply(higherSplit, function(x)
		   {
			   paste(c("Network ", x[1], " is higher than network ", x[2],
						".\n"), sep="")
		  })
		message(report)
		cat("This will be respected in the simulations. ")
		cat("If this is not desired, change attribute 'higher'\n")
		cat("by function sienaDataConstraint.\n")
	}
	if (any(disjoint))
	{
		disjointSplit <- strsplit(names(disjoint)[disjoint],',')
		report <- sapply(disjointSplit, function(x)
		   {
			   paste(c("Network ", x[1], " is disjoint from network ",
						x[2], ".\n"), sep="")
		  })
		message(report)
		cat("This will be respected in the simulations.\n")
		cat("If this is not desired, change attribute 'disjoint'\n")
		cat("by function sienaDataConstraint.\n")

	}
	if (any(atLeastOne))
	{
		atLeastOneSplit <- strsplit(names(atLeastOne)[atLeastOne],',')
		report <- sapply(atLeastOneSplit, function(x)
		   {
			   paste(c("A link in at least one of networks ",
						x[1], " and", x[2],
					   " always exists.\n"), sep="")
		  })
		message(report)
		cat("This will be respected in the simulations.")
		cat("If this is not desired, change attribute 'atLeastOne'\n")
		cat("by function sienaDataConstraint.\n")
	}
	z
}

##@rangeAndSimilarity DataCreate
rangeAndSimilarity <- function(vals, rvals=NULL)
{
	zeroOrNA <- function(x){ifelse(is.na(x), TRUE, (x==0))}
	if (is.null(rvals))
	{
		rvals <- range(vals, na.rm=TRUE)
	}
	if (zeroOrNA(var(as.vector(vals), na.rm=TRUE)))
	{
		simTotal <- 0
		simCnt <- sum(!is.na(vals))^2
		simMean <- 0
	}
	else
	{
		vals <- as.matrix(vals)
		rvals1 <- rvals[2] - rvals[1]
		tmp <- apply(vals, 2, function(v)
			 {
				sapply(1: length(v), function(x, y, r)
					{
						z <- y
						z[x] <- NA
						##browser()
						tmp1 <- 1 - abs(y[x] - z) / r
						list(sum(tmp1, na.rm=TRUE), sum(!is.na(tmp1)))
					},
						y=v, r=rvals1)
			 }
					)
		tmp <- unlist(tmp)
		raw <- tmp[seq(1, length(tmp), by=2)]
		cnts <- tmp[seq(2, length(tmp), by=2)]
		simTotal <- sum(raw)
		simCnt <- sum(cnts)
		simMean <- ifelse(simCnt==0, 0, simTotal/simCnt)
    }
    
    # and variance
    sum <- sum(c(vals), na.rm = TRUE)
    sumSq <- sum(c(vals)^2, na.rm = TRUE)
    nonmis <- sum(!is.na(c(vals)))
    variance <- ifelse(nonmis==0, 0, (sumSq/nonmis) - (sum/nonmis)^2)

	list(simTotal=simTotal, simMean=simMean, range=rvals, simCnt=simCnt,
        sum=sum, sumSq=sumSq, variance=variance, varCnt=nonmis)
}
##@groupRangeAndSimilarityAndMean DataCreate
## calculates attributes at group level and re-centers actor covariates
groupRangeAndSimilarityAndMean <- function(group)
{
	atts <- attributes(group)
	##behavs <- atts$types == "behavior"
	netnames <- atts$netnames
	bRange <- namedVector(NA, netnames)
	behRange <- matrix(NA, ncol=length(netnames), nrow=2)
	colnames(behRange) <- netnames
	bSim <- namedVector(NA, netnames)
    bVar <- namedVector(NA, netnames)
	bPoszvar <- namedVector(NA, netnames)
	bMoreThan2 <- namedVector(NA, netnames)
	bAnyMissing <- namedVector(FALSE, netnames)
	for (net in which(atts$types %in% c("behavior", "continuous")))
	{
		simTotal <- 0
		simCnt <- 0
        sumTotal <- 0
        sumSqTotal <- 0
        varCnt <- 0
		anyMissing <- FALSE
		bPoszvar[net] <- TRUE
		thisrange <- matrix(NA, ncol=length(group), nrow=2)
		for (i in 1:length(group))
		{
			j <- match(netnames[net], names(group[[i]]$depvars))
			if (is.na(j))
				stop("network names not consistent")
			depvar <- group[[i]]$depvars[[j]][, 1, ]
			## this should be a matrix
			thisrange[, i] <- range(depvar, na.rm=TRUE)
			if (any(is.na(depvar)))
			{
				anyMissing <- TRUE
			}
		}
		behRange[, net] <- round(range(thisrange, na.rm=TRUE))
		bRange[net] <- behRange[, net][2] - behRange[, net][1]
		if (behRange[, net][2] == behRange[, net][1] && !anyMissing)
		{
			bPoszvar[net] <- FALSE
		}
		values <- NULL
		for (i in 1:length(group))
		{
			j <- match(netnames[net], names(group[[i]]$depvars))
			depvar <- group[[i]]$depvars[[j]][, 1, ]
			## this should be a matrix

			tmp <- rangeAndSimilarity(depvar[, -ncol(depvar)],
									  behRange[, net])
			simTotal <- simTotal + tmp$simTotal
			simCnt <- simCnt + tmp$simCnt
			values <- c(values, unique(depvar))
            
            sumTotal <- sumTotal + tmp$sum
            sumSqTotal <- sumSqTotal + tmp$sumSq
            varCnt <- varCnt + tmp$varCnt
		}
		simMean <- ifelse(simCnt==0, 0, simTotal/simCnt)
		bSim[net] <- simMean
		bMoreThan2[net] <- length(unique(values)) > 2
		if (anyMissing)
			bAnyMissing[net] <- TRUE
        variance <- ifelse(varCnt==0, 0, sumSqTotal/varCnt - (sumTotal/varCnt)^2)
        bVar[net] <- variance
	}
	## constant ones will not exist unless there is only one data object
	cCovarRange <- namedVector(NA, atts$cCovars)
	cCovarSim <- namedVector(NA, atts$cCovars)
	cCovarPoszvar <- namedVector(TRUE, atts$cCovars)
	cCovarMoreThan2 <- namedVector(FALSE, atts$cCovars)
	cCovarMean <- namedVector(NA, atts$cCovars)
	cCovarRange2 <- matrix(NA, ncol=length(atts$cCovars), nrow=2)
	colnames(cCovarRange2) <- atts$cCovars
	for (covar in seq(along=atts$cCovars))
	{
		if (length(group) > 1)
			stop("group create constant covariate error")
		simTotal <- 0
		simCnt <- 0
		anyMissing <- FALSE
		## first find the range
		thisrange <- matrix(NA, ncol=length(group),nrow=2)
		for (i in 1:length(group))
		{
			j <- match(atts$cCovars[covar], names(group[[i]]$cCovars))
			if (is.na(j))
			{
				stop("inconsistent actor covariate names")
			}
			thisrange[, i] <- range(group[[i]]$cCovars[[j]],
								   na.rm=TRUE)
			if (any(is.na(group[[i]]$cCovars[[j]])))
			{
				anyMissing <- TRUE
			}
		}
		rr <- range(thisrange, na.rm=TRUE)
		cCovarRange[covar] <- rr[2] - rr[1]
		if (rr[2] == rr[1] && !anyMissing)
				cCovarPoszvar[covar] <- FALSE
		##then calculate similarity
		for (i in 1:length(group))
		{
			j <- match(atts$cCovars[covar], names(group[[i]]$cCovars))
			tmp <- rangeAndSimilarity(group[[i]]$cCovars[[covar]],
									  rr)
			simTotal <- simTotal + tmp$simTotal
			simCnt <- simCnt + tmp$simCnt
		}
		simMean <- simTotal/simCnt
		cCovarSim[covar] <- simMean
		cCovarMoreThan2[covar] <- attr(group[[i]]$cCovars[[covar]], "moreThan2")
		cCovarMean[covar] <- attr(group[[i]]$cCovars[[covar]], "mean")
		cCovarRange2[, covar] <- attr(group[[i]]$cCovars[[covar]], "range2")
   }

	vCovarRange <- namedVector(NA, atts$vCovars)
	vCovarSim <- namedVector(NA, atts$vCovars)
	vCovarPoszvar <- namedVector(TRUE, atts$vCovars)
	vCovarMoreThan2 <- namedVector(FALSE, atts$vCovars)
	vCovarMean <- namedVector(NA, atts$vCovars)
	for (covar in seq(along=atts$vCovars))
	{
		vartotal <- 0
		nonMissingCount <- 0
		## need to re-centre these values
		for (i in 1:length(group))
		{
			j <- match(atts$vCovars[covar], names(group[[i]]$vCovars))
			j1 <- match(atts$vCovars[covar], names(group[[1]]$vCovars))
			if (is.na(j))
			{
				stop("inconsistent actor covariate names")
			}
			if (attr(group[[i]]$vCovars[[j]],"centered") != attr(group[[1]]$vCovars[[j]],"centered"))
			{
				stop(paste("Inconsistent centering for covariate", names(group[[i]]$vCovars)[j]))
			}
			vartotal <- vartotal + attr(group[[i]]$vCovars[[j]], "vartotal")
			nonMissingCount <- nonMissingCount +
					attr(group[[i]]$vCovars[[j]], "nonMissingCount")
			if (attr(group[[i]]$vCovars[[j]],"centered"))
			{
				group[[i]]$vCovars[[j]] <- group[[i]]$vCovars[[j]] +
					attr(group[[i]]$vCovars[[j]], "vartotal") /
						attr(group[[i]]$vCovars[[j]], "nonMissingCount")
			}
		}
		varmean <- vartotal / nonMissingCount
		j <- match(atts$vCovars[covar], names(group[[1]]$vCovars))
		if (attr(group[[1]]$vCovars[[j]],"centered"))
		{
			for (i in 1:length(group))
			{
				j <- match(atts$vCovars[covar], names(group[[i]]$vCovars))
				if (is.na(j))
				{
					stop("inconsistent actor covariate names")
				}
				group[[i]]$vCovars[[j]] <- group[[i]]$vCovars[[j]] -
					varmean
			}
		}
		simTotal <- 0
		simCnt <- 0
		anyMissing <- FALSE
		## first find the range
		thisrange <- matrix(NA, ncol=length(group), nrow=2)
		values <- NULL
		for (i in 1:length(group))
		{
			j <- match(atts$vCovars[covar], names(group[[i]]$vCovars))
			if (is.na(j))
			{
				stop("inconsistent actor covariate names")
			}
			thisrange[, i] <- range(group[[i]]$vCovars[[j]],
									na.rm=TRUE)
			if (any(is.na(group[[i]]$vCovars[[j]])))
			{
				anyMissing <- TRUE
			}
			values <- c(values, unique(group[[i]]$vCovars[[j]]))
		}
		rr <- range(thisrange, na.rm=TRUE)
		if (rr[2] == rr[1] && !anyMissing)
		{
			vCovarPoszvar[covar] <-	 FALSE
		}
		vCovarRange[covar] <- rr[2] - rr[1]
		##then calculate similarity. Note ignore final observation
		for (i in 1:length(group))
		{
			j <- match(atts$vCovars[covar], names(group[[i]]$vCovars))
			tmpmat <- group[[i]]$vCovars[[covar]]
			tmp <- rangeAndSimilarity(tmpmat, rr)
			simTotal <- simTotal + tmp$simTotal
			simCnt <- simCnt + tmp$simCnt
		}
		simMean <- simTotal/simCnt
		vCovarSim[covar] <- simMean
		## storage.mode(attr(vCovars[[i]], 'range')) <- 'double'
		vCovarMean[covar] <- varmean
		vCovarMoreThan2[covar] <- length(unique(values)) > 2
	}
	dycCovarMean <- namedVector(NA, atts$dycCovars)
	dycCovarRange <- namedVector(NA, atts$dycCovars)
	dycCovarRange2 <- matrix(NA, 2, length(atts$dycCovars))
	colnames(dycCovarRange2) <- atts$dycCovars
	for (covar in seq(along=atts$dycCovars))
	{
		if (length(group) > 1)
			stop("error in dyadic constant covariate, group create")
		j <- match(atts$dycCovars[covar], names(group[[1]]$dycCovars))
		if (is.na(j))
		{
			stop("inconsistent dyadic covariate names")
		}
		dycCovarMean[covar] <- attr(group[[1]]$dycCovars[[j]], "mean")
		dycCovarRange[covar] <- attr(group[[1]]$dycCovars[[j]], "range")
		dycCovarRange2[, covar] <- attr(group[[1]]$dycCovars[[j]], "range2")
	}
	dyvCovarMean <- namedVector(NA, atts$dyvCovars)
	dyvCovarRange <- namedVector(NA, atts$dyvCovars)
	dyvCovarRange2 <- matrix(NA, 2, length(atts$dyvCovars))
	colnames(dyvCovarRange2) <- atts$dyvCovars
	for (covar in seq(along=atts$dyvCovars))
	{
		vartotal <- 0
		nonMissingCount <- 0
		thisrange <- matrix(NA, ncol=length(group), nrow=2)
		for (i in 1:length(group))
		{
			j <- match(atts$dyvCovars[covar], names(group[[i]]$dyvCovars))
			if (is.na(j))
			{
				stop("inconsistent dyadic covariate names")
			}
			sparse <- attr(group[[i]]$dyvCovars[[j]], "sparse")
			vardims <- attr(group[[i]]$dyvCovars[[j]], "vardims")
			if (sparse)
			{
				for (obs in 1:vardims[3])
				{
					vartotal <- vartotal + sum(group[[i]]$dyvCovars[[j]][[obs]],
											   na.rm=TRUE)
					nonMissingCount <- nonMissingCount +
						sum(!is.na(group[[i]]$dyvCovars[[j]][[obs]]))
				}
				thisrange[, i] <- range(sapply(group[[i]]$dyvCovars[[j]], range,
											   na.rm=TRUE),
										na.rm=TRUE)
			}
			else
			{
				vartotal <- vartotal + sum(group[[i]]$dyvCovars[[j]], na.rm=TRUE)
				nonMissingCount <- nonMissingCount +
					sum(!is.na(group[[i]]$dyvCovars[[j]]))
				thisrange[, i] <- range(group[[i]]$dyvCovars[[j]],
										na.rm=TRUE)
			}
			centered <- attr(group[[i]]$dyvCovars[[j]], "centered")
			if (i <= 1) {centered1 <- centered}
			if (centered != centered1)
			{
				stop("inconsistent centering of dyadic covariates")
			}
		}
		dyvCovarMean[covar] <- ifelse(centered, vartotal / nonMissingCount, 0)
		rr <- range(thisrange, na.rm=TRUE)
		dyvCovarRange2[,covar] <- rr
		dyvCovarRange[covar] <- rr[2] - rr[1]
   }
	attr(group, "bRange") <- bRange
	attr(group, "behRange") <- behRange
	attr(group, "bSim") <- bSim
	attr(group, "bPoszvar") <- bPoszvar
	attr(group, "bMoreThan2") <- bMoreThan2
	attr(group, "bAnyMissing") <- bAnyMissing
	attr(group, "cCovarPoszvar") <- cCovarPoszvar
	attr(group, "cCovarMoreThan2") <- cCovarMoreThan2
	attr(group, "cCovarRange") <- cCovarRange
	attr(group, "cCovarRange2") <- cCovarRange2
	attr(group, "cCovarSim") <- cCovarSim
	attr(group, "cCovarMean") <- cCovarMean
	attr(group, "vCovarRange") <- vCovarRange
	attr(group, "vCovarSim") <- vCovarSim
	attr(group, "vCovarMoreThan2") <- vCovarMoreThan2
	attr(group, "vCovarPoszvar") <- vCovarPoszvar
	attr(group, "vCovarMean") <- vCovarMean
	attr(group, "dycCovarMean") <- dycCovarMean
	attr(group, "dycCovarRange") <- dycCovarRange
	attr(group, "dycCovarRange2") <- dycCovarRange2
	attr(group, "dyvCovarRange") <- dyvCovarRange
	attr(group, "dyvCovarRange2") <- dyvCovarRange2
	attr(group, "dyvCovarMean") <- dyvCovarMean
	group
}
##@namedVector DataCreate Utility to create a vector with appropriate names
namedVector <- function(vectorValue, vectorNames, listType=FALSE)
{
	if (listType)
	{
		tmp <- vector("list", length(vectorNames))
	}
	else
	{
		tmp <- rep(vectorValue, length(vectorNames))
	}
	names(tmp) <- vectorNames
	tmp
}


##@createSettings DataCreate
# create a settings structure for sienaData object x
createSettings <- function(x, varName=1, model=TRUE)
##
{
	if (!inherits(x, 'siena'))
	{
		stop('x should be a siena data object')
	}
	if (is.null(x$depvars[[varName]]))
	{
		stop('varName should refer to a dependent network variable in x')
	}
	if ((attr(x$depvars[[varName]], 'type')) != 'oneMode')
	{
		stop('varName should refer to a dependent network variable in x')
	}
	universalOnly <- ifelse(model, "up", "none")
	attr(x$depvars[[varName]], 'settingsinfo') <- list(
		list(id="universal", type="universal", only=universalOnly, covariate=""),
		list(id="primary", type="primary", only="both", covariate=""))
# universal MUST be mentioned first here, and primary second.
	x
}

##@hasSettings DataCreate
# do the dependent variables in sienaData object x have a settings structure
# if varName is specified (number or name),
# this refers only to the indicated variable.
hasSettings <- function(x, varName=NULL)
{
	if (!inherits(x, 'siena'))
	{
		stop('x should be a siena data object')
	}
	if (is.null(varName))
	{
		varNames <- 1:length(x$depvars)
	}
	else
	{
		varNames <- varName
		if (is.numeric(varName))
		{
			if ((varName < 1) | (varName > length(x$depvars)))
			{
				stop('varName should indicate one among the dependent variables')
			}
		}
		else
		{
			if (is.null(x$depvars[[varName]]))
			{
				stop('varName should indicate one among the dependent variables')
			}
		}
	}
	hasSettingsi <- rep(FALSE, length(varNames))
	for (i in seq(along = varNames))
	{
		hasSettingsi[i] <- (!is.null(attr(x$depvars[[varNames[i]]], 'settingsinfo')))
	}
	as.logical(hasSettingsi)
}

##@describeTheSetting DataCreate
# gives a description of the settings structure of depvar
describeTheSetting <- function(depvar)
{
	if (!inherits(depvar, 'sienaDependent'))
	{
		stop('depvar should be a dependent variable')
	}
	settingsinfo <- attr(depvar, 'settingsinfo')
	if (is.null(settingsinfo))
	{
		dts <- ""
	}
	else
	{
		dts <- list()
		for (i in seq_along(settingsinfo))
		{
			dts[[i]] <- unlist(c(name=attr(depvar, 'name'), settingsinfo[[i]]))
		}
		dtsn <- sapply(dts, function(x){x['name']})
		dtsid <- sapply(dts, function(x){x['id']})
		dtstype <- sapply(dts, function(x){x['type']})
		dtscovar <- sapply(dts, function(x){x['covar']})
		dtscovar <- sapply(dtscovar, function(x){ifelse(is.na(x),'none',x)})
		dtsonly <- sapply(dts, function(x){x['only']})
		dtsonly  <- sapply(dtsonly, function(x){ifelse(is.na(x),'both',x)})
		dts <- cbind(dtsn, dtsid, dtstype, dtscovar, dtsonly)
		colnames(dts) <- c('dependent variable', 'setting', 'type', 'covariate', 'direction')
		rownames(dts) <- 1:seq_len(dtsn)
	}
	dts
}



##@sienaGroupCreate DataCreate
sienaGroupCreate <- function(objlist, singleOK=FALSE, getDocumentation=FALSE)
{
	## suppressPackageStartupMessages(require(network))
	##@copyAttributes internal sienaGroupCreate
	copyAttributes <- function(x, y)
	{
		atts <- attributes(y)
		attr(x, "rangep") <- matrix(atts$range2, ncol=ncol(x), nrow=2)
		attr(x, "meanp") <- rep(atts$mean, ncol(x))
		attr(x, "range") <- atts$range
		attr(x, 'mean') <- atts$mean
		attr(x, 'centered') <- atts$centered
		attr(x, 'vartotal') <- atts$vartotal
		attr(x, 'nonMissingCount') <- atts$nonMissingCount
		attr(x, 'simMeans') <- atts$simMeans
		rr <- rangeAndSimilarity(x, atts$range2)
		attr(x, 'poszvar') <- atts$poszvar
		attr(x, 'similarity') <- rr$simValues
		attr(x, 'simMean') <- rr$simMean
		attr(x, 'moreThan2') <- atts$moreThan2
		attr(x, 'name') <- atts$name
		## storage.mode(attr(vCovars[[i]], 'range')) <- 'double'
	   x
	}
	if (getDocumentation)
	{
		return(getInternals())
	}
	if (!is.list(objlist))
	{
		stop('Need a list of objects')
	}
	if (any (sapply(objlist, function(x) !inherits(x, 'siena'))))
	{
		stop('Not a list of valid siena objects')
	}
	if (length(objlist) == 1 && !singleOK)
	{
		stop('Need more than one siena object')
	}
	## get hold of the attributes from the networks and lists of periods etc.
	group <- objlist
	## get the names of the nets and other objects

	netnames <-	 names(objlist[[1]]$depvars)
	pairsnames <- names(attr(objlist[[1]], "higher"))
	cCovars <- names(objlist[[1]]$cCovars)
	if (is.null(cCovars))
	{
		cCovars <- character(0)
	}
	vCovars <- names(objlist[[1]]$vCovars)
	if (is.null(vCovars))
	{
		vCovars <- character(0)
	}
	dycCovars <-  names(objlist[[1]]$dycCovars)
	if (is.null(dycCovars))
	{
		dycCovars <- character(0)
	}
	dyvCovars <- names(objlist[[1]]$dyvCovars)
	if (is.null(dyvCovars))
	{
		dyvCovars <- character(0)
	}
	## nNetnames <- length(netnames)
	anyMissing <- namedVector(FALSE, netnames)
	symmetric <- namedVector(TRUE, netnames)
	structural <- namedVector(FALSE, netnames)
	allUpOnly <- namedVector(TRUE, netnames)
	allDownOnly <- namedVector(TRUE, netnames)
	anyUpOnly <- namedVector(FALSE, netnames)
	anyDownOnly <- namedVector(FALSE, netnames)
	allHigher <- namedVector(TRUE,pairsnames )
	allDisjoint <- namedVector(TRUE, pairsnames)
	allAtLeastOne <- namedVector(TRUE, pairsnames)
	anyHigher <- namedVector(FALSE, pairsnames)
	anyDisjoint <- namedVector(FALSE, pairsnames)
	anyAtLeastOne <- namedVector(FALSE, pairsnames)
	types <- namedVector(NA, netnames)
	nodeSets <- namedVector(NA, netnames, listType=TRUE)
	ccnodeSets <- namedVector(NA, cCovars)
	cvnodeSets <- namedVector(NA, vCovars)
	dycnodeSets <- namedVector(NA, dycCovars, listType=TRUE)
	dyvnodeSets <- namedVector(NA, dyvCovars, listType=TRUE)
	dyctype <- namedVector(NA, dycCovars)
	dyvtype <- namedVector(NA, dyvCovars)
  #	 totalMissings <- namedVector(0, netnames, listType=TRUE)
  #	 nonMissingCount <- namedVector(0, netnames, listType=TRUE)
	observations <- 0
	periodNos <- rep(NA, 2)
	numberMissingNetwork <- rep(0, 2)
	numberMissingBehavior <- rep(0, 2)
	numberNonMissingNetwork <- rep(0, 2)
	numberNonMissingBehavior <- rep(0, 2)
	groupPeriods <- namedVector(NA, names(objlist))
	for (i in 1:length(objlist))
	{
		newobs <- objlist[[i]]$observations
		periodNos[observations + (1 : (newobs - 1))] <-
			observations + i - 1 + (1 : (newobs - 1))
		numberMissingBehavior[observations + (1 : (newobs - 1))] <- 0
		numberNonMissingBehavior[observations + (1 : (newobs - 1))] <- 0
		numberMissingNetwork[observations + (1 : (newobs - 1))] <- 0
		numberNonMissingNetwork[observations + (1 : (newobs - 1))] <- 0
		for (j in 1:length(objlist[[i]]$depvars))
		{
			varname <- names(objlist[[i]]$depvars)[j]
			netnamesub <- match(varname, netnames)
			if (is.na(netnamesub))
			{
				stop('object names inconsistent')
			}
			attribs <- attributes(objlist[[i]]$depvars[[j]])
			if (is.na(types[netnamesub]))
			{
				types[netnamesub] <- attribs[['type']]
			}
			else if (types[netnamesub] != attribs[['type']])
			{
				stop('Inconsistent network types')
			}
			if (is.null(nodeSets[[netnamesub]]))
			{
				  nodeSets[[netnamesub]] <- attribs[['nodeSet']]
			}
			else if (any(nodeSets[[netnamesub]] != attribs[['nodeSet']]))
			{
				stop('Inconsistent node Sets')
			}
			if (attribs[['type']] == 'oneMode')
			{
				if (!attribs[['symmetric']])
					symmetric[netnamesub] <- FALSE
			}
			if (any(attribs[['uponly']]))
			{
				anyUpOnly[netnamesub] <- TRUE
			}
			if (any(attribs[['downonly']]))
			{
				anyDownOnly[netnamesub] <- TRUE
			}
			if (!all(attribs[['uponly']]))
			{
				allUpOnly[netnamesub] <- FALSE
			}
			if (!all(attribs[['downonly']]))
			{
				allDownOnly[netnamesub] <- FALSE
			}
			if (attribs[["missing"]])
			{
				anyMissing[netnamesub] <- TRUE
			}
			if (attribs[["structural"]])
			{
				structural[netnamesub] <- TRUE
			}
			if (attribs[["type"]] == "behavior")
			{
				numberMissingBehavior[observations + (1 : (newobs - 1))] <-
					numberMissingBehavior[observations + (1 : (newobs - 1))] +
						attribs[["noMissingEither"]]
				numberNonMissingBehavior[observations + (1 : (newobs - 1))] <-
					numberNonMissingBehavior[observations +
											 (1 : (newobs - 1))] +
												 attribs[["nonMissingEither"]]
			}
			else
			{
				numberMissingNetwork[observations + (1 : (newobs - 1))] <-
					numberMissingNetwork[observations + (1 : (newobs - 1))] +
						attribs[["noMissingEither"]]
				numberNonMissingNetwork[observations + (1 : (newobs - 1))] <-
					numberNonMissingNetwork[observations + (1 : (newobs - 1))] +
						attribs[["nonMissingEither"]]
			}
		}
		thisHigher <- attr(objlist[[i]], "higher")
		thisDisjoint <- attr(objlist[[i]], "disjoint")
		thisAtLeastOne <- attr(objlist[[i]], "atLeastOne")

		anyHigher[names(thisHigher)[thisHigher]] <- TRUE
		anyDisjoint[names(thisDisjoint)[thisDisjoint]] <- TRUE
		anyAtLeastOne[names(thisAtLeastOne)[thisAtLeastOne]] <- TRUE
		allHigher[names(thisHigher)[!thisHigher]] <- FALSE
		allDisjoint[names(thisDisjoint)[!thisDisjoint]] <- FALSE
		allAtLeastOne[names(thisAtLeastOne)[!thisAtLeastOne]] <- FALSE
		for (j in seq(along=objlist[[i]]$cCovars))
		{
			varname <- names(objlist[[i]]$cCovars)[j]
			covarsub <- match(varname, cCovars)
			if (is.na(covarsub))
			{
				stop('actor covariate names inconsistent')
			}
			attribs <- attributes(objlist[[i]]$cCovars[[j]])
			if (is.na(ccnodeSets[covarsub]))
			{
				ccnodeSets[covarsub] <- attribs[['nodeSet']]
			}
			else if (ccnodeSets[covarsub] != attribs[['nodeSet']])
			{
				stop('Inconsistent node Sets')
			}
		}
		for (j in seq(along=objlist[[i]]$vCovars))
		{
			varname <- names(objlist[[i]]$vCovars)[j]
			covarsub <- match(varname, vCovars)
			if (is.na(covarsub))
			{
				stop('covariate names inconsistent')
			}
			attribs <- attributes(objlist[[i]]$vCovars[[j]])
			if (is.na(cvnodeSets[covarsub]))
			{
				cvnodeSets[covarsub] <- attribs[['nodeSet']]
			}
			else if (cvnodeSets[covarsub] != attribs[['nodeSet']])
			{
				stop('Inconsistent node Sets')
			}
		}

		 for (j in seq(along=objlist[[i]]$dycCovars))
		{
			varname <- names(objlist[[i]]$dycCovars)[j]
			covarsub <- match(varname, dycCovars)
			if (is.na(covarsub))
			{
				stop('dyadic covariate names inconsistent')
			}
			attribs <- attributes(objlist[[i]]$dycCovars[[j]])
			if (is.null(dycnodeSets[[covarsub]]))
			{
				dycnodeSets[[covarsub]] <- attribs[['nodeSet']]
				dyctype[[covarsub]] <- attribs[['type']]
			}
			else
			{
				if (any(dycnodeSets[[covarsub]] != attribs[['nodeSet']]))
				{
					stop('Inconsistent node Sets')
				}
				if (dyctype[[covarsub]] != attribs[['type']])
				{
					stop("Inconsistent constant dyadic covariate types")
				}
			}
		}
		for (j in seq(along=objlist[[i]]$dyvCovars))
		{
			varname <- names(objlist[[i]]$dyvCovars)[j]
			covarsub <- match(varname, dyvCovars)
			if (is.na(covarsub))
			{
				stop('dyadic covariate names inconsistent')
			}
			attribs <- attributes(objlist[[i]]$dyvCovars[[j]])
			if (is.null(dyvnodeSets[[covarsub]]))
			{
				dyvnodeSets[[covarsub]] <- attribs[['nodeSet']]
				dyvtype[[covarsub]] <- attribs[['type']]
			}
			else
			{
				if (any(dyvnodeSets[[covarsub]] != attribs[['nodeSet']]))
				{
					stop('Inconsistent node Sets')
				}
				if (dyvtype[[covarsub]] != attribs[['type']])
				{
					stop("Inconsistent changing dyadic covariate types")
				}
			}
		}
		observations <- observations + objlist[[i]]$observations - 1
		groupPeriods[i] <- newobs
	}
	## if more than one object, now create the group proper
	if (length(objlist) > 1)
	{
		for (i in 1:length(objlist))
		{
			## constant covars turn into changing ones, probably because the
			## actors change. Change on the individual data object.
			const <- objlist[[i]]$cCovars
			vars <- objlist[[i]]$vCovars
			oldnames <- names(vars)
			nVCovar <- length(vars)
			for (j in seq(along=const))
			{
				newcovar <-
					varCovar(matrix(const[[j]],
									ncol = objlist[[i]]$observations - 1,
									nrow=length(const[[j]])),
							 nodeSet=attr(const[[j]], "nodeSet"), warn=FALSE)
				newcovar <- copyAttributes(newcovar, const[[j]])
				nVCovar <- nVCovar + 1
				vars[[nVCovar]] <- newcovar
			}
			names(vars) <- c(oldnames, names(const))
			objlist[[i]]$vCovars <- vars
			objlist[[i]]$cCovars <- list()

			## now the dyadic ones similarly.

			const <- objlist[[i]]$dycCovars
			vars <- objlist[[i]]$dyvCovars
			oldnames <- names(vars)
			nVCovar <- length(vars)
			for (j in seq(along=const))
			{
				oneCentered <- FALSE
				oneNonCentered <- FALSE
				dim3 <- objlist[[i]]$observations - 1
				newcovar <-
					varDyadCovar(array(const[[j]], dim=c(dim(const[[j]]),
												   dim3)),
								 nodeSets=attr(const[[j]], "nodeSet"), warn=FALSE)
				attr(newcovar, "vartotal") <- attr(const[[j]], "vartotal")
				attr(newcovar, "nonMissingCount") <-
					 attr(const[[j]], "nonMissingCount")
				attr(newcovar, "mean") <- attr(const[[j]], "mean")
				attr(newcovar, "meanp") <- rep(attr(const[[j]], "mean"), dim3)
				attr(newcovar, "centered") <- attr(const[[j]], "centered")
				attr(newcovar, "range") <- attr(const[[j]], "range")
				attr(newcovar, "rangep") <- rep(attr(const[[j]], "range"),
												dim3)
				attr(newcovar, "range2") <- matrix(attr(const[[j]], "range2"),
												   ncol=dim3, nrow=2)
				attr(newcovar, 'name') <- attr(const[[j]], "name")
				nVCovar <- nVCovar + 1
				vars[[nVCovar]] <- newcovar
			}
			names(vars) <- c(oldnames, names(const))
			objlist[[i]]$dyvCovars <- vars
			objlist[[i]]$dycCovars <- list()
		}
		## update the working lists as well
		cCovars <- names(objlist[[1]]$cCovars)
		dycCovars <- names(objlist[[1]]$dycCovars)
		vCovars <- names(objlist[[1]]$vCovars)
		dyvCovars <- names(objlist[[1]]$dyvCovars)
	}
	symmetric[types=='behavior'] <- NA
	symmetric[types=='bipartite'] <- FALSE
	group <-  objlist
	attr(group, 'netnames') <- netnames
	attr(group, 'symmetric') <- symmetric
	attr(group, 'structural') <- structural
 #	 attr(group, "totalMissings") <- totalMissings
	attr(group, "numberNonMissingNetwork") <-
		numberNonMissingNetwork[1:observations]
	attr(group, "numberMissingNetwork") <-
		numberMissingNetwork[1:observations]
	attr(group, "numberNonMissingBehavior") <-
		numberNonMissingBehavior[1:observations]
	attr(group, "numberMissingBehavior") <-
		numberMissingBehavior[1:observations]
	attr(group, 'allUpOnly') <- allUpOnly
	attr(group, 'allDownOnly') <- allDownOnly
	attr(group, 'anyUpOnly') <- anyUpOnly
	attr(group, 'anyDownOnly') <- anyDownOnly
	attr(group, 'allHigher') <- allHigher
	attr(group, 'allDisjoint') <- allDisjoint
	attr(group, 'allAtLeastOne') <- allAtLeastOne
	attr(group, 'anyHigher') <- anyHigher
	attr(group, 'anyDisjoint') <- anyDisjoint
	attr(group, 'anyAtLeastOne') <- anyAtLeastOne
	attr(group, 'types') <- types
	attr(group, 'observations') <- observations
	attr(group, 'periodNos') <- periodNos[1:observations]
	attr(group, 'groupPeriods') <- groupPeriods
	attr(group, 'netnodeSets') <- nodeSets
	attr(group, 'cCovars') <- cCovars
	attr(group, 'vCovars') <- vCovars
	attr(group, 'dycCovars') <- dycCovars
	attr(group, 'dyvCovars') <- dyvCovars
	attr(group, 'ccnodeSets') <- ccnodeSets
	attr(group, 'cvnodeSets') <- cvnodeSets
	attr(group, 'dycnodeSets') <- dycnodeSets
	attr(group, 'dyvnodeSets') <- dyvnodeSets
	attr(group, 'compositionChange') <-
		any( sapply(objlist, function(x) length(x$compositionChange)>0))
	## in fact assume all the same - validate earlier!
	if (attr(group, 'compositionChange'))
	{
		exooptions <- sapply(objlist[[1]]$compositionChange, function(x)
							 attr(x, "ccOption"))
		names(exooptions) <- sapply(objlist[[1]]$compositionChange, function(x)
									attr(x, "nodeSet"))
		attr(group, 'exooptions') <- exooptions
	}
	else
	{
		attr(group, 'exooptions') <- character(0)
	}
	if (is.null(names(group)))
	{
		names(group) <- paste('Data', 1:length(group), sep="")
	}
# This is where the names Data1 etc. are created.
	class(group)<- c("sienaGroup", "siena")
	attr(group, "version") <- packageDescription(pkgname, fields = "Version")
	balmeans <- calcBalmeanGroup (group)
	names(balmeans) <- netnames
	attr(group, "balmean") <- balmeans
	structmeans <- calcStructmeanGroup (group)
	names(structmeans) <- netnames
	attr(group, "structmean") <- structmeans
	## calculate overall degree averages
	atts <- attributes(group)
	netnames <- atts$netnames
	types <- atts$types
	## cat(types,'\n')
	degrees <- namedVector(NA, netnames)
	for (net in seq(along=netnames))
	{
		if (types[net] != "behavior")
		{
			degree <- 0
			nDegree <- 0
			for (i in 1: length(group))
			{
				j <- match(netnames[net], names(group[[i]]$depvars))
				if (is.na(j))
					stop("network names not consistent")
				depvar <- group[[i]]$depvars[[j]]
				degs <- attr(depvar, "degree")
				degree <- degree + sum(degs)
				nDegree <- nDegree + length(degs)
			}
			degrees[net] <- degree / nDegree
		}
	}
	attr(group, "averageOutDegree") <- degrees
	attr(group, "averageInDegree") <- degrees
	group <- groupRangeAndSimilarityAndMean(group)
	bAnyMissing <- attr(group, "bAnyMissing")
	attr(group, "anyMissing") <- anyMissing | bAnyMissing
	attr(group, "bAnyMissing") <- NULL
	tmp <- getGroupNetRanges(group)
	colnames(tmp) <- netnames
	attr(group, "netRanges") <- tmp
	group <- copyGroupAttributes(group, "depvars", "behRange", "behRange")
	if (length(group) >= 2)
	{
		##copy the global attributes down to individual level where appropriate
		##group <- copyGroupAttributes(group, "depvars", "balmean", "balmean")
		##group <- copyGroupAttributes(group, "depvars", "structmean",
		## "structmean")
		group <- copyGroupAttributes(group, "depvars", "symmetric", "symmetric")
		##group <- copyGroupAttributes(group, "depvars", "averageInDegree",
		##							   "averageInDegree")
		##group <- copyGroupAttributes(group, "depvars", "averageOutDegree",
		##							   "averageOutDegree")
		##group <- copyGroupAttributes(group, "depvars", "bSim", "simMean")
		group <- copyGroupAttributes(group, "depvars", "bposzvar", "poszvar")
		group <- copyGroupAttributes(group, "depvars", "bRange", "range")
		group <- copyGroupAttributes(group, "depvars", "bMoreThan2",
								 "moreThan2")
		group <- copyGroupAttributes(group, "depvars", "anyMissing", "missing")
		group <- copyGroupAttributes(group, "depvars", "structural",
								 "structural")
		##group <- copyGroupAttributes(group, "vCovars", "vCovarSim", "simMean")
		group <- copyGroupAttributes(group, "vCovars", "vCovarRange", "range",
								 TRUE)
		##group <- copyGroupAttributes(group, "vCovars", "vCovarMean", "mean",
		## TRUE)
		group <- copyGroupAttributes(group, "vCovars", "vCovarPoszvar",
								 "poszvar")
		group <- copyGroupAttributes(group, "vCovars", "vCovarMoreThan2",
								 "moreThan2")
		##group <- copyGroupAttributes(group, "dycCovars", "dycCovarMean",
		##	"mean")
		group <- copyGroupAttributes(group, "dycCovars", "dycCovarRange2",
								 "range2", TRUE)
		##group <- copyGroupAttributes(group, "dyvCovars", "dyvCovarMean",
		## "mean")
		group <- copyGroupAttributes(group, "dyvCovars", "dyvCovarRange",
								 "range", TRUE)
		group <- copyGroupAttributes(group, "dyvCovars", "dyvCovarRange2",
								 "range2", TRUE)
	}
	group
}
##@copyGroupAttributes DataCreate
## copies attribute groupattrname from the top level to the attribute
## attrname of the corresponding objects of vartype in the subsidiary data
## objects.
copyGroupAttributes <- function(group, vartype, groupattrname, attrname,
								storage.mode=FALSE)
{
	attributeValue <- attr(group, groupattrname)
	for (i in names(group))
	{
		for (j in names(group[[i]][[vartype]]))
		{
			if (is.matrix(attributeValue))
			{
				attr(group[[i]][[vartype]][[j]], attrname) <-
					as.vector(attributeValue[, j])
			}
			else
			{
				attr(group[[i]][[vartype]][[j]], attrname) <-
					as.vector(attributeValue[j])
			}
			if (storage.mode)
				 storage.mode(attr(group[[i]][[vartype]][[j]], attrname)) <-
					 "double"
		}
	}
	group
}

##@calcBalmeanGroup DataCreate
calcBalmeanGroup <- function(data)
{
	atts <- attributes(data)
	netnames <- atts$netnames
	types <- atts$types
   ## cat(types,'\n')
	balmeans <- rep(NA, length(netnames))
	for (net in seq(along=netnames))
	{
		if (types[net] == "oneMode")
		{
			tempra <- 0
			temprb <- 0
			for (i in 1: length(data))
			{
				j <- match(netnames[net], names(data[[i]]$depvars))
				if (is.na(j))
					stop("network names not consistent")
				depvar <- data[[i]]$depvars[[j]]
				dims <- attr(depvar, "netdims")
				for (k in 1 : (dims[3] - 1))
				{
					## remove structural values here?
					if (attr(depvar, "sparse"))
					{
						tmp <- depvar[[k]]
						diag(tmp)  <-  NA ## just in case
						x1 <- tmp@x
						struct <- !is.na(x1) & x1 %in% c(10, 11)
						x1[struct] <- x1[struct] - 10
						tmp@x <- x1
						tmp1 <- colSums(is.na(tmp))
						tempra <- tempra +
							sum(2 * colSums(tmp, na.rm=TRUE) *
								(dims[1] - tmp1 - colSums(tmp, na.rm=TRUE)),
								na.rm=TRUE)
						tmp2 <- colSums(!is.na(tmp))
						temprb <- temprb + sum(tmp2 * (tmp2 - 1))
					}
					else
					{
						tmp <- depvar[, , k]
						diag(tmp) <- NA ## just in case
						struct <- !is.na(tmp) & tmp %in% c(10, 11)
						tmp[struct] <- tmp[struct] - 10
						tempra <- tempra +
							sum(2 * colSums(tmp, na.rm=TRUE) *
								colSums(1 - tmp, na.rm=TRUE),
								na.rm=TRUE)
						tmp1 <- colSums(is.na(tmp))
						tmp2 <- colSums(!is.na(tmp))
						temprb <- temprb + sum(tmp2 * (tmp2 - 1))
					}
				}
			}
			balmeans[net] <- tempra / temprb
		}
	}
	balmeans
}
##@calcStructmeanGroup DataCreate
calcStructmeanGroup <- function(data)
{
	atts <- attributes(data)
	netnames <- atts$netnames
	types <- atts$types
   ## cat(types,'\n')
	structmeans <- rep(NA, length(netnames))
	for (net in seq(along=netnames))
	{
		if (types[net] == "oneMode")
		{
			tempra <- 0
			temprb <- 0
			for (i in 1: length(data))
			{
				j <- match(netnames[net], names(data[[i]]$depvars))
				if (is.na(j))
					stop("network names not consistent")
				depvar <- data[[i]]$depvars[[j]]
				dims <- attr(depvar, "netdims")
				for (k in 1 : (dims[3] - 1))
				{
					## remove structural values here?
					if (attr(depvar, "sparse"))
					{
						tmp <- depvar[[k]]
						diag(tmp)  <-  NA ## just in case
						x1 <- tmp@x
						struct <- !is.na(x1) & x1 %in% c(10, 11)
						x1[struct] <- x1[struct] - 10
						tmp@x <- x1
						tmp1 <- rowSums(is.na(tmp))
						tempra <- tempra +
							sum(2 * rowSums(tmp, na.rm=TRUE) *
								(dims[1] - tmp1 - rowSums(tmp, na.rm=TRUE)),
								na.rm=TRUE)
						tmp2 <- rowSums(!is.na(tmp))
						temprb <- temprb + sum(tmp2 * (tmp2 - 1))
					}
					else
					{
						tmp <- depvar[, , k]
						diag(tmp) <- NA ## just in case
						struct <- !is.na(tmp) & tmp %in% c(10, 11)
						tmp[struct] <- tmp[struct] - 10
						tempra <- tempra +
							sum(2 * rowSums(tmp, na.rm=TRUE) *
								rowSums(1 - tmp, na.rm=TRUE),
								na.rm=TRUE)
						tmp1 <- rowSums(is.na(tmp))
						tmp2 <- rowSums(!is.na(tmp))
						temprb <- temprb + sum(tmp2 * (tmp2 - 1))
					}
				}
			}
			structmeans[net] <- tempra / temprb
		}
	}
	structmeans
}
##@calcStructmean DataCreate
calcStructmean <- function(depvar)
{
	tempra <- 0
	temprb <- 0
	dims <- attr(depvar, "netdims")
	for (k in 1 : (dims[3] - 1))
	{
		if (attr(depvar, "sparse"))
		{
			tmp <- depvar[[k]]
			diag(tmp)  <-  NA ## just in case
			x1 <- tmp@x
			struct <- !is.na(x1) & x1 %in% c(10, 11)
			x1[struct] <- x1[struct] - 10
			tmp@x <- x1
			tmp1 <- rowSums(is.na(tmp))
			tempra <- tempra +
				sum(2 * rowSums(tmp, na.rm=TRUE) *
					(dims[1] - tmp1 - rowSums(tmp, na.rm=TRUE)),
					na.rm=TRUE)
			tmp2 <- rowSums(!is.na(tmp))
			temprb <- temprb + sum(tmp2 * (tmp2 - 1))
	 }
		else
		{
			##extract this observation
			tmp <- depvar[, , k]
			##remove diagonal
			diag(tmp) <- NA ## just in case
			##subtract 10 from structurals
			struct <- !is.na(tmp) & tmp %in% c(10, 11)
			tmp[struct] <- tmp[struct] - 10
			## rowSums here counts how many 1's or 0's in tmp
			tempra <- tempra +
				sum(2 * rowSums(tmp, na.rm=TRUE) *
					rowSums(1 - tmp, na.rm=TRUE),
					na.rm=TRUE)
			tmp2 <- rowSums(!is.na(tmp)) ## counts non-missings in row
			temprb <- temprb + sum(tmp2 * (tmp2 - 1))
		}
	}
	structmean <- tempra / temprb
	#  cat(tempra, temprb, structmean,'\n')
	structmean
}
##@calcBalmean DataCreate
calcBalmean <- function(depvar)
{
	tempra <- 0
	temprb <- 0
	dims <- attr(depvar, "netdims")
	for (k in 1 : (dims[3] - 1))
	{
		if (attr(depvar, "sparse"))
		{
			tmp <- depvar[[k]]
			diag(tmp)  <-  NA ## just in case
			x1 <- tmp@x
			struct <- !is.na(x1) & x1 %in% c(10, 11)
			x1[struct] <- x1[struct] - 10
			tmp@x <- x1
			tmp1 <- colSums(is.na(tmp))
			tempra <- tempra +
				sum(2 * colSums(tmp, na.rm=TRUE) *
					(dims[1] - tmp1 - colSums(tmp, na.rm=TRUE)),
					na.rm=TRUE)
			tmp2 <- colSums(!is.na(tmp))
			temprb <- temprb + sum(tmp2 * (tmp2 - 1))
	 }
		else
		{
			##extract this observation
			tmp <- depvar[, , k]
			##remove diagonal
			diag(tmp) <- NA ## just in case
			##subtract 10 from structurals
			struct <- !is.na(tmp) & tmp %in% c(10, 11)
			tmp[struct] <- tmp[struct] - 10
			## colSums here counts how many 1's or 0's in tmp
			tempra <- tempra +
				sum(2 * colSums(tmp, na.rm=TRUE) *
					colSums(1 - tmp, na.rm=TRUE),
					na.rm=TRUE)
			tmp2 <- colSums(!is.na(tmp)) ## counts non-missings in column
			temprb <- temprb + sum(tmp2 * (tmp2 - 1))
		}
	}
	balmean <- tempra / temprb
	#  cat(tempra, temprb,balmean,'\n')
	balmean
}

##@getGroupNetRanges DataCreate
getGroupNetRanges <- function(data)
{
	atts <- attributes(data)
	netnames <- atts$netnames
	types <- atts$types
   ## cat(types,'\n')
	ranges <- matrix(NA, ncol=length(netnames), nrow=2)
	varmin <- NA
	varmax <- NA
	for (net in seq(along=netnames))
	{
		if (types[net] %in% c("oneMode", "bipartite"))
		{
			varmin <- NA
			varmax <- NA
			for (i in 1: length(data))
			{
				j <- match(netnames[net], names(data[[i]]$depvars))
				if (is.na(j))
					stop("network names not consistent")
				depvar <- data[[i]]$depvars[[j]]
				### need to exclude structural values
				if (attr(depvar, "sparse"))
				{
			  ##   varmin <- min(varmin, sapply(depvar, function(x)
			  ##							   {
			  ##								   tmp <- x@x
			  ##								   min(tmp[!(is.na(tmp) |
			  ##										 tmp %in% c(10, 11))])
			  ##							   }), na.rm=TRUE)
					varmin <- 0 ## in sparse matrices 0's are not there
					varmax <- max(varmax, sapply(depvar, function(x)
											 {
												 tmp <- x@x
												 ifelse(length(tmp)==0,0,
												 max(tmp[!(is.na(tmp) |
														   tmp %in% c(10, 11))]))
											 }), na.rm=TRUE)
				}
				else
				{
					varmin <- min(varmin, depvar[!(is.na(depvar) |
												   depvar %in% c(10,11))],
								  na.rm=TRUE)
					varmax <- max(varmax, depvar[!(is.na(depvar) |
												   depvar %in% c(10,11))],
								  na.rm=TRUE)
				}
			}
			ranges[, net] <- c(varmin, varmax)
		}
	}
	ranges
}

##@covarDist2 DataCreate calculate average alter values for covariate wrt each net
covarDist2 <- function(z)
{
	netNames <- names(z$depvars)
	netTypes <- sapply(z$depvars, function(x)attr(x, "type"))
	netActorSet <- sapply(z$depvars, function(x)
					  {
						  if (attr(x, "type") == "bipartite")
						  {
							  attr(x, "nodeSet")[2]
						  }
						  else
						  {
							  attr(x, "nodeSet")
						  }
					  }
						  )
	for (i in seq(along=z$cCovars))
	{
		nodeSet <- attr(z$cCovars[[i]], "nodeSet")
		use <- (!(netTypes %in% c("behavior", "continuous")) & (netActorSet == nodeSet))
		simMeans <- namedVector(NA, netNames[use])
		for (j in which(use))
		{
			simMeans[netNames[j]] <-
				calcCovarDist2(matrix(z$cCovars[[i]], ncol=1),
							   z$depvars[[j]])
		}
		attr(z$cCovars[[i]], "simMeans") <- simMeans
	}
	for (i in seq(along=z$vCovars))
	{
		nodeSet <- attr(z$vCovars[[i]], "nodeSet")
		use <- (!(netTypes %in% c("behavior", "continuous")) & (netActorSet == nodeSet))
		simMeans <- namedVector(NA, netNames[use])
		for (j in which(use))
		{
			simMeans[netNames[j]] <-
				calcCovarDist2(z$vCovars[[i]], z$depvars[[j]])
		}
		attr(z$vCovars[[i]], "simMeans") <- simMeans

	}
	for (i in seq(along=z$depvars))
	{
		if (netTypes[i] %in% c("behavior", "continuous"))
		{
			beh <- z$depvars[[i]][, 1, ]
			## take off the mean NB no structurals yet!
			beh <- beh - mean(beh, na.rm=TRUE)
			nodeSet <- netActorSet[i]
			use <- (!(netTypes %in% c("behavior", "continuous"))
			        & netActorSet == nodeSet)
			simMeans <- namedVector(NA, netNames[use])
			for (j in which(use))
			{
				simMeans[netNames[j]] <-
					calcCovarDist2(beh[, -ncol(beh), drop=FALSE],
								   z$depvars[[j]], rval=range(beh, na.rm=TRUE))
			}
			attr(z$depvars[[i]], "simMeans") <- simMeans
		}
	}
	z
}
##@calcCovarDist2 DataCreate similarity mean for alter of this depvar
calcCovarDist2 <- function(covar, depvar, rval=NULL)
{
	## remove final obs from depvars
	observations <- attr(depvar, "netdims")[3] - 1

	simTotal <- rep(NA, observations)
	simCnt <- rep(NA, observations)

	for (i in 1:observations)
	{
		if (ncol(covar) == 1) # maybe constant covariate used for each obs
		{
			xx <- covar[, 1]
		}
		else
		{
			xx <- covar[, i]
		}
		if (attr(depvar, "sparse"))
		{
			dep <- depvar[[i]]
			dep@x[dep@x %in% c(10, 11)] <- dep@x[dep@x %in% c(10, 11)] - 10
		}
		else
		{
			dep <- depvar[, , i, drop=FALSE]
			dep[dep %in% c(10, 11)] <- dep[dep %in% c(10, 11)] - 10
		}
		vi <- apply(dep, 1, function(x)
				{
					if (sum(x, na.rm=TRUE) > 0)
					{
						nonMissing <- sum(!is.na(xx[!is.na(x) & x > 0]))
						if (nonMissing > 0)
						{
							sum(x * xx, na.rm=TRUE)/ sum(x, na.rm=TRUE)
						}
						else
						{
							NA
						}
					}
					else
					{
						0
					}

				}
					)

		tmp <- rangeAndSimilarity(vi, rval)
		simTotal[i] <- tmp$simTotal
		simCnt[i] <- tmp$simCnt
	}
	ifelse(sum(simCnt)==0, 0, sum(simTotal)/sum(simCnt))
}
