#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: effects.r
# *
# * Description: This module contains the code for the creation of the
# * effects object to go with a Siena data object or group object.
# *****************************************************************************/
##@substituteNames DataCreate replace xxxxxx, yyyyyy, zzzzzz
substituteNames <- function(nameVectors, xName = NULL, yName = NULL, zName = NULL) {
	effects <- nameVectors[, c("effectName", "functionName", "interaction1",
		"interaction2")]
	if (!is.null(xName)) {
		effects <- sapply(effects, function(x) gsub("xxxxxx", xName, x))
	}
	if (!is.null(yName)) {
		effects <- sapply(effects, function(x) gsub("yyyyyy", yName, x))
	}
	if (!is.null(zName)) {
		effects <- sapply(effects, function(x) gsub("zzzzzz", zName, x))
	}
	nameVectors[, c("effectName", "functionName", "interaction1", "interaction2")] <- effects
	nameVectors
}

##@createEffects DataCreate Extract required rows and change text
createEffects <- function(effectGroup, xName=NULL, yName=NULL, zName = NULL,
	name, groupName, group, netType)
{
	effects <- RSiena::allEffects[RSiena::allEffects$effectGroup == effectGroup, ]
	if (nrow(effects) == 0)
	{
		stop("empty effect group", effectGroup)
	}
	if (any(is.na(effects$effectName)))
	{
		stop("missing effect name")
	}
	effects <- substituteNames(effects, xName, yName, zName)
	effects$effectGroup <- NULL
	nn <- nrow(effects)
	# ignore the endowment field of 'gmm' type effects
	effects$endowment[which(effects$type == 'gmm')] <- FALSE
	# If we have some with a valid endowment field in the selection.
	if (!all(is.na(effects$endowment))) {
		# create a new vector with the right dimension
		neweffects <- effects[rep(1:nn, times=(1 + 2 * as.numeric(effects$endowment))), ]
		# fill with eval, endow and creation effects
		neweffects$type <- unlist(by(effects, 1:nn,
				function(e) {
					if (e$endowment) {
						c('eval', 'endow', 'creation')
					} else {
						e$type
					}
				}))
		effects <- neweffects
		nn <- nrow(effects)
	}
	effects$endowment <-  NULL
	effectFn <- vector('list', nn)
	statisticFn <- vector('list', nn)
	effects$effectFn <- effectFn
	effects$statisticFn <- statisticFn
	effects$netType <- netType
	effects$groupName <- groupName
	effects$group <- group
	effectsname <- rep(name, nn)
	effects <- data.frame(name=effectsname, effects, stringsAsFactors=FALSE)
	effects
}

##@getEffects DataCreate create effects object
getEffects<- function(x, nintn = 10, behNintn=4, getDocumentation=FALSE, onePeriodSde=FALSE)
{
	##@duplicateDataFrameRow internal getEffects Put period numbers in
	duplicateDataFrameRow <- function(x, n)
	{
		tmp <- NULL
		for (i in 1:n)
		{
			xx <- x
			xx[, c("effectName", "functionName", "period")] <-
				sub("nnnnnn", periodNos[i], xx[, c("effectName", "functionName",
						"period")])
			tmp <- rbind(tmp, xx)
		}
		tmp
	}

	##@networkRateEffects internal getEffects create a set of rate effects
	networkRateEffects <- function(depvar, varname, symmetric, bipartite)
	{
		if (symmetric)
		{
			rateEffects <- createEffects("symmetricRate", varname, name=varname,
				groupName=groupName, group=group,
				netType=netType)
		}
		else if (bipartite)
		{
			rateEffects <- createEffects("bipartiteRate", varname, name=varname,
				groupName=groupName, group=group,
				netType=netType)
		}
		else
		{
			rateEffects <- createEffects("nonSymmetricRate", varname,
				name=varname,
				groupName=groupName, group=group,
				netType=netType)
		}
		if (observations == 1)
		{
			rateEffects <- rateEffects[-2, ] ## remove the extra period
		}
		else
		{
			## get correct number of rows
			rateEffects <- rbind(duplicateDataFrameRow(rateEffects[2, ],
					observations),
				rateEffects[-c(1, 2), ])
		}
		rateEffects
	}
	##@oneModeNet internal getEffects process a one mode net
	oneModeNet <- function(depvar, varname)
	{
		symmetric <- attr(depvar, "symmetric")
		nodeSet <- attr(depvar, 'nodeSet')

		rateEffects <- networkRateEffects(depvar, varname, symmetric=symmetric,
			bipartite=FALSE)

		if (symmetric)
		{
			objEffects <- createEffects("symmetricObjective", varname,
				name=varname,
				groupName=groupName, group=group,
				netType=netType)
		}
		else
		{
			objEffects <- createEffects("nonSymmetricObjective", varname,
				name=varname,
				groupName=groupName, group=group,
				netType=netType)
		}
		for (j in seq(along = xx$dycCovars))
		{
			if (attr(xx$dycCovars[[j]], "type") == "oneMode" &&
				attr(xx$dycCovars[[j]], 'nodeSet')[1] == nodeSet)
			{
				objEffects <- rbind(objEffects,
					createEffects("dyadObjective",
						names(xx$dycCovars)[j],
						name=varname,
						groupName=groupName,
						group=group,
						netType=netType))
			}
		}
		for (j in seq(along = xx$dyvCovars))
		{
			if (attr(xx$dyvCovars[[j]], "type") == "oneMode" &&
				attr(xx$dyvCovars[[j]], 'nodeSet')[1] == nodeSet)
			{
				objEffects <- rbind(objEffects,
					createEffects("dyadObjective",
						names(xx$dyvCovars)[j],
						name=varname,
						groupName=groupName,
						group=group,
						netType=netType))
			}
		}
		for (j in seq(along = xx$cCovars))
		{
			if (attr(xx$cCovars[[j]], 'nodeSet') == nodeSet)
			{
				tmp <- covarOneModeEff(names(xx$cCovars)[j],
					attr(xx$cCovars[[j]], 'poszvar'),
					attr(xx$cCovars[[j]], 'moreThan2'),
					symmetric, constant=TRUE, name=varname)
				objEffects <-  rbind(objEffects, tmp$objEff)
				rateEffects <- rbind(rateEffects, tmp$rateEff)
			}
		}
		for (j in seq(along=xx$depvars))
		{
			if (types[j] %in% c('behavior', 'continuous') &&
				attr(xx$depvars[[j]], 'nodeSet') == nodeSet)
			{
				tmp <- covarOneModeEff(names(xx$depvars)[j],
					poszvar=TRUE,
					attr(xx$depvars[[j]], 'moreThan2'),
					symmetric, constant=FALSE, name=varname)
				objEffects <- rbind(objEffects, tmp$objEff)
				rateEffects <- rbind(rateEffects, tmp$rateEff)
			}
		}
		for (j in seq(along=xx$vCovars))
		{
			if (attr(xx$vCovars[[j]], 'nodeSet') == nodeSet)
			{
				tmp <- covarOneModeEff(names(xx$vCovars)[j],
					attr(xx$vCovars[[j]], 'poszvar'),
					attr(xx$vCovars[[j]], 'moreThan2'),
					symmetric, constant=FALSE, name=varname)
				objEffects <- rbind(objEffects,tmp$objEff)
				rateEffects<- rbind(rateEffects,tmp$rateEff)
			}
		}

		if (length(xx$cCovars) + length(xx$vCovars) +
			length(xx$dycCovars) + length(xx$dyvCovars) +
			length(types=='behavior') + length(types=='continuous') > 0)
		{
			interaction <- createEffects("unspecifiedNetInteraction",
				name=varname,
				groupName=groupName, group=group,
				netType=netType)
			objEffects <-  rbind(objEffects, interaction[rep(1:nrow(interaction), nintn), ])
		}

		for (j in seq(along=xx$depvars))
		{
			otherName <- names(xx$depvars)[j]
			if (types[j] == 'oneMode' &&
				attr(xx$depvars[[j]], 'nodeSet') == nodeSet &&
				varname != otherName)
			{
				if (attr(xx$depvars[[j]], "symmetric"))
				{
					objEffects <-
						rbind(objEffects,
							createEffects("nonSymmetricSymmetricSObjective",
								otherName, name=varname,
								groupName=groupName, group=group,
								netType=netType))
					objEffects <-
						rbind(objEffects,
							createEffects("nonSymmetricSymmetricObjective",
								otherName, name=varname,
								groupName=groupName, group=group,
								netType=netType))
				}
				else
				{
					objEffects <-
						rbind(objEffects,
							createEffects("nonSymmetricNonSymmetricObjective",
								otherName, name=varname,
								groupName=groupName, group=group,
								netType=netType))
				}
			}

			for (k in seq(along=xx$depvars))
			{
				# This is an effect group that can be used
				# for any three dependent networks, provided the first is one-mode;
				# the receiver sets of the second and third networks
				# must be the same.
				# Note that some or all of the three networks may be the same.
				if ((types[j] == "bipartite") && (types[k] == "bipartite"))
				{
					suitable <- (attr(xx$depvars[[j]], 'nodeSet')[2] ==
						attr(xx$depvars[[k]], 'nodeSet')[2])
				}
				else
				{
					suitable <- ((types[j] == "oneMode") && (types[k] == "oneMode"))
				}
				if (suitable)
				{
					thirdName <- names(xx$depvars)[k]
					objEffects <-
						rbind(objEffects,
							createEffects("tripleNetworkObjective",
								otherName, thirdName, name=varname,
								groupName=groupName, group=group,
								netType=netType))
				}
			}
			if (types[j] == 'bipartite' &&
				(attr(xx$depvars[[j]], 'nodeSet')[1] == nodeSet))
			{
				objEffects <-
					rbind(objEffects,
						createEffects("nonSymmetricBipartiteObjective",
							otherName, name=varname,
							groupName=groupName, group=group,
							netType=netType))
			}
			if ((!(types[j] %in% c('behavior', 'continuous'))) && (varname != otherName))
			{
				for (k in seq(along=xx$cCovars))
				{
					if (attr(xx$cCovars[[k]], 'nodeSet') == nodeSet)
					{
						objEffects <-
							rbind(objEffects,
								covarNetNetEff(otherName,
									covarname=names(xx$cCovars)[k],
									attr(xx$depvars[[j]], 'nodeSet'),
									attr(xx$cCovars[[k]], 'nodeSet'),
									attr(xx$cCovars[[k]], 'poszvar'),
									name=varname))
					}
				}
				for (k in seq(along=xx$vCovars))
				{
					if (attr(xx$vCovars[[k]], 'nodeSet') == nodeSet)
					{
						objEffects <-
							rbind(objEffects,
								covarNetNetEff(otherName,
									covarname=names(xx$vCovars)[k],
									attr(xx$depvars[[j]], 'nodeSet'),
									attr(xx$vCovars[[k]], 'nodeSet'),
									attr(xx$vCovars[[k]], 'poszvar'),
									name=varname))
					}
				}
				for (k in seq(along = xx$dycCovars))
				{
					if (attr(xx$dycCovars[[k]], "type") == "oneMode" &&
						attr(xx$dycCovars[[k]], 'nodeSet')[1] == nodeSet)
					{
						othervarname <- names(xx$dycCovars)[k]
						objEffects <- rbind(objEffects,
							createEffects("dyadANetNetObjective",
								varname, otherName, othervarname,
								name=varname, groupName=groupName,
								group=group, netType=netType))
					}
				}
				for (k in seq(along = xx$dyvCovars))
				{
					if (attr(xx$dyvCovars[[k]], "type") == "oneMode" &&
						attr(xx$dyvCovars[[k]], 'nodeSet')[1] == nodeSet)
					{
						othervarname <- names(xx$dyvCovars)[k]
						objEffects <- rbind(objEffects,
							createEffects("dyadANetNetObjective",
								varname, otherName, othervarname,
								name=varname, groupName=groupName,
								group=group, netType=netType))
					}
				}
				for (k in seq(along=xx$depvars))
				{
					if ((types[k] %in% c('behavior', 'continuous')) &&
						(attr(xx$depvars[[k]], 'nodeSet') == nodeSet))
					{
						objEffects <-
							rbind(objEffects,
								covarNetNetEff(otherName,
									covarname=names(xx$depvars)[k],
									attr(xx$depvars[[j]], 'nodeSet'),
									attr(xx$depvars[[k]], 'nodeSet'),
									poszvar=TRUE,
									name=varname))
					}
				}

			}
		}
		if ((nOneModes + nBipartites) > 1) ## add the network name
		{
			objEffects$functionName <- paste(varname, ': ',
				objEffects$functionName, sep = '')
			objEffects$effectName <- paste(varname, ': ',
				objEffects$effectName, sep = '')
		}

		## replace the text for endowment and creation effects
		tmp <- objEffects$functionName[objEffects$type =='endow']
		tmp <- paste('Lost ties:', tmp)
		objEffects$functionName[objEffects$type == 'endow'] <- tmp
		tmp <- objEffects$functionName[objEffects$type =='creation']
		tmp <- paste('New ties:', tmp)
		objEffects$functionName[objEffects$type == 'creation'] <- tmp
		tmp <- objEffects$functionName[objEffects$type =='gmm']
		tmp <- paste('Gmm statistic:', tmp)
		objEffects$functionName[objEffects$type == 'gmm'] <- tmp

		## get starting values
		starts <- getNetworkStartingVals(depvar)
		##set defaults
		rateEffects[1:noPeriods, "include"] <- TRUE
		rateEffects[1:noPeriods, "initialValue"] <- starts$startRate
		rateEffects$basicRate[1:observations] <- TRUE

		objEffects[objEffects$shortName == "density" &
			objEffects$type == "eval",'randomEffects'] <- TRUE
		objEffects[objEffects$shortName == "linear" &
			objEffects$type == "eval",'randomEffects'] <- TRUE

		objEffects$untrimmedValue <- rep(0, nrow(objEffects))
		if (attr(depvar,'symmetric'))
		{
			objEffects[objEffects$shortName == "density" &
				objEffects$type == "eval",
			c('include', "initialValue", "untrimmedValue")] <-
				list(TRUE, starts$degree, starts$untrimmed)
#			objEffects[objEffects$shortName=='transTriads' &
#					   objEffects$type=='eval','include'] <- TRUE
		}
		else
		{
			if (!(attr(depvar,'allUpOnly') || attr(depvar, 'allDownOnly')))
			{
				objEffects[objEffects$shortName == "density" &
					objEffects$type == 'eval',
				c('include', "initialValue", "untrimmedValue")] <-
					list(TRUE, starts$degree, starts$untrimmed)
			}
			else
			{
				objEffects <-
					objEffects[!objEffects$shortName == "density", ]
			}
			objEffects[objEffects$shortName == 'recip'&
				objEffects$type == 'eval', 'include'] <- TRUE
		}
		rateEffects$basicRate[1:observations] <- TRUE

		if (!is.null(attr(depvar,"settingsinfo")))
		{
			settingIds <- sapply(attr(depvar,"settingsinfo"), function(s) s$id)
			nbrSettings <- length(settingIds)

			if ("primary" %in% settingIds) {
				objEffects <- rbind(objEffects, createEffects(
					"settingsObjective", varname, name=varname,
					groupName=groupName, group=group, netType=netType))
				# append effects with an interaction on the primary settings network of `varname`
				objEffects <- rbind(objEffects, createEffects(
					"nonSymmetricSymmetricSObjective", paste0("primary(", varname, ")") , name=varname,
					groupName=groupName, group=group, netType=netType))
			}

			# duplicate rate effects, split by periods
			setRate <- rateEffects[1:observations, ]
			setRate <- setRate[rep(1:nrow(setRate), each=nbrSettings), ]
			setRateByPeriod <- split(setRate, list(setRate$group, setRate$period))

			# for each period, set "setting" and modify "effectName"
			setRateByPeriod <- lapply(setRateByPeriod, function(dd) {
					dd$setting <- settingIds
					dd$interaction1 <- settingIds
# also see below with 0.75
					dd$initialValue[dd$setting != 'primary'] <- 0.5
					dd$initialValue[dd$setting == 'primary'] <-
										0.75*dd$initialValue[dd$setting == 'primary']
					dd$effectName <- paste(dd$setting, ' setting rate ',
										varname,' (period ', dd$period, ')', sep='')
					dd$functionName <- paste('distance in ',dd$setting,
										' setting (period ', dd$period, ')', sep='')
					dd
				})
			setRate <- do.call(rbind, setRateByPeriod)
			## add the extra column also to the other effects
			rateEffects$setting <- rep("", nrow(rateEffects))
			objEffects$setting <- rep("", nrow(objEffects))
			# get the settings description
			settingsDescription <- describeTheSetting(depvar)
			rateEffects <- rbind(setRate, rateEffects[!rateEffects$basicRate, ])
		}
		else
		{
			settingsDescription <- ""
		}
		list(effects=rbind(rateEffects = rateEffects, objEffects = objEffects),
			starts=starts, settingsDescription=settingsDescription)
	}

	##@behaviornet internal getEffects
	behaviorNet <- function(depvar, varname)
	{
		nodeSet <- attr(depvar,'nodeSet')

		rateEffects <- createEffects("behaviorRate", varname, name=varname,
			groupName=groupName, group=group,
			netType=netType)
		if (observations == 1)
		{
			rateEffects <- rateEffects[-2, ] ## remove the extra period
		}
		else
		{
			## get correct number of rows
			rateEffects <- rbind(duplicateDataFrameRow(rateEffects[2, ],
					observations),
				rateEffects[-c(1, 2), ])
		}

		objEffects <- createEffects("behaviorObjective", varname, name=varname,
			groupName=groupName, group=group,
			netType=netType)

		for (j in seq(along=xx$depvars))
		{
			if (types[j] == "oneMode" &&
				attr(xx$depvars[[j]], "nodeSet") == nodeSet)
			{
				depvarname <- names(xx$depvars)[j]
				if (attr(xx$depvars[[j]], "symmetric"))
				{
					tmpObjEffects <-
						createEffects("behaviorSymmetricObjective",
							varname, depvarname, name=varname,
							groupName=groupName, group=group,
							netType=netType)
					tmpObjEffects2 <- NULL
					for (k in seq(along=xx$depvars))
					{
						if ((j != k) && (types[k] == "oneMode") &&
							(attr(xx$depvars[[k]], "symmetric")) &&
							(attr(xx$depvars[[k]], "nodeSet") == nodeSet))
						{
							othervarname <- names(xx$depvars)[k]
							tmpObjEffects3 <-
								createEffects("behaviorSymSymObjective",
									varname, depvarname, othervarname,
									name=varname,
									groupName=groupName, group=group,
									netType=netType)
							tmpObjEffects3$functionName <-
								paste(tmpObjEffects3$functionName,
									" (", depvarname, ",", othervarname, ")", sep="")
							tmpObjEffects3$effectName <-
								paste(tmpObjEffects3$effectName,
									" (", depvarname, ",", othervarname, ")", sep = "")
							tmpObjEffects2 <- rbind(tmpObjEffects2, tmpObjEffects3)
						}
					}
					tmpRateEffects <-
						createEffects("behaviorSymmetricRate",
							varname, depvarname, name=varname,
							groupName=groupName, group=group,
							netType=netType)
				}
				else # non-symmetric
				{
					tmpObjEffects <-
						createEffects("behaviorOneModeObjective",
							varname, depvarname, name=varname,
							groupName=groupName, group=group,
							netType=netType)
					tmpObjEffects2 <- NULL
					for (k in seq(along=xx$depvars))
					{
						if ((k != j) && (types[k] == "oneMode") &&
							(!attr(xx$depvars[[k]], "symmetric")) &&
							(attr(xx$depvars[[k]], "nodeSet") == nodeSet))
						{
							othervarname <- names(xx$depvars)[k]
							tmpObjEffects3 <-
								createEffects("behaviorOneOneModeObjective",
									varname, depvarname, othervarname,
									name=varname,
									groupName=groupName, group=group,
									netType=netType)
							tmpObjEffects3$functionName <-
								paste(tmpObjEffects3$functionName,
									" (", depvarname, ",", othervarname, ")", sep="")
							tmpObjEffects3$effectName <-
								paste(tmpObjEffects3$effectName,
									" (", depvarname, ",", othervarname, ")", sep = "")
							tmpObjEffects2 <- rbind(tmpObjEffects2, tmpObjEffects3)
						}
						if ((types[k] == "oneMode") &&
							(attr(xx$depvars[[k]], "symmetric")) &&
							(attr(xx$depvars[[k]], "nodeSet") == nodeSet))
						{
							othervarname <- names(xx$depvars)[k]
							tmpObjEffects3 <-
								createEffects("behaviorOneModeSymObjective",
									varname, depvarname, othervarname,
									name=varname,
									groupName=groupName, group=group,
									netType=netType)
							tmpObjEffects3$functionName <-
								paste(tmpObjEffects3$functionName,
									" (", depvarname, ",", othervarname, ")", sep="")
							tmpObjEffects3$effectName <-
								paste(tmpObjEffects3$effectName,
									" (", depvarname, ",", othervarname, ")", sep = "")
							tmpObjEffects2 <- rbind(tmpObjEffects2, tmpObjEffects3)
						}
					}
					tmpRateEffects <-
						createEffects("behaviorOneModeRate",
							varname, depvarname, name=varname,
							groupName=groupName, group=group,
							netType=netType)
				}
				# Now irrespective of (attr(xx$depvars[[j]], "symmetric"))
				for (k in seq(along = xx$dycCovars))
				{
					if (attr(xx$dycCovars[[k]], "type") == "oneMode" &&
						attr(xx$dycCovars[[k]], 'nodeSet')[1] == nodeSet)
					{
						othervarname <- names(xx$dycCovars)[k]
						tmpObjEffects3 <- createEffects("dyadBehaviorNetObjective",
							varname, depvarname, othervarname,
							name=varname, groupName=groupName,
							group=group, netType=netType)
						tmpObjEffects2 <- rbind(tmpObjEffects2, tmpObjEffects3)
					}
				}
				for (k in seq(along = xx$dyvCovars))
				{
					if (attr(xx$dyvCovars[[k]], "type") == "oneMode" &&
						attr(xx$dyvCovars[[k]], 'nodeSet')[1] == nodeSet)
					{
						othervarname <- names(xx$dyvCovars)[k]
						tmpObjEffects3 <- createEffects("dyadBehaviorNetObjective",
							varname, depvarname, othervarname,
							name=varname, groupName=groupName,
							group=group, netType=netType)
						tmpObjEffects2 <- rbind(tmpObjEffects2, tmpObjEffects3)
					}
				}
				if ((nOneModes + nBipartites) > 1) ## add the network name
				{
					tmpObjEffects$functionName <-
						paste(tmpObjEffects$functionName,
							" (", depvarname, ")", sep="")
					tmpObjEffects$effectName <-
						paste(tmpObjEffects$effectName,
							" (", depvarname, ")", sep = "")
					tmpRateEffects$functionName <-
						paste(tmpRateEffects$functionName,
							" (", depvarname, ")", sep="")
					tmpRateEffects$effectName <-
						paste(tmpRateEffects$effectName,
							" (", depvarname, ")", sep = "")
				}
				objEffects <- rbind(objEffects, tmpObjEffects, tmpObjEffects2)
				rateEffects <- rbind(rateEffects, tmpRateEffects)
			}
			if (types[j] == 'bipartite' &&
				(attr(xx$depvars[[j]], 'nodeSet')[1] == nodeSet))
			{
				depvarname <- names(xx$depvars)[j]
				tmpObjEffects <-
					createEffects("behaviorBipartiteObjective",
						varname, depvarname, name=varname,
						groupName=groupName, group=group,
						netType=netType)
				tmpObjEffects2 <- NULL
				for (k in seq(along=xx$depvars))
				{
					if ((k != j) && (types[k] == "bipartite") &&
						(attr(xx$depvars[[k]], "nodeSet")[1] == nodeSet) &&
						(attr(xx$depvars[[j]], "nodeSet")[2] ==
							attr(xx$depvars[[k]], "nodeSet")[2]))
					{
						othervarname <- names(xx$depvars)[k]
						tmpObjEffects3 <-
							createEffects("behaviorBipBipObjective",
								varname, depvarname, othervarname,
								name=varname,
								groupName=groupName, group=group,
								netType=netType)
						tmpObjEffects3$functionName <-
							paste(tmpObjEffects3$functionName,
								" (", depvarname, ",", othervarname, ")", sep="")
						tmpObjEffects3$effectName <-
							paste(tmpObjEffects3$effectName,
								" (", depvarname, ",", othervarname, ")", sep = "")
						tmpObjEffects2 <- rbind(tmpObjEffects2, tmpObjEffects3)
					}
				}
				tmpRateEffects <-
					createEffects("behaviorBipartiteRate",
						varname, depvarname, name=varname,
						groupName=groupName, group=group,
						netType=netType)
				if ((nOneModes + nBipartites) > 1) ## add the network name
				{
					tmpObjEffects$functionName <-
						paste(tmpObjEffects$functionName,
							" (", depvarname, ")", sep="")
					tmpObjEffects$effectName <-
						paste(tmpObjEffects$effectName,
							" (", depvarname, ")", sep = "")
					tmpRateEffects$functionName <-
						paste(tmpRateEffects$functionName,
							" (", depvarname, ")", sep="")
					tmpRateEffects$effectName <-
						paste(tmpRateEffects$effectName,
							" (", depvarname, ")", sep = "")
				}

				objEffects <- rbind(objEffects, tmpObjEffects, tmpObjEffects2)
				rateEffects <- rbind(rateEffects, tmpRateEffects)
			}
		}
		for (j in seq(along = xx$cCovars))
		{
			if (attr(xx$cCovars[[j]], 'nodeSet') == nodeSet)
			{
				tmp <- covBehEff(varname, names(xx$cCovars)[j], nodeSet,
					type='', name=varname)
				objEffects<- rbind(objEffects, tmp$objEff)
				rateEffects<- rbind(rateEffects, tmp$rateEff)
			}
			else
			{
				nodeSet2 <- attr(xx$cCovars[[j]], 'nodeSet')
				tmp <- covBBehEff(varname, names(xx$cCovars)[j], nodeSet2,
					name=varname)
				objEffects<- rbind(objEffects, tmp$objEff)
			}
		}
		for (j in seq(along=xx$depvars))
		{
			if (types[j] %in% c('behavior', 'continuous') &&
				attr(xx$depvars[[j]], 'nodeSet') == nodeSet)
			{
				tmp <- covBehEff(varname, names(xx$depvars)[j], nodeSet, j==i,
					type='Beh', name=varname)
				objEffects<- rbind(objEffects, tmp$objEff)
				rateEffects<- rbind(rateEffects, tmp$rateEff)
			}
		}
		for (j in seq(along=xx$vCovars))
		{
			if (attr(xx$vCovars[[j]], 'nodeSet') == nodeSet)
			{
				tmp <- covBehEff(varname, names(xx$vCovars)[j], nodeSet,
					type='Var', name=varname)
				objEffects<- rbind(objEffects, tmp$objEff)
				rateEffects<- rbind(rateEffects, tmp$rateEff)
			}
			else
			{
				nodeSet2 <- attr(xx$vCovars[[j]], 'nodeSet')
				tmp <- covBBehEff(varname, names(xx$vCovars)[j], nodeSet2,
					name=varname)
				objEffects<- rbind(objEffects, tmp$objEff)
			}
		}
		interaction <- createEffects("unspecifiedBehaviorInteraction",
			varname, name=varname,
			groupName=groupName, group=group,
			netType=netType)
		objEffects <- rbind(objEffects, interaction[rep(1:nrow(interaction), behNintn),])

		## get starting values
		starts <- getBehaviorStartingVals(depvar)
		## set defaults
		if (!(attr(depvar,'allUpOnly') || attr(depvar, 'allDownOnly')))
		{
			objEffects[grepl("linear shape", objEffects$effectName) &
				objEffects$type == 'eval',
			c('include', 'initialValue','untrimmedValue')]  <-
				list(TRUE, starts$tendency, starts$untrimmed)
		}
		else
		{
			objEffects <- objEffects[objEffects$shortName != "linear", ]
		}
		if (attr(depvar, "range") >= 2)
		{
			objEffects[grepl("quadratic shape", objEffects$effectName) &
				objEffects$type == 'eval','include'] <- TRUE
			## no starting value for quadratic effect
		}


		rateEffects[1:observations, 'include'] <- TRUE
		rateEffects[1:noPeriods, 'initialValue'] <- starts$startRate
		rateEffects$basicRate[1:observations] <- TRUE

		## alter the text for endowment and creation names
		objEffects$effectName[objEffects$type == 'endow'] <-
			sub('behavior', 'dec. beh.',
				objEffects$effectName[objEffects$type == 'endow'])
		objEffects$effectName[objEffects$type == 'creation'] <-
			sub('behavior', 'inc. beh.',
				objEffects$effectName[objEffects$type == 'creation'])

		list(effects = rbind(rateEffects = rateEffects,
				objEffects = objEffects), starts=starts)
	}

	##@continuousNet internal getEffects
    continuousNet <- function(depvars, varnames)
    {
        nodeSet <- attr(depvars[[1]],'nodeSet')	## NN: nodeset should be the same for
						##     all continuous depvars

        rateEffects <- createEffects("continuousRate", name="sde",
                                         groupName=groupName, group=group,
                                         netType=netType)
        if (observations == 1)
        {
            rateEffects <- rateEffects[-2, ] ## remove the extra period
        }
        else
        {
            ## get correct number of rows
            rateEffects <- rbind(duplicateDataFrameRow(rateEffects[2, ],
                                                       observations),
                                 rateEffects[-c(1, 2), ])
        }

        objEffects <- fbicEffects <- wEffects <- NULL # general effects, feedback
                                                      # and intercept, wiener
        for (j in seq(along=varnames)) # for all continuous variables
        {
			for (k in seq(along=varnames))
            {
				fbicEffects <- rbind(fbicEffects, createEffects("continuousFeedback",
                                xName = varnames[j], yName = varnames[k],
							    name=varnames[j], groupName=groupName, group=group,
                                netType=netType))
				if (j <= k)
				    wEffects <- rbind(wEffects, createEffects("continuousWiener",
                                xName = varnames[k], yName = varnames[j],
							    name=varnames[k], groupName=groupName, group=group,
                                netType=netType))
			}
            fbicEffects <- rbind(fbicEffects, createEffects("continuousIntercept",
							xName = varnames[j], name = varnames[j],
							groupName=groupName, group=group, netType=netType))

            for (k in seq(along=depvars))
            {
                if (types[k] == "oneMode" &&
                    attr(xx$depvars[[k]], "nodeSet") == nodeSet)
                {
                    depvarname <- names(xx$depvars)[k]

                    tmpObjEffects <-
                            createEffects("continuousOneModeObjective",
                                          varnames[j], depvarname, name=varnames[j],
                                          groupName=groupName, group=group,
                                          netType=netType)
                    }
                    if ((nOneModes) > 1) # add the network name, TODO: same for nBipartites
                    {
                        tmpObjEffects$functionName <-
                            paste(tmpObjEffects$functionName,
                                  " (", depvarname, ")", sep="")
                        tmpObjEffects$effectName <-
                            paste(tmpObjEffects$effectName,
                                  " (", depvarname, ")", sep = "")
                    }

                objEffects <- rbind(objEffects, tmpObjEffects)
            }
            for (k in seq(along = xx$cCovars))
            {
                if (attr(xx$cCovars[[k]], 'nodeSet') == nodeSet)
                {
                    tmp <- covContEff(varnames[j], names(xx$cCovars)[k], nodeSet,
                                     type='', name=varnames[j])
                    objEffects <- rbind(objEffects, tmp$objEff)
                }
            }
            for (k in seq(along=xx$depvars))
            {
                if (types[k] %in% c('behavior', 'continuous') &&
                    attr(xx$depvars[[k]], 'nodeSet') == nodeSet)
                {
                    tmp <- covContEff(varnames[j], names(xx$depvars)[k], nodeSet,
                                     varnames[j] == names(xx$depvars)[k],
                                     type='Beh', name=varnames[j])
                    objEffects <- rbind(objEffects, tmp$objEff)
                }
            }
            for (k in seq(along=xx$vCovars))
            {
                if (attr(xx$vCovars[[k]], 'nodeSet') == nodeSet)
                {
                    tmp <- covContEff(varnames[j], names(xx$vCovars)[k], nodeSet,
                                     type='Var', name=varnames[j])
                    objEffects <- rbind(objEffects, tmp$objEff)
                }
            }
			interaction <- createEffects("unspecifiedContinuousInteraction",
                                     varnames[j], name=varnames[j],
									 groupName=groupName, group=group,
									 netType=netType)

			objEffects <- rbind(objEffects, interaction[rep(1, behNintn),])
		}

		fbicEffects$include <- TRUE

        if (onePeriodSde)
        {
            wEffects$include <- TRUE
            rateEffects$fix[1] <- TRUE
            rateEffects$include[1] <- TRUE
        }
        else
        {
            wEffects$fix[1] <- TRUE
            rateEffects$include[1:observations] <- TRUE
            rateEffects$basicRate[1:observations] <- TRUE
        }

		if (nContinuous > 1)
			wEffects$include <- TRUE

        starts <- getContinuousStartingVals(depvars, onePeriodSde)
        wEffects$initialValue <- starts$startWiener
		rateEffects$initialValue[1:noPeriods] <- starts$startScale
		fbicEffects$initialValue <- starts$startFbic

        list(effects = rbind(rateEffects, wEffects, fbicEffects, objEffects),
             starts = starts)
    }

	##@bipartiteNet internal getEffects
	bipartiteNet <- function(depvar, varname)
	{
		nodeSets <- attr(depvar, 'nodeSet')

		rateEffects <- networkRateEffects(depvar, varname, symmetric=FALSE,
			bipartite=TRUE)

		objEffects <- createEffects("bipartiteObjective", varname,
			name=varname,
			groupName=groupName, group=group,
			netType=netType)

		for (j in seq(along = xx$dycCovars))
		{
			if (attr(xx$dycCovars[[j]], "type") == "bipartite" &&
				all(nodeSets == attr(xx$dycCovars[[j]], 'nodeSet')))
			{
				objEffects <- rbind(objEffects,
					createEffects("dyadBipartiteObjective",
						names(xx$dycCovars)[j],
						name=varname,
						groupName=groupName, group=group,
						netType=netType))
			}
			if (attr(xx$dycCovars[[j]], "type") == "oneMode" &&
				(nodeSets[2] == attr(xx$dycCovars[[j]], 'nodeSet')[1]) &&
				(nodeSets[2] == attr(xx$dycCovars[[j]], 'nodeSet')[2]) )
			{
				objEffects <- rbind(objEffects,
					createEffects("dyadSecondBipartiteObjective",
						names(xx$dycCovars)[j],
						name=varname,
						groupName=groupName, group=group,
						netType=netType))
			}
		}
		for (j in seq(along = xx$dyvCovars))
		{
			if (attr(xx$dyvCovars[[j]], "type") == "bipartite" &&
				all(nodeSets == attr(xx$dyvCovars[[j]], 'nodeSet')))
			{
				objEffects <- rbind(objEffects,
					createEffects("dyadBipartiteObjective",
						names(xx$dyvCovars)[j],
						name=varname,
						groupName=groupName, group=group,
						netType=netType))
			}
			if (attr(xx$dyvCovars[[j]], "type") == "oneMode" &&
				(nodeSets[2] == attr(xx$dyvCovars[[j]], 'nodeSet')[1]) &&
				(nodeSets[2] == attr(xx$dyvCovars[[j]], 'nodeSet')[2]) )
			{
				objEffects <- rbind(objEffects,
					createEffects("dyadSecondBipartiteObjective",
						names(xx$dycCovars)[j],
						name=varname,
						groupName=groupName, group=group,
						netType=netType))
			}
		}
		for (j in seq(along = xx$cCovars))
		{
			covNodeset <- match(attr(xx$cCovars[[j]], "nodeSet"),
				nodeSets)
			if (!is.na(covNodeset))
			{
				tmp <- covarBipartiteEff(names(xx$cCovars)[j],
					attr(xx$cCovars[[j]],
						'poszvar'),
					attr(xx$cCovars[[j]],
						'moreThan2'),
					covNodeset, name=varname)
				objEffects <- rbind(objEffects, tmp$objEff)
				rateEffects <- rbind(rateEffects, tmp$rateEff)
			}
		}
		for (j in seq(along=xx$depvars))
		{
			if (types[j] %in% c('behavior', 'continuous'))
			{
				covNodeset <- match(attr(xx$depvars[[j]], "nodeSet"),
					nodeSets)
				if (!is.na(covNodeset))
				{
					tmp <- covarBipartiteEff(names(xx$depvars)[j],
						poszvar=TRUE,
						attr(xx$depvars[[j]],
							'moreThan2'),
						covNodeset, name=varname)
					objEffects <- rbind(objEffects,tmp$objEff)
					rateEffects <- rbind(rateEffects,tmp$rateEff)
				}
			}
		}
		for (j in seq(along=xx$vCovars))
		{
			covNodeset <- match(attr(xx$vCovars[[j]], "nodeSet"),
				nodeSets)
			if (!is.na(covNodeset))
			{
				tmp <- covarBipartiteEff(names(xx$vCovars)[j],
					attr(xx$vCovars[[j]],
						'poszvar'),
					attr(xx$vCovars[[j]],
						'moreThan2'),
					covNodeset, name=varname)
				objEffects <- rbind(objEffects, tmp$objEff)
				rateEffects <- rbind(rateEffects, tmp$rateEff)
			}
		}
		if (length(xx$cCovars) + length(xx$vCovars) +
			length(xx$dycCovars) + length(xx$dyvCovars) +
			length(types=='behavior') + length(types=='continuous') > 0)
		{
			interaction <- createEffects("unspecifiedNetInteraction",
				name=varname,
				groupName=groupName, group=group,
				netType=netType)
			objEffects <-  rbind(objEffects, interaction[rep(1:nrow(interaction), nintn), ])
		}

		for (j in seq(along=xx$depvars))
		{
			otherName <- names(xx$depvars)[j]
			if (types[j] == 'oneMode' &&
				attr(xx$depvars[[j]], 'nodeSet') == nodeSets[1] )
			{
				if (attr(xx$depvars[[j]], "symmetric"))
				{
					objEffects <-
						rbind(objEffects,
							createEffects("bipartiteSymmetricObjective",
								names(xx$depvars)[[j]],
								name=varname,
								groupName=groupName, group=group,
								netType=netType))
				}
				else
				{
					objEffects <-
						rbind(objEffects,
							createEffects("bipartiteNonSymmetricObjective",
								names(xx$depvars)[[j]],
								name=varname,
								groupName=groupName, group=group,
								netType=netType))
				}
			}
			if (types[j] == 'bipartite' &&
				(attr(xx$depvars[[j]], 'nodeSet')[1] == nodeSets[1]) &&
				varname != otherName)
			{
				objEffects <-
					rbind(objEffects,
						createEffects("bipartiteBipartiteObjective",
							names(xx$depvars)[[j]],
							name=varname,
							groupName=groupName, group=group,
							netType=netType))
			}
			if (!(types[j] %in% c('behavior', 'continuous'))) # other networks, irrespective of onemode/symmetric/twomode
			{
				for (k in seq(along=xx$cCovars))
				{
					if (attr(xx$cCovars[[k]], 'nodeSet') == nodeSets[1])
					{
						objEffects <-
							rbind(objEffects,
								createEffects("covarABipNetObjective", otherName,
									names(xx$cCovars)[k], name=varname,
									groupName=groupName, group=group,
									netType=netType))
					}
				}
				for (k in seq(along=xx$vCovars))
				{
					if (attr(xx$vCovars[[k]], 'nodeSet') == nodeSets[1])
					{
						objEffects <-
							rbind(objEffects,
								createEffects("covarABipNetObjective", otherName,
									names(xx$vCovars)[k], name=varname,
									groupName=groupName, group=group,
									netType=netType))
					}
				}
				for (k in seq(along=xx$depvars))
				{
					if (types[k] %in% c('behavior', 'continuous') &&
						attr(xx$depvars[[k]], 'nodeSet') == nodeSets[1])
					{
						objEffects <-
							rbind(objEffects,
								createEffects("covarABipNetObjective", otherName,
									names(xx$depvars)[k], name=varname,
									groupName=groupName, group=group,
									netType=netType))
					}
				}
			}
		}
		if ((nOneModes + nBipartites) > 1) ## add the network name
		{
			objEffects$functionName <- paste(varname, ': ',
				objEffects$functionName, sep = '')
			objEffects$effectName <- paste(varname, ': ',
				objEffects$effectName, sep = '')
		}

		objEffects$functionName[objEffects$type == 'endow'] <-
			paste('Lost ties:',
				objEffects$functionName[objEffects$type =='endow'])
		objEffects$functionName[objEffects$type == 'creation'] <-
			paste('New ties:',
				objEffects$functionName[objEffects$type =='creation'])

		## get starting values
		starts <- getBipartiteStartingVals(depvar)
		##set defaults
		rateEffects[1:observations, 'include'] <- TRUE
		rateEffects[1:noPeriods, 'initialValue'] <- starts$startRate
		rateEffects$basicRate[1:observations] <- TRUE

		if (!(attr(depvar,'allUpOnly') || attr(depvar, 'allDownOnly')))
		{
			objEffects[objEffects$shortName =='density' &
				objEffects$type == 'eval',
			c('include', 'initialValue', 'untrimmedValue')] <-
				list(TRUE, starts$degree, starts$untrimmed)
		}
		else
		{
			objEffects <-
				objEffects[!objEffects$shortName == "density", ]
		}

		rateEffects$basicRate[1:observations] <- TRUE

		list(effects=rbind(rateEffects = rateEffects, objEffects = objEffects),
			starts=starts)
	}

	##@covarOneModeEff internal getEffects
	covarOneModeEff<- function(covarname, poszvar, moreThan2, symmetric, constant,
		name)
	{
		if (symmetric)
		{
			covObjEffects <- createEffects("covarSymmetricObjective", covarname,
				name=name,
				groupName=groupName, group=group,
				netType=netType)
			covRateEffects <- createEffects("covarSymmetricRate", covarname,
				name=name,
				groupName=groupName, group=group,
				netType=netType)
		}
		else
		{
			covObjEffects <- createEffects("covarNonSymmetricObjective",
				covarname, name=name,
				groupName=groupName, group=group,
				netType=netType)
			covRateEffects <- createEffects("covarNonSymmetricRate", covarname,
				name=name,
				groupName=groupName, group=group,
				netType=netType)
		}

		if (constant)
		{
			covObjEffects <-
				covObjEffects[!(covObjEffects$shortName %in% c("avGroupEgoX")),]
		}

# these lines tentatively dropped version 1.2-5
#		if (!tr & (!poszvar))  # not (positive variance of z, or any z missing)
#		{
#			if (symmetric)
#			{
#				covObjEffects <-
#					covObjEffects[covObjEffects$shortName %in%
#					c("altX", "altSqX"), ]
#			}
#			else
#			{
#				covObjEffects <-
#					covObjEffects[covObjEffects$shortName %in%
#					c("egoX", "egoSqX"), ]
#			}
#		}
#		if (!tr & (!moreThan2))
#		{
#			covObjEffects <-
#				covObjEffects[!covObjEffects$shortName %in%
#				c("altSqX", "egoPlusAltSqX"), ]
#		}

		list(objEff=covObjEffects, rateEff=covRateEffects)
	}
	##@covarBipartiteEff internal getEffects
	covarBipartiteEff<- function(covarname, poszvar, moreThan2, nodesetNbr,
		name)
	{
		covRateEffects <- NULL
		if (nodesetNbr == 1)
		{
			covObjEffects <-
				createEffects("covarBipartiteObjective", covarname,
					name=varname,
					groupName=groupName, group=group,
					netType=netType)
			# restrict to covariates on first node set
			covObjEffects <-
				covObjEffects[covObjEffects$shortName %in%
				c("egoX", "egoSqX", "egoLThresholdX", "egoRThresholdX",
					"degAbsDiffX", "degPosDiffX", "degNegDiffX",
					"altInDist2", "totInDist2",
					"simEgoInDist2", "sameXInPop", "diffXInPop",
					"sameXCycle4", "inPopX", "inActX", "avGroupEgoX"), ]
			covRateEffects <- createEffects("covarBipartiteRate", covarname,
				name=varname,
				groupName=groupName, group=group,
				netType=netType)
		}
		else if (tr | poszvar) # positive variance of z, or any z missing
		{
			covObjEffects <- createEffects("covarBipartiteObjective", covarname,
				name=varname,
				groupName=groupName, group=group,
				netType=netType)
			# restrict to covariates on second node set
			covObjEffects <-
				covObjEffects[covObjEffects$shortName %in%
				c("altX", "altSqX",  "altLThresholdX", "altRThresholdX",
					"homXOutAct", "altXOutAct",
					"inActX", "outActX"), ]
			if (!tr & (!moreThan2))
			{
				covObjEffects <-
					covObjEffects[!covObjEffects$shortName %in% c("altSqX"), ]
			}
		}
		else
		{
			covObjEffects <- NULL
		}

		list(objEff=covObjEffects, rateEff=covRateEffects)
	}
	##@covBehEff internal getEffects
	## beware: this function does lots of things!
	covBehEff<- function(varname, covarname, nodeSet, same=FALSE,
		## same indicates that varname and covarname are
		## the same: just one rate effect required
		## type is no longer used
		type=c('', 'Var', 'Beh'), name)
	{
		covObjEffects <-  NULL
		if (!same)
		{
			covObjEffects <- createEffects("covarBehaviorObjective", varname,
				covarname, name=name,
				groupName=groupName, group=group,
				netType=netType)
		}

		covRateEffects <- createEffects("covarBehaviorRate", varname, covarname,
			name=name,
			groupName=groupName, group=group,
			netType=netType)
		## if we have a real covariate, need to add other effects for every
		## appropriate network
		if (!same)
		{
			for (j in seq(along=xx$depvars))
			{
				if (types[j] == 'oneMode' &&
					attr(xx$depvars[[j]], 'nodeSet') == nodeSet)
				{
					newEffects <-
						createEffects("covarBehaviorNetObjective", varname,
							covarname, names(xx$depvars)[j],
							groupName=groupName, group=group,
							netType=netType, name=name)

					covObjEffects <- rbind(covObjEffects, newEffects)
					if (attr(xx$depvars[[j]], "symmetric"))
					{
						covBehRateEffects <-
							createEffects("covarBehaviorSymmetricRate", varname,
								yName=names(xx$depvars)[j],
								zName=covarname,
								groupName=groupName, group=group,
								netType=netType, name=name)
					}
					else
					{
						covBehRateEffects <-
							createEffects("covarBehaviorOneModeRate", varname,
								yName=names(xx$depvars)[j],
								zName=covarname,
								groupName=groupName, group=group,
								netType=netType, name=name)
					}
					covRateEffects <- rbind(covRateEffects, covBehRateEffects)
				}
				if ((types[j] == "bipartite" &&
						attr(xx$depvars[[j]], 'nodeSet')[1] == nodeSet))
				{
					newEffects <-
						createEffects("covarABehaviorBipartiteObjective", varname,
							covarname, names(xx$depvars)[j],
							groupName=groupName, group=group,
							netType=netType, name=name)
					covObjEffects <- rbind(covObjEffects, newEffects)
				}
			}
		}

		list(objEff=covObjEffects, rateEff=covRateEffects)
	}
	##@covBBehEff internal getEffects
	covBBehEff<- function(varname, covarname, nodeSetB, name)
	{
		covObjEffects <-  NULL
		for (j in seq(along=xx$depvars))
		{
			if ((types[j] == "bipartite" &&
					attr(xx$depvars[[j]], 'nodeSet')[2] == nodeSetB))
			{
				newEffects <-
					createEffects("covarBBehaviorBipartiteObjective", varname,
						covarname, names(xx$depvars)[j],
						groupName=groupName, group=group,
						netType=netType, name=name)
				covObjEffects <- rbind(covObjEffects, newEffects)
			}
		}
		list(objEff=covObjEffects)
	}

	##@covarNetNetEff internal getEffects
	covarNetNetEff<- function(othernetname, covarname, nodeSetsj, nodeSetk, poszvar, name)
	{
		if (tr | poszvar) # positive variance of z, or any z missing
		{
			if (length(nodeSetsj) <= 1) ## second network onemode
			{
				processA <- TRUE
				processB <- TRUE
				objEffectsOneMode <- createEffects("covarNetNetObjective", othernetname,
					covarname, name=name,
					groupName=groupName, group=group,
					netType=netType)
			}
			else ## second network twomode
			{
				processA <- (nodeSetsj[[1]] == nodeSetk)
				processB <- (nodeSetsj[[2]] == nodeSetk)
				objEffectsOneMode <- NULL
			}
			if (processA)
			{
				objEffectsA <- createEffects("covarANetNetObjective", othernetname,
					covarname, name=name,
					groupName=groupName, group=group,
					netType=netType)
			}
			else
			{
				objEffectsA <- NULL
			}
			if (processB)
			{
				objEffectsB <- createEffects("covarBNetNetObjective", othernetname,
					covarname, name=name,
					groupName=groupName, group=group,
					netType=netType)
			}
			else
			{
				objEffectsB <- NULL
			}
			objEffectsAB <- createEffects("covarABNetNetObjective", othernetname,
				covarname, name=name,
				groupName=groupName, group=group,
				netType=netType)
			objEffects <- rbind(objEffectsOneMode, objEffectsAB, objEffectsA, objEffectsB)
		}
		else
		{
			objEffects <- NULL
		}
		objEffects
	}

	##@covContEff internal getEffects
	covContEff <- function(varname, covarname, nodeSet, same=FALSE,
                         ## same indicates that varname and covarname are
                         ## the same: just one rate effect required
                         ## type is no longer used
                           type=c('', 'Var', 'Beh'), name)
	{
        	covObjEffects <-  NULL
		if (!same)
		{
			covObjEffects <- createEffects("covarContinuousObjective", varname,
                                    covarname, name=name,
                                    groupName=groupName, group=group,
                                         netType=netType)
		}
		list(objEff=covObjEffects)
	}

	###################################
	## start of function getEffects
	##################################
	if (getDocumentation)
	{
		tt <- getInternals()
		return(tt)
	}
	if (!inherits(x, 'sienaGroup') && !inherits(x, 'siena'))
	{
		stop('Not a valid siena data object or group')
	}
	if (inherits(x, 'sienaGroup'))
	{
		groupx <- TRUE
	}
	else
	{
		groupx <- FALSE
	}
	tr <- TRUE; # supersedes restrictions by poszvar and moreThan2; 1.1-306
	## validate the object?
	## find the total number of periods to be processed = local var observations
	## then process the first or only data object. Fill in starting values
	## for other periods from the other objects, if any.
	if (groupx)
	{
		groupNames <- names(x)
		observations <- attr(x, 'observations')
		xx <- x[[1]]
		periodNos <- attr(x, "periodNos")
	}
	else
	{
		groupNames <- 'Group1'
		observations <- x$observations - 1
		xx <- x
		periodNos <- 1:observations
	}
	n <- length(xx$depvars)
	types <- sapply(xx$depvars, function(x)attr(x, 'type'))
	#sparses <- sapply(xx$depvars, function(x)attr(x, 'sparse'))
	nOneModes <- sum(types == 'oneMode')
	#nBehaviors <- sum(types == 'behavior')
	nBipartites <- sum(types =='bipartite')
	nContinuous <- sum(types == 'continuous')
	effects <- vector('list',n+1) 			# n+1 th place for all sde parameters
	#nodeSetNames <- sapply(xx$nodeSets, function(x)attr(x, 'nodeSetName'))
	names(effects) <- names(xx$depvars) 	# n+1 th place has no name

	if (onePeriodSde && xx$observations > 2)
		stop('onePeriodSde only possible in case of 2 observations')

	if (onePeriodSde && groupx)
		stop('onePeriodSde not possible in combination with multi-group')

	for (i in 1:(n-nContinuous))
	{
		varname<- names(xx$depvars)[i]
		groupName <- groupNames[1]
		group <- 1
		noPeriods <- xx$observations - 1 ## xx is a single object
		depvar <- xx$depvars[[i]]
		if (groupx)
		{
			netnamesub <- match(varname, attr(x, 'netnames'))
			attr(depvar, 'allUpOnly') <- attr(x, 'allUpOnly')[netnamesub]
			attr(depvar, 'anyUpOnly') <- attr(x, 'anyUpOnly')[netnamesub]
			attr(depvar, 'anyDownOnly') <- attr(x, 'anyDownOnly')[netnamesub]
			attr(depvar, 'allDownOnly') <- attr(x, 'allDownOnly')[netnamesub]
			if (types[i] == 'oneMode')
				attr(depvar, 'symmetric') <- attr(x, 'symmetric')[netnamesub]
		}
		else
		{
			attr(depvar, 'allUpOnly') <- all(attr(depvar, 'uponly'))
			attr(depvar, 'anyUpOnly') <- any(attr(depvar, 'uponly'))
			attr(depvar, 'anyDownOnly') <- any(attr(depvar, 'downonly'))
			attr(depvar, 'allDownOnly') <- all(attr(depvar, 'downonly'))
		}
		switch(types[i],
			behavior =
			{
				netType <- "behavior"
				tmp <- behaviorNet(depvar, varname)
				effects[[i]] <- tmp$effects
				attr(effects[[i]], 'starts') <- tmp$starts
				attr(effects[[i]], 'settings') <- ''
			},
			oneMode =
			{
				netType <- "oneMode"
				tmp <- oneModeNet(depvar, varname)
				effects[[i]] <- tmp$effects
				attr(effects[[i]], 'starts') <- tmp$starts
				attr(effects[[i]], 'settings') <- tmp$settingsDescription
			},
			bipartite =
			{
				netType <- "bipartite"
				tmp <- bipartiteNet(depvar, varname)
				effects[[i]] <- tmp$effects
				attr(effects[[i]], 'starts') <- tmp$starts
				attr(effects[[i]], 'settings') <- ''
			},
			stop('error type'))
	}
	settingsList <- lapply(effects, function(ef){attr(ef,'settings')})
	if (nContinuous > 0)
	{
		groupName <- groupNames[1]
        	group <- 1
	        noPeriods <- xx$observations - 1
		netType <- "continuous"
		contIndices <- (n-nContinuous+1):n ## indicates continuous depvars
		varnames <- names(xx$depvars)[contIndices]
        	depvars <- xx$depvars[contIndices]
		tmp <- continuousNet(depvars,varnames)
		effects[[n+1]] <- tmp$effects
		attr(effects[[n+1]], 'starts') <- tmp$starts
		for (i in contIndices)
		{
			# all the continuous variable specific effects are currently
			# also part of effects[[n+1]]
		}
	}
	## add starting values for the other objects
	if (groupx && length(x) > 1)
	{
		period <-  xx$observations ##periods used so far

		for (group in 2:length(x))
		{
			xx <- x[[group]]
			n <- length(xx$depvars)
			types <- sapply(xx$depvars, function(x)attr(x, 'type'))
			noPeriods <- xx$observations - 1
			nContinuous <- sum(types == 'continuous')

			for (i in 1:(n-nContinuous))
			{
				varname<- names(xx$depvars)[i]
				depvar <- xx$depvars[[i]]
				netnamesub <- match(varname, attr(x, 'netnames'))
				if (types[i] == 'oneMode')
				{
					attr(depvar, 'symmetric') <- attr(x, 'symmetric')[netnamesub]
				}
					switch(types[i],
						behavior =
						{
							starts <-  getBehaviorStartingVals(depvar)
						## first for the rate parameters
							## find the appropriate set of effects
							eff <- match(varname, names(effects))
							if (is.na(eff))
						{
								stop("depvars don't match")
						}
							effectname <- paste('rate ', varname,' (period ',
							period + 1:noPeriods, ')',sep='')
						use <- effects[[eff]]$effectName %in% effectname
							effects[[eff]][use, c('include','initialValue',
								'groupName', 'group', 'period')] <-
									list(TRUE, starts$startRate,
									groupNames[group], group, 1:noPeriods)
									## now sort out the tendency and update the
									## attribute on the effects list:
						newdif <- c(starts$dif, attr(effects[[eff]], "starts")$dif)
									meandif <- mean(newdif, na.rm=TRUE)
									vardif <- var(as.vector(newdif), na.rm=TRUE)
									if (meandif < 0.9 * vardif)
									{
										tendency <- 0.5 * log((meandif + vardif)/
											(vardif - meandif))
									}
									else
									{
										tendency <- meandif / (vardif + 1)
									}
									untrimmed <- tendency
									tendency <- ifelse(tendency < -3.0, -3.0,
										ifelse(tendency > 3/0, 3.0, tendency))
									use <- (effects[[eff]]$shortName == "linear" &
										effects[[eff]]$type == "eval")
									effects[[eff]][use, c("include", "initialValue",
							"untrimmedValue")] <- list(TRUE, tendency, untrimmed)
											attr(effects[[eff]], 'starts')$dif <- newdif
						},
						oneMode =
						{
							starts <- getNetworkStartingVals(depvar)
						## first for the rate parameters
							## find the appropriate set of effects
							eff <- match(varname, names(effects))
							if (is.na(eff))
							{
								stop("depvars don't match")
							}
						thisSettingDescription <- paste(unlist(describeTheSetting(depvar)), collapse=" ")
						if (thisSettingDescription !=
									paste(unlist(settingsList[[varname]]), collapse=" "))
						{
							stop(paste('setting definitions do not match (group ',
								group,', variable ', varname,')', sep=''))

						}
						if (thisSettingDescription == "")
						{
							effectname <- paste('constant ', varname, ' rate (period ',
								period + 1:noPeriods,')', sep='')
							use <- effects[[eff]]$effectName %in% effectname
							effects[[eff]][use, c('include', 'initialValue',
								"groupName", "group", "period")] <-
									list(TRUE, starts$startRate,
									groupNames[group], group, 1:noPeriods)
						}
						else
						{
							effectname <- paste('setting rate ', varname, ' (period ',
								period + 1:noPeriods,')', sep='')
							use <- grep(effectname, effects[[eff]]$effectName, fixed=TRUE)
							effects[[eff]][use, c('include', 'initialValue',
								"groupName", "group", "period")] <-
									list(TRUE, starts$startRateSett,
										groupNames[group], group, 1:noPeriods)
						}
									## now sort out the degree and
									## update the attribute on the effects list
									oldstarts <- attr(effects[[eff]], "starts")
									alpha <- c(oldstarts$alpha, starts$alpha)
									prec <- c(oldstarts$prec, starts$prec)
									degree <- sum(alpha * prec) / sum(prec)
									untrimmed <- degree
									degree <- ifelse (degree < -3, -3,
										ifelse(degree > 3, 3, degree))
									attr(effects[[eff]], "starts")$alpha <- alpha
									attr(effects[[eff]], "starts")$prec <-  prec
									if (attr(depvar, 'symmetric'))
									{
										effects[[eff]][effects[[eff]]$shortName ==
											'density' &
											effects[[eff]]$type == 'eval',
										c('initialValue','untrimmedValue')] <-
											list(degree, untrimmed)
									}
									else
									{
										if (!(attr(x,'anyUpOnly') || attr(x, 'anyDownOnly')))
										{
											effects[[eff]][effects[[eff]]$shortName ==
												'density' &
												effects[[eff]]$type == 'eval',
											c('initialValue',
												"untrimmedValue")] <-
													list(degree, untrimmed)
										}
									}
						},
						bipartite =
						{
							starts <- getBipartiteStartingVals(depvar)
						## first for the rate parameters
							## find the appropriate set of effects
							eff <- match(varname, names(effects))
							if (is.na(eff))
							{
								stop("depvars don't match")
							}
							effectname <- paste('constant ', varname,
								' rate (period ',
								period + 1:noPeriods,')', sep='')
							use <- effects[[eff]]$effectName %in% effectname
							effects[[eff]][use, c('include', 'initialValue',
								'groupName', 'group',
								'period')] <-
									list(TRUE, starts$startRate,
										groupNames[group],
										group, 1:noPeriods)
									## now sort out the degree and
									## update the attribute on the effects list
									oldstarts <- attr(effects[[eff]], "starts")
									alpha <- c(oldstarts$alpha, starts$alpha)
									prec <- c(oldstarts$prec, starts$prec)
									degree <- sum(alpha * prec) / sum(prec)
									untrimmed <- degree
									degree <- ifelse (degree < -3, -3,
										ifelse(degree > 3, 3, degree))
									attr(effects[[eff]], "starts")$alpha <- alpha
									attr(effects[[eff]], "starts")$prec <-  prec
									if (!(attr(x,'anyUpOnly') || attr(x, 'anyDownOnly')))
									{
										effects[[eff]][effects[[eff]]$shortName ==
											'density' &
											effects[[eff]]$type == 'eval',
										c('initialValue',
											"untrimmedValue")] <-
												list(degree, untrimmed)
									}
						},
						stop('error type'))
			}
			if (nContinuous > 0)
			{
				contIndices <- (n-nContinuous+1):n ## indicates continuous depvars
				varnames <- names(xx$depvars)[contIndices]
				depvars <- xx$depvars[contIndices]
				starts <- getContinuousStartingVals(depvars, onePeriodSde = FALSE)
				## At this point, onePeriodSde is always FALSE, as the combination
				## of multi-group and onePeriodSde = TRUE is not feasible
				for (i in 1:nContinuous) ## check whether all continous vars match
				{
					eff <- match(varnames[i], names(effects))
					if (is.na(eff))
					{
						stop("depvars don't match")
					}
				}
				effectname <- paste('scale parameter (period ',
									period + 1:noPeriods, ')', sep='')
				use <- effects[[n+1]]$effectName %in% effectname
				effects[[n+1]][use, c('include','initialValue',
										'groupName', 'group', 'period')] <-
											 list(TRUE, starts$startScale,
												  groupNames[group], group,
												  1:noPeriods)
			}

			period <-  period + xx$observations ##periods used so far
		}
	}
	effects <- do.call(rbind, effects)
	attr(effects, "starts") <- NULL
	effects <- cbind(effects, effectNumber=1:nrow(effects))
	attr(effects, "settings") <- settingsList
	cl <- class(effects)
	if (groupx)
		class(effects) <- c('sienaGroupEffects','sienaEffects', cl)
	else
		class(effects) <- c('sienaEffects', cl)
	myrownames <- paste(sapply(strsplit(row.names(effects), ".", fixed=TRUE),
			function(x)paste(x[1:2], collapse='.')),
		effects$type, sep='.')
	mytab <- table(myrownames)
	for (i in 1:length(mytab))
	{
		myrownames[myrownames == names(mytab)[i]] <-
			paste(myrownames[myrownames == names(mytab)[i]], 1:mytab[i], sep=".")
	}
	myrownames <- sub("Effects", "", myrownames)
	myrownames <- sub("rate.rate", "rate", myrownames)
	rownames(effects) <- myrownames
	effects
}
##@getBehaviorStartingVals DataCreate
getBehaviorStartingVals <- function(depvar)
{
	drange <- attr(depvar, 'range')
	depvar <- depvar[,1,]
	nactors <- nrow(depvar)
	## rate
	if (round(drange) == 2)
	{
		rr <- range(depvar, na.rm=TRUE)
		dif <- t(diff(t(depvar))) ##calculate column differences
		tmp <- sapply(1:ncol(dif), function(x, y, z){
			mintab <- c("FALSE"=0, "TRUE" = 0)
			mintab[1] <- 1 + sum(dif[z[, x] == rr[1], x] == 0, na.rm=TRUE)
			mintab[2] <- 1 + sum(dif[z[, x] == rr[1], x] > 0, na.rm=TRUE)
			maxtab <- c("FALSE"=0, "TRUE" = 0)
			maxtab[1] <- 1 + sum(dif[z[, x] == rr[2], x] == 0, na.rm=TRUE)
			maxtab[2] <- 1 + sum(dif[z[, x] == rr[2], x] < 0, na.rm=TRUE)
			val <- mintab[2] / sum(mintab) + maxtab[2] / sum(maxtab)
			if (val > 0.9) val <- 0.5
			c(-log(1 - val), mintab, maxtab)
		}, z = depvar, y = dif)
		startRate <- tmp[1, ]
		##tendency
		tmp <- rowSums(tmp[-1, , drop=FALSE])
		tendency <- log((tmp[2] * (tmp[3] + tmp[4])) /
			(tmp[4] * (tmp[1] + tmp[2])))
		untrimmed <- tendency
		tendency <- ifelse(tendency < -2, -2, ifelse (tendency > 2, 2,
				tendency))
	}
	else
	{
		dif <- t(round(diff(t(depvar)))) ## get column differences
		absdif <- colSums(abs(dif), na.rm=TRUE) ##sum abs over columns
		# sqrdif <- colSums(dif*dif, na.rm=TRUE) ##sum sqr over columns
		# dif <- colSums(dif,na.rm=TRUE) ## sum over columns
		# startRate <- nactors / (nactors-1) * (sqrdif/nactors -
		#   dif * dif/nactors/nactors)
		startRate <- apply(dif, 2, var, na.rm=TRUE)
		startRate <- pmax(startRate, 0.1 + absdif/nactors)
		## tendency
		##  dif<- sum(dif/nactors,na.rm=TRUE)
		##  dif2<- nactors/(nactors-1)*(sum(sqrdif,na.rm=TRUE)/nactors-dif^2)
		dif1 <- mean(dif, na.rm=TRUE)
		dif2 <- var(as.vector(dif), na.rm=TRUE)
		if (is.na(dif2))
		{
			stop('Too little information in the dependent variable')
		}
		if (abs(dif1) < 0.9 * dif2)
		{
			tendency <- 0.5 * sign(dif1) * log((abs(dif1)+dif2)/(dif2-abs(dif1)))
		}
		else
		{
			tendency <- dif1/(dif2+1)
		}
	}
	untrimmed <- tendency
	tendency <- ifelse(tendency < -3, -3,
		ifelse (tendency > 3, 3, tendency))
	list(startRate=startRate, tendency=tendency, untrimmed = untrimmed, dif=dif)
}
##@getContinuousStartingVals DataCreate
getContinuousStartingVals <- function(depvars, onePeriodSde)
{
	depvar <- depvars[[1]][,1,]       # first continuous dependent variable
	nActors <- nrow(depvar)           # no. of actors
	nPeriods <- ncol(depvar) - 1	  # no. of periods
	nCont <- length(depvars)          # no. of continuous variables

	if(nCont > 1)
		stop("getContinuousStartingVals not defined for more than one continuous variable")

	# determine SDE parameters by regressing each observation on the
	# preceding observation, using the Bergstrom formula
	# note: wiener process parameter G is set to 1 for identifiability
	LL <- function(theta)
	{
		a <- theta[1]
		b0 <- theta[2]
		tau <- theta[3:(2+nPeriods)]
		minlogliks <- rep(0, nPeriods)

		for(i in 1:nPeriods)
		{	# actors present at time i and i+1
			act <- which(!is.na(depvar[,i] + depvar[,i+1]))
			R <- depvar[act,i+1] - exp(a * tau[i]) * depvar[act,i] -
				(exp(a*tau[i]) - 1) * b0 / a
			R <- suppressWarnings(dnorm(R, 0, sqrt((exp(2 * a * tau[i]) - 1) / (2*a)),
				log = TRUE))
			minlogliks[i] <- -sum(R)
		}
		sum(minlogliks)
    	}

	LLonePeriodSde <- function(theta)
	{
		a <- theta[1]
		b0 <- theta[2]
		g <- theta[3]
		# actors present at time 1 and 2
		act <- which(!is.na(depvar[,1] + depvar[,2]))
		R <- depvar[act,2] - exp(a) * depvar[act,1] -
			(exp(a) - 1) * b0 / a
		R <- suppressWarnings(dnorm(R, 0, sqrt((exp(2 * a) - 1) * g^2 / (2*a)),
			log = TRUE))
		-sum(R)
	}

	if (onePeriodSde)
		fit <- optim(theta <- c(-0.5, 3, 1), LLonePeriodSde, hessian = TRUE)
	else
		fit <- optim(theta <- c(-0.5, 3, rep(1, nPeriods)), LL, hessian = TRUE,
				method = "BFGS", control = list(maxit = 10000))

	if (fit$convergence != 0)
		stop("Initial SDE parameter values could not be determined.")

	## NN: still prints the SDE initial parameters and standard errors when getEffects
	##     is called
	cat("SDE init parameters: ", fit$par, '\n')
	cat("SDE par stand errors:", sqrt(diag(solve(fit$hessian))), '\n')

	if (onePeriodSde) # return: tau, g, a, b0
		list(startScale = 1, startWiener = fit$par[3], startFbic = fit$par[1:2])
	else
		list(startScale = fit$par[3:(2+nPeriods)], startWiener = 1,
			 startFbic = fit$par[1:2])
}
##@getNetworkStartingVals DataCreate
getNetworkStartingVals <- function(depvar)
{
	noPeriods <- attr(depvar, "netdims")[3] - 1
	##rate
	##get distance and number of valid links
	if (!attr(depvar,'sparse'))
	{
		nactors <- nrow(depvar)
		use <- !is.na(depvar) & (depvar == 10 | depvar == 11)
		depvar[use] <- depvar[use] - 10  ## remove structural values
		tmp <- sapply(1:noPeriods, function(x, z){
			diag(z[ , , x]) <- NA
			diag(z[, , x + 1]) <- NA
			matdiff <- sum(z[, , x + 1] != z[, , x], na.rm=TRUE)
			# matchange0 <- table(z[, , x + 1], z[, , x])
			# Changed to protect against zero rows or columns
			mc00 <- sum((1 - z[ , , x+1])*(1 - z[ , , x]), na.rm=TRUE)
			mc01 <- sum(z[ , , x+1]*(1 - z[ , , x]), na.rm=TRUE)
			mc10 <- sum((1 - z[ , , x+1])*z[ , , x], na.rm=TRUE)
			mc11 <- sum(z[ , , x+1]*z[ , , x], na.rm=TRUE)
			matchange <- matrix(c(mc00, mc01, mc10, mc11), 2, 2)
			#cat(matchange0,'\n',matchange,'\n')
			matcnt <- nactors * nactors -
				sum(is.na(z[, , x + 1]) | is.na(z[, , x]))
			tmp <- c(matcnt=matcnt, matdiff=matdiff, matchange=matchange)
			names(tmp) <- c("matcnt", "matdiff", "matchangeFrom0To0",
				"matchangeFrom0To1",
				"matchangeFrom1To0", "matchangeFrom1To1")
			tmp
		}, z=depvar)
	}
	else
	{
		nactors <- nrow(depvar[[1]])
		matdiff<- rep(NA, noPeriods)
		matcnt<- rep(NA, noPeriods)
		matchange<- matrix(NA, nrow=4, ncol=noPeriods)
		for (i in 1: noPeriods)
		{
			mymat1 <- depvar[[i]]
			mymat2 <- depvar[[i+1]]
			use <- mymat1@x %in% c(10, 11)
			mymat1@x[use] <- mymat1@x[use] - 10
			use <- mymat2@x %in% c(10, 11)
			mymat2@x[use] <- mymat2@x[use] - 10
			mymat1 <- drop0(mymat1)
			mymat2 <- drop0(mymat2)
			diag(mymat1) <- NA
			diag(mymat2) <- NA
			mydif <- mymat2 - mymat1
			matdiff[i] <- sum(abs(mydif), na.rm=TRUE)
			tmp <- table(mydif@x)
			dummy <- factor(NA, levels=c(-1,0,1))
			dummy <- table(dummy)
			dummy[names(tmp)] <- tmp
			tmp <- dummy
			tmp00 <- nactors * nactors - length(mydif@x)
			tmp <- c(tmp00, tmp[c(3, 1, 2)])
			matchange[, i] <- tmp
			matcnt[i] <- sum(tmp)
		}
		matchange <- data.frame(matchange)
		tmp <- as.matrix(rbind(matcnt=matcnt, matdiff=matdiff,
				matchange=matchange))
		row.names(tmp) <- c("matcnt", "matdiff", "matchangeFrom0To0",
			"matchangeFrom0To1",
			"matchangeFrom1To0", "matchangeFrom1To1")
	}
	distance <- attr(depvar, "distance" )
	if (attr(depvar,'symmetric'))
		startRate <- nactors * (0.2 + distance)/(tmp['matcnt',] %/% 2 +1)
	else
		startRate <- nactors * (0.2 + 2 * distance)/(tmp['matcnt',] + 1)
	startRate <- pmax(0.1, startRate)
	startRate <- pmin(100, startRate)
	if (length(attr(depvar,'settings')) >= 2) # also see above with 0.75
	{
		startR <- rep(startRate, length(attr(depvar,'settings')))
#		setR <- rep(sapply(attr(depvar,'settings'), function(x){x$id}), length(startRate))
		setR <- rep(attr(depvar,'settings'), length(startRate))
		startR[setR=='primary'] <- 0.75 * startR[setR=='primary']
		startR[setR!='primary'] <- 0.5
		startRateSett <- startR
	}
	else
	{
		startRateSett <- NA
	}
	##degree
	matchange<- as.matrix(tmp[grep("matchange", rownames(tmp)),,drop=FALSE])
	if (attr(depvar,'symmetric'))
	{
		matchange <- matchange %/% 2
	}
	p01 <- ifelse (matchange["matchangeFrom0To0", ] +
		matchange["matchangeFrom0To1", ] >=1,
	matchange["matchangeFrom0To1", ] /
		(matchange["matchangeFrom0To0", ] +
			matchange["matchangeFrom0To1", ]), 0.5)
	p10 <- ifelse (matchange["matchangeFrom1To0", ]
		+ matchange["matchangeFrom1To1", ] >=1,
		matchange["matchangeFrom1To0", ] /
			(matchange["matchangeFrom1To0", ] +
				matchange["matchangeFrom1To1", ]), 0.5)
	p01 <- pmax(0.02, p01)
	p10 <- pmax(0.02, p10)
	p01 <- pmin(0.98, p01)
	p10 <- pmin(0.98, p10)
	alpha <- 0.5 * log(p01 / p10)
	p00 <- ifelse (matchange["matchangeFrom0To0", ] +
		matchange["matchangeFrom0To1", ] >=1,
	matchange["matchangeFrom0To0", ] /
		(matchange["matchangeFrom0To0", ] +
			matchange["matchangeFrom0To1", ]), 0.0)
	p11 <- ifelse (matchange["matchangeFrom1To0", ]
		+ matchange["matchangeFrom1To1", ] >=1,
		matchange["matchangeFrom1To1", ] /
			(matchange["matchangeFrom1To0", ] +
				matchange["matchangeFrom1To1", ]), 0.0)
	p00 <- pmax(0.02, p00)
	p11 <- pmax(0.02, p11)
	p00 <- pmin(0.98, p00)
	p11 <- pmin(0.98, p11)
	prec <- ifelse(matchange["matchangeFrom0To1", ] *
		matchange["matchangeFrom1To0", ] >= 1,
	4 / ((p00 / matchange["matchangeFrom0To1", ]) +
		(p11 / matchange["matchangeFrom1To0", ])), 1e-6)
	alphaf1 <- sum(alpha * prec / sum(prec))
	untrimmed <- alphaf1
	alphaf1 <- ifelse(alphaf1 < -3, -3, ifelse(alphaf1 > 3, 3, alphaf1))
	list(startRate=startRate, degree=alphaf1, alpha=alpha, prec=prec, tmp=tmp,
		untrimmed = untrimmed, startRateSett=startRateSett)
}
##@getBipartiteStartingVals DataCreate
getBipartiteStartingVals <- function(depvar)
{
	noPeriods <- attr(depvar, "netdims")[3] - 1
	##rate
	##get distance and number of valid links
	if (!attr(depvar,'sparse'))
	{
		nsenders<- nrow(depvar)
		nreceivers <- ncol(depvar)
		use <- !is.na(depvar) & (depvar == 10 | depvar == 11)
		depvar[use] <- depvar[use] - 10  ## remove structural values
		tmp <- sapply(1:noPeriods, function(x, z){
			matdiff <- sum(z[, , x + 1] != z[, , x], na.rm=TRUE)
			# matchange0 <- table(z[, , x + 1], z[, , x])
			# Changed to protect against zero rows or columns
			mc00 <- sum((1 - z[ , , x+1])*(1 - z[ , , x]), na.rm=TRUE)
			mc01 <- sum(z[ , , x+1]*(1 - z[ , , x]), na.rm=TRUE)
			mc10 <- sum((1 - z[ , , x+1])*z[ , , x], na.rm=TRUE)
			mc11 <- sum(z[ , , x+1]*z[ , , x], na.rm=TRUE)
			matchange <- matrix(c(mc00, mc01, mc10, mc11), 2, 2)
			#cat(matchange0,'\n',matchange,'\n')
			matcnt <- nsenders * nreceivers -
				sum(is.na(z[, , x + 1]) | is.na(z[, , x]))
			tmp <- c(matcnt=matcnt, matdiff=matdiff, matchange=matchange)
			names(tmp) <- c("matcnt", "matdiff", "matchangeFrom0To0",
				"matchangeFrom0To1",
				"matchangeFrom1To0", "matchangeFrom1To1")
			tmp
		}, z=depvar)
	}
	else
	{
		nsenders <- nrow(depvar[[1]])
		nreceivers <- ncol(depvar[[1]]) # CS: Was 2, but why?
		matdiff<- rep(NA, noPeriods)
		matcnt<- rep(NA, noPeriods)
		matchange<- matrix(NA, nrow=4, ncol=noPeriods)
		for (i in 1: noPeriods)
		{
			mymat1 <- depvar[[i]]
			mymat2 <- depvar[[i+1]]
			use <- mymat1@x %in% c(10, 11)
			mymat1@x[use] <- mymat1@x[use] - 10
			use <- mymat2@x %in% c(10, 11)
			mymat2@x[use] <- mymat2@x[use] - 10
			mymat1 <- drop0(mymat1)
			mymat2 <- drop0(mymat2)
			mydif <- mymat2 - mymat1
			matdiff[i] <- sum(abs(mydif), na.rm=TRUE)
			tmp <- table(mydif@x)
			tmp00 <- nsenders * nreceivers - length(mydif@x)
			tmp <- c(tmp00, tmp[c(3, 1, 2)])
			matchange[,i] <- tmp
			matcnt[i] <- sum(tmp)
		}
		matchange <- data.frame(matchange)
		tmp <-as.matrix(rbind(matcnt=matcnt, matdiff=matdiff,
				matchange=matchange))
		row.names(tmp) <- c("matcnt", "matdiff", "matchangeFrom0To0",
			"matchangeFrom0To1",
			"matchangeFrom1To0", "matchangeFrom1To1")
	}
	distance <- attr(depvar, "distance" )
	startRate <- nreceivers * (0.2 + 2 * distance)/(tmp['matcnt',] + 1)
	# CS: the above used to be 'nsenders' instead of 'nreceivers';
	#     this was a wrong calculation and led to extremely high
	#     rate parameters for comparatively small receiver nodesets
	#     slowing down estimation and prohibiting identification
	#     of parameters.
	startRate <- pmax(0.1, startRate)
	startRate <- pmin(100, startRate)
	##degree
	matchange<- as.matrix(tmp[grep("matchange", rownames(tmp)),,drop=FALSE])
	p01 <- ifelse (matchange["matchangeFrom0To0", ] +
		matchange["matchangeFrom0To1", ] >=1,
	matchange["matchangeFrom0To1", ] /
		(matchange["matchangeFrom0To0", ] +
			matchange["matchangeFrom0To1", ]), 0.5)
	p10 <- ifelse (matchange["matchangeFrom1To0", ]
		+ matchange["matchangeFrom1To1", ] >=1,
		matchange["matchangeFrom1To0", ] /
			(matchange["matchangeFrom1To0", ] +
				matchange["matchangeFrom1To1", ]), 0.5)
	p01 <- pmax(0.02,p01)
	p10 <- pmax(0.02,p10)
	p01 <- pmin(0.98,p01)
	p10 <- pmin(0.98,p10)
	alpha <- 0.5 * log(p01/p10)
	p00 <- ifelse (matchange["matchangeFrom0To0", ] +
		matchange["matchangeFrom0To1", ] >=1,
	matchange["matchangeFrom0To0", ] /
		(matchange["matchangeFrom0To0", ] +
			matchange["matchangeFrom0To1", ]), 0.0)
	p11 <- ifelse (matchange["matchangeFrom1To0", ]
		+ matchange["matchangeFrom1To1", ] >=1,
		matchange["matchangeFrom1To1", ] /
			(matchange["matchangeFrom1To0", ] +
				matchange["matchangeFrom1To1", ]), 0.0)
	p00 <- pmax(0.02,p00)
	p11 <- pmax(0.02,p11)
	p00 <- pmin(0.98,p00)
	p11 <- pmin(0.98,p11)
	prec <- ifelse(matchange["matchangeFrom0To1", ] *
		matchange["matchangeFrom1To0", ] >= 1,
	4 / ((p00 / matchange["matchangeFrom0To1", ]) +
		(p11 / matchange["matchangeFrom1To0", ])), 1e-6)
	alphaf1 <- sum(alpha*prec/sum(prec))
	## }
	untrimmed <- alphaf1
	alphaf1 <- ifelse(alphaf1 < -3, -3, ifelse(alphaf1 > 3, 3, alphaf1))
	list(startRate=startRate, degree=alphaf1, alpha=alpha, prec=prec, tmp=tmp,
		untrimmed = untrimmed)
}
