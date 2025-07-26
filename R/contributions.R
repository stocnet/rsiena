#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena/
# *
# * File: contributions.r
# *
# * Description: This module contains the functions to extract contributions
# * at observation moments or over simulated sequences of ministeps
# *****************************************************************************/

##@calculateContribution NOT USED YET IN sienaRI (this is more transforming then calculating)
calculateContribution <- function(ans, data,
                                  depvar = NULL,
                                  algorithm = NULL){
  ego <- choice <- period <- NULL # To resolve R CMD checks not understanding data.table syntax
  ## NOT NECESSARY?!
  effects <- ans$effects
  noRate <- effects$type != "rate"
  effects <- effects[noRate, ]
  effectNames <- effects[effects$include == TRUE,"shortName"]
  
  if(is.null(depvar)){
    depvar <- names(data$depvars)[1]
  }
  
  if(is.null(algorithm)){
    algorithm <- ans$x
  }
  
  ## should be possible to use counterfactual data here
  
  contributions <- getChangeContributions(
    algorithm,
    data,
    effects
  )
  
  n_periods <- length(contributions[[1]])
  n_effects <- length(effectNames)
  n_egos <- dim(data[["depvars"]][[depvar]])[1]
  n_choices <- dim(data[["depvars"]][[depvar]])[2]

  # Preallocate a big array: [period, effect, ego, choice]
  ## this looks super inefficient
  result <- array(NA_real_, dim = c(n_periods, n_effects, n_egos, n_choices),
               dimnames = list(period = seq_len(n_periods),
                               effect = effectNames,
                               ego = seq_len(n_egos),
                               choice = seq_len(n_choices)))
  for (p in seq_len(n_periods)) {
    for (e in seq_len(n_effects)) {
      # contributions[[1]][[p]][[e]] is a list of ego -> choice vector
      mat <- do.call(rbind, contributions[[1]][[p]][[e]]) # [n_egos, n_choices]
      result[p, e, , ] <- mat
    }
  }
  
  # indices <- expand.grid(
  #   choice = 1:n_choices,
  #   ego = 1:n_egos,
  #   effect = effectNames,
  #   period = seq_along(contributions[[1]])
  # )
  # 
  # df <- data.frame(
  #   period = indices$period,
  #   ego = indices$ego,
  #   choice = indices$choice,
  #   effect = indices$effect,
  #   # quite slow
  #   cont = mapply(function(x, y, z, w) contributions[[1]][[x]][[y]][[z]][w], 
  #                 indices$period, indices$effect, indices$ego, indices$choice)
  # )
  # # very inefficient! especially "interaction"?
  # df <- reshape(df, idvar = c("period", "ego", "choice"), 
  #               timevar = "effect", direction = "wide")
  # df <- setNames(df,
  #                gsub("(.*)\\.","",names(df)))
  
  ## base R solution is much more inefficient
  result <- data.table::as.data.table(as.table(result))
  data.table::setnames(result, c("period", "effect", "ego", "choice", "cont"))
  ## only works when [ from data.dable is imported?
  result[, ego := as.integer(ego)]
  result[, choice := as.integer(choice)]
  result[, period := as.integer(period)]
  # Now dcast to wide, as in your function:
  result <- data.table::dcast(
    result,
    period + ego + choice ~ effect, 
    value.var = "cont"
  )
  data.table::setDF(result)
  result
  # df <- as.data.frame.table(arr, responseName = "cont")
  # 
  # # Reshape to wide
  # df_wide <- reshape(df,
  #                    idvar = c("period", "ego", "choice"),
  #                    timevar = "effect",
  #                    direction = "wide")
  # # Optionally fix column names
  # colnames(df_wide) <- sub("^cont\\.", "", colnames(df_wide))
  # 
  # df_wide
}

##@calculateContributionsDynamic sienaRIDynamics (this is more transforming then calculating)
calculateContributionsDynamic <- function(ans,
                                         data,
                                         theta,
                                         algorithm,
                                         effects,
                                         depvar,
                                         n3 = NULL, 
                                         useChangeContributions = FALSE) {
  changeContributions <- getChangeContributionsDynamic(
        data = data,
        ans = ans,
        theta = theta,
        algorithm = algorithm,
        effects = effects,
        n3 = n3, 
        useChangeContributions = useChangeContributions)                                     
  chains <- algorithm$n3
	periods <- data$observation-1
	effects <- effects[effects$include==TRUE,]
    # really restrict to noRate here?
	noRate <- effects$type != "rate"
	effectNames <- effects[noRate,"shortName"]
	depvar <- depvar
#	effectTypes <- effects$type[noRate]
#	networkName <- effects$name[noRate]
#	networkInteraction <- effects$interaction1[noRate]
#	effectIds <- paste(effectNames,effectTypes,networkInteraction, sep = ".")
	currentNetObjEffs <- effects$name[noRate] == depvar
    all_rows <- list()
    row_idx <- 1L
	for (chain in (1:chains))
	{
		for(period in 1:periods)
		{
			ministeps <- length(changeContributions[[chain]][[1]][[period]])
			for(ministep in 1:ministeps)
			{
				#depvar <- attr(changeContributions[[1]][[period]][[ministep]],"networkName")
				if(attr(changeContributions[[chain]][[1]][[period]][[ministep]],
												"networkName")==depvar)
				{
					cdec <- ans$changeContributions[[1]][[period]][[ministep]]
                    # cdec is a contributions x choices matrix
                    n_choices <- ncol(cdec)
                    for(choice in seq_len(n_choices)) {
                        all_rows[[row_idx]] <- c(chain, 
                                                period, 
                                                ministep, 
                                                choice,
                                                cdec[, choice, drop = TRUE])
                    row_idx <- row_idx + 1L
				    }
			    }
		    }
    	}
    }
  # combine rows into one matrix (inefficient?)
  outmat <- do.call(rbind, all_rows)
  colnames(outmat) <- c("chain","period","ministep","choice", effectNames)
  outmat
}

##@getChangeContributions. Use as RSiena:::getChangeContributions - copied from sienaRI
getChangeContributions <- function(algorithm, data, effects)
{
    ## Daniel: Probably should be changed similarily to the dynamic contributions

	## Gets the simulated statistics.
	## The following initializations data, effects, and model
	## for calling "getTargets" in "siena07.setup.h"
	## is more or less copied from "getTargets" in "getTargets.r".
	## However, some modifications have been necessary to get it to work.
	f <- unpackData(data,algorithm)

	effects <- effects[effects$include,]
	if (!is.null(algorithm$settings))
	{
		stop('not implemented: RI together with settings')
		# effects <- addSettingsEffects(effects, algorithm)
	}
	else
	{
		effects$setting <- rep("", nrow(effects))
	}
	pData <- .Call(C_setupData, PACKAGE=pkgname,
		list(as.integer(f$observations)),
		list(f$nodeSets))
	## register a finalizer
	ans <- reg.finalizer(pData, clearData, onexit = FALSE)
	ans<- .Call(C_OneMode, PACKAGE=pkgname,
		pData, list(f$nets))
	ans <- .Call(C_Bipartite, PACKAGE=pkgname, # added 1.1-299
		pData, list(f$bipartites))
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

	effects$effectPtr <- rep(NA, nrow(effects))
	depvarnames <- names(data$depvars)
	tmpeffects <- split(effects, effects$name)
	myeffectsOrder <- match(depvarnames, names(tmpeffects))
	ans <- .Call(C_effects, PACKAGE=pkgname, pData, tmpeffects)
	pModel <- ans[[1]][[1]]
	for (i in 1:length(ans[[2]]))
	{
		effectPtr <- ans[[2]][[i]]
		tmpeffects[[i]]$effectPtr <- effectPtr
	}
	myeffects <- tmpeffects
	for(i in 1:length(myeffectsOrder)){
		myeffects[[i]]<-tmpeffects[[myeffectsOrder[i]]]
	}
#cat("e\n") #Hier gaat hij fout. 
# Voor returnStaticChangeContributions=TRUE gaat het wel goed.
#browser()
	ans <- .Call(C_getTargets, PACKAGE=pkgname, pData, pModel, myeffects,
		parallelrun=TRUE, returnActorStatistics=FALSE,
		returnStaticChangeContributions=TRUE)
	# See getTargets in siena07setup.cpp; also see rTargets in StatisticsSimulation.cpp
	ans # why call this ans?
}

## extracts dynamic contributions from ans or simulates sequences ministeps to generate them. Changed and extracted from sienaRIDynamics
getChangeContributionsDynamic <- function(data=NULL, ans = NULL, theta=NULL, algorithm, effects, n3 = NULL, useChangeContributions = FALSE){
	if (!is.null(n3)) 
  {
   	algorithm$n3 <- as.integer(n3)
	}
	if (useChangeContributions) 
	{
    if (is.null(ans) || is.null(ans$changeContributions)) 
		{
      warning("useChangeContributions=TRUE, but 'ans' does not contain 'changeContributions'. Will rerun siena07 to generate them.")
      useChangeContributions <- FALSE
    }
		if (!is.null(n3))
		{
	    warning("n3 cannot be changed when useChangeContributions=TRUE; using the number of chains stored in 'ans'.")
   		n3 <- NULL
		}
  }
  if (!useChangeContributions)
	{
		algorithm$nsub <- 0
    if (!is.null(ans))
		{
		  prevAns <- ans
		} else {
		  prevAns <- NULL
		  effects$initialValue[effects$include] <- theta # somewhat unsafe, it theta is not in the correct order
		}
    ans <- siena07(algorithm, data=data, effects=effects,
                       prevAns=prevAns, initC=FALSE, returnDeps=FALSE, 
                       returnChangeContributions=TRUE)
  }
	changeContributions <- ans$changeContributions
}