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
  
  staticChangeContributions <- getChangeContributions(
    algorithm,
    data,
    effects
  )
  
  n_periods <- length(staticChangeContributions[[1]])
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
      # staticChangeContributions[[1]][[p]][[e]] is a list of ego -> choice vector
      mat <- do.call(rbind, staticChangeContributions[[1]][[p]][[e]]) # [n_egos, n_choices]
      result[p, e, , ] <- mat
    }
  }
  
  if (requireNamespace("data.table", quietly = TRUE)) {
    DT <- data.table::as.data.table(as.table(result))
    data.table::setnames(DT, c("period", "effect", "ego", "choice", "cont"))
    DT[, ego := as.integer(ego)]
    DT[, choice := as.integer(choice)]
    DT[, period := as.integer(period)]
    DT <- data.table::dcast(DT, period + ego + choice ~ effect, value.var = "cont")
    # data.table::setDF(DT) only if we want to use pure data.frame later!
    DT
  } else {
    df <- as.data.frame(as.table(result))
    names(df) <- c("period","effect","ego","choice","cont")
    df$period <- as.integer(df$period)
    df$ego <- as.integer(df$ego)
    df$choice <- as.integer(df$choice)

    df_wide <- reshape(df, idvar = c("period", "ego", "choice"), timevar = "effect",
                       direction = "wide")
    colnames(df_wide) <- sub("^cont\\.", "", colnames(df_wide))
    df_wide
  }
}

# ##@calculateContributionsDynamic sienaRIDynamics (this is more transforming then calculating)
# calculateContributionsDynamic <- function(ans,
#                                          data,
#                                          theta,
#                                          algorithm,
#                                          effects, #effects object from ans$effects can not be used?!
#                                          depvar,
#                                          n3 = NULL, 
#                                          useChangeContributions = TRUE) {

#   ## should add some checks, that effects is a real effects object
#   changeContributions <- getChangeContributionsDynamic(
#         data = data,
#         ans = ans,
#         theta = theta,
#         algorithm = algorithm,
#         effects = effects,
#         n3 = n3, 
#         useChangeContributions = useChangeContributions,
#         returnDataFrame = TRUE)                         
#   chains <- algorithm$n3
# 	periods <- data$observation-1
# 	effects <- effects[effects$include==TRUE,]
#     # really restrict to noRate here?
# 	noRate <- effects$type != "rate"
# 	effectNames <- effects[noRate,"shortName"]
# 	depvar <- depvar

#   # changeContributions <- contrib_to_matrix(changeContributions, chains, periods, depvar, effectNames)
  
#   if(requireNamespace("data.table", quietly = TRUE)) {
#     DT <- rbindlist(
#       lapply(seq_along(changeContributions), function(i) {
#         DT <- as.data.table(changeContributions[[i]])
#         DT[, chain := i]
#         DT
#       })
#     )
#     DT <- dcast(DT, chain + group + period + ministep + choice ~ effectname, value.var="contribution")
#     return(DT)
#   } else {
#     outmat <- changeContributions
#     return(outmat)
#   }
# }

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

	effects <- effects[effects$include,] # unnecessary?
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
  ## is this really all the same ans?
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
	staticChangeContributions <- ans
  staticChangeContributions # why call this ans?
}

## extracts dynamic contributions from ans or simulates sequences ministeps to generate them. Changed and extracted from sienaRIDynamics
getChangeContributionsDynamic <- function(ans = NULL, theta=NULL, data=NULL, 
                                          algorithm, 
                                          effects,
                                          depvar = NULL,
                                          n3 = NULL, 
                                          useChangeContributions = FALSE, 
                                          returnDataFrame = FALSE){
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
      stopifnot(!is.null(theta))
		  prevAns <- NULL
		  effects$initialValue[effects$include] <- theta # somewhat unsafe, it theta is not in the correct order
      # thetaValues <- t(theta)
		}
    ans <- siena07(algorithm, data=data, effects=effects,
                       prevAns=prevAns, initC=FALSE, returnDeps=FALSE, 
                       returnChangeContributions=TRUE, returnDataFrame = returnDataFrame)
  }
  changeContributions <- ans$changeContributions

  if (returnDataFrame) {
    if (is.list(changeContributions) && !is.data.frame(changeContributions[[1]])) 
    {
      warning("SIENA backend returned nested lists, not data.frames; auto-flattening each chain.")
      changeContributions <- lapply(changeContributions, flattenChangeContributionsList)
    }
  # filtering before binding proofed to be quite slow
    if (requireNamespace("data.table", quietly = TRUE)) {
      changeContributions <- data.table::rbindlist(
        lapply(seq_along(changeContributions), function(i) {
          df <- changeContributions[[i]]
          df$chain <- i
          df
        }),
        use.names = TRUE, fill = TRUE
      )
      if (!is.null(depvar))
      {
        changeContributions <- changeContributions[changeContributions[["networkName"]] %in% depvar, , drop=FALSE]
      }
    } else {
      changeContributions <- do.call(rbind, lapply(seq_along(changeContributions), function(i) {
        df <- changeContributions[[i]]
        df$chain <- i
        df
      }))
      rownames(changeContributions) <- NULL
      if (!is.null(depvar))
      {
        changeContributions <- changeContributions[changeContributions[["networkName"]] %in% depvar, , drop=FALSE]
      }
    }
    return(changeContributions)
  }
  changeContributions
}

## Not doing this might reduce run time quite a bit! -> code would need to be adjusted though
widenContribution <- function(changeContributions){
  ## currently only works for dynamic case and data has to be pre filtered to only one depvar
  if (all(c("effectname", "contribution") %in% names(changeContributions))) {
    if (requireNamespace("data.table", quietly = TRUE)) {
        changeContributions <- data.table::dcast(
          changeContributions, chain + group + period + ministep + choice ~ effectname,
          value.var = "contribution"
        )
      }
    } else {
      if (all(c("effectname", "contribution") %in% names(changeContributions))) {
        changeContributions <- reshape(
          changeContributions,
          idvar = c("chain", "group", "period", "ministep", "choice"),
          timevar = "effectname", v.names = "contribution",
          direction = "wide"
        )
        colnames(changeContributions) <- sub("^contribution\\.", "", colnames(changeContributions))
      }
    }
  changeContributions
}

# Wrapper to use in R code (RCPP style)
flattenChangeContributionsList <- function(x) {
  .Call("C_flattenChangeContributionsList", x)
}