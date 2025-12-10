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

getStaticChangeContributions <- function(ans = NULL,
                                         data,
                                         algorithm = NULL,
                                         effects = NULL,
                                         depvar = NULL,
                                         returnDataFrame = FALSE) {
  # Flexible argument handling
  if (!is.null(ans)) {
    if (is.null(algorithm)) algorithm <- ans$x
    if (is.null(effects)) effects <- ans$requestedEffects 
    # does not work for interactions without base effects?
  }
  if (is.null(data) || is.null(algorithm) || is.null(effects)) {
    stop("Must provide either 'ans' or all of 'algorithm' and 'effects'.")
  }
  # Remove rate effects
  noRate <- effects$type != "rate"
  effects <- effects[noRate, ] # not sure if correct
  effectNames <- effects[effects$include == TRUE, "shortName"]
  effectDepvars <- effects[effects$include == TRUE, "name"]
  effectNetTypes <- effects[effects$include == TRUE, "netType"]

  depvars_available <- names(data$depvars)
  if (is.null(depvar)) {
    depvar <- depvars_available
  }

  # Prepare for C++ call
  if (!is.null(algorithm$settings)) {
    stop('not implemented: RI together with settings')
  } else {
    pData <- sienaSetupDataForCpp(algorithm,
                                  data,
                                  includeBehavior = TRUE,
                                  includeBipartite = FALSE)
    setup <- sienaSetupEffectsForCpp(pData, data, effects)
  }

  staticChangeContributions <- .Call(C_getStaticChangeContributions, PACKAGE = pkgname,
                                     pData, setup$pModel, setup$myeffects,
                                     parallelrun = FALSE)
  #attr(staticChangeContributions, "effectNames") <- effects$shortName

  if (!returnDataFrame) {
    return(staticChangeContributions)
  }

  staticChangeContributions_effectNames <- attr(staticChangeContributions, "effectNames")
  results <- list()
  for (dv in depvar) {
    effect_idx <- which(effectDepvars == dv)
    effectNames_depvar <- effectNames[effect_idx]
    match_idx <- match(effectNames_depvar, staticChangeContributions_effectNames)
    n_effects <- length(match_idx)
    n_egos <- dim(data[["depvars"]][[dv]])[1]
    netType <- unique(effectNetTypes[effectDepvars == dv])
    if (netType == "oneMode") {
      n_choices <- dim(data[["depvars"]][[dv]])[2]
    } else if (netType == "behavior") {
      n_choices <- 3
    } else {
      stop("Unsupported netType: ", netType)
    }
    for (g in seq_along(staticChangeContributions)) {
      group_contrib <- staticChangeContributions[[g]]
      n_periods <- length(group_contrib)
      result <- array(NA_real_, dim = c(n_periods, n_effects, n_egos, n_choices),
                      dimnames = list(
                        period = seq_len(n_periods),
                        effect = effectNames_depvar,
                        ego = seq_len(n_egos),
                        choice = seq_len(n_choices)
                      ))
      for (p in seq_len(n_periods)) {
        for (e in seq_len(n_effects)) {
          mat <- do.call(rbind, group_contrib[[p]][[match_idx[e]]])
          result[p, e, , ] <- mat
        }
      }
      if (requireNamespace("data.table", quietly = TRUE)) {
        # To resolve R CMD checks not understanding data.table syntax
        ego <- choice <- period <- group <- networkName <- NULL 
        DT <- data.table::as.data.table(as.table(result))
        data.table::setnames(DT, c("period", "effectname", "ego", "choice", "contribution"))
        DT[, ego := as.integer(ego)]
        DT[, choice := as.integer(choice)]
        DT[, period := as.integer(period)]
        DT[, group := g]
        DT[, networkName := dv]
        results[[paste0("g", g, "_", dv)]] <- DT
      } else {
        df <- as.data.frame(as.table(result))
        names(df) <- c("period", "effectname", "ego", "choice", "contribution")
        df$period <- as.integer(df$period)
        df$ego <- as.integer(df$ego)
        df$choice <- as.integer(df$choice)
        df$group <- g
        df$networkName <- dv
        results[[paste0("g", g, "_", dv)]] <- df
      }
    }
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
  } else {
    do.call(rbind, results)
  }
}

##@sienaSetupDataForCpp Use as RSiena:::sienaSetupDataForCpp
sienaSetupDataForCpp <- function(algorithm, data, 
                             includeBehavior = TRUE, 
                             includeBipartite = TRUE) {
    f <- unpackData(data, algorithm)
    pData <- .Call(C_setupData, PACKAGE=pkgname,
                   list(as.integer(f$observations)),
                   list(f$nodeSets))
    reg.finalizer(pData, clearData, onexit = FALSE)
    .Call(C_OneMode, PACKAGE=pkgname, pData, list(f$nets))
    if (includeBipartite && !is.null(f$bipartites)) {
        .Call(C_Bipartite, PACKAGE=pkgname, pData, list(f$bipartites))
    }
    if (includeBehavior && !is.null(f$behavs)) {
        .Call(C_Behavior, PACKAGE=pkgname, pData, list(f$behavs))
    }
    .Call(C_ConstantCovariates, PACKAGE=pkgname, pData, list(f$cCovars))
    .Call(C_ChangingCovariates, PACKAGE=pkgname, pData, list(f$vCovars))
    .Call(C_DyadicCovariates, PACKAGE=pkgname, pData, list(f$dycCovars))
    .Call(C_ChangingDyadicCovariates, PACKAGE=pkgname, pData, list(f$dyvCovars))
    # if (is.null(f$exog) || !is.list(f$exog)) {
	# 	f$exog <- vector("list", length(f$depvars))
	# }
	# ## split the names of the constraints
	# higher <- attr(f, "allHigher")
	# disjoint <- attr(f, "allDisjoint")
	# atLeastOne <- attr(f, "allAtLeastOne")
	# froms <- sapply(strsplit(names(higher), ","), function(x)x[1])
	# tos <- sapply(strsplit(names(higher), ","), function(x)x[2])
	# .Call(C_Constraints, PACKAGE=pkgname,
	# 	pData, froms[higher], tos[higher],
	# 	froms[disjoint], tos[disjoint],
	# 	froms[atLeastOne], tos[atLeastOne])
	# siena07 only does this with !initC
    return(pData)
}

##@sienaSetupEffectsForCpp Use as RSiena:::sienaSetupEffectsForCpp
sienaSetupEffectsForCpp <- function(pData, data, effects) {
    effects$setting <- rep("", nrow(effects))
    storage.mode(effects$parm) <- 'integer'
    storage.mode(effects$group) <- 'integer'
    storage.mode(effects$period) <- 'integer'
	## find any effects not included which are needed for interactionst
	tmpEffects <- effects[effects$include, ]
	interactionNos <- unique(c(tmpEffects$effect1, tmpEffects$effect2,
			tmpEffects$effect3))
	interactionNos <- interactionNos[interactionNos > 0]
	interactions <- effects$effectNumber %in%
		interactionNos
	effects$requested <- effects$include
	# requestedEffects <- effects[effects$include, ] 

	effects$include[interactions] <- TRUE
	effects <- effects[effects$include, ]

    effects$effectPtr <- rep(NA, nrow(effects))

    myeffects <- split(effects, effects$name)
	myCompleteEffects <- myeffects
	## remove interaction effects and save till later
	basicEffects <-
		lapply(myeffects, function(x)
			{
				x[!x$shortName %in% c("unspInt", "behUnspInt", "contUnspInt"), ]
			}
			)
	basicEffectsl <-
		lapply(myeffects, function(x)
			{
				!x$shortName %in% c("unspInt", "behUnspInt", "contUnspInt")
			}
			)

	interactionEffects <-
		lapply(myeffects, function(x)
			{
				x[x$shortName %in% c("unspInt", "behUnspInt", "contUnspInt"), ]
			}
			)
	interactionEffectsl <-
		lapply(myeffects, function(x)
			{
				x$shortName %in% c("unspInt", "behUnspInt", "contUnspInt")
			}
			)
    ans <- .Call(C_effects, PACKAGE=pkgname, pData, basicEffects)
    pModel <- ans[[1]][[1]]
    reg.finalizer(pModel, clearModel, onexit = FALSE)
    for (i in 1:length(ans[[2]])) ## ans[[2]] is a list of lists of
        ## pointers to effects. Each list corresponds to one
        ## dependent variable
    {
		effectPtr <- ans[[2]][[i]]
		basicEffects[[i]]$effectPtr <- effectPtr
		interactionEffects[[i]]$effect1 <-
			basicEffects[[i]]$effectPtr[match(interactionEffects[[i]]$effect1,
				basicEffects[[i]]$effectNumber)]
		interactionEffects[[i]]$effect2 <-
			basicEffects[[i]]$effectPtr[match(interactionEffects[[i]]$effect2,
				basicEffects[[i]]$effectNumber)]
		interactionEffects[[i]]$effect3 <-
			basicEffects[[i]]$effectPtr[match(interactionEffects[[i]]$effect3,
				basicEffects[[i]]$effectNumber)]
    }
	ans <- .Call(C_interactionEffects, PACKAGE=pkgname,
		pModel, interactionEffects)
		## copy these pointers to the interaction effects and then insert in
	## effects object in the same rows for later use
	for (i in 1:length(ans[[1]])) ## ans is a list of lists of
		## pointers to effects. Each list corresponds to one
		## dependent variable
	{
		if (nrow(interactionEffects[[i]]) > 0)
		{
			effectPtr <- ans[[1]][[i]]
			interactionEffects[[i]]$effectPtr <- effectPtr
		}
		myCompleteEffects[[i]][basicEffectsl[[i]], ] <- basicEffects[[i]]
		myCompleteEffects[[i]][interactionEffectsl[[i]],] <-
			interactionEffects[[i]]
		##myeffects[[i]] <- myeffects[[i]][order(myeffects[[i]]$effectNumber),]
	}
	## remove the effects only created as underlying effects
	## for interaction effects. first store the original for use next time
	myeffects <- lapply(myCompleteEffects, function(x)
		{
			x[x$requested, ]
		}
	)
    list(pModel = pModel, myeffects = myeffects)
}

widenStaticContribution <- function(changeContributions){
  ## currently only works for dynamic case and data has to be pre filtered to only one depvar
  if (all(c("effectname", "contribution") %in% names(changeContributions))) {
    if (requireNamespace("data.table", quietly = TRUE)) {
        changeContributions <- data.table::dcast(
          changeContributions, group + period + ego + choice ~ effectname,
          value.var = "contribution"
        )
      } else {
        ## in the dynamic case, this might lose some information that we should keep
        needed <- c("group", 
          "period", 
          "ego", 
          "choice", 
          "effectname", 
          "contribution")
        changeContributions <- changeContributions[, needed]
        changeContributions <- reshape(
          changeContributions,
          idvar = c("group", "period", "ego", "choice"),
          timevar = "effectname", v.names = "contribution",
          direction = "wide"
        )
        colnames(changeContributions) <- sub("^contribution\\.", "", 
          colnames(changeContributions))
      }
    }
  changeContributions
}

## extracts dynamic contributions from ans or simulates sequences ministeps to 
## generate them. Changed and extracted from sienaRIDynamics
getDynamicChangeContributions <- function(ans = NULL, 
                                          theta=NULL, 
                                          data = NULL, 
                                          algorithm = NULL, 
                                          effects = NULL,
                                          depvar = NULL,
                                          n3 = NULL, 
                                          useChangeContributions = FALSE, 
                                          returnDataFrame = FALSE,
                                          silent = TRUE,
                                          seed = NULL) {
  # Flexible argument handling allows either extracted pre
  # calculated contributions from ans or simulating them

  if (useChangeContributions && (is.null(ans) || is.null(ans$changeContributions))) 
  {
            warning("useChangeContributions=TRUE, but 'ans' does not 
			      contain 'changeContributions'. Will rerun siena07 to generate them.")
            useChangeContributions <- FALSE
  }
  if (useChangeContributions && !returnDataFrame && 
    (is.data.frame(ans$changeContributions) || 
      is.data.frame(ans$changeContributions[[1]]))) 
  {
            stop("useChangeContributions=TRUE & returnDataFrame=FALSE; but
              'ans$changeContributions' contains data.frames")
  }
  if (useChangeContributions && is.null(n3)) 
  {
    warning("n3 cannot be changed when useChangeContributions=TRUE; 
    using the number of chains stored in 'ans'.")
  }
  if(!useChangeContributions) # could be extracted to a separate function
  {
    if (is.null(algorithm))
    {
      stop("Must provide 'algorithm' when useChangeContributions=FALSE.")
    } 
    else
    {
      if (!inherits(algorithm, "sienaAlgorithm")) 
      {
        stop("algorithm is not a legitimate Siena algorithm specification")
      }
    }
    if(is.null(data)) 
    {
      stop("Must provide 'data' when useChangeContributions=FALSE ")
    } else 
    {
      if (!inherits(data, "siena")) {
        stop("'data' is not a legitimate Siena data specification")
      }
    }
    if(is.null(effects))
    {
      if (is.null(ans)) 
      {
        stop("Must provide 'ans' or 'effects' when useChangeContributions=FALSE.")
      }
      if (ans$cconditional) 
      {
        stop("Must provide 'effects' for conditional estimation")
      }
      effects <- ans$requestedEffects
      if (!inherits(effects, "sienaEffects")) 
      {
        stop("effects is not a legitimate Siena effects object")
      }
    }
    if(is.null(theta)) {
      if (is.null(ans)) {
        stop("Must provide 'ans' or 'theta' when useChangeContributions=FALSE.")
      }
      theta <- ans$theta
      if (!is.numeric(theta)) 
      {
        stop("theta is not a legitimate parameter vector")
      }
    }
    # Prepare algorithm for simulating only phase 3
    algorithm$nsub <- 0
    if (!is.null(n3)) 
    {
      algorithm$n3 <- as.integer(n3)
    }
    if(!is.null(seed))
    {
      if(is.numeric(seed))
      {
        algorithm$seed <- as.integer(seed)
      } else {
        warning("'seed' has to be of type 'numeric' \n used default settings")
      }
    }
    # Update effects initial values with provided theta to sample chains from
    # specificed parameter values in phase 3 of siena07
    include <- effects[["include"]]
    includedEffects <- effects[include, ]
    noRateIncluded <- includedEffects[["type"]] != "rate"
    if(length(theta) > length(effects[include, "initialValue"])) {
        theta <- theta[noRateIncluded] # if conditional estimation?
    } else if(length(theta) < length(effects[include, "initialValue"])) {
        theta <- c(ans$rate, theta) # if unconditional estimation?
    }
    effects$initialValue[effects$include] <- theta
    ans <- siena07(
      algorithm, 
      data=data, 
      effects=effects,
      initC=FALSE, # why?
      returnDeps=FALSE, 
      returnChangeContributions=TRUE, 
      returnDataFrame = returnDataFrame,
      silent = silent)
  }

  changeContributions <- ans$changeContributions

  if (!returnDataFrame) { # could be extracted to a separate function
    return(changeContributions) # what exactly is the list format?
  } else {
    if (is.list(changeContributions) && !is.data.frame(changeContributions[[1]])) 
    {
      warning("SIENA backend returned nested lists, not data.frames; 
      auto-flattening each chain.")
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
        changeContributions <- changeContributions[changeContributions[["networkName"]] 
          %in% depvar, , drop=FALSE]
      }
    } else {
      changeContributions <- do.call(rbind, 
        lapply(seq_along(changeContributions), function(i) {
          df <- changeContributions[[i]]
          df$chain <- i
          df
      }))
      rownames(changeContributions) <- NULL
      if (!is.null(depvar))
      {
        changeContributions <- changeContributions[
          changeContributions[["networkName"]] %in% depvar, , drop=FALSE]
      }
    }
    return(changeContributions)
  }
}

## Not doing this might reduce run time quite a bit! -> code would need to be adjusted though
widenDynamicContribution <- function(changeContributions){
  ## currently only works for dynamic case and data has to be pre filtered to only one depvar
  if (all(c("effectname", "contribution") %in% names(changeContributions))) {
    if (requireNamespace("data.table", quietly = TRUE)) {
        changeContributions <- data.table::dcast(
          changeContributions, chain + group + period + ministep + choice ~ effectname,
          value.var = "contribution"
        )
      } else {
        ## in the dynamic case, this might lose some information that we should keep
        needed <- c("chain", 
          "group", 
          "period", 
          "ministep", 
          "choice", 
          "effectname", 
          "contribution")
        changeContributions <- changeContributions[, needed]
        changeContributions <- reshape(
          changeContributions,
          idvar = c("chain", "group", "period", "ministep", "choice"),
          timevar = "effectname", v.names = "contribution",
          direction = "wide"
        )
        colnames(changeContributions) <- sub("^contribution\\.", "", 
          colnames(changeContributions))
      }
    }
  changeContributions
}

# Wrapper to use in R code (RCPP style)
flattenChangeContributionsList <- function(x) {
  .Call("C_flattenChangeContributionsList", x)
}
