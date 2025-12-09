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
  }
  if (is.null(data) || is.null(algorithm) || is.null(effects)) {
    stop("Must provide either 'ans' or all of 'algorithm' and 'effects'.")
  }

  # Remove rate effects
  noRate <- effects$type != "rate"
  effects <- effects[noRate, ]
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
    effects <- effects[effects$include, ]
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

# ##@calculateContribution NOT USED YET IN sienaRI (this is more transforming then calculating)
# calculateContribution <- function(ans, data,
#                                   depvar = NULL,
#                                   algorithm = NULL){
#   ego <- choice <- period <- NULL # To resolve R CMD checks not understanding data.table syntax
#   ## NOT NECESSARY?!
#   effects <- ans$requestedEffects
#   noRate <- effects$type != "rate"
#   effects <- effects[noRate, ]
#   effectNames <- effects[effects$include == TRUE,"shortName"]
#   effectDepvars <- effects[effects$include == TRUE, "name"]
#   effectNetTypes <- effects[effects$include == TRUE, "netType"]


#   depvars_available <- names(data$depvars)
#   if (is.null(depvar)) {
#     depvar <- depvars_available
#   }
 
#   if(is.null(algorithm)){
#     algorithm <- ans$x
#   }
  
#   ## should be possible to use counterfactual data here
#   ## potentially restrict to only one depvar here already?
#   staticChangeContributions <- getStaticChangeContributions(
#     algorithm,
#     data,
#     effects
#   )
  
#   staticChangeContributions_effectNames <- attr(staticChangeContributions, "effectNames")

#   # Filter effects for the selected depvar
#   effect_idx <- which(effectDepvars == depvar)
#   effectNames_depvar <- effectNames[effect_idx]

#   # Match effectNames_depvar to staticChangeContributions effect names
#   match_idx <- match(effectNames_depvar, staticChangeContributions_effectNames)
#   n_effects <- length(match_idx)

#   results <- list()
#   for (dv in depvar) {
#     # Filter effects for this depvar
#     effect_idx <- which(effectDepvars == dv)
#     effectNames_depvar <- effectNames[effect_idx]
#     match_idx <- match(effectNames_depvar, staticChangeContributions_effectNames)
#     n_effects <- length(match_idx)
#     n_egos <- dim(data[["depvars"]][[dv]])[1]
#     netType <- unique(effectNetTypes[effectDepvars == dv])
#     if (netType == "oneMode") {
#       n_choices <- dim(data[["depvars"]][[dv]])[2]
#     } else if (netType == "behavior") {
#       beh <- data[["depvars"]][[dv]]
#       n_choices <- 3
#     } else {
#       stop("Unsupported netType: ", netType)
#     }
#     for (g in seq_along(staticChangeContributions)) {
#       group_contrib <- staticChangeContributions[[g]]
#       n_periods <- length(group_contrib)
#       # Preallocate array for this group
#       result <- array(NA_real_, dim = c(n_periods, n_effects, n_egos, n_choices),
#                       dimnames = list(
#                         period = seq_len(n_periods),
#                         effect = effectNames_depvar,
#                         ego = seq_len(n_egos),
#                         choice = seq_len(n_choices)
#                       ))
#       for (p in seq_len(n_periods)) {
#         for (e in seq_len(n_effects)) {
#           mat <- do.call(rbind, group_contrib[[p]][[match_idx[e]]])
#           result[p, e, , ] <- mat
#         }
#       }

#       if (requireNamespace("data.table", quietly = TRUE)) {
#         DT <- data.table::as.data.table(as.table(result))
#         data.table::setnames(DT, c("period", "effectname", "ego", "choice", "contribution"))
#         DT[, ego := as.integer(ego)]
#         DT[, choice := as.integer(choice)]
#         DT[, period := as.integer(period)]
#         DT[, group := g]
#         DT[, networkName := dv]
#         # DT <- data.table::dcast(DT, networkName + group + period + ego + choice ~ effect, value.var = "cont")
#         results[[paste0("g", g, "_", dv)]] <- DT
#       } else {
#         df <- as.data.frame(as.table(result))
#         names(df) <- c("period","effectname","ego","choice","contribution")
#         df$period <- as.integer(df$period)
#         df$ego <- as.integer(df$ego)
#         df$choice <- as.integer(df$choice)
#         df$group <- g
#         df$networkName <- dv
#         # df_wide <- reshape(df, idvar = c("networkName", "group","period", "ego", "choice"), timevar = "effect",
#         #                   direction = "wide")
#         # colnames(df_wide) <- sub("^cont\\.", "", colnames(df_wide))
#         # results[[dv]] <- df_wide
#         results[[paste0("g", g, "_", dv)]] <- df
#       }
#     }
#   }
#   # Combine all depvars
#   # could be kept as a list of dt/dfs instead
#   if (requireNamespace("data.table", quietly = TRUE)) {
#     data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
#   } else {
#     do.call(rbind, results)
#   }
# }

# ##@calculateContributionsDynamic sienaRIDynamics (this is more transforming then calculating)
# calculateContributionsDynamic <- function(ans,
#                                          data,
#                                          theta,
#                                          algorithm,
#                                          effects, #effects object from ans$requestedEffects can not be used?!
#                                          depvar,
#                                          n3 = NULL, 
#                                          useChangeContributions = TRUE) {

#   ## should add some checks, that effects is a real effects object
#   changeContributions <- getDynamicChangeContributions(
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

##@getStaticChangeContributions. Use as RSiena:::getStaticChangeContributions
# getStaticChangeContributions <- function(algorithm, data, effects) {
# 	if (!is.null(algorithm$settings)) {
#         stop('not implemented: RI together with settings')
#         # effects <- addSettingsEffects(effects, algorithm)
#     } else {
#         pData <- sienaSetupDataForCpp(algorithm,
#                                       data,
#                                       includeBehavior = TRUE,
#                                       includeBipartite = FALSE)
#         effects <- effects[effects$include, ]
#         setup <- sienaSetupEffectsForCpp(pData,
#                                        data, 
#                                        effects)
#   }

#     ans <- .Call(C_getStaticChangeContributions, PACKAGE=pkgname, 
#                  pData, setup$pModel, setup$myeffects,
#                  parallelrun = FALSE)
#     ans
# }

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
      effects <- ans$effects
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

# Wrapper to use in R code (RCPP style)
flattenChangeContributionsList <- function(x) {
  .Call("C_flattenChangeContributionsList", x)
}