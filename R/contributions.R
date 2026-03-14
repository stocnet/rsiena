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
                                         effects = NULL,
                                         depvar = NULL,
                                         returnDataFrame = FALSE,
                                         returnWide = FALSE) {
  # Flexible argument handling
  if (!is.null(ans)) {
    if (is.null(effects)) effects <- ans$requestedEffects 
  }
  if (is.null(data) || is.null(effects)) {
    stop("Must provide either 'ans' or all of 'data' and 'effects'.")
  }
  noRate <- effects$type != "rate"
  effects <- effects[noRate, ]
  effectNames <- effects[effects$include == TRUE, "shortName"]
  effectDepvars <- effects[effects$include == TRUE, "name"]
  effectNetTypes <- effects[effects$include == TRUE, "netType"]
  effectTypes <- effects[effects$include == TRUE, "type"]


  depvars_available <- names(data$depvars)
  if (is.null(depvar)) {
    depvar <- depvars_available
  }

  # Prepare for C++ call
  # Settings can not be checked without algorithm?
  # if (!is.null(algorithm$settings)) {
  #   stop('not implemented: RI together with settings')
  # } else {
  pData <- sienaSetupDataForCpp(data,
                                includeBehavior = TRUE,
                                includeBipartite = FALSE)
  setup <- sienaSetupEffectsForCpp(pData, data, effects)
  #}

  staticChangeContributions <- .Call(
    C_getStaticChangeContributions, 
    PACKAGE = pkgname,
    pData, setup$pModel, setup$myeffects,
    parallelrun = FALSE
  )
  #attr(staticChangeContributions, "effectNames") <- effects$shortName

  if (!returnDataFrame && !returnWide) {
    return(staticChangeContributions)
  }

  staticChangeContributions_effectNames <- attr(staticChangeContributions,
    "effectNames")

  if (returnWide) {
    # Accumulate across all groups for one depvar into a single unified struct.
    if (length(depvar) != 1L)
      stop("returnWide = TRUE requires exactly one depvar; got ",
           length(depvar), ". Call once per depvar.")
    all_contribMat <- list()
    all_period     <- list()
    all_ego        <- list()
    all_choice     <- list()
    all_group_id   <- list()
    all_group      <- list()
    effectNames_out <- NULL
    thetaIdx_out    <- NULL
    group_id_offset <- 0L

    for (dv in depvar) {
      effectIdx      <- which(effectDepvars == dv)
      effectNamesDv  <- effectNames[effectIdx]
      compositeNames <- paste(effectNamesDv, effectTypes[effectIdx], sep = "_")
      nEffects       <- length(effectIdx)
      nEgos          <- dim(data[["depvars"]][[dv]])[1]
      netType        <- unique(effectNetTypes[effectDepvars == dv])
      nChoices       <- if (netType == "oneMode") dim(data[["depvars"]][[dv]])[2L]
                        else if (netType == "behavior") 3L
                        else stop("Unsupported netType: ", netType)

      if (is.null(effectNames_out)) {
        effectNames_out <- compositeNames
        thetaIdx_out    <- effectIdx
      } else {
        effectNames_out <- c(effectNames_out, compositeNames)
        thetaIdx_out    <- c(thetaIdx_out, effectIdx)
      }

      for (g in seq_along(staticChangeContributions)) {
        nPeriods <- length(staticChangeContributions[[g]])
        arr <- array(NA_real_, dim = c(nPeriods, nEffects, nEgos, nChoices))
        for (p in seq_len(nPeriods)) {
          for (e in seq_len(nEffects)) {
            arr[p, e, , ] <- do.call(rbind,
              staticChangeContributions[[g]][[p]][[effectIdx[e]]])
          }
        }
        arr <- aperm(arr, c(4L, 3L, 1L, 2L))
        N          <- nPeriods * nEgos * nChoices
        contribMat <- matrix(arr, nrow = N, ncol = nEffects)
        colnames(contribMat) <- compositeNames
        period_vec <- rep(seq_len(nPeriods), each = nEgos * nChoices)
        ego_vec    <- rep(rep(seq_len(nEgos), each = nChoices), times = nPeriods)
        choice_vec <- rep(seq_len(nChoices), times = nPeriods * nEgos)
        group_id   <- group_id_offset + rep(seq_len(nPeriods * nEgos), each = nChoices)

        all_contribMat[[length(all_contribMat) + 1L]] <- contribMat
        all_period[[length(all_period) + 1L]]         <- period_vec
        all_ego[[length(all_ego) + 1L]]               <- ego_vec
        all_choice[[length(all_choice) + 1L]]         <- choice_vec
        all_group_id[[length(all_group_id) + 1L]]     <- group_id
        all_group[[length(all_group) + 1L]]            <- rep(as.integer(g), N)
        group_id_offset <- group_id_offset + nPeriods * nEgos
      }
    }

    return(list(
      contribMat  = do.call(rbind, all_contribMat),
      period      = unlist(all_period,   use.names = FALSE),
      ego         = unlist(all_ego,      use.names = FALSE),
      choice      = unlist(all_choice,   use.names = FALSE),
      group_id    = unlist(all_group_id, use.names = FALSE),
      group       = unlist(all_group,    use.names = FALSE),
      effectNames = effectNames_out,
      thetaIdx    = thetaIdx_out
    ))
  }

  results <- list()
  for (dv in depvar) {
    effectIdx <- which(effectDepvars == dv)
    effectNamesDepvar <- effectNames[effectIdx]
    compositeNames    <- paste(effectNamesDepvar, effectTypes[effectIdx], sep = "_")
    typeMap           <- setNames(effectTypes[effectIdx], compositeNames)
    nameMap           <- setNames(effectNamesDepvar, compositeNames)
    matchIdx <- match(effectNamesDepvar, staticChangeContributions_effectNames)
    nEffects <- length(matchIdx)
    nEgos <- dim(data[["depvars"]][[dv]])[1]
    netType <- unique(effectNetTypes[effectDepvars == dv])
    if (netType == "oneMode") {
      nChoices <- dim(data[["depvars"]][[dv]])[2]
    } else if (netType == "behavior") {
      nChoices <- 3
    } else {
      stop("Unsupported netType: ", netType)
    }
    for (g in seq_along(staticChangeContributions)) {
      groupContributions <- staticChangeContributions[[g]]
      nPeriods <- length(groupContributions)
      result <- array(NA_real_, dim = c(nPeriods, nEffects, nEgos, nChoices),
                      dimnames = list(
                        period = seq_len(nPeriods),
                        effect = compositeNames,
                        ego = seq_len(nEgos),
                        choice = seq_len(nChoices)
                      ))
      for (p in seq_len(nPeriods)) {
        for (e in seq_len(nEffects)) {
          mat <- do.call(rbind, groupContributions[[p]][[matchIdx[e]]])
          result[p, e, , ] <- mat
        }
      }
      if (requireNamespace("data.table", quietly = TRUE)) {
        # To resolve R CMD checks not understanding data.table syntax
        ego <- choice <- period <- group <- networkName <- effecttype <- NULL
        DT <- data.table::as.data.table(as.table(result))
        data.table::setnames(DT, c("period", "effectname", "ego", "choice", "contribution"))
        DT[, effecttype := typeMap[effectname]]
        DT[, effectname := nameMap[effectname]]
        DT[, ego := as.integer(ego)]
        DT[, choice := as.integer(choice)]
        DT[, period := as.integer(period)]
        DT[, group := as.integer(g)]
        DT[, networkName := dv]
        data.table::setcolorder(DT, c("group", "period", "networkName",
                                      "ego", "choice", "effectname",
                                      "effecttype", "contribution"))
        results[[paste0("g", g, "_", dv)]] <- DT
      } else {
        df <- as.data.frame(as.table(result))
        names(df) <- c("period", "effectname", "ego", "choice", "contribution")
        df$effecttype <- typeMap[df$effectname]
        df$effectname <- nameMap[df$effectname]
        df$period <- as.integer(df$period)
        df$ego <- as.integer(df$ego)
        df$choice <- as.integer(df$choice)
        df$group <- as.integer(g)
        df$networkName <- dv
        df <- df[, c("group", "period", "networkName", "ego", "choice",
                     "effectname", "effecttype", "contribution"), drop = FALSE]
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

# widenStaticContribution <- function(changeContributions){
#   ## currently only works for dynamic case and data has to be pre filtered to only one depvar
#   if (all(c("effectname", "contribution") %in% names(changeContributions))) {
#     if (requireNamespace("data.table", quietly = TRUE)) {
#         changeContributions <- data.table::dcast(
#           changeContributions, group + period + ego + choice ~ effectname,
#           value.var = "contribution"
#         )
#       } else {
#         ## in the dynamic case, this might lose some information that we should keep
#         needed <- c("group", 
#           "period", 
#           "ego", 
#           "choice", 
#           "effectname", 
#           "contribution")
#         changeContributions <- as.data.frame(changeContributions)[, needed, drop = FALSE]
#         changeContributions <- reshape(
#           changeContributions,
#           idvar = c("group", "period", "ego", "choice"),
#           timevar = "effectname", v.names = "contribution",
#           direction = "wide"
#         )
#         colnames(changeContributions) <- sub("^contribution\\.", "", 
#           colnames(changeContributions))
#       }
#     }
#   changeContributions
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
                                          returnWide = FALSE,
                                          batch = TRUE,
                                          silent = TRUE,
                                          seed = NULL) {
  if(returnWide && returnDataFrame) stop("Can not return wide data.frame")
  # Set silent according to batch if not explicitly set

  # Flexible argument handling allows either extracted pre
  # calculated contributions from ans or simulating them
  if (useChangeContributions && (is.null(ans) || is.null(ans$changeContributions))) 
  {
            warning("useChangeContributions=TRUE, but 'ans' does not 
			      contain 'changeContributions'. Will rerun siena07 to generate them.")
            useChangeContributions <- FALSE
  }
  if (useChangeContributions && (!returnDataFrame || returnWide) && 
    (is.data.frame(ans$changeContributions) || 
      is.data.frame(ans$changeContributions[[1]]))) 
  {
            stop("useChangeContributions=TRUE but 'ans$changeContributions' contains 
              data.frames; cannot use for returnDataFrame=FALSE or returnWide=TRUE")
  }
  if (useChangeContributions && is.null(n3)) 
  {
    warning("n3 cannot be changed when useChangeContributions=TRUE; 
    using the number of chains stored in 'ans'.")
  }
  if(!useChangeContributions) # could be extracted to a separate function
  {    ## algorithm should not be necessary anymore when using new siena() function
    if (is.null(algorithm))
    {
      stop("Must provide 'algorithm' when useChangeContributions=FALSE.")
    } 
    else
    {
      if (!inherits(algorithm, "sienaAlgorithm") && 
        !inherits(algorithm, "sienaAlgorithmSettings")) 
      {
        stop("algorithm is not a legitimate Siena algorithm specification")
      }
    }
    if(is.null(data)) 
    {
      stop("Must provide 'data' when useChangeContributions=FALSE ")
    } else 
    {
      if (!inherits(data, "sienadata") && !inherits(data, "siena")) {
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

    if (is.null(silent)) silent <- batch
    ans <- siena07(
        algorithm, 
        data=data, 
        effects=effects,
        initC=FALSE, # why?
        returnDeps=FALSE, 
        returnChangeContributions=TRUE, 
        returnDataFrame = returnDataFrame,
        batch = batch,
        silent = silent)
    }

  changeContributions <- ans$changeContributions

if (!returnDataFrame && !returnWide) {
  if (!is.null(depvar) && is.list(changeContributions) && length(changeContributions) > 0) {
    # Structure: [[nit]][[group]][[period]][[ministep_matrix]]
    # Must descend to the ministep level to filter by networkName attribute
    changeContributions <- lapply(changeContributions, function(chain_nit)
      lapply(chain_nit, function(chain_group)
        lapply(chain_group, function(chainPeriod)
          Filter(function(ms) {
            nm <- attr(ms, "networkName")
            is.null(nm) || as.character(nm) %in% depvar
          }, chainPeriod)
        )
      )
    )
  }
  return(changeContributions)
}
  if (returnWide) {
    if (is.null(effects))
      stop("'effects' must be provided when returnWide = TRUE")
    if (is.null(depvar))
      stop("'depvar' must be provided when returnWide = TRUE (use a single depvar name)")
      effsNoRate     <- effects[effects$include & effects$type != "rate", ]
      thetaIdx       <- which(effsNoRate$name == depvar)
      bareNames      <- effsNoRate$shortName[thetaIdx]
      compositeNames <- paste(bareNames, effsNoRate$type[thetaIdx], sep = "_")
      firstMs        <- findFirstMinistepForDepvar(changeContributions, depvar)
      matchIdx       <- matchRawEffectIdx(compositeNames, firstMs)
      result         <- flattenContributionsWide(changeContributions, compositeNames, matchIdx, depvar)
      result$thetaIdx    <- thetaIdx
      return(result)
  }
  { # returnDataFrame = TRUE
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
        # in current implementation, each chain contains contributions for only one depvar, 
        # but we keep the filtering here in case this changes in the future
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

# Not used currently
# ## Not doing this might reduce run time quite a bit! -> code would need to be adjusted though
# widenDynamicContribution <- function(changeContributions){
#   ## currently only works for dynamic case and data has to be pre filtered to only one depvar
#   if (all(c("effectname", "contribution") %in% names(changeContributions))) {
#     if (requireNamespace("data.table", quietly = TRUE)) {
#         changeContributions <- data.table::dcast(
#           changeContributions, chain + group + period + ministep + choice ~ effectname,
#           value.var = "contribution"
#         )
#       } else {
#         ## in the dynamic case, this might lose some information that we should keep
#         needed <- c("chain", 
#           "group", 
#           "period", 
#           "ministep", 
#           "choice", 
#           "effectname", 
#           "contribution")
#         changeContributions <- as.data.frame(changeContributions)[, needed, drop = FALSE]
#         changeContributions <- reshape(
#           changeContributions,
#           idvar = c("chain", "group", "period", "ministep", "choice"),
#           timevar = "effectname", v.names = "contribution",
#           direction = "wide"
#         )
#         colnames(changeContributions) <- sub("^contribution\\.", "", 
#           colnames(changeContributions))
#       }
#     }
#   changeContributions
# }

#' R wrapper for C++ flattenChangeContributionsWide.
#' Converts raw nested ministep list → wide list_wide struct:
#'   contribMat [N × nEffects], chain, group, period, ministep, choice, group_id
#' matchIdx: 0-based integer indices into raw ministep matrix rows.
flattenContributionsWide <- function(changeContributions, effectNames,
                                     matchIdx, depvar = NULL) {
  changeContributions <- .Call(C_flattenChangeContributionsWide,
              changeContributions,
              as.character(effectNames),
              as.integer(matchIdx),
              if (is.null(depvar)) NULL else as.character(depvar))
  changeContributions$effectNames <- as.character(effectNames)
  changeContributions
}

# really better in R but not in c++?

# Find the first ministep matrix for a given depvar in the nested list.
findFirstMinistepForDepvar <- function(changeContributions, depvar) {
  for (ch in changeContributions)
    for (grp in ch)
      for (per in grp)
        for (ms in per) {
          nn <- attr(ms, "networkName")
          if (is.null(nn) || as.character(nn) %in% depvar) return(ms)
        }
  NULL
}

# Match compositeNames (shortName_type) against attrs of first ministep.
# Returns 0-based integer indices for C++.
matchRawEffectIdx <- function(compositeNames, firstMs) {
  rawNames <- attr(firstMs, "effectNames")
  rawTypes <- attr(firstMs, "effectTypes")
  rawComp  <- paste(rawNames, rawTypes, sep = "_")
  idx      <- match(compositeNames, rawComp) - 1L
  if (anyNA(idx))
    stop("Effects not found in ministep matrix: ",
         paste(compositeNames[is.na(idx)], collapse = ", "))
  idx
}

# Wrapper to use in R code (RCPP style)
flattenChangeContributionsList <- function(x) {
  .Call(C_flattenChangeContributionsList, x)
}
