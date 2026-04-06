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

# ---- User-facing change statistics ------------------------------------------

# Extract change statistics (sign-adjusted contributions) as a data.frame.
#
# Contributions from the SAOM represent the difference in the statistic
# when switching from "no change" to a specific choice. For the density
# effect the "no change" row has contribution -1; for all other effects
# the sign is flipped on those rows so that the values represent genuine
# change statistics.
#
# ans:   A sienaFit object.
# data:  A siena / sienadata object.
# effects: Optional sienaEffects; defaults to ans$requestedEffects.
# depvar:  Character: dependent variable name(s). NULL = all.
# dynamic: Logical: if TRUE, simulate chains via getDynamicChangeContributions.
# algorithm: Required when dynamic = TRUE.
# n3:    Number of chains for dynamic simulation.
# useChangeContributions: Logical: reuse stored contributions in ans.
# flip:  Logical: if TRUE (default), apply sign adjustment so values
#        represent change statistics rather than raw contributions.
# Returns a data.frame (or data.table) with columns for
#   period, ego, choice (static) or chain, ministep, choice (dynamic),
#   plus one column per effect containing the (possibly sign-adjusted) values.

getChangeStatistics <- function(ans, data, effects = NULL, depvar = NULL,
                                 dynamic = FALSE, algorithm = NULL,
                                 n3 = 500, useChangeContributions = FALSE,
                                 flip = TRUE) {
    if (is.null(effects)) effects <- ans$requestedEffects
    if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]

    if (dynamic) {
        if (is.null(algorithm))
            stop("'algorithm' must be provided when dynamic = TRUE")
        changeStats <- getDynamicChangeContributions(
            ans = ans, theta = ans$theta, data = data,
            algorithm = algorithm, effects = effects,
            depvar = depvar, n3 = n3,
            useChangeContributions = useChangeContributions,
            returnWide = TRUE
        )
    } else {
        changeStats <- getStaticChangeContributions(
            ans = ans, data = data, effects = effects,
            depvar = depvar, returnWide = TRUE
        )
    }

    out <- groupColsList(changeStats)
    out <- attachContributions(out, changeStats$effectNames,
                                changeStats$contribMat, flip = flip)
    attr(out, "row.names") <- .set_row_names(length(out[[1L]]))
    class(out) <- "data.frame"
    out
}

sienaSetupModelOptionsForCpp <- function(ans = NULL, algorithm = NULL, data, effects, pData, pModel,
                                         profileData = FALSE, parallelrun = FALSE) {
  if (is.null(algorithm) && !is.null(ans)) {
    algorithm <- ans$x
  }
  if (!is.null(algorithm)) {
    MAXDEGREE <- if (is.null(algorithm$MaxDegree) || all(algorithm$MaxDegree == 0)) {
      NULL
    } else {
      tmp <- algorithm$MaxDegree
      storage.mode(tmp) <- "integer"
      tmp
    }
    UNIVERSALOFFSET <- if (is.null(algorithm$UniversalOffset) || all(algorithm$UniversalOffset == 0)) {
      NULL
    } else {
      tmp <- algorithm$UniversalOffset
      storage.mode(tmp) <- "double"
      tmp
    }
    MODELTYPE <- if (is.null(algorithm$modelType) || all(algorithm$modelType == 0)) NULL else {
      tmp <- algorithm$modelType
      storage.mode(tmp) <- "integer"
      tmp
    }
    BEHMODELTYPE <- if (is.null(algorithm$behModelType) || all(algorithm$behModelType == 0)) NULL else {
      tmp <- algorithm$behModelType
      storage.mode(tmp) <- "integer"
      tmp
    }
    if (!is.null(algorithm$cconditional) && algorithm$cconditional) {
      if (!is.null(algorithm$condname) && nzchar(algorithm$condname)) {
        CONDVAR <- algorithm$condname
      } else if (!is.null(algorithm$condvarno) && length(algorithm$condvarno) > 0 &&
                 algorithm$condvarno >= 1 && algorithm$condvarno <= length(data$depvars)) {
        CONDVAR <- names(data$depvars)[algorithm$condvarno]
      } else {
        CONDVAR <- NULL
      }

      if (!is.null(CONDVAR)) {
        if (inherits(data, "sienaGroup")) {
          CONDTARGET <- as.integer(sapply(data, function(xx)
            attr(xx$depvars[[CONDVAR]], "distance")))
        } else {
          CONDTARGET <- as.integer(attr(data$depvars[[CONDVAR]], "distance"))
        }
      } else {
        CONDTARGET <- NULL
      }
    } else {
      CONDVAR <- NULL
      CONDTARGET <- NULL
    }
    normSetRates <- if (!is.null(algorithm$normSetRates)) as.logical(algorithm$normSetRates) else FALSE
  } else {
    MAXDEGREE <- NULL
    UNIVERSALOFFSET <- NULL
    MODELTYPE <- NULL
    BEHMODELTYPE <- NULL
    CONDVAR <- NULL
    CONDTARGET <- NULL
    normSetRates <- FALSE
  }

  simpleRates <- TRUE
  if (!is.null(effects) && any(!effects$basicRate & effects$type == "rate")) {
    simpleRates <- FALSE
  }

  .Call(C_setupModelOptions, PACKAGE = pkgname,
        pData, pModel,
        MAXDEGREE, UNIVERSALOFFSET,
        CONDVAR, CONDTARGET,
        profileData, parallelrun,
        MODELTYPE, BEHMODELTYPE,
        simpleRates, normSetRates)
}

getStaticChangeContributions <- function(ans = NULL,
                                         data,
                                         effects = NULL,
                                         depvar = NULL,
                                         algorithm = NULL,
                                         includePermitted = FALSE,
                                         returnDataFrame = FALSE,
                                         returnWide = FALSE) {
  # Flexible argument handling
  if (!is.null(ans)) {
    if (is.null(effects)) effects <- ans$requestedEffects 
  }
  if (is.null(ans) && is.null(algorithm)) {
    warning("Neither 'ans' nor 'algorithm' provided to getStaticChangeContributions; ",
            "\ model options (maxDegree, conditional, etc.) are not set and defaults are used.")
  }
  if (is.null(data) || is.null(effects)) {
    stop("Must provide either 'ans' or all of 'data' and 'effects'.")
  }
  noRate <- effects$type != "rate"
  effects <- effects[noRate, ]

  # Interaction effects (unspInt, behUnspInt, contUnspInt) reference component
  # effects via effect1/effect2/effect3 (effectNumber).  The full effects object
  # (from getEffects / includeEffects / includeInteraction) contains these as
  # rows with include=FALSE.  ans$requestedEffects does not — it only keeps
  # included rows — so sienaSetupEffectsForCpp cannot resolve the components.
  intShort <- c("unspInt", "behUnspInt", "contUnspInt")
  intIdx   <- effects$shortName %in% intShort
  if (any(intIdx)) {
    compNos <- unique(c(effects$effect1[intIdx],
                        effects$effect2[intIdx],
                        effects$effect3[intIdx]))
    compNos <- compNos[compNos > 0]
    missing <- compNos[!compNos %in% effects$effectNumber]
    if (length(missing) > 0) {
      stop("The effects object is missing component effects (effectNumber ",
           paste(missing, collapse = ", "), ") needed by interaction effects.\n",
           "This happens when using ans$requestedEffects (the default) ",
           "instead of the full effects object.\n",
           "Please pass the full effects object.")
    }
  }

  depvars_available <- names(data$depvars)
  if (is.null(depvar)) {
    depvar <- depvars_available
  }

  # Prepare for C++ call
  pData <- sienaSetupDataForCpp(data,
                                includeBehavior = TRUE,
                                includeBipartite = FALSE)
  setup <- sienaSetupEffectsForCpp(pData, data, effects)

  # Propagate model options to C++ (maxDegree, universalOffset, conditional, etc.)
  sienaSetupModelOptionsForCpp(ans = ans, algorithm = algorithm, data = data, effects = effects,
                               pData = pData, pModel = setup$pModel,
                               profileData = FALSE, parallelrun = FALSE)

  staticChangeContributions <- .Call(
    C_getStaticChangeContributions, 
    PACKAGE = pkgname,
    pData, setup$pModel, setup$myeffects,
    parallelrun = FALSE
  )

  if (!returnDataFrame && !returnWide) {
    if (includePermitted) {
      return(list(contributions = staticChangeContributions,
                  permitted = attr(staticChangeContributions, "permitted")))
    }
    return(staticChangeContributions)
  }

  # Derive effect metadata in C++ output order.
  # sienaSetupEffectsForCpp uses split(effects, effects$name) which sorts
  # dep var names alphabetically.  The C++ iterates over this split list,
  # so its output order is: effects for the alphabetically-first dep var,
  # then the next, etc.  We must use the same ordering here.
  cppEffects <- do.call(rbind, setup$myeffects)
  effectNames            <- numberIntShortNames(cppEffects$shortName)
  effectDepvars          <- cppEffects$name
  effectNetTypes         <- cppEffects$netType
  effectTypes            <- cppEffects$type
  effectInteractionTypes <- cppEffects$interactionType
  # Enrich shortNames with covariate identifiers so e.g. two egoX effects
  # (for different covariates) get unique names: egoX_gender, egoX_age.
  effectCovarSuffixes <- effectCovarSuffix(
    if (!is.null(cppEffects$interaction1)) cppEffects$interaction1 else rep("", length(effectNames)),
    if (!is.null(cppEffects$interaction2)) cppEffects$interaction2 else rep("", length(effectNames))
  )
  effectNamesEnriched <- ifelse(effectCovarSuffixes == "", effectNames,
                                paste(effectNames, effectCovarSuffixes, sep = "_"))

  staticChangeContributions_effectNames <- attr(staticChangeContributions,
    "effectNames")

  if (returnWide) {
    # Accumulate across all groups for one depvar into a single unified struct.
    if (length(depvar) != 1L)
      stop("returnWide = TRUE requires exactly one depvar; got ",
           length(depvar), ". Call once per depvar.")
    all_contribMat <- list()
    all_permitted <- list()
    all_period <- list()
    all_ego <- list()
    all_choice <- list()
    all_group_id <- list()
    all_group <- list()
    effectNames_out <- NULL
    group_id_offset <- 0L

    permData <- attr(staticChangeContributions, "permitted")
    for (dv in depvar) {
      effectIdx <- which(effectDepvars == dv)
      effectNamesDv <- effectNamesEnriched[effectIdx]
      compositeNames <- paste(dv, effectNamesDv, effectTypes[effectIdx], sep = "_")
      nEffects <- length(effectIdx)
      nEgos <- dim(data[["depvars"]][[dv]])[1]
      netType <- unique(effectNetTypes[effectDepvars == dv])
      nChoices <- if (netType == "oneMode") dim(data[["depvars"]][[dv]])[2L]
                        else if (netType == "behavior") 3L
                        else stop("Unsupported netType: ", netType)

      if (is.null(effectNames_out)) {
        effectNames_out <- compositeNames
      } else {
        effectNames_out <- c(effectNames_out, compositeNames)
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
        N <- nPeriods * nEgos * nChoices
        contribMat <- matrix(arr, nrow = N, ncol = nEffects)
        colnames(contribMat) <- compositeNames
        period_vec <- rep(seq_len(nPeriods), each = nEgos * nChoices)
        ego_vec <- rep(rep(seq_len(nEgos), each = nChoices), times = nPeriods)
        choice_vec <- rep(seq_len(nChoices), times = nPeriods * nEgos)
        group_id <- group_id_offset + rep(seq_len(nPeriods * nEgos), each = nChoices)

        all_contribMat[[length(all_contribMat) + 1L]] <- contribMat
        if (includePermitted) {
          # Permitted is per-alter, not per-effect; take first effect (all identical)
          permArr <- array(FALSE, dim = c(nPeriods, nEgos, nChoices))
          for (p in seq_len(nPeriods)) {
            permArr[p, , ] <- do.call(rbind,
              permData[[g]][[p]][[effectIdx[1]]])
          }
          permArr <- aperm(permArr, c(3L, 2L, 1L))
          all_permitted[[length(all_permitted) + 1L]] <- as.logical(permArr)
        }
        all_period[[length(all_period) + 1L]] <- period_vec
        all_ego[[length(all_ego) + 1L]]  <- ego_vec
        all_choice[[length(all_choice) + 1L]] <- choice_vec
        all_group_id[[length(all_group_id) + 1L]] <- group_id
        all_group[[length(all_group) + 1L]] <- rep(as.integer(g), N)
        group_id_offset <- group_id_offset + nPeriods * nEgos
      }
    }

    out <- list(
      contribMat  = do.call(rbind, all_contribMat),
      period = unlist(all_period, use.names = FALSE),
      ego = unlist(all_ego, use.names = FALSE),
      choice = unlist(all_choice, use.names = FALSE),
      group_id = unlist(all_group_id, use.names = FALSE),
      group = unlist(all_group, use.names = FALSE),
      effectNames = effectNames_out,
      effectInteractionTypes = setNames(effectInteractionTypes[
          match(effectNames_out,
                paste(effectDepvars, effectNamesEnriched, effectTypes, sep = "_"))],
          effectNames_out)
    )
    if (includePermitted) {
      out$permitted <- unlist(all_permitted, use.names = FALSE)
    }
    return(out)
  }

  results <- list()
  for (dv in depvar) {
    effectIdx <- which(effectDepvars == dv)
    effectNamesDepvar <- effectNamesEnriched[effectIdx]
    compositeNames <- paste(dv, effectNamesDepvar, effectTypes[effectIdx], sep = "_")
    typeMap <- setNames(effectTypes[effectIdx], compositeNames)
    nameMap <- setNames(effectNamesDepvar, compositeNames)
    nEffects <- length(effectIdx)
    nEgos <- dim(data[["depvars"]][[dv]])[1]
    netType <- unique(effectNetTypes[effectDepvars == dv])
    if (netType == "oneMode") {
      nChoices <- dim(data[["depvars"]][[dv]])[2]
    } else if (netType == "behavior") {
      nChoices <- 3
    } else {
      stop("Unsupported netType: ", netType)
    }
    permData <- attr(staticChangeContributions, "permitted")
    for (g in seq_along(staticChangeContributions)) {
      groupContributions <- staticChangeContributions[[g]]
      groupPermitted <- permData[[g]]
      nPeriods <- length(groupContributions)
      result <- array(NA_real_, dim = c(nPeriods, nEffects, nEgos, nChoices),
                      dimnames = list(
                        period = seq_len(nPeriods),
                        effect = compositeNames,
                        ego = seq_len(nEgos),
                        choice = seq_len(nChoices)
                      ))
      permResult <- array(FALSE, dim = c(nPeriods, nEffects, nEgos, nChoices),
                      dimnames = list(
                        period = seq_len(nPeriods),
                        effect = compositeNames,
                        ego = seq_len(nEgos),
                        choice = seq_len(nChoices)
                      ))
      for (p in seq_len(nPeriods)) {
        for (e in seq_len(nEffects)) {
          mat <- do.call(rbind, groupContributions[[p]][[effectIdx[e]]])
          result[p, e, , ] <- mat
          permResult[p, e, , ] <- do.call(rbind, groupPermitted[[p]][[effectIdx[e]]])
        }
      }
      # data.table path removed — using base R path below
      # if (requireNamespace("data.table", quietly = TRUE)) {
      #   # To resolve R CMD checks not understanding data.table syntax
      #   ego <- choice <- period <- group <- networkName <- effecttype <- NULL
      #   DT <- data.table::as.data.table(as.table(result))
      #   data.table::setnames(DT, c("period", "effectname", "ego", "choice", "contribution"))
      #   DT[, effecttype := typeMap[effectname]]
      #   DT[, effectname := nameMap[effectname]]
      #   DT[, ego := as.integer(ego)]
      #   DT[, choice := as.integer(choice)]
      #   DT[, period := as.integer(period)]
      #   DT[, group := as.integer(g)]
      #   DT[, networkName := dv]
      #   data.table::setcolorder(DT, c("group", "period", "networkName",
      #                                 "ego", "choice", "effectname",
      #                                 "effecttype", "contribution"))
      #   results[[paste0("g", g, "_", dv)]] <- DT
      # } else {
        df <- as.data.frame(as.table(result))
        names(df) <- c("period", "effectname", "ego", "choice", "contribution")
        df$effecttype <- typeMap[df$effectname]
        df$effectname <- nameMap[df$effectname]
        df$period <- as.integer(df$period)
        df$ego <- as.integer(df$ego)
        df$choice <- as.integer(df$choice)
        df$group <- as.integer(g)
        df$networkName <- dv
        df$permitted <- c(permResult)
        df <- df[, c("group", "period", "networkName", "ego", "choice",
                     "effectname", "effecttype", "contribution", "permitted"), drop = FALSE]
        results[[paste0("g", g, "_", dv)]] <- df
      # }
    }
  }
  # data.table removed — use do.call(rbind, ...)
  # if (requireNamespace("data.table", quietly = TRUE)) {
  #   data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
  # } else {
    do.call(rbind, results)
  # }
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
getDynamicChangeContributions <- function(
  ans = NULL, 
  theta = NULL, 
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
  if (useChangeContributions && (is.null(ans) || 
    is.null(ans$changeContributions))) {
      warning("useChangeContributions=TRUE, but 'ans' does not 
      contain 'changeContributions'. Will rerun siena07 to generate them.")
            useChangeContributions <- FALSE
  }
  if (useChangeContributions && (!returnDataFrame || returnWide) && 
    (is.data.frame(ans$changeContributions) || 
      is.data.frame(ans$changeContributions[[1]]))) {
        stop("useChangeContributions=TRUE but 'ans$changeContributions' contains 
          data.frames; cannot use for returnDataFrame=FALSE or returnWide=TRUE")
  }
  if (useChangeContributions && is.null(n3)) {
    warning("n3 cannot be changed when useChangeContributions=TRUE; 
    using the number of chains stored in 'ans'.")
  }
  if(!useChangeContributions) {   
    # could be extracted to a separate function
    if (is.null(algorithm)) {
       ## algorithm should not be necessary anymore when using new siena() function
      stop("Must provide 'algorithm' when useChangeContributions=FALSE.")
    } else {
      if (!inherits(algorithm, "sienaAlgorithm") && 
        !inherits(algorithm, "sienaAlgorithmSettings")) {
        stop("algorithm is not a legitimate Siena algorithm specification")
      }
    }
    if(is.null(data)) {
      stop("Must provide 'data' when useChangeContributions=FALSE ")
    } else {
      if (!inherits(data, "sienadata") && !inherits(data, "siena")) {
        stop("'data' is not a legitimate Siena data specification")
      }
    }
    if(is.null(effects)) {
      if (is.null(ans)) {
        stop("Must provide 'ans' or 'effects' when useChangeContributions=FALSE.")
      }
      if (ans$cconditional) {
        stop("Must provide 'effects' for conditional estimation")
      }
      effects <- ans$requestedEffects
      if (!inherits(effects, "sienaEffects")) {
        stop("effects is not a legitimate Siena effects object")
      }
    }
    # siena07 -> initializeFRAN -> sienaTimeFix always computes groupPeriods from
    # effects$period. Rate effects carry the period values; without them every
    # period is NA, max(..., na.rm=TRUE) returns -Inf, and the subsequent rep()
    # crashes. This can happen when effects was supplied (or defaulted) from
    # ans$requestedEffects of a conditionally-estimated model, where rate effects
    # are stripped out. Check here — after effects is resolved — so the guard
    # fires whether effects was NULL-defaulted or explicitly passed in.
    # later fix sienaTimeFix to not crash on this and add a more 
    # specific error message about missing rate effects.
    if ("basicRate" %in% names(effects) && !any(effects$basicRate)) {
      stop(
        "The 'effects' object contains no rate effects. This is required when ",
        "useChangeContributions = FALSE or estimating dynamic uncertainty,
        because siena07 is rerun internally.\n",
        "Please pass the full effects object, e.g.:\n",
        "  marginalEffects(..., dynamic = TRUE, effects = mymodel, ...)"
      )
    }
    if(is.null(theta)) {
      if (is.null(ans)) {
        stop("Must provide 'ans' or 'theta' when useChangeContributions=FALSE.")
      }
      theta <- ans$theta
      if (!is.numeric(theta)) {
        stop("theta is not a legitimate parameter vector")
      }
    }
    if(algorithm$maxlike) {
      # improve: restore all other settings from the provided algorithm object instead of using defaults
      warning("maxLike=TRUE is not compatible with dynamic change contributions; 
       creating algorithm with maxLike=FALSE; using default settings for all other parameters.")
      algorithm <- sienaAlgorithmCreate(projname = NULL, maxlike = FALSE)
    }
    algorithm$nsub <- 0
    algorithm$simOnly <- TRUE
    if (!is.null(n3)) {
      # should maybe be done differently
      algorithm$n3 <- as.integer(n3)
    }
    if(!is.null(seed)) {
      if(is.numeric(seed)) {
        algorithm$seed <- as.integer(seed)
      } else {
        warning("'seed' has to be of type 'numeric' \n used default settings")
      }
    }
    # Update effects initial values with provided theta to sample chains from
    # specificed parameter values in phase 3 of siena07
    allEffNames <- getNamesFromEffects(effects[effects$include, ])
    # Defensive backfill: if theta is named but shorter than allEffNames
    # (e.g. caller used coef() with dropRates=TRUE, or conditional model
    # where rate params live in ans$rate), fill in missing effects from
    # ans$theta or — for conditional models — from ans$rate spliced into
    # the condEffects rows that marginalEffects added to effects.
    # The preferred path is for callers to supply full theta, making
    # this a no-op safety net.
    if (!is.null(names(theta)) && length(theta) < length(allEffNames)) {
      if (!is.null(ans) && !is.null(ans$theta) &&
          length(ans$theta) == length(allEffNames)) {
        ansNamed <- setNames(ans$theta, allEffNames)
        missingIdx <- which(!allEffNames %in% names(theta))
        theta <- c(theta, ansNamed[allEffNames[missingIdx]])
      } else if (!is.null(ans) && isTRUE(ans$cconditional) &&
                 !is.null(ans$rate)) {
        # Conditional model: ans$theta covers non-rate effects;
        # ans$rate covers rate effects.  effects may already include
        # rate rows (e.g. via condEffects splice in marginalEffects).
        rateIdx <- which(effects$basicRate[effects$include])
        nonRateIdx <- which(!effects$basicRate[effects$include])
        if (length(rateIdx) > 0L && length(ans$rate) == length(rateIdx)) {
          rateNamed <- setNames(ans$rate, allEffNames[rateIdx])
          theta <- c(theta, rateNamed[setdiff(names(rateNamed), names(theta))])
        }
      }
    }
    inTheta <- allEffNames %in% names(theta)
    effects$initialValue[effects$include][inTheta] <- 
      theta[allEffNames[inTheta]]

    if (is.null(silent)) silent <- batch
    ans <- siena07(
        algorithm, 
        data=data, 
        effects=effects,
        nbrNodes=1,  # enforce single-threaded: outer mclapply already parallelises over draws
        initC=FALSE, # parent already called initializeFRAN (really?); workers always get initC=TRUE in robmon
        returnDeps=FALSE, 
        returnChangeContributions=TRUE, 
        returnDataFrame = returnDataFrame,
        batch = batch,
        silent = silent)
  }

  changeContributions <- ans$changeContributions

  if (!returnDataFrame && !returnWide) {
    if (!is.null(depvar) && 
      is.list(changeContributions) && length(changeContributions) > 0) {
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
    if (is.null(depvar))
      stop("'depvar' must be provided when returnWide = TRUE (use a single depvar name)")
    wide <- flattenChangeContributionsWide(changeContributions, depvar)
    # Enrich column names with covariate identifiers from the effects data frame,
    # matching positionally (C++ processes non-rate included effects sorted
    # alphabetically by name, then in effects data frame order within each depvar).
    if (!is.null(effects) && length(depvar) == 1L) {
      noRate <- effects[["type"]] != "rate"
      inc    <- effects[noRate & effects[["include"]], ]
      inc    <- inc[order(inc[["name"]]), ]   # alphabetical by depvar name
      inc_dv <- inc[inc[["name"]] == depvar, ]
      if (nrow(inc_dv) == ncol(wide$contribMat)) {
        sn <- numberIntShortNames(inc_dv[["shortName"]])
        i1 <- if (!is.null(inc_dv[["interaction1"]])) inc_dv[["interaction1"]] else rep("", nrow(inc_dv))
        i2 <- if (!is.null(inc_dv[["interaction2"]])) inc_dv[["interaction2"]] else rep("", nrow(inc_dv))
        cs <- effectCovarSuffix(i1, i2)
        snWithCovar    <- ifelse(cs == "", sn, paste(sn, cs, sep = "_"))
        enrichedNames  <- paste(depvar, snWithCovar, inc_dv[["type"]], sep = "_")
        colnames(wide$contribMat) <- enrichedNames
        wide$effectNames          <- enrichedNames
      }
    }
    # Attach effectInteractionTypes from the effects object if available
    if (!is.null(effects) && "interactionType" %in% names(effects)) {
      noRate <- effects$type != "rate"
      eff <- effects[noRate & effects$include, ]
      iTypes <- eff[eff$name == depvar, "interactionType"]
      names(iTypes) <- wide$effectNames
      wide$effectInteractionTypes <- iTypes
    }
    return(wide)
  }
  if(returnDataFrame == TRUE) {
    if (is.list(changeContributions) && 
      !is.data.frame(changeContributions[[1]])) {
      warning("SIENA backend returned nested lists, not data.frames; 
      auto-flattening each chain.")
      changeContributions <- lapply(changeContributions, 
        flattenChangeContributionsList)
    }
  # filtering before binding proofed to be quite slow
    # data.table path removed — using base R path below
    # if (requireNamespace("data.table", quietly = TRUE)) {
    #   changeContributions <- data.table::rbindlist(
    #     lapply(seq_along(changeContributions), function(i) {
    #       df <- changeContributions[[i]]
    #       df$chain <- i
    #       df
    #     }),
    #     use.names = TRUE, fill = TRUE
    #   )
    #   if (!is.null(depvar)) {
    #     # in current implementation, each chain contains contributions for only one depvar, 
    #     # but we keep the filtering here in case this changes in the future
    #     changeContributions <- changeContributions[changeContributions[["networkName"]] 
    #       %in% depvar, , drop=FALSE]
    #   }
    # } else {
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
    # }
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

# R wrapper for C++ flattenChangeContributionsWide.
# Converts raw nested ministep list -> wide list_wide struct:
#   contribMat [N x nEffects], chain, group, period, ministep, choice, group_id
flattenChangeContributionsWide <- function(changeContributions, depvar = NULL) {
  result <- .Call(C_flattenChangeContributionsWide,
              changeContributions,
              if (is.null(depvar)) NULL else as.character(depvar))
  if (!is.null(depvar) && length(depvar) == 1L) {
    cols <- numberIntColNames(colnames(result$contribMat))
    colnames(result$contribMat) <- paste(depvar, cols, sep = "_")
  }
  result$effectNames <- colnames(result$contribMat)
  result
}

# Wrapper to use in R code (RCPP style)
flattenChangeContributionsList <- function(x) {
  .Call(C_flattenChangeContributionsList, x)
}

