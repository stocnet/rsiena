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

  staticChangeContributions <- .Call(C_getStaticChangeContributions, PACKAGE = pkgname,
                                     pData, setup$pModel, setup$myeffects,
                                     parallelrun = FALSE)
  #attr(staticChangeContributions, "effectNames") <- effects$shortName

  if (!returnDataFrame && !returnWide) {
    return(staticChangeContributions)
  }

  staticChangeContributions_effectNames <- attr(staticChangeContributions, "effectNames")

  if (returnWide) {
    results <- list()
    for (dv in depvar) {
      effect_idx     <- which(effectDepvars == dv)
      effectNames_dv <- effectNames[effect_idx]
      match_idx      <- match(effectNames_dv, staticChangeContributions_effectNames)
      n_effects      <- length(match_idx)
      n_egos  <- dim(data[["depvars"]][[dv]])[1]
      netType <- unique(effectNetTypes[effectDepvars == dv])
      n_choices <- if (netType == "oneMode") dim(data[["depvars"]][[dv]])[2L]
                   else if (netType == "behavior") 3L
                   else stop("Unsupported netType: ", netType)
      for (g in seq_along(staticChangeContributions)) {
        n_periods <- length(staticChangeContributions[[g]])
        arr <- array(NA_real_, dim = c(n_periods, n_effects, n_egos, n_choices))
        for (p in seq_len(n_periods)) {
          for (e in seq_len(n_effects)) {
            arr[p, e, , ] <- do.call(rbind,
              staticChangeContributions[[g]][[p]][[match_idx[e]]])
          }
        }
        arr <- aperm(arr, c(4L, 3L, 1L, 2L))
        N   <- n_periods * n_egos * n_choices
        contrib_mat <- matrix(arr, nrow = N, ncol = n_effects)
        colnames(contrib_mat) <- effectNames_dv
        period_vec <- rep(seq_len(n_periods), each = n_egos * n_choices)
        ego_vec    <- rep(rep(seq_len(n_egos), each = n_choices), times = n_periods)
        choice_vec <- rep(seq_len(n_choices), times = n_periods * n_egos)
        group_id   <- rep(seq_len(n_periods * n_egos), each = n_choices)
        results[[paste0("g", g, "_", dv)]] <- list(
          contrib_mat = contrib_mat,
          period      = period_vec,
          ego         = ego_vec,
          choice      = choice_vec,
          group_id    = group_id,
          n_choices   = as.integer(n_choices),
          group       = as.integer(g),
          networkName = dv,
          effectNames = effectNames_dv
        )
      }
    }
    return(results)
  }

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
    return(changeContributions)
  }
  if (returnWide) {
    if (is.null(effects))
      stop("'effects' must be provided when returnWide = TRUE")
    if (is.null(depvar))
      stop("'depvar' must be provided when returnWide = TRUE (use a single depvar name)")
    effectNames <- getEffectNamesNoRate(effects, depvar)
    return(flattenContributionsWide(changeContributions, effectNames, depvar))
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

# Wrapper to use in R code (RCPP style)
flattenChangeContributionsList <- function(x) {
  .Call("C_flattenChangeContributionsList", x)
}
