##@predict.sienaFit
predict.sienaFit <- function(
    object,
    newdata,
    type = c("changeProb", "tieProb"),
    newParams = NULL,
    effects = NULL,
    depvar = NULL,
    dynamic = FALSE,
    algorithm = NULL,
    n3 = 1000,
    useChangeContributions = FALSE,
    level = "period",
    condition = NULL,
    sum.fun = mean,
    na.rm = TRUE,
    uncertainty = TRUE,
    uncertainty_mode = c("batch", "stream"),
    useCluster = FALSE,
    nbrNodes = 1,
    nsim = 1000,
    uncertainty_sd = TRUE,
    uncertainty_ci = TRUE,
    uncertainty_probs = c(0.025, 0.5, 0.975),
    uncertainty_mcse = FALSE,
    uncertainty_mcse_batches = NULL,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batch_dir = "temp",
    prefix = "simBatch_b",
    combine_batch = TRUE,
    batch = TRUE,
    silent = NULL,
    batch_size = NULL,
    keep_batch = FALSE,
    verbose = TRUE,
    memory_scale = NULL,
    batch_unit_budget = 5e8,
    dynamic_ministep_factor = 10,
    ...
) {
    if (inherits(newdata, "sienaGroup"))
      stop("predict.sienaFit does not support multi-group data (sienaGroup).")
    uncertainty_mode <- match.arg(uncertainty_mode)
    type             <- match.arg(type)
    if (is.null(depvar)) depvar <- names(newdata[["depvars"]])[1]
    if (dynamic && is.null(algorithm)) stop("'algorithm' must be provided when dynamic = TRUE")
    if (dynamic && is.null(silent)) silent <- batch
    if (is.null(batch_size)) {
        batch_size <- plan_batch(
          data = newdata, 
          depvar = depvar, 
          nsim = nsim,
          nbrNodes = nbrNodes, 
          useCluster = useCluster,
          dynamic = dynamic, 
          n3 = n3,
          unit_budget = batch_unit_budget,
          dynamic_ministep_factor = dynamic_ministep_factor,
          memory_scale = memory_scale
        )
    }
    if (dynamic) {
    predictFun  <- predictProbabilityDynamic
    predictArgs <- list(
        ans                    = object,
        data                   = newdata,
        algorithm              = algorithm,
        effects                = effects,
        type                   = type,
        depvar                 = depvar,
        n3                     = n3,
        useChangeContributions = useChangeContributions,
        batch                  = batch,
        silent                 = silent
    )
    } else {
    # Build the contribution matrix ONCE — reused across all theta draws
    staticContributions <- getStaticChangeContributions(
        ans     = object,
        data    = newdata,
        effects = effects,
        depvar  = depvar,
        returnWide = TRUE
    )
    predictFun  <- predictProbabilityStatic
    predictArgs <- list(
        ans      = object,
        staticContributions = staticContributions,
        type     = type
    )
    }

    sienaPostestimate(
    predictFun               = predictFun,
    predictArgs              = predictArgs,
    outcome                  = type,
    level                    = level,
    condition                = condition,
    sum.fun                  = sum.fun,
    na.rm                    = na.rm,
    theta_hat                = object[["theta"]],
    cov_theta                = object[["covtheta"]],
    uncertainty              = uncertainty,
    uncertainty_mode         = uncertainty_mode,
    nsim                     = nsim,
    uncertainty_sd           = uncertainty_sd,
    uncertainty_ci           = uncertainty_ci,
    uncertainty_probs        = uncertainty_probs,
    uncertainty_mcse         = uncertainty_mcse,
    uncertainty_mcse_batches = uncertainty_mcse_batches,
    useCluster               = useCluster,
    nbrNodes                 = nbrNodes,
    clusterType              = clusterType,
    cluster                  = cluster,
    batch_dir                = batch_dir,
    prefix                   = prefix,
    combine_batch            = combine_batch,
    batch_size               = batch_size,
    keep_batch               = keep_batch,
    verbose                  = verbose,
    useChangeContributions   = if (dynamic) useChangeContributions else NULL
    )
}

plan_batch <- function(data, depvar, nsim,
                       nbrNodes = 1L, useCluster = FALSE,
                       dynamic = FALSE, n3 = NULL,
                       unit_budget = 2.5e8,
                       dynamic_ministep_factor = 10,
                       memory_scale = NULL) {
  dv     <- data$depvars[[depvar]]
  dv_dim <- dim(dv)
  n_ego  <- if (!is.null(dv_dim) && length(dv_dim) >= 1L) as.integer(dv_dim[1]) else as.integer(length(dv))
  n_choice   <- if (!is.null(dv_dim) && length(dv_dim) >= 2L) as.integer(dv_dim[2]) else 1L
  n_per  <- if (!is.null(dv_dim) && length(dv_dim) >= 3L) max(1L, as.integer(dv_dim[3] - 1L)) else 1L
  units  <- as.numeric(n_ego) * as.numeric(n_choice) * as.numeric(n_per)

  effective_workers <- if (isTRUE(useCluster)) max(1L, as.integer(nbrNodes)) else 1L

  if (dynamic) {
    n3_val <- if (is.null(n3)) 1L else max(1L, as.integer(n3))
    units_per_call <- units * as.numeric(dynamic_ministep_factor) * as.numeric(n3_val)
  } else {
    n3_val <- 1L
    units_per_call <- units
  }


  units_per_agg <- max(1.0, units * as.numeric(n3_val) *
                              if (dynamic) as.numeric(dynamic_ministep_factor) else 1.0)

  budget_for_agg <- as.numeric(unit_budget) - effective_workers * units_per_call

  if (budget_for_agg <= 0) {
    if (effective_workers > 1L) {
      warning(sprintf(
        "Memory budget (%.0f units) may be insufficient for %d parallel worker(s) at %.0f units each. Consider reducing nbrNodes or increasing batch_unit_budget.",
        unit_budget, effective_workers, units_per_call
      ))
    }

    max_batch <- max(1L, as.integer(floor(as.numeric(unit_budget) / units_per_agg)))
  } else {
    max_batch <- max(1L, as.integer(floor(budget_for_agg / units_per_agg)))
  }

  batch_raw <- min(as.integer(nsim), max_batch)

  if (!is.null(memory_scale) && as.integer(memory_scale) > 1L)
    batch_raw <- max(1L, as.integer(floor(batch_raw / as.integer(memory_scale))))

  k <- if (isTRUE(useCluster) && as.integer(nbrNodes) > 1L) as.integer(nbrNodes) else 1L
  if (k > 1L && as.integer(nsim) >= k) {
    b2 <- (batch_raw %/% k) * k
    batch_raw <- if (b2 < k) k else b2
  }
  min(max(1L, batch_raw), as.integer(nsim))
}

#' Per-theta static prediction using staticContributions matrices.
predictProbabilityStatic <- function(ans,
                                        staticContributions,
                                        theta,
                                        type = "changeProb") {
  effectNames <- staticContributions[[1]]$effectNames
  thetaNoRate <- alignThetaNoRate(theta, effectNames, ans)

  results <- vector("list", length(staticContributions))
  for (i in seq_along(staticContributions)) {
    staticCont        <- staticContributions[[i]]
    theta_use <- thetaNoRate[staticCont$effectNames]
    utility <- calculateUtility(staticCont$contrib_mat, theta_use)
    prob    <- as.numeric(softmax_arma_by_group(utility, staticCont$group_id))
    out <- data.frame(.groupColsList(staticCont),
                      changeUtil = utility, changeProb = prob,
                      stringsAsFactors = FALSE)
    out <- .attachContribColumns(out, staticCont$effectNames, staticCont$contrib_mat,
                                flip = TRUE)
    tp <- .computeTieProb(prob, staticCont$contrib_mat[, "density"], type)
    if (!is.null(tp)) out[["tieProb"]] <- tp

    if (requireNamespace("data.table", quietly = TRUE)) data.table::setDT(out)
    results[[i]] <- out
  }
  if (requireNamespace("data.table", quietly = TRUE))
    data.table::rbindlist(results, use.names = TRUE)
  else {
    out <- do.call(rbind, results)
    rownames(out) <- NULL
    out
  }
}


# ---- Dynamic -------------------------------------------------------

#' Per-theta dynamic prediction using matrix-internal pipeline.
predictProbabilityDynamic <- function(ans, data, theta, algorithm, effects,
                                         type = "changeProb", depvar,
                                         n3 = NULL,
                                         useChangeContributions = FALSE,
                                         batch = TRUE, silent = TRUE) {
  dynamicChangeContributions <- getDynamicChangeContributions(
    ans = ans, theta = theta, data = data, algorithm = algorithm,
    effects = effects, depvar = depvar, n3 = n3,
    useChangeContributions = useChangeContributions,
    returnWide = TRUE, batch = batch, silent = silent
  )
  thetaNoRate <- alignThetaNoRate(theta, getEffectNamesNoRate(effects, depvar), ans)
  theta_use   <- thetaNoRate[dynamicChangeContributions$effectNames]
  utility <- calculateUtility(dynamicChangeContributions$contrib_mat, theta_use)
  prob    <- as.numeric(softmax_arma_by_group(utility, dynamicChangeContributions$group_id))

  # Build output table — all computation done in vectors above
  out <- data.frame(.groupColsList(dynamicChangeContributions),
                    changeUtil = utility, changeProb = prob,
                    stringsAsFactors = FALSE)
  tp <- .computeTieProb(prob, dynamicChangeContributions$contrib_mat[, "density"], type)
  if (!is.null(tp)) out[["tieProb"]] <- tp
  out <- .attachContribColumns(out, dynamicChangeContributions$effectNames, dynamicChangeContributions$contrib_mat, flip = TRUE)
  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDT(out)
  out
}

#' R wrapper for C++ flattenChangeContributionsWide.
#' Converts raw nested ministep list → wide list_wide struct:
#'   contrib_mat [N × n_effects], chain, group, period, ministep, choice, group_id
flattenContributionsWide <- function(changeContributions, effectNames,
                                     depvar = NULL) {
  changeContributions <- .Call(C_flattenChangeContributionsWide,
              changeContributions,
              as.character(effectNames),
              if (is.null(depvar)) NULL else as.character(depvar))
  changeContributions$effectNames <- as.character(effectNames)
  changeContributions
}

# ---- Shared helpers -------------------------------------------------

calculateUtility <- function(mat, theta) {
  stopifnot(is.matrix(mat))
  as.numeric(mat %*% theta)
}

.computeTieProb <- function(prob, density_col, type) {
  if (type != "tieProb") return(NULL)
  prob <- prob
  neg1 <- density_col == -1L
  prob[neg1] <- 1 - prob[neg1]
  prob[density_col == 0L] <- NA
  prob
}


# ===========================================================================
# === V1 BACKUP: predict.sienaFit + static/dynamic per-theta functions    ===
# === (data.frame-internal pipeline, superseded by V2 above)              ===
# ===========================================================================

# ##@predict.sienaFit
# predict.sienaFit <- function(
#     object,
#     newdata,
#     type = c("changeProb", "tieProb"), #behavior not implemented *here*
#     newParams = NULL,
#     effects = NULL,
#     depvar = NULL,
#     dynamic = FALSE,
#     algorithm = NULL,
#     n3 = 1000,
#     useChangeContributions = FALSE,
#     level = "period",
#     condition = NULL,
#     sum.fun = mean,
#     na.rm = TRUE,
#     uncertainty = TRUE,
#     uncertainty_mode = c("batch", "stream"),
#     useCluster = FALSE,
#     nbrNodes = 1,
#     nsim = 1000,
#     uncertainty_sd = TRUE,
#     uncertainty_ci = TRUE,
#     uncertainty_probs = c(0.025, 0.5, 0.975),
#     uncertainty_mcse = FALSE,
#     uncertainty_mcse_batches = NULL,
#     clusterType = c("PSOCK", "FORK"),
#     cluster = NULL,
#     batch_dir = "temp",
#     prefix = "simBatch_b",
#     combine_batch = TRUE,
#     batch = TRUE,
#     silent = NULL,
#     batch_size = NULL,
#     keep_batch = FALSE,
#     verbose = TRUE,
#     memory_scale = NULL,
#     batch_unit_budget = 5e6,
#     dynamic_ministep_factor = 10,
#     ...
#   ){
#       uncertainty_mode <- match.arg(uncertainty_mode)
#       type <- match.arg(type)
#       if (is.null(depvar)) depvar <- names(newdata[["depvars"]])[1]
#       if (dynamic && is.null(algorithm)) {
#         stop("'algorithm' must be provided when dynamic = TRUE")
#         # could default now instead
#       }
#       if (dynamic && is.null(silent)) silent <- batch
#       if (is.null(batch_size)) {
#         batch_size <- plan_batch(
#           data = newdata, 
#           depvar = depvar, 
#           nsim = nsim,
#           nbrNodes = nbrNodes, 
#           useCluster = useCluster,
#           dynamic = dynamic, 
#           n3 = n3,
#           unit_budget = batch_unit_budget,
#           dynamic_ministep_factor = dynamic_ministep_factor,
#           memory_scale = memory_scale
#         )
#       }
#       if (dynamic) {
#         predictFun <- predictProbabilityDynamic
#         predictArgs <- list(
#           ans = object,
#           data = newdata,
#           algorithm = algorithm,
#           effects = effects,
#           type = type,
#           depvar = depvar,
#           n3 = n3,
#           useChangeContributions = useChangeContributions,
#           batch = batch,
#           silent = silent
#         )
#       } else {
#           staticContributions <- getStaticChangeContributions(
#             ans = object,
#             data = newdata,
#             effects = effects,
#             depvar = depvar,
#             returnDataFrame = TRUE
#           )
#         predictFun <- predictProbabilityStatic
#         predictArgs <- list(
#           ans = object,
#           staticContributions = staticContributions,
#           type = type
#         )
#       }
# 
#       sienaPostestimate(
#         predictFun = predictFun,
#         predictArgs = predictArgs,
#         outcome = type,
#         level = level,
#         condition = condition,
#         sum.fun = sum.fun,
#         na.rm = na.rm,
#         theta_hat = object[["theta"]], # change into coef
#         cov_theta = object[["covtheta"]], # change into cov
#         uncertainty = uncertainty,
#         uncertainty_mode = uncertainty_mode,
#         nsim = nsim,
#         uncertainty_sd = uncertainty_sd,
#         uncertainty_ci = uncertainty_ci,
#         uncertainty_probs = uncertainty_probs,
#         uncertainty_mcse = uncertainty_mcse,
#         uncertainty_mcse_batches = uncertainty_mcse_batches,
#         useCluster = useCluster,
#         nbrNodes = nbrNodes,
#         clusterType = clusterType,
#         cluster = cluster,
#         batch_dir = batch_dir,
#         prefix = prefix,
#         combine_batch = combine_batch,
#         batch_size = batch_size,
#         keep_batch = keep_batch,
#         verbose = verbose,
#         useChangeContributions = if (dynamic) useChangeContributions else NULL
#       )
# }
# 
# 
# 
# predictProbabilityStatic <- function(
#   ans, 
#   staticContributions, 
#   theta, 
#   type = "changeProb") {
#     effectNames <- getEffectNamesNoRate(ans[["requestedEffects"]])
#     thetaNoRate <- alignThetaNoRate(theta, effectNames, ans)
#     staticContributions <- widenStaticContribution(staticContributions)
#     df <- predictProbability(staticContributions, effectNames, thetaNoRate,
#         group_vars = c("period", "ego"), type = type)
#     df <- conditionalReplace(df, df[["density"]] == -1,
#         setdiff(effectNames, "density"), function(x) x * -1)
#     df
# }
# 
# predictDynamic <- function(
#     ans,
#     newdata,
#     effects,
#     algorithm,
#     type = c("changeProb", "tieProb"),
#     depvar = NULL,
#     level = "period",
#     condition = NULL,
#     sum.fun = mean,
#     na.rm = TRUE,
#     n3 = 1000,
#     useChangeContributions = FALSE,
#     uncertainty = TRUE,
#     uncertainty_mode = c("batch", "stream"),
#     useCluster = FALSE,
#     nbrNodes = 1,
#     nsim = 100,
#     uncertainty_sd = TRUE,
#     uncertainty_ci = TRUE,
#     uncertainty_probs = c(0.025, 0.5, 0.975),
#     uncertainty_mcse = FALSE,
#     uncertainty_mcse_batches = NULL,
#     clusterType = c("PSOCK", "FORK"),
#     cluster = NULL,
#     batch_dir = "temp",
#     prefix = "simBatch_b",
#     combine_batch = TRUE,
#     batch_size = NULL,
#     keep_batch = FALSE,
#     verbose = TRUE,
#     batch = TRUE,
#     silent = NULL,
#     memory_scale = NULL,
#     batch_unit_budget = 5e6,
#     dynamic_ministep_factor = 10
# ){
#   predict.sienaFit(
#     object = ans,
#     newdata = newdata,
#     type = type,
#     effects = effects,
#     depvar = depvar,
#     dynamic = TRUE,
#     algorithm = algorithm,
#     n3 = n3,
#     useChangeContributions = useChangeContributions,
#     level = level,
#     condition = condition,
#     sum.fun = sum.fun,
#     na.rm = na.rm,
#     uncertainty = uncertainty,
#     uncertainty_mode = uncertainty_mode,
#     useCluster = useCluster,
#     nbrNodes = nbrNodes,
#     nsim = nsim,
#     uncertainty_sd = uncertainty_sd,
#     uncertainty_ci = uncertainty_ci,
#     uncertainty_probs = uncertainty_probs,
#     uncertainty_mcse = uncertainty_mcse,
#     uncertainty_mcse_batches = uncertainty_mcse_batches,
#     clusterType = clusterType,
#     cluster = cluster,
#     batch_dir = batch_dir,
#     prefix = prefix,
#     combine_batch = combine_batch,
#     batch_size = batch_size,
#     keep_batch = keep_batch,
#     verbose = verbose,
#     batch = batch,
#     silent = silent,
#     memory_scale = memory_scale,
#     batch_unit_budget = batch_unit_budget,
#     dynamic_ministep_factor = dynamic_ministep_factor
#   )
# }
# 
# predictProbabilityDynamic <- function(ans, data, theta, algorithm, effects,
#     type = "changeProb", depvar, n3 = NULL, useChangeContributions = FALSE,   
#     batch = TRUE, silent = TRUE) {
#     effectNames <- getEffectNamesNoRate(effects)
#     thetaNoRate <- alignThetaNoRate(theta, effectNames, ans)
# 
#     df <- getDynamicChangeContributions(
#         ans = ans, 
#         theta = theta, 
#         data = data, 
#         algorithm = algorithm,
#         effects = effects, 
#         depvar = depvar, 
#         n3 = n3, 
#         useChangeContributions = useChangeContributions, 
#         returnDataFrame = TRUE,
#         batch = batch,
#         silent = silent
#     )
#     df <- widenDynamicContribution(df)
#     df <- predictProbability(df, effectNames, thetaNoRate,
#         group_vars = c("chain", "period", "ministep"), type = type)
#     conditionalReplace(df, df[["density"]] == -1,
#         setdiff(effectNames, "density"), function(x) x * -1)
# }

# ===========================================================================
# === V1 BACKUP: predictProbability, addUtilityColumn, addProbabilityColumn ===
# === (V1-only per-row helpers; calculateUtility is now a live function)    ===
# ===========================================================================

# 
# #' Shared core: compute utility and probability columns
# #' Used by predictProbabilityStatic, predictProbabilityDynamic,
# #' predictFirstDiff, and predictSecondDiff.
# predictProbability <- function(df, effectNames, thetaNoRate, group_vars, type) {
#     df <- addUtilityColumn(df, effectNames, thetaNoRate)
#     df <- addProbabilityColumn(df, group_vars = group_vars, type = type)
#     df
# }
# 
# addUtilityColumn <- function(df, effectNames, theta) {
#   available_effects <- effectNames[effectNames %in% names(df)]
#   if (!length(available_effects)) {
#     stop("No requested effect columns found in data for utility calculation")
#   }
# 
#   # checked already?
#   if (!is.null(names(theta))) {
#     theta_use <- theta[available_effects]
#   } else {
#     theta_use <- theta[seq_along(available_effects)]
#   }
# 
#   if (requireNamespace("data.table", quietly = TRUE) && 
#       data.table::is.data.table(df)) {
#     df[, ("changeUtil") := calculateUtility(
#       as.matrix(.SD), 
#       theta_use), 
#       .SDcols = available_effects]
#   } else {
#     df[["changeUtil"]] <- calculateUtility(
#       as.matrix(df[, available_effects, drop = FALSE]), 
#       theta_use)
#     return(df)
#   }
# }
# 
# # for density == 0: can we filter density 0 out?
# addProbabilityColumn <- function(
#     df,
#     group_vars,
#     type = "changeProb"
# ) {
#   stopifnot(is.data.frame(df))
#   # Zero-copy list of grouping column pointers -> single C++ call
#   group_cols <- lapply(group_vars, function(v) as.integer(df[[v]]))
#   df[["changeProb"]] <- as.vector(softmax_rcpp_grouped_lst(df[["changeUtil"]], group_cols))
#   if (type == "tieProb") {
#     df[["tieProb"]] <- df[["changeProb"]]
#     if ("density" %in% names(df)) {
#       idx_neg1 <- which(df[["density"]] == -1)
#       df[["tieProb"]][idx_neg1] <- 1 - df[["tieProb"]][idx_neg1]
#       df[["tieProb"]][df[["density"]] == 0] <- NA
#     }
#   }
#   df
# }
