sienaPostestimate <- function(
    predictFun,
    predictArgs,
    outcomeName,
    level = "period",
    condition = NULL,
    sum_fun = mean,
    na.rm = TRUE,
    thetaHat,
    covTheta,
    uncertainty = TRUE,
    uncertaintyMode = c("batch", "stream"),
    nsim = 1000,
    uncertaintySd = TRUE,
    uncertaintyCi = TRUE,
    uncertaintyProbs = c(0.025, 0.5, 0.975),
    uncertaintyMcse = FALSE,
    uncertaintymcseBatches = NULL,
    useCluster = FALSE,
    nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batchDir = "temp",
    prefix = "simBatch_b",
    combineBatch = TRUE,
    batchSize = NULL,
    keepBatch = FALSE,
    verbose = TRUE,
    useChangeContributions = NULL,
    gcEachBatch = TRUE,
    gcEachSim = FALSE,
    metadata = NULL,
    egoNormalize = TRUE,
    returnDecisionDetails = FALSE
) {
  uncertaintyMode <- match.arg(uncertaintyMode)
    if (length(outcomeName) != 1L) {
        stop("'outcomeName' must be a single column name.")
    }

    if(!is.null(condition)) condition <- resolveCondition(condition)

    estimator <- makeEstimator(
        predictFun = predictFun,
        predictArgs = predictArgs,
        outcomeName = outcomeName,
        level = level,
        condition = condition,
        sum_fun = sum_fun,
        na.rm = na.rm,
        thetaNames = names(thetaHat),
        egoNormalize = egoNormalize
    )

    if (!is.null(useChangeContributions)) {
        expect <- estimator(thetaHat, 
          useChangeContributions = useChangeContributions,
          returnDecisionDetails = returnDecisionDetails)
    } else {
        expect <- estimator(thetaHat,
          returnDecisionDetails = returnDecisionDetails)
    }

    # Unpack decision details from point estimate if requested
    decisionDetails <- NULL
    if (returnDecisionDetails && is.list(expect) &&
        !is.null(expect$decisionDetails)) {
      decisionDetails <- expect$decisionDetails
      expect <- expect$summary
    }

    if (!uncertainty) {
      result <- stampPostestimate(expect, metadata)
      if (!is.null(decisionDetails)) attr(result, "decisionDetails") <- decisionDetails
      return(result)
    }


    # Resolve condition to actual column name now that we have the output shape.
    # agg() resolves lazily per call, but drawSimStream needs it up-front for
    # group_vars.  Doing it here keeps all downstream paths consistent.
    if (!is.null(condition) && !is.null(colnames(expect))) {
        condition <- tryCatch(
            resolveEffectName(condition, colnames(expect)),
            error = function(e) condition
        )
    }
    
    uncertainty_summary_fun <- makeUncertaintySummarizer(
      return_sd = uncertaintySd,
      return_ci = uncertaintyCi,
      probs = uncertaintyProbs,
      return_mcse = uncertaintyMcse,
      mcseBatches = uncertaintymcseBatches
    )

    if (identical(uncertaintyMode, "stream")) {
      uncert <- drawSimStream(
        estimator = estimator,
        outcomeName = outcomeName,
        level = level,
        condition = condition,
        uncertainty_summary_fun = uncertainty_summary_fun,
        thetaHat = thetaHat,
        covTheta = covTheta,
        nbrNodes = if (useCluster) nbrNodes else 1L,
        nsim = nsim,
        clusterType = clusterType,
        cluster = cluster,
        batchSize = batchSize,
        verbose = verbose,
        gcEachBatch = gcEachBatch,
        gcEachSim = gcEachSim
      )
    } else {
      uncert <- drawSim(
        estimator = estimator,
        thetaHat = thetaHat,
        covTheta = covTheta,
        nbrNodes = if (useCluster) nbrNodes else 1L,
        nsim = nsim,
        clusterType = clusterType,
        cluster = cluster,
        batchDir = batchDir,
        prefix = prefix,
        combineBatch = combineBatch,
        batchSize = batchSize,
        keepBatch = keepBatch,
        verbose = verbose,
        gcEachBatch = gcEachBatch,
        gcEachSim = gcEachSim
      )
      uncert <- agg(
        outcomeName, 
        uncert, 
        level = level, 
        condition = condition, 
        sum_fun = uncertainty_summary_fun
      )
    }
    result <- mergeEstimates(expect, uncert, level = level, condition = condition)
    result <- stampPostestimate(result, metadata)
    if (!is.null(decisionDetails)) attr(result, "decisionDetails") <- decisionDetails
    result
}

makeUncertaintySummarizer <- function(
    return_sd = TRUE,
    return_ci = TRUE,
    probs = c(0.025, 0.5, 0.975),
    return_mcse = FALSE,
    mcseBatches = NULL
) {
  if (length(probs) != 3L) {
    stop("'uncertaintyProbs' must be length 3: lower, median, upper.")
  }
  probs <- as.numeric(probs)
  if (any(!is.finite(probs)) || any(probs <= 0) || any(probs >= 1)) {
    stop("'uncertaintyProbs' must be strictly between 0 and 1.")
  }

  function(x, na.rm = TRUE) {
    if (na.rm) x <- x[!is.na(x)]
    n <- length(x)

    out <- list(
      Mean = if (n) mean(x) else NA_real_
    )

    if (isTRUE(return_sd)) {
      out[["SE"]] <- if (n >= 2L) stats::sd(x) else NA_real_
    }

    if (isTRUE(return_ci)) {
      if (n) {
        qu <- unname(stats::quantile(x, probs = probs, na.rm = FALSE))
        out[["Median"]] <- qu[2]
        out[["q_025"]] <- qu[1]
        out[["q_975"]] <- qu[3]
      } else {
        out[["Median"]] <- NA_real_
        out[["q_025"]] <- NA_real_
        out[["q_975"]] <- NA_real_
      }
    }

    out[["cases"]] <- n

    if (isTRUE(return_mcse)) {
      out <- c(out, computeMcse(
        x = x,
        return_sd = return_sd,
        return_ci = return_ci,
        probs = probs,
        batches = mcseBatches
      ))
    }

    out
  }
}

computeMcse <- function(x,
                                   return_sd = TRUE,
                                   return_ci = TRUE,
                                   probs = c(0.025, 0.5, 0.975),
                                   batches = NULL) {
  n <- length(x)
  if (n < 2L) {
    out <- list(mcse_Mean = NA_real_)
    if (isTRUE(return_sd)) out[["mcse_SE"]] <- NA_real_
    if (isTRUE(return_ci)) {
      out[["mcse_Median"]] <- NA_real_
      out[["mcse_q_025"]] <- NA_real_
      out[["mcse_q_975"]] <- NA_real_
    }
    return(out)
  }

  # Batch MCSE: split into B batches, compute the statistic per batch,
  # then MCSE(stat) = sd(stat_b) / sqrt(B).
  if (is.null(batches)) {
    B <- max(2L, min(as.integer(floor(sqrt(n))), n))
  } else {
    B <- as.integer(batches)
    if (!is.finite(B) || B < 2L) {
      stop("'uncertaintymcseBatches' must be NULL or an integer >= 2.")
    }
    B <- min(B, n)
  }

  m <- as.integer(floor(n / B))
  if (m < 1L) {
    # Should be unreachable given B <= n, but keep defensive.
    B <- n
    m <- 1L
  }
  n_use <- as.integer(B * m)
  x_use <- x[seq_len(n_use)]

  batch_stats_mean <- numeric(B)
  batch_stats_sd <- if (isTRUE(return_sd)) numeric(B) else NULL
  batch_stats_q <- if (isTRUE(return_ci)) matrix(NA_real_, nrow = B, ncol = 3L) else NULL

  for (b in seq_len(B)) {
    idx <- ((b - 1L) * m + 1L):(b * m)
    xb <- x_use[idx]
    batch_stats_mean[b] <- mean(xb)
    if (isTRUE(return_sd)) {
      batch_stats_sd[b] <- if (m >= 2L) stats::sd(xb) else NA_real_
    }
    if (isTRUE(return_ci)) {
      batch_stats_q[b, ] <- unname(stats::quantile(xb, probs = probs, na.rm = FALSE))
    }
  }

  out <- list(
    mcse_Mean = stats::sd(batch_stats_mean) / sqrt(B)
  )

  if (isTRUE(return_sd)) {
    out[["mcse_SE"]] <- stats::sd(batch_stats_sd) / sqrt(B)
  }

  if (isTRUE(return_ci)) {
    out[["mcse_q_025"]] <- stats::sd(batch_stats_q[, 1]) / sqrt(B)
    out[["mcse_Median"]] <- stats::sd(batch_stats_q[, 2]) / sqrt(B)
    out[["mcse_q_975"]] <- stats::sd(batch_stats_q[, 3]) / sqrt(B)
  }

  out
}

# thetaNames should not be necessary anymore?
makeEstimator <- function(predictFun, predictArgs, outcomeName,
    level = "period", condition = NULL, sum_fun = mean, na.rm = TRUE,
    thetaNames = NULL, egoNormalize = TRUE) {
    function(theta, useChangeContributions = FALSE,
             returnDecisionDetails = FALSE) {
        if (!is.null(thetaNames) && is.null(names(theta)))
            names(theta) <- thetaNames
        if(!is.null(predictArgs[["useChangeContributions"]])){
          predictArgs[["useChangeContributions"]] <- useChangeContributions
        }
        callArgs <- c(list(theta = theta), predictArgs)
        unit_pred <- do.call(predictFun, callArgs)
        main_result <- agg(
            outcomeName = outcomeName,
            data = unit_pred,
            level = level,
            condition = condition,
            sum_fun = sum_fun,
            na.rm = na.rm,
            egoNormalize = egoNormalize
        )

        # Aggregate mass contrast columns alongside the main outcome.
        # Mass contrasts are ego-level quantities (already summed across the
        # density dimension internally)
        massCols <- intersect(c("massCreation", "massDissolution"), names(unit_pred))
        for (mc in massCols) {
            mc_agg <- agg(
                outcomeName = mc,
                data = unit_pred,
                level = level,
                condition = NULL,
                sum_fun = sum_fun,
                na.rm = na.rm,
                egoNormalize = egoNormalize
            )
            main_result[[mc]] <- mc_agg[[mc]]
        }
        if (isTRUE(predictArgs$details)) {
          details <- unit_pred
          return(list(summary = main_result, details = details))
        }
        if (returnDecisionDetails) {
          return(list(summary = main_result,
                      decisionDetails = unit_pred))
        }
        return(main_result)
    }
}

drawSim <- function(
    estimator,
    thetaHat,
    covTheta,
    nbrNodes = 1,
    nsim = 1000,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batchDir = "temp",
    prefix = "simBatch_b",
    batchSize = NULL,
    combineBatch = TRUE,
    keepBatch = FALSE,
    verbose = TRUE,
    gcEachBatch = TRUE,
    gcEachSim = FALSE
) {

    cli <- openCluster(nbrNodes, clusterType, cluster,
      export_vars  = c("estimator", "thetaHat", "covTheta"),
      export_envir = environment(),
      verbose      = verbose
    )
    if (!dir.exists(batchDir)) dir.create(batchDir, recursive = TRUE)
    nbatches <- ceiling(nsim / batchSize) 
    if (verbose) {
      message(
        "Starting simulations: nsim=", nsim,
        ", batchSize=", batchSize,
        ", batches=", nbatches,
        ". (Last batch may be smaller.)"
      )
    }
    
    for (b in seq_len(nbatches)) {
      first_sim <- 1 + (b - 1) * batchSize
      last_sim  <- min(b * batchSize, nsim)
      batch     <- first_sim:last_sim
      if (verbose) {
        message(
          sprintf(
            "Batch %d/%d: simulations %d-%d (%d sims)",
            b, nbatches, first_sim, last_sim, length(batch)
          )
        )
      }

      if (cli$usePSOCK) {
        batch_result <- parallel::parLapply(cli$cl, batch, function(i) {
          theta_sim <- MASS::mvrnorm(1, mu = thetaHat, Sigma = covTheta)
          sim <- do.call(estimator, list(theta = theta_sim))
          sim[, "sim"] <- i
          sim
        })
      } else if (cli$useFORK) {
        batch_result <- parallel::mclapply(batch, function(i) {
          # if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1) # data.table removed
          theta_sim <- MASS::mvrnorm(1, mu = thetaHat, Sigma = covTheta)
          sim <- do.call(estimator, list(theta = theta_sim))
          sim[, "sim"] <- i
          sim
        }, mc.cores = nbrNodes)
      } else {
        batch_result <- vector("list", length(batch))
        for (k in seq_along(batch)) {
          i <- batch[k]
          theta_sim <- MASS::mvrnorm(1, mu = thetaHat, Sigma = covTheta)
          sim <- do.call(estimator, list(theta = theta_sim))
          sim[, "sim"] <- i
          batch_result[[k]] <- sim
          if (gcEachSim) gc(verbose = FALSE)
        }
      }

      file_name <- file.path(batchDir, sprintf("%s%03d.rds", prefix, b))
      saveRDS(batch_result, file = file_name)
      rm(batch_result)
      if (gcEachBatch) gc(verbose = FALSE)
      if (verbose) {
        done <- last_sim
        pct <- round(100 * done / nsim)
        message(sprintf("Completed batch %d/%d (%d/%d sims, %d%%)", b, nbatches, done, nsim, pct))
      }
    }
    
    closteCluster(cli, verbose)
    
    if (verbose) message("All batches complete. Now combining and cleaning up...")
    
    # if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads() # data.table removed
    if(combineBatch){
      # Combine all batches
      batch_files <- list.files(batchDir, 
        pattern = sprintf("^%s\\d{3}\\.rds$", prefix), full.names = TRUE)
      batch_files <- sort(batch_files)
      combined_results <- vector("list", nsim)
      
      for (b in seq_along(batch_files)) {
        file_name <- batch_files[b]
        results <- readRDS(file_name)
        first <- 1 + (b - 1) * batchSize
        last  <- min(b * batchSize, nsim)
        combined_results[first:last] <- results
        if(!keepBatch){
          file.remove(file_name)
          if (verbose) message(sprintf("Deleted batch %d (%d-%d)", b, first, last))
        }
      }
      if (verbose) message("Returning combined results.")
      # data.table removed — use do.call(rbind, ...)
      # if (requireNamespace("data.table", quietly = TRUE)) {
      #   return(data.table::rbindlist(combined_results, use.names = TRUE, fill = TRUE))
      # } else {
        out <- do.call(rbind, combined_results)
        rownames(out) <- NULL
        return(out)
      # }
    }
}

drawSimStream <- function(
    estimator,
    outcomeName,
    level,
    condition,
    uncertainty_summary_fun,
    thetaHat,
    covTheta,
    nbrNodes = 1,
    nsim = 1000,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batchSize = NULL,
    verbose = TRUE,
    gcEachBatch = TRUE,
    gcEachSim = FALSE
) {

    group_vars <- getGroupVars(level = level, condition = condition)


    cli <- openCluster(nbrNodes, clusterType, cluster,
      export_vars  = c("estimator", "thetaHat", "covTheta"),
      export_envir = environment(),
      verbose      = verbose
    )
    nbatches <- ceiling(nsim / batchSize)
    if (verbose) {
      message(
        "Starting streaming simulations: nsim=", nsim,
        ", batchSize=", batchSize,
        ", batches=", nbatches,
        ". (Last batch may be smaller.)"
      )
    }

    stream_state <- new.env(parent = emptyenv(), hash = TRUE)
    group_state <- new.env(parent = emptyenv(), hash = TRUE)

    for (b in seq_len(nbatches)) {
      first_sim <- 1 + (b - 1) * batchSize
      last_sim  <- min(b * batchSize, nsim)
      batch     <- first_sim:last_sim

      if (verbose) {
        message(
          sprintf(
            "Batch %d/%d: simulations %d-%d (%d sims)",
            b, nbatches, first_sim, last_sim, length(batch)
          )
        )
      }

      if (cli$usePSOCK) {
        batch_result <- parallel::parLapply(cli$cl, batch, function(i) {
          theta_sim <- MASS::mvrnorm(1, mu = thetaHat, Sigma = covTheta)
          do.call(estimator, list(theta = theta_sim))
        })
      } else if (cli$useFORK) {
        batch_result <- parallel::mclapply(batch, function(i) {
          # if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1) # data.table removed
          theta_sim <- MASS::mvrnorm(1, mu = thetaHat, Sigma = covTheta)
          do.call(estimator, list(theta = theta_sim))
        }, mc.cores = nbrNodes)
      } else {
        batch_result <- vector("list", length(batch))
        for (k in seq_along(batch)) {
          theta_sim <- MASS::mvrnorm(1, mu = thetaHat, Sigma = covTheta)
          batch_result[[k]] <- do.call(estimator, list(theta = theta_sim))
          if (gcEachSim) gc(verbose = FALSE)
        }
      }

      for (sim_df in batch_result) {
        updateStream(
          stream_state = stream_state,
          group_state = group_state,
          sim_df = sim_df,
          outcomeName = outcomeName,
          group_vars = group_vars
        )
      }

      rm(batch_result)
      if (gcEachBatch) gc(verbose = FALSE)
      if (verbose) {
        done <- last_sim
        pct <- round(100 * done / nsim)
        message(sprintf("Completed batch %d/%d (%d/%d sims, %d%%)", b, nbatches, done, nsim, pct))
      }
    }

    closteCluster(cli, verbose)

    finalizeStream(
      stream_state = stream_state,
      group_state = group_state,
      group_vars = group_vars,
      uncertainty_summary_fun = uncertainty_summary_fun
    )
}

openCluster <- function(nbrNodes, clusterType, cluster, export_vars, export_envir, verbose) {
  clusterType <- if (nbrNodes > 1) match.arg(clusterType, c("PSOCK", "FORK")) else "FORK"
  usePSOCK <- (clusterType == "PSOCK" && nbrNodes > 1)
  useFORK  <- (clusterType == "FORK"  && nbrNodes > 1 && .Platform$OS.type != "windows")

  if (clusterType == "FORK" && capabilities("X11") && dev.cur() == 1)
    warning("An X11 graphics device is open. This can cause parallel processing to hang. ...")

  cluster_created <- FALSE
  cl <- cluster
  if (usePSOCK && is.null(cl)) {
    cl <- parallel::makeCluster(nbrNodes, type = "PSOCK")
    cluster_created <- TRUE
    if (verbose) message("Created new parallel cluster with ", nbrNodes, " cores.")
  }
  if (usePSOCK) {
    parallel::clusterExport(cl, export_vars, envir = export_envir)
    # data.table removed — no need for setDTthreads
    # parallel::clusterEvalQ(cl, {
    #   if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
    # })
  }
  list(cl = cl, cluster_created = cluster_created, usePSOCK = usePSOCK, useFORK = useFORK)
}

closteCluster <- function(cli, verbose) {
  if (cli$cluster_created) {
    parallel::stopCluster(cli$cl)
    if (verbose) message("Stopped internal cluster.")
  }
}

updateStream <- function(stream_state,
                                group_state,
                                sim_df,
                                outcomeName,
                                group_vars) {
  if (length(outcomeName) != 1L) stop("'outcomeName' must be a single column name.")

  # data.table removed — sim_df is always a data.frame now

  keys <- makeGroupKey(sim_df, group_vars)
  vals <- sim_df[[outcomeName]]

  for (i in seq_along(vals)) {
    key <- keys[[i]]
    old <- if (exists(key, envir = stream_state, inherits = FALSE)) {
      get(key, envir = stream_state, inherits = FALSE)
    } else {
      numeric(0)
    }
    assign(key, c(old, vals[[i]]), envir = stream_state)

    if (!exists(key, envir = group_state, inherits = FALSE) && length(group_vars)) {
      assign(key, sim_df[i, group_vars, drop = FALSE], envir = group_state)
    }
  }
}

finalizeStream <- function(stream_state,
                                  group_state,
                                  group_vars,
                                  uncertainty_summary_fun) {
  # should be able to handle data.table as well
  keys <- ls(stream_state, all.names = TRUE)
  if (!length(keys)) {
    return(data.frame())
  }

  rows <- lapply(keys, function(key) {
    draws <- get(key, envir = stream_state, inherits = FALSE)
    summary_vals <- as.data.frame(as.list(uncertainty_summary_fun(draws, na.rm = TRUE)), stringsAsFactors = FALSE)
    if (length(group_vars)) {
      grp <- get(key, envir = group_state, inherits = FALSE)
      cbind(grp, summary_vals)
    } else {
      summary_vals
    }
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (length(group_vars) > 0L && nrow(out) > 0L) {
    enc <- encodeGroupKeys(out, group_vars)
    ord <- do.call(order, as.data.frame(enc$G))
    out <- out[ord, , drop = FALSE]
    rownames(out) <- NULL
  }
  out
}

getGroupVars <- function(level = "none", condition = NULL) {
  levels <- list(
    none = character(0),
    period = "period",
    ego = c("period", "ego"),
    egoChoice = c("period", "ego", "choice"),
    chain = c("period", "chain"),
    chainEgo = c("period", "chain", "ego"),
    ministep = c("period", "chain", "ego","ministep"),
    ministepChoice = c("period", "chain", "ego","ministep", "choice")
  )
  c(levels[[level]], condition)
}

makeGroupKey <- function(df, group_vars) {
  # data.table removed — df is always a data.frame now
  if (!length(group_vars)) {
    rep("__all__", nrow(df))
  } else {
    do.call(paste, c(df[group_vars], sep = "\r"))
  }
}

getEffectNamesNoRate <- function(effects, depvar) {
    include <- effects[["include"]]
    includedEffects <- effects[include, ]
    keep <- includedEffects[["type"]] != "rate" & includedEffects[["name"]] == depvar
    inc <- includedEffects[keep, ]
    paste(inc[["name"]], inc[["shortName"]], inc[["type"]], sep = "_")
}



# Resolve the perturbation type for a named effect.
#
# Used by marginalEffects (and potentially predict) to
# decide whether a counterfactual shift uses the one-alternative
# ("alter") or the ego-wide ("ego") mlogit update.
#
# Auto-detection relies on the interactionType metadata stored
# in the wide struct's effectInteractionTypes vector:
#   "ego"    -> perturbType "ego"
#   "dyadic" -> perturbType "alter"
#   anything else (empty, "OK", structural effects) ->
#             perturbType "alter" (safe default for network effects).
#
# effectName: Character: (resolved) composite effect name.
# effectInteractionTypes: Named character vector; may be NULL.
# override: Optional user override, one of "alter", "ego", or NULL (auto-detect).
# Returns character scalar: "alter" or "ego".
resolvePerturbType <- function(effectName,
                               effectInteractionTypes = NULL,
                               override = NULL) {
    if (!is.null(override)) {
        override <- match.arg(override, c("alter", "ego"))
        return(override)
    }
    if (is.null(effectInteractionTypes) ||
        !(effectName %in% names(effectInteractionTypes))) {
        return("alter")
    }
    iType <- effectInteractionTypes[[effectName]]
    if (identical(iType, "ego")) "ego" else "alter"
}





# Convert perturbType string to the integer code expected by
# the mlogit_update Rcpp function.
# perturbType: "alter" or "ego".
# Returns integer: 0L for "alter", 1L for "ego".
perturbTypeToInt <- function(perturbType) {
    switch(perturbType, alter = 0L, ego = 1L,
           stop("Invalid perturbType: ", perturbType))
}

# Multinomial-logit probability update with string-based perturbation type.
#
# Thin wrapper around the Rcpp mlogit_update that accepts
# perturbType as "alter" or "ego" (instead of 0L/1L)
# and returns a plain numeric vector.
#
# p:           Numeric vector of baseline probabilities.
# delta_u:     Numeric vector of utility shifts (same length as p).
# group_id:    Integer vector of group identifiers.
# perturbType: Character: "alter" (one-alternative update)
#              or "ego" (ego-wide renormalization).
# Returns numeric vector of updated probabilities.
mlogit_update_r <- function(p, delta_u, group_id, perturbType) {
    if (is.null(group_id)) group_id <- integer(length(p))
    as.vector(mlogit_update(p, delta_u, group_id,
                            perturbTypeToInt(perturbType)))
}

# Compute ego-level probability-mass contrasts.
#
# For ego-wide perturbations, the natural actor-level QOIs are the
# creation and dissolution probability-mass contrasts
# Delta P_i+ = sum_{j in C_i} (p'_ij - p_ij) and
# Delta P_i- = sum_{j in D_i} (p'_ij - p_ij),
# where the sums run over the creation and dissolution risk sets
# respectively.
#
# firstDiff: Numeric vector of dyad-level first differences
#            (on whichever scale: changeProb or tieProb).
# density:   Integer vector: 1 = creation, -1 = dissolution.
#            Rows with density = 0 should already be removed.
# ego:       Integer vector of ego identifiers.
# period:    Integer vector of period identifiers.
# group:     Integer/character vector of group identifiers.
# type:      Character: "changeProb" or "tieProb".
#            Needed to recover change-probability-scale diffs for tieProb mode.
# Returns a data.frame with columns massCreation and
#   massDissolution, one value per row (broadcast from ego-level aggregate).
computeMassContrasts <- function(firstDiff, density, ego, period, group,
                                 type = "changeProb") {
    n <- length(firstDiff)
    massCreation <- rep(NA_real_, n)
    massDissolution <- rep(NA_real_, n)

    # Convert firstDiff to change-probability scale if needed.
    # For tieProb, dissolution rows have firstDiff = -(changeProbDiff),
    # so we flip them back.
    changeProbDiff <- firstDiff
    if (type == "tieProb") {
        diss <- density == -1L
        changeProbDiff[diss] <- -changeProbDiff[diss]
    }

    # Build a grouping key for ego × period × group
    grpKey <- paste(group, period, ego, sep = "\r")
    ukeys <- unique(grpKey)

    for (k in ukeys) {
        idx <- which(grpKey == k)
        cre <- idx[density[idx] ==  1L]
        dis <- idx[density[idx] == -1L]
        mc <- if (length(cre)) sum(changeProbDiff[cre], na.rm = TRUE) else 0
        md <- if (length(dis)) sum(changeProbDiff[dis], na.rm = TRUE) else 0
        massCreation[idx] <- mc
        massDissolution[idx] <- md
    }
    data.frame(massCreation = massCreation, massDissolution = massDissolution)
}

conditionalReplace <- function(df, row_ids, cols, fun) {
  # data.table path removed (this function is not called in current code)
  # if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
  #   df[eval(substitute(row_ids)), (cols) := lapply(.SD, fun), .SDcols = cols]
  # } else {
    df[row_ids, cols] <- fun(df[row_ids, cols])
  # }
  df
}

# Detect columns that identify a single ego-decision in the data.
# Dynamic data has chain/ministep; static data does not.
detectEgoUnit <- function(data) {
  cols <- names(data)
  if ("chain" %in% cols) {
    intersect(c("chain", "group", "period", "ego", "ministep"), cols)
  } else {
    intersect(c("group", "period", "ego"), cols)
  }
}

# Encode group columns as a contiguous integer matrix for grouped_agg_cpp.
# Non-integer/non-numeric columns are factor-encoded.
# Returns list(G = integer matrix, decode = list of level vectors per column).
encodeGroupKeys <- function(data, group_vars) {
  n <- nrow(data)
  ncols <- length(group_vars)
  G <- matrix(0L, nrow = n, ncol = ncols)
  decode <- vector("list", ncols)
  names(decode) <- group_vars
  for (j in seq_len(ncols)) {
    col <- data[[group_vars[j]]]
    if (is.integer(col)) {
      G[, j] <- col
      decode[j] <- list(NULL)  # identity mapping (note: [[j]]<-NULL would delete)
    } else if (is.numeric(col)) {
      icol <- as.integer(col)

      if (identical(as.numeric(icol), as.numeric(col))) {
        G[, j] <- icol
        decode[j] <- list(NULL)
      } else {
        f <- factor(col)
        G[, j] <- as.integer(f)
        decode[[j]] <- levels(f)
      }
    } else {
      vals_num <- suppressWarnings(as.numeric(as.character(col)))
      f <- if (all(!is.na(vals_num))) {
        factor(col, levels = as.character(sort(unique(vals_num))))
      } else {
        factor(col)
      }
      G[, j] <- as.integer(f)
      decode[[j]] <- levels(f)
    }
  }
  storage.mode(G) <- "integer"
  list(G = G, decode = decode)
}

# Decode integer group key matrix back to original values using decode info.
decodeGroupKeys <- function(key_mat, group_vars, decode) {
  out <- data.frame(key_mat)
  names(out) <- group_vars
  for (j in seq_along(group_vars)) {
    if (!is.null(decode[[j]])) {
      out[[group_vars[j]]] <- decode[[j]][out[[group_vars[j]]]]
    }
  }
  out
}

# Works but could be nicer - push more to rcpp?
# Pre-aggregate outcomeName within each ego-unit: within-ego mean.
# Returns one row per ego-unit with pre_agg_vars + outcomeName columns.
preAggEgo <- function(data, outcomeName, group_vars, ego_id_cols, na.rm) {
  pre_agg_vars <- unique(c(group_vars, ego_id_cols))

  # data.table path removed — using Rcpp grouped_agg_cpp for speed
  # if (requireNamespace("data.table", quietly = TRUE) &&
  #     data.table::is.data.table(data)) {
  #   return(data[, setNames(list(mean(get(outcomeName), na.rm = na.rm)),
  #                          outcomeName),
  #               by = pre_agg_vars])
  # }

  # Rcpp fast path: use grouped_agg_cpp for contiguous integer keys
  # grouped_agg_cpp requires contiguous groups — sort by encoded keys
  enc <- encodeGroupKeys(data, pre_agg_vars)
  ord <- do.call(order, as.data.frame(enc$G))
  G_sorted <- enc$G[ord, , drop = FALSE]
  x_sorted <- data[[outcomeName]][ord]
  res <- grouped_agg_cpp(x_sorted, G_sorted, na_rm = na.rm, do_mean = TRUE)
  out <- decodeGroupKeys(res$key, pre_agg_vars, enc$decode)
  out[[outcomeName]] <- res$value
  out
  # old base R fallback (kept for reference)
  # vals <- data[[outcomeName]]
  # ego_key <- do.call(paste, c(data[pre_agg_vars], sep = "\r"))
  # ukeys <- unique(ego_key)
  # first_idx <- match(ukeys, ego_key)
  # split_vals <- split(vals, ego_key)
  # ego_means <- vapply(split_vals[ukeys],
  #                     function(x) mean(x, na.rm = na.rm),
  #                     numeric(1))
  # out <- data[first_idx, pre_agg_vars, drop = FALSE]
  # out[[outcomeName]] <- unname(ego_means)
  # rownames(out) <- NULL
  # out
}


agg <- function(outcomeName,
                data,
                level = "none",
                condition = NULL,
                sum_fun = mean,
                na.rm = TRUE,
                egoNormalize = TRUE) {
  if (length(outcomeName) != 1L) stop("'outcomeName' must be a single column name.")
  if (!is.null(condition)) condition <- resolveEffectName(condition, colnames(data))
  group_vars <- getGroupVars(level = level, condition = condition)
  # Helper for complex output (vector/list output from sum_fun)
  expand_summary <- function(val) {
    if (length(val) == 1 && (is.null(names(val)) || names(val) == "")) {
      setNames(list(val), outcomeName)
    } else {
      as.list(val)
    }
  }

  # Ego-first pre-aggregation: compute the mean within each realized
  # ego-decision first, then average those decision-level means.
  if (egoNormalize) {
    ego_id_cols <- detectEgoUnit(data)
    extra <- setdiff(ego_id_cols, group_vars)
    if (length(extra) > 0) {
      data <- preAggEgo(data, outcomeName, group_vars, ego_id_cols, na.rm)
    }
  }

  # data.table path removed — use Rcpp fast path for mean/sum, base R for custom sum_fun
  # if (requireNamespace("data.table", quietly = TRUE)) {
  #   if (!data.table::is.data.table(data)) {
  #     data <- data.table::as.data.table(data)
  #   }
  #   result <- data[, {
  #     res <- expand_summary(sum_fun(get(outcomeName), na.rm = na.rm))
  #     res
  #   }, by = group_vars]
  #   return(result)
  # }

  # ---- Rcpp fast path for simple mean/sum on grouped integer keys ----
  is_simple_mean <- identical(sum_fun, mean) || identical(sum_fun, base::mean)
  if (is_simple_mean && length(group_vars) > 0) {
    enc <- encodeGroupKeys(data, group_vars)
    ord <- do.call(order, as.data.frame(enc$G))
    G_sorted <- enc$G[ord, , drop = FALSE]
    x_sorted <- data[[outcomeName]][ord]
    res <- grouped_agg_cpp(x_sorted, G_sorted, na_rm = na.rm, do_mean = TRUE)
    out <- decodeGroupKeys(res$key, group_vars, enc$decode)
    out[[outcomeName]] <- res$value
    return(out)
  }

  # ---- Base R path ----
  if (length(group_vars) == 0) {
    output <- expand_summary(sum_fun(data[[outcomeName]], na.rm = na.rm))
    output <- as.data.frame(output)
    return(output)
  }

  # General base R path (custom sum_fun or multi-column output)
  grouping <- interaction(data[, group_vars, drop = FALSE], drop = TRUE)
  split_data <- split(data, grouping)
  agg_list <- lapply(split_data, function(subdf) {
    vals <- subdf[[outcomeName]]
    res <- expand_summary(sum_fun(vals, na.rm = na.rm))
    group_vals <- subdf[1, group_vars, drop = FALSE]
    cbind(group_vals, as.data.frame(res, stringsAsFactors = FALSE))
  })
  out <- do.call(rbind, agg_list)
  rownames(out) <- NULL
  return(out)
}

summarizeValue <- function(x, na.rm = TRUE){
  if(na.rm){
    x <- x[! is.na(x)]
  }
  qu <- unname(quantile(x, probs = c(.025,0.5,.975)))
  mn <- mean(x)
  se <- sd(x)
  n <- length(x)
  as.list(c("Mean" = mn, "SE" = se, "Median" = qu[2], "q_025" = qu[1], "q_975" = qu[3], "cases" = n))
}

mergeEstimates <- function(df1, df2, level = "none", condition = NULL) {
  levels_list <- list(
    none = character(0),
    period = "period",
    ego = c("period", "ego"),
    egoChoice = c("period", "ego", "choice"),
    chain = c("period", "chain"),
    ministep = c("period", "chain", "ministep"),
    ministepChoice = c("period", "chain", "ministep", "choice")
  )
  if (!is.null(condition) && is.data.frame(df1) && ncol(df1) > 0L)
    condition <- resolveEffectName(condition, colnames(df1))
  group_vars <- c(levels_list[[level]], condition)
  # inefficient for data.table? Even for data.frame this could be more efficient
  # if we avoid the merge and just do a match on the group vars, but for now we
  # can keep it simple and use merge, and optimize later if needed
  if (length(group_vars) > 0) {
    merge(df1, df2, by = group_vars, all = TRUE)
  } else {
    cbind(df1, df2)
  }
}

# Ensure theta has names. If theta is already named, return it unchanged.
# If unnamed (e.g. user-supplied), warn and generate names from effects using
# the same convention as getNamesFromEffects() in robmon.r.
nameThetaFromEffects <- function(theta, effects) {
    if (!is.null(names(theta))) return(theta)
    warning("theta has no names; generating names from effects")
    setNames(theta, getNamesFromEffects(effects))
}

# Align theta to a contributions effectNames vector via name-based subsetting.
# theta must already be named (e.g. via nameThetaFromEffects).
alignThetaNoRate <- function(theta, effectNames) {
    theta[effectNames]
}

# Attach effect contribution columns from a matrix to a data.frame.
# Builds all columns at once via cbind to avoid repeated data.frame copies.
# If flip = TRUE, negates non-density columns where density == -1.
attachContribColumns <- function(out, effectNames, contrib, flip = TRUE) {
  density_j <- grep("density", effectNames, fixed = TRUE)
  if (flip && length(density_j) > 0) {
    neg1 <- contrib[, density_j[1L]] == -1L
    neg1[is.na(neg1)] <- FALSE
    if (any(neg1)) {
      flip_idx <- grep("density", effectNames, fixed = TRUE, invert = TRUE)
      for (j in flip_idx) contrib[neg1, j] <- -contrib[neg1, j]
    }
  }
  # Build data.frame from list of column vectors — avoids slow
  # as.data.frame.matrix which copies via as.vector per column.
  nc <- ncol(contrib)
  extra <- vector("list", nc)
  names(extra) <- effectNames
  for (j in seq_len(nc)) extra[[j]] <- contrib[, j]
  attr(extra, "row.names") <- .set_row_names(nrow(contrib))
  class(extra) <- "data.frame"
  cbind(out, extra)
}

# Extract group-column vectors from a Contributions struct as a named list.
# keep is an optional logical/integer subset vector.
# For static, pb$group is scalar (recycled by data.frame()).
# For dynamic, ego is included when present in the struct.
groupColsList <- function(pb, keep = NULL) {
  if (!is.null(pb$chain)) {
    if (is.null(keep)) {
      out <- list(chain = pb$chain, group = pb$group, period = pb$period,
           ego = pb$ego, ministep = pb$ministep, choice = pb$choice)
    } else {
      out <- list(chain = pb$chain[keep], group = pb$group[keep],
           period = pb$period[keep], ego = pb$ego[keep],
           ministep = pb$ministep[keep], choice = pb$choice[keep])
    }
    # Drop ego if not present (backward compat with old structs)
    if (is.null(out$ego)) out$ego <- NULL
    out
  } else {
    if (is.null(keep)) {
      list(group = pb$group, period = pb$period, ego = pb$ego,
           choice = pb$choice)
    } else {
      list(group = pb$group[keep], period = pb$period[keep],
           ego = pb$ego[keep], choice = pb$choice[keep])
    }
  }
}

resolveCondition <- function(condition) {
  ifelse(grepl("_", condition, fixed = TRUE), condition, paste0(condition, "_eval"))
}

# ---- Class stamping and S3 methods for postestimate output -------------------

# Convert postestimate result to the appropriate S3 class with metadata.
stampPostestimate <- function(result, metadata = NULL) {
  # data.table removed — result is always a data.frame now
  # if (!returnDataTable && inherits(result, "data.table")) {
  #   result <- as.data.frame(result)
  # }
  if (!is.null(metadata)) {
    # Store resolved condition column names from the output
    cond_cols <- grep("_eval$", names(result), value = TRUE)
    if (length(cond_cols) > 0) metadata$conditionCols <- cond_cols
    for (nm in names(metadata)) {
      attr(result, nm) <- metadata[[nm]]
    }
    cls <- switch(metadata$method,
      predict         = "sienaPrediction",
      marginalEffects = "sienaMarginalEffect",
      NULL
    )
    if (!is.null(cls)) {
      class(result) <- c(cls, class(result))
    }
  }
  result
}

##@print.sienaPrediction S3 print
print.sienaPrediction <- function(x, ...) {
  cat("SAOM Prediction\n")
  cat("  Type:      ", if (!is.null(attr(x, "type"))) attr(x, "type") else "unknown", "\n")
  cat("  Dep. var.: ", if (!is.null(attr(x, "depvar"))) attr(x, "depvar") else "unknown", "\n")
  cat("  Level:     ", if (!is.null(attr(x, "level"))) attr(x, "level") else "unknown", "\n")
  cat("  Dynamic:   ", if (!is.null(attr(x, "dynamic"))) attr(x, "dynamic") else FALSE, "\n")
  cat("  nsim:      ", if (!is.null(attr(x, "nsim"))) attr(x, "nsim") else NA, "\n")
  cat("\n")
  print.data.frame(x, ...)
  invisible(x)
}

##@summary.sienaPrediction S3 summary
summary.sienaPrediction <- function(object, ...) {
  cat("SAOM Prediction Summary\n")
  cat("  Type:      ", if (!is.null(attr(object, "type"))) attr(object, "type") else "unknown", "\n")
  cat("  Dep. var.: ", if (!is.null(attr(object, "depvar"))) attr(object, "depvar") else "unknown", "\n")
  cat("  Level:     ", if (!is.null(attr(object, "level"))) attr(object, "level") else "unknown", "\n")
  cat("  Dynamic:   ", if (!is.null(attr(object, "dynamic"))) attr(object, "dynamic") else FALSE, "\n")
  cat("  nsim:      ", if (!is.null(attr(object, "nsim"))) attr(object, "nsim") else NA, "\n")
  if (!is.null(attr(object, "condition")))
    cat("  Condition: ", paste(attr(object, "condition"), collapse = ", "), "\n")
  cat("\n")
  cat("  Rows:      ", nrow(object), "\n")
  cat("  Columns:   ", paste(names(object), collapse = ", "), "\n\n")
  summary.data.frame(object, ...)
}

##@print.sienaMarginalEffect S3 print
print.sienaMarginalEffect <- function(x, ...) {
  cat("SAOM Marginal Effect\n")
  cat("  Effect:    ", if (!is.null(attr(x, "effectName1"))) attr(x, "effectName1") else "unknown", "\n")
  if (isTRUE(attr(x, "second"))) {
    cat("  Effect 2:  ", if (!is.null(attr(x, "effectName2"))) attr(x, "effectName2") else "unknown", "\n")
  }
  me <- if (!is.null(attr(x, "mainEffect"))) attr(x, "mainEffect") else "riskDifference"
  cat("  Scale:     ", me, "\n")
  cat("  Type:      ", if (!is.null(attr(x, "type"))) attr(x, "type") else "unknown", "\n")
  cat("  Dep. var.: ", if (!is.null(attr(x, "depvar"))) attr(x, "depvar") else "unknown", "\n")
  cat("  Level:     ", if (!is.null(attr(x, "level"))) attr(x, "level") else "unknown", "\n")
  cat("  Dynamic:   ", if (!is.null(attr(x, "dynamic"))) attr(x, "dynamic") else FALSE, "\n")
  cat("  nsim:      ", if (!is.null(attr(x, "nsim"))) attr(x, "nsim") else NA, "\n")
  cat("\n")
  print.data.frame(x, ...)
  invisible(x)
}

##@summary.sienaMarginalEffect S3 summary
summary.sienaMarginalEffect <- function(object, ...) {
  cat("SAOM Marginal Effect Summary\n")
  cat("  Effect:    ", if (!is.null(attr(object, "effectName1"))) attr(object, "effectName1") else "unknown", "\n")
  if (isTRUE(attr(object, "second"))) {
    cat("  Effect 2:  ", if (!is.null(attr(object, "effectName2"))) attr(object, "effectName2") else "unknown", "\n")
  }
  me <- if (!is.null(attr(object, "mainEffect"))) attr(object, "mainEffect") else "riskDifference"
  cat("  Scale:     ", me, "\n")
  cat("  Type:      ", if (!is.null(attr(object, "type"))) attr(object, "type") else "unknown", "\n")
  cat("  Dep. var.: ", if (!is.null(attr(object, "depvar"))) attr(object, "depvar") else "unknown", "\n")
  cat("  Level:     ", if (!is.null(attr(object, "level"))) attr(object, "level") else "unknown", "\n")
  cat("  Dynamic:   ", if (!is.null(attr(object, "dynamic"))) attr(object, "dynamic") else FALSE, "\n")
  cat("  nsim:      ", if (!is.null(attr(object, "nsim"))) attr(object, "nsim") else NA, "\n")
  if (!is.null(attr(object, "condition")))
    cat("  Condition: ", paste(attr(object, "condition"), collapse = ", "), "\n")
  cat("\n")
  cat("  Rows:      ", nrow(object), "\n")
  cat("  Columns:   ", paste(names(object), collapse = ", "), "\n\n")
  summary.data.frame(object, ...)
}