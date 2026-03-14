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
    gcEachSim = FALSE
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
        na.rm = na.rm
    )

    if (!is.null(useChangeContributions)) {
        expect <- estimator(thetaHat, 
          useChangeContributions = useChangeContributions)
    } else {
        expect <- estimator(thetaHat)
    }

    if (!uncertainty) return(expect)

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
    mergeEstimates(expect, uncert, level = level, condition = condition)
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

makeEstimator <- function(predictFun, predictArgs, outcomeName,
    level = "period", condition = NULL, sum_fun = mean, na.rm = TRUE) {
    function(theta, useChangeContributions = FALSE) {
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
            na.rm = na.rm
        )
        if (isTRUE(predictArgs$details)) {
          details <- unit_pred
          return(list(summary = main_result, details = details))
        } else {
          return(main_result)
        }
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
          if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
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
    
    if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads()
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
      if (requireNamespace("data.table", quietly = TRUE)) {
        return(data.table::rbindlist(combined_results, use.names = TRUE, fill = TRUE))
      } else {
        out <- do.call(rbind, combined_results)
        rownames(out) <- NULL
        return(out)
      }
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
          if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
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
    parallel::clusterEvalQ(cl, {
      if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
    })
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

#we should not convert to data.frame, adjust code to handle data.table if needed, but for now we convert to data.frame to avoid issues with list columns in data.table
  if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(sim_df)) {
    sim_df <- as.data.frame(sim_df)
  }

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
  out
}

getGroupVars <- function(level = "none", condition = NULL) {
  levels <- list(
    none = character(0),
    period = "period",
    ego = c("period", "ego"),
    egoChoice = c("period", "ego", "choice"),
    chain = c("period", "chain"),
    ministep = c("period", "chain", "ministep"),
    ministepChoice = c("period", "chain", "ministep", "choice")
  )
  c(levels[[level]], condition)
}

makeGroupKey <- function(df, group_vars) {
  # should not convert to data.frame, adjust code to handle data.table if needed, but for now we convert to data.frame to avoid issues with list columns in data.table
  if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
    df <- as.data.frame(df)
  }
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
    paste(inc[["shortName"]], inc[["type"]], sep = "_")
}

conditionalReplace <- function(df, row_ids, cols, fun) {
  if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
    df[eval(substitute(row_ids)), (cols) := lapply(.SD, fun), .SDcols = cols]
  } else {
    df[row_ids, cols] <- fun(df[row_ids, cols])
  }
  df
}


agg <- function(outcomeName,
                data,
                level = "none",
                condition = NULL,
                sum_fun = mean,
                na.rm = TRUE) {
  if (length(outcomeName) != 1L) stop("'outcomeName' must be a single column name.")
  if (!is.null(condition)) condition  <- resolveCondition(condition)
  group_vars <- getGroupVars(level = level, condition = condition)
  # Helper for complex output (vector/list output from sum_fun)
  # extract from agg?
  expand_summary <- function(val) {
    if (length(val) == 1 && (is.null(names(val)) || names(val) == "")) {
      setNames(list(val), outcomeName)
    } else {
      as.list(val)
    }
  }

  # Use data.table if available
  if (requireNamespace("data.table", quietly = TRUE)) {
    if (!data.table::is.data.table(data)) {
      data <- data.table::as.data.table(data)
    }
    result <- data[, expand_summary(sum_fun(get(outcomeName), na.rm = na.rm)), 
      by = group_vars]
    return(result)
  }
  # ---- Fallback base R: ----
  if (length(group_vars) == 0) {
    output <- expand_summary(sum_fun(data[[outcomeName]], na.rm = na.rm))
    output <- as.data.frame(output)
    return(output)
  }

  grouping <- interaction(data[, group_vars, drop = FALSE], drop = TRUE)
  split_data <- split(data, grouping)
  agg_list <- lapply(split_data, function(subdf) {
    vals <- subdf[[outcomeName]]
    res <- expand_summary(sum_fun(vals, na.rm = na.rm))
    # Add group columns
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

#' Align theta by name and subset to non-rate effects.
#' effectNames come from contributions structs (shortName_type, no depvar prefix).
#' theta may be named with depvar-prefixed names (from nameThetaFromEffects).
#' Positional matching via req handles both formats.
alignThetaNoRate <- function(theta, effectNames, ans = NULL) {
    # Fast path: theta is already named with the same format as effectNames
    if (!is.null(names(theta)) && all(effectNames %in% names(theta))) {
        return(theta[effectNames])
    }
    if (!is.null(ans)) {
        req <- ans$requestedEffects
        if (is.data.frame(req)) {
            noRate   <- req$type != "rate"
            # Full names (may include depvar prefix): match theta's naming convention
            fullNm   <- getNamesFromEffects(req)
            # Short names (shortName_type, no depvar): match contributions effectNames format
            shortNm  <- paste(req$shortName, req$type, sep = "_")
            useNoRate <- sum(noRate) == length(theta)
            src_full  <- if (useNoRate) fullNm[noRate]  else fullNm
            src_short <- if (useNoRate) shortNm[noRate] else shortNm
            pos <- match(effectNames, src_full)
            if (anyNA(pos)) pos <- match(effectNames, src_short)
            if (!anyNA(pos)) {
                result        <- theta[pos]
                names(result) <- effectNames
                return(result)
            }
        }
    }
    theta[seq_along(effectNames)]
}

#' Attach effect contribution columns from a matrix to a data.frame.
#' Extracts column-by-column (avoids full-matrix copy).
#' If \code{flip = TRUE}, negates non-density columns where density == -1.
attachContribColumns <- function(out, effectNames, contrib, flip = TRUE) {
  density_j <- grep("density", effectNames, fixed = TRUE)
  if (flip && length(density_j) > 0) {
    neg1 <- contrib[, density_j[1L]] == -1L
    neg1[is.na(neg1)] <- FALSE          # guard against NA in contrib matrix
    has_neg1 <- any(neg1)
  } else {
    has_neg1 <- FALSE
  }
  flip_idx <- grep("density", effectNames, fixed = TRUE, invert = TRUE)
  for (j in seq_along(effectNames)) {
    col <- contrib[, j]
    if (has_neg1 && j %in% flip_idx) col[neg1] <- -col[neg1]
    out[[effectNames[j]]] <- col
  }
  out
}

#' Extract group-column vectors from a Contributions struct as a named list.
#' \code{keep} is an optional logical/integer subset vector.
#' For static, pb$group is scalar (recycled by data.frame()).
groupColsList <- function(pb, keep = NULL) {
  if (!is.null(pb$chain)) {
    if (is.null(keep)) {
      list(chain = pb$chain, group = pb$group, period = pb$period,
           ministep = pb$ministep, choice = pb$choice)
    } else {
      list(chain = pb$chain[keep], group = pb$group[keep],
           period = pb$period[keep], ministep = pb$ministep[keep],
           choice = pb$choice[keep])
    }
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