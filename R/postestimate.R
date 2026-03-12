sienaPostestimate <- function(
    predictFun,
    predictArgs,
    outcome,
    level = "period",
    condition = NULL,
    sum.fun = mean,
    na.rm = TRUE,
    theta_hat,
    cov_theta,
    uncertainty = TRUE,
    uncertainty_mode = c("batch", "stream"),
    nsim = 1000,
    uncertainty_sd = TRUE,
    uncertainty_ci = TRUE,
    uncertainty_probs = c(0.025, 0.5, 0.975),
    uncertainty_mcse = FALSE,
    uncertainty_mcse_batches = NULL,
    useCluster = FALSE,
    nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batch_dir = "temp",
    prefix = "simBatch_b",
    combine_batch = TRUE,
    batch_size = NULL,
    keep_batch = FALSE,
    verbose = TRUE,
    useChangeContributions = NULL,
    gc_each_batch = TRUE,
    gc_each_sim = FALSE
) {
  uncertainty_mode <- match.arg(uncertainty_mode)
    if (length(outcome) != 1L) {
        stop("'outcome' must be a single column name.")
    }

    estimator <- makeEstimator(
        predictFun = predictFun,
        predictArgs = predictArgs,
        outcome = outcome,
        level = level,
        condition = condition,
        sum.fun = sum.fun,
        na.rm = na.rm
    )

    if (!is.null(useChangeContributions)) {
        expect <- estimator(theta_hat, 
          useChangeContributions = useChangeContributions)
    } else {
        expect <- estimator(theta_hat)
    }

    if (!uncertainty) return(expect)

    uncertainty_summary_fun <- make_uncertainty_summarizer(
      return_sd = uncertainty_sd,
      return_ci = uncertainty_ci,
      probs = uncertainty_probs,
      return_mcse = uncertainty_mcse,
      mcse_batches = uncertainty_mcse_batches
    )

    if (identical(uncertainty_mode, "stream")) {
      uncert <- drawSim_stream(
        estimator = estimator,
        outcome = outcome,
        level = level,
        condition = condition,
        uncertainty_summary_fun = uncertainty_summary_fun,
        theta_hat = theta_hat,
        cov_theta = cov_theta,
        nbrNodes = if (useCluster) nbrNodes else 1L,
        nsim = nsim,
        clusterType = clusterType,
        cluster = cluster,
        batch_size = batch_size,
        verbose = verbose,
        gc_each_batch = gc_each_batch,
        gc_each_sim = gc_each_sim
      )
    } else {
      uncert <- drawSim(
        estimator = estimator,
        theta_hat = theta_hat,
        cov_theta = cov_theta,
        nbrNodes = if (useCluster) nbrNodes else 1L,
        nsim = nsim,
        clusterType = clusterType,
        cluster = cluster,
        batch_dir = batch_dir,
        prefix = prefix,
        combine_batch = combine_batch,
        batch_size = batch_size,
        keep_batch = keep_batch,
        verbose = verbose,
        gc_each_batch = gc_each_batch,
        gc_each_sim = gc_each_sim
      )
      uncert <- agg(outcome, uncert, level = level, condition = condition, sum.fun = uncertainty_summary_fun)
    }
    mergeEstimates(expect, uncert, level = level, condition = condition)
}

make_uncertainty_summarizer <- function(
    return_sd = TRUE,
    return_ci = TRUE,
    probs = c(0.025, 0.5, 0.975),
    return_mcse = FALSE,
    mcse_batches = NULL
) {
  if (length(probs) != 3L) {
    stop("'uncertainty_probs' must be length 3: lower, median, upper.")
  }
  probs <- as.numeric(probs)
  if (any(!is.finite(probs)) || any(probs <= 0) || any(probs >= 1)) {
    stop("'uncertainty_probs' must be strictly between 0 and 1.")
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
      out <- c(out, compute_mcse_from_draws(
        x = x,
        return_sd = return_sd,
        return_ci = return_ci,
        probs = probs,
        batches = mcse_batches
      ))
    }

    out
  }
}

compute_mcse_from_draws <- function(x,
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
      stop("'uncertainty_mcse_batches' must be NULL or an integer >= 2.")
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

makeEstimator <- function(predictFun, predictArgs, outcome,
    level = "period", condition = NULL, sum.fun = mean, na.rm = TRUE) {
    function(theta, useChangeContributions = FALSE) {
        if(!is.null(predictArgs[["useChangeContributions"]])){
          predictArgs[["useChangeContributions"]] <- useChangeContributions
        }
        callArgs <- c(list(theta = theta), predictArgs)
        unit_pred <- do.call(predictFun, callArgs)
        main_result <- agg(
            ME = outcome,
            data = unit_pred,
            level = level,
            condition = condition,
            sum.fun = sum.fun,
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
    theta_hat,
    cov_theta,
    nbrNodes = 1,
    nsim = 1000,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batch_dir = "temp",
    prefix = "simBatch_b",
    batch_size = NULL,
    combine_batch = TRUE,
    keep_batch = FALSE,
    verbose = TRUE,
    gc_each_batch = TRUE,
    gc_each_sim = FALSE
) {

    cli <- open_cluster(nbrNodes, clusterType, cluster,
      export_vars  = c("estimator", "theta_hat", "cov_theta"),
      export_envir = environment(),
      verbose      = verbose
    )
    if (!dir.exists(batch_dir)) dir.create(batch_dir, recursive = TRUE)
    nbatches <- ceiling(nsim / batch_size) 
    if (verbose) {
      message(
        "Starting simulations: nsim=", nsim,
        ", batch_size=", batch_size,
        ", batches=", nbatches,
        ". (Last batch may be smaller.)"
      )
    }
    
    for (b in seq_len(nbatches)) {
      first_sim <- 1 + (b - 1) * batch_size
      last_sim  <- min(b * batch_size, nsim)
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
          theta_sim <- MASS::mvrnorm(1, mu = theta_hat, Sigma = cov_theta)
          sim <- do.call(estimator, list(theta = theta_sim))
          sim[, "sim"] <- i
          sim
        })
      } else if (cli$useFORK) {
        batch_result <- parallel::mclapply(batch, function(i) {
          if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
          theta_sim <- MASS::mvrnorm(1, mu = theta_hat, Sigma = cov_theta)
          sim <- do.call(estimator, list(theta = theta_sim))
          sim[, "sim"] <- i
          sim
        }, mc.cores = nbrNodes)
      } else {
        batch_result <- vector("list", length(batch))
        for (k in seq_along(batch)) {
          i <- batch[k]
          theta_sim <- MASS::mvrnorm(1, mu = theta_hat, Sigma = cov_theta)
          sim <- do.call(estimator, list(theta = theta_sim))
          sim[, "sim"] <- i
          batch_result[[k]] <- sim
          if (gc_each_sim) gc(verbose = FALSE)
        }
      }

      file_name <- file.path(batch_dir, sprintf("%s%03d.rds", prefix, b))
      saveRDS(batch_result, file = file_name)
      rm(batch_result)
      if (gc_each_batch) gc(verbose = FALSE)
      if (verbose) {
        done <- last_sim
        pct <- round(100 * done / nsim)
        message(sprintf("Completed batch %d/%d (%d/%d sims, %d%%)", b, nbatches, done, nsim, pct))
      }
    }
    
    close_cluster(cli, verbose)
    
    if (verbose) message("All batches complete. Now combining and cleaning up...")
    
    if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads()
    if(combine_batch){
      # Combine all batches
      batch_files <- list.files(batch_dir, 
        pattern = sprintf("^%s\\d{3}\\.rds$", prefix), full.names = TRUE)
      batch_files <- sort(batch_files)
      combined_results <- vector("list", nsim)
      
      for (b in seq_along(batch_files)) {
        file_name <- batch_files[b]
        results <- readRDS(file_name)
        first <- 1 + (b - 1) * batch_size
        last  <- min(b * batch_size, nsim)
        combined_results[first:last] <- results
        if(!keep_batch){
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

drawSim_stream <- function(
    estimator,
    outcome,
    level,
    condition,
    uncertainty_summary_fun,
    theta_hat,
    cov_theta,
    nbrNodes = 1,
    nsim = 1000,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batch_size = NULL,
    verbose = TRUE,
    gc_each_batch = TRUE,
    gc_each_sim = FALSE
) {
    group_vars <- get_group_vars(level = level, condition = condition)

    cli <- open_cluster(nbrNodes, clusterType, cluster,
      export_vars  = c("estimator", "theta_hat", "cov_theta"),
      export_envir = environment(),
      verbose      = verbose
    )
    nbatches <- ceiling(nsim / batch_size)
    if (verbose) {
      message(
        "Starting streaming simulations: nsim=", nsim,
        ", batch_size=", batch_size,
        ", batches=", nbatches,
        ". (Last batch may be smaller.)"
      )
    }

    stream_state <- new.env(parent = emptyenv(), hash = TRUE)
    group_state <- new.env(parent = emptyenv(), hash = TRUE)

    for (b in seq_len(nbatches)) {
      first_sim <- 1 + (b - 1) * batch_size
      last_sim  <- min(b * batch_size, nsim)
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
          theta_sim <- MASS::mvrnorm(1, mu = theta_hat, Sigma = cov_theta)
          do.call(estimator, list(theta = theta_sim))
        })
      } else if (cli$useFORK) {
        batch_result <- parallel::mclapply(batch, function(i) {
          if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
          theta_sim <- MASS::mvrnorm(1, mu = theta_hat, Sigma = cov_theta)
          do.call(estimator, list(theta = theta_sim))
        }, mc.cores = nbrNodes)
      } else {
        batch_result <- vector("list", length(batch))
        for (k in seq_along(batch)) {
          theta_sim <- MASS::mvrnorm(1, mu = theta_hat, Sigma = cov_theta)
          batch_result[[k]] <- do.call(estimator, list(theta = theta_sim))
          if (gc_each_sim) gc(verbose = FALSE)
        }
      }

      for (sim_df in batch_result) {
        update_stream_state(
          stream_state = stream_state,
          group_state = group_state,
          sim_df = sim_df,
          outcome = outcome,
          group_vars = group_vars
        )
      }

      rm(batch_result)
      if (gc_each_batch) gc(verbose = FALSE)
      if (verbose) {
        done <- last_sim
        pct <- round(100 * done / nsim)
        message(sprintf("Completed batch %d/%d (%d/%d sims, %d%%)", b, nbatches, done, nsim, pct))
      }
    }

    close_cluster(cli, verbose)

    finalize_stream_state(
      stream_state = stream_state,
      group_state = group_state,
      group_vars = group_vars,
      uncertainty_summary_fun = uncertainty_summary_fun
    )
}

open_cluster <- function(nbrNodes, clusterType, cluster, export_vars, export_envir, verbose) {
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

close_cluster <- function(cli, verbose) {
  if (cli$cluster_created) {
    parallel::stopCluster(cli$cl)
    if (verbose) message("Stopped internal cluster.")
  }
}

update_stream_state <- function(stream_state,
                                group_state,
                                sim_df,
                                outcome,
                                group_vars) {
  if (length(outcome) != 1L) stop("'outcome' must be a single column name.")

#we should not convert to data.frame, adjust code to handle data.table if needed, but for now we convert to data.frame to avoid issues with list columns in data.table
  if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(sim_df)) {
    sim_df <- as.data.frame(sim_df)
  }

  keys <- make_group_key(sim_df, group_vars)
  vals <- sim_df[[outcome]]

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

finalize_stream_state <- function(stream_state,
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

get_group_vars <- function(level = "none", condition = NULL) {
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

make_group_key <- function(df, group_vars) {
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
    includedEffects[["shortName"]][keep]
}

conditionalReplace <- function(df, row_ids, cols, fun) {
  if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
    df[eval(substitute(row_ids)), (cols) := lapply(.SD, fun), .SDcols = cols]
  } else {
    df[row_ids, cols] <- fun(df[row_ids, cols])
  }
  df
}


agg <- function(ME,
                data,
                level = "none",
                condition = NULL,
                sum.fun = mean,
                na.rm = TRUE) {
  if (length(ME) != 1L) stop("'ME' must be a single column name.")
  group_vars <- get_group_vars(level = level, condition = condition)

  # Helper for complex output (vector/list output from sum.fun)
  # extract from agg?
  expand_summary <- function(val) {
    if (length(val) == 1 && (is.null(names(val)) || names(val) == "")) {
      setNames(list(val), ME)
    } else {
      as.list(val)
    }
  }

  # Use data.table if available
  if (requireNamespace("data.table", quietly = TRUE)) {
    if (!data.table::is.data.table(data)) {
      data <- data.table::as.data.table(data)
    }
    result <- data[, expand_summary(sum.fun(get(ME), na.rm = na.rm)), 
      by = group_vars]
    return(result)
  }
  # ---- Fallback base R: ----
  if (length(group_vars) == 0) {
    output <- expand_summary(sum.fun(data[[ME]], na.rm = na.rm))
    output <- as.data.frame(output)
    return(output)
  }

  grouping <- interaction(data[, group_vars, drop = FALSE], drop = TRUE)
  split_data <- split(data, grouping)
  agg_list <- lapply(split_data, function(subdf) {
    vals <- subdf[[ME]]
    res <- expand_summary(sum.fun(vals, na.rm = na.rm))
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

#' Align theta by name and subset to non-rate effects
alignThetaNoRate <- function(theta, effectNames, ans = NULL) {
    if (is.null(names(theta)) && !is.null(ans)) {
        # Use requestedEffects (not effects): the latter may include extra base
        # effects injected for interactions, causing a length mismatch with theta.
        allNames <- ans$requestedEffects$shortName
        if (length(allNames) == length(theta)) {
            names(theta) <- allNames
        }
    }
    if (!is.null(names(theta))) {
        theta[effectNames]
    } else {
        # positional fallback when names cannot be determined
        theta[seq_along(effectNames)]
    }
}

#' Attach effect contribution columns from a matrix to a data.frame.
#' Extracts column-by-column (avoids full-matrix copy).
#' If \code{flip = TRUE}, negates non-density columns where density == -1.
.attachContribColumns <- function(out, effectNames, contrib, flip = TRUE) {
  if (flip && "density" %in% effectNames) {
    neg1 <- contrib[, "density"] == -1L
    has_neg1 <- any(neg1)
  } else {
    has_neg1 <- FALSE
  }
  flip_idx <- which(effectNames != "density")
  for (j in seq_along(effectNames)) {
    col <- contrib[, j]
    if (has_neg1 && j %in% flip_idx) col[neg1] <- -col[neg1]
    out[[effectNames[j]]] <- col
  }
  out
}

#' Extract group-column vectors from a staticContributions struct as a named list.
#' \code{keep} is an optional logical/integer subset vector.
#' For static, pb$group is scalar (recycled by data.frame()).
.groupColsList <- function(pb, keep = NULL) {
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
      list(group = pb$group, period = pb$period[keep],
           ego = pb$ego[keep], choice = pb$choice[keep])
    }
  }
}