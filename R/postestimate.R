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
    nsim = 1000,
    uncertaintySd = TRUE,
    uncertaintyCi = TRUE,
    uncertaintyMean = FALSE,
    uncertaintyMedian = FALSE,
    uncertaintyProbs = c(0.025, 0.5, 0.975),
    uncertaintyMcse = FALSE,
    uncertaintymcseBatches = NULL,
    useCluster = FALSE,
    nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    cl = NULL,
    batchDir = "temp",
    prefix = "simBatch_b",
    combineBatch = TRUE,
    batchSize = NULL,
    keepBatch = FALSE,
    verbose = TRUE,
    useChangeContributions = NULL,
    gcEachBatch = FALSE,
    gcEachSim = FALSE,
    metadata = NULL,
    egoNormalize = TRUE,
    returnDecisionDetails = FALSE
) {
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


    uncertainty_summary_fun <- makeUncertaintySummarizer(
      return_sd = uncertaintySd,
      return_ci = uncertaintyCi,
      return_mean = uncertaintyMean,
      return_median = uncertaintyMedian,
      probs = uncertaintyProbs,
      return_mcse = uncertaintyMcse,
      mcseBatches = uncertaintymcseBatches
    )

    raw_sims <- drawSim(
      estimator = estimator,
      thetaHat = thetaHat,
      covTheta = covTheta,
      nbrNodes = if (useCluster) nbrNodes else 1L,
      nsim = nsim,
      clusterType = clusterType,
      cl = cl,
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
      raw_sims, 
      level = level, 
      condition = condition, 
      sum_fun = uncertainty_summary_fun
    )
    result <- mergeEstimates(expect, uncert, level = level, condition = condition)

    # Mass contrast uncertainty
    massCols <- intersect(c("massCreation", "massDissolution"), names(raw_sims))
    for (mc in massCols) {
        mc_uncert <- agg(mc, raw_sims, level = level, condition = NULL,
                         sum_fun = uncertainty_summary_fun)
        level_by <- intersect(
            getGroupVars(level = level, condition = NULL),
            names(mc_uncert)
        )
        uc_cols <- setdiff(names(mc_uncert), level_by)
        for (uc in uc_cols) {
            names(mc_uncert)[names(mc_uncert) == uc] <- paste0(mc, "_", uc)
        }
        if (length(level_by) > 0L) {
            result <- merge(result, mc_uncert, by = level_by,
                            all.x = TRUE, sort = FALSE)
        } else {
            for (col in setdiff(names(mc_uncert), level_by)) {
                result[[col]] <- mc_uncert[[col]]
            }
        }
    }

    result <- stampPostestimate(result, metadata)
    if (!is.null(decisionDetails)) attr(result, "decisionDetails") <- decisionDetails
    result
}

makeUncertaintySummarizer <- function(
    return_sd = TRUE,
    return_ci = TRUE,
    return_mean = FALSE,
    return_median = FALSE,
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

    out <- list()

    if (isTRUE(return_mean)) {
      out[["Mean"]] <- if (n) mean(x) else NA_real_
    }

    if (isTRUE(return_sd)) {
      out[["SE"]] <- if (n >= 2L) stats::sd(x) else NA_real_
    }

    if (isTRUE(return_ci) || isTRUE(return_median)) {
      if (n) {
        qu <- unname(stats::quantile(x, probs = probs, na.rm = FALSE))
      } else {
        qu <- rep(NA_real_, 3L)
      }
      if (isTRUE(return_median)) out[["Median"]] <- qu[2]
      if (isTRUE(return_ci)) {
        out[["q_025"]] <- qu[1]
        out[["q_975"]] <- qu[3]
      }
    }

    out[["cases"]] <- n

    if (isTRUE(return_mcse)) {
      out <- c(out, computeMcse(
        x = x,
        return_sd = return_sd,
        return_ci = return_ci,
        return_mean = return_mean,
        return_median = return_median,
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
                                   return_mean = FALSE,
                                   return_median = FALSE,
                                   probs = c(0.025, 0.5, 0.975),
                                   batches = NULL) {
  n <- length(x)
  need_q <- isTRUE(return_ci) || isTRUE(return_median)
  if (n < 2L) {
    out <- list()
    if (isTRUE(return_mean)) out[["mcse_Mean"]] <- NA_real_
    if (isTRUE(return_sd)) out[["mcse_SE"]] <- NA_real_
    if (isTRUE(return_median)) out[["mcse_Median"]] <- NA_real_
    if (isTRUE(return_ci)) {
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

  batch_stats_mean <- if (isTRUE(return_mean)) numeric(B) else NULL
  batch_stats_sd <- if (isTRUE(return_sd)) numeric(B) else NULL
  batch_stats_q <- if (need_q) matrix(NA_real_, nrow = B, ncol = 3L) else NULL

  for (b in seq_len(B)) {
    idx <- ((b - 1L) * m + 1L):(b * m)
    xb <- x_use[idx]
    if (isTRUE(return_mean)) batch_stats_mean[b] <- mean(xb)
    if (isTRUE(return_sd)) {
      batch_stats_sd[b] <- if (m >= 2L) stats::sd(xb) else NA_real_
    }
    if (need_q) {
      batch_stats_q[b, ] <- unname(stats::quantile(xb, probs = probs, na.rm = FALSE))
    }
  }

  out <- list()

  if (isTRUE(return_mean)) {
    out[["mcse_Mean"]] <- stats::sd(batch_stats_mean) / sqrt(B)
  }

  if (isTRUE(return_sd)) {
    out[["mcse_SE"]] <- stats::sd(batch_stats_sd) / sqrt(B)
  }

  if (isTRUE(return_median)) {
    out[["mcse_Median"]] <- stats::sd(batch_stats_q[, 2]) / sqrt(B)
  }

  if (isTRUE(return_ci)) {
    out[["mcse_q_025"]] <- stats::sd(batch_stats_q[, 1]) / sqrt(B)
    out[["mcse_q_975"]] <- stats::sd(batch_stats_q[, 3]) / sqrt(B)
  }

  out
}

# ---- Shared-cache infrastructure (Step 3 of the refactoring plan) ----------
#
# A shared cache allows multiple estimators (effects) to share the expensive
# contribution-matrix and baseline-probability computations within a single
# theta draw.  The cache is a private environment, never exposed to the user.
#
# Static path:  getContribFun(theta) returns the same pre-built struct
#               regardless of theta, so the cache effectively stores it once.
# Dynamic path: getContribFun(theta) re-simulates chains; the cache stores
#               the result for the current theta and invalidates on change.
#
# Memory: only the *last* theta's results are kept.  No accumulation.

makeSharedCache <- function() {
    env <- new.env(parent = emptyenv())
    env$contribTheta  <- NULL
    env$contrib       <- NULL
    env$baselineTheta <- NULL
    env$baseline      <- NULL
    env
}

getCachedContrib <- function(cache, getContribFun, theta) {
    if (!identical(cache$contribTheta, theta)) {
        cache$contrib      <- getContribFun(theta)
        cache$contribTheta <- theta
        # Invalidate baseline when contributions change
        cache$baselineTheta <- NULL
        cache$baseline      <- NULL
    }
    cache$contrib
}

getCachedBaseline <- function(cache, cc, theta, type) {
    if (!identical(cache$baselineTheta, theta)) {
        cache$baseline      <- predictProbability(cc, theta, type,
                                                   returnComponents = TRUE)
        cache$baselineTheta <- theta
    }
    cache$baseline
}

makeEstimator <- function(predictFun, predictArgs, outcomeName,
    level = "period", condition = NULL, sum_fun = mean, na.rm = TRUE,
    egoNormalize = TRUE) {
    function(theta, useChangeContributions = FALSE,
             returnDecisionDetails = FALSE) {
        if(!is.null(predictArgs[["useChangeContributions"]])){
          predictArgs[["useChangeContributions"]] <- useChangeContributions
        }
        # Only attach contribution columns when the caller needs decision details.
        # Skipping this saves materializing the full contribMat into a data.frame
        # for every theta draw, which is expensive for large networks.
        if (!is.null(predictArgs[["attachContribs"]])) {
          predictArgs[["attachContribs"]] <- returnDecisionDetails ||
              !is.null(condition)
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
            # Join by the level-based group variables so that the mass contrast
            # value for each period (and group, if present) is matched to all
            # condition × period rows in main_result rather than recycled.
            mc_by <- intersect(
                getGroupVars(level = level, condition = NULL),
                intersect(names(main_result), names(mc_agg))
            )
            if (length(mc_by) > 0L) {
                main_result <- merge(main_result, mc_agg, by = mc_by,
                                     all.x = TRUE, sort = FALSE)
            } else {
                # level = "none": single scalar result — assign directly.
                main_result[[mc]] <- mc_agg[[mc]]
            }
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
    cl = NULL,
    batchDir = "temp",
    prefix = "simBatch_b",
    batchSize = NULL,
    combineBatch = TRUE,
    keepBatch = FALSE,
    verbose = TRUE,
    gcEachBatch = FALSE,
    gcEachSim = FALSE
) {

    cli <- openCluster(nbrNodes, clusterType, cl,
      export_vars  = c("estimator", "thetaHat", "covTheta"),
      export_envir = environment(),
      verbose      = verbose
    )
    nbatches <- ceiling(nsim / batchSize)
    skip_disk <- (nbatches == 1L && !keepBatch)
    if (!skip_disk && !dir.exists(batchDir)) dir.create(batchDir, recursive = TRUE)
    in_memory_batch <- NULL
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
          Sys.setenv(OMP_NUM_THREADS = "1",
                     OPENBLAS_NUM_THREADS = "1",
                     MKL_NUM_THREADS = "1")
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
      if (skip_disk) {
        in_memory_batch <- batch_result
      } else {
        saveRDS(batch_result, file = file_name)
      }
      rm(batch_result)
      if (gcEachBatch) gc(verbose = FALSE)
      if (verbose) {
        done <- last_sim
        pct <- round(100 * done / nsim)
        message(sprintf("Completed batch %d/%d (%d/%d sims, %d%%)", b, nbatches, done, nsim, pct))
      }
    }
    
    closeCluster(cli, verbose)
    
    if (verbose) message("All batches complete. Now combining and cleaning up...")
    
    # if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads() # data.table removed
    if(combineBatch){
      # Single-batch in-memory path: no disk I/O
      if (skip_disk) {
        out <- do.call(rbind, in_memory_batch)
        rownames(out) <- NULL
        return(out)
      }
      # Combine all batches from disk
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


# ---- Batch-shared simulation loop ----------------------------------------
#
# Runs one set of theta draws and evaluates N effects per draw, sharing the
# expensive contribution-matrix computation step.
#
# Static path:  getStaticChangeContributions() computed once (theta-independent),
#               then N effects × nsim draws evaluated cheaply from the same matrix.
# Dynamic path: getDynamicChangeContributions() runs n3 forward sims ONCE per
#               theta draw (instead of N × n3 per draw), all N effects evaluated
#               from that one set of simulated ministeps.
#
# getContribFun(theta) -> wide changeContributions struct
# effectSpecList: named list of per-effect specs (see marginalEffects effectList)
#
# Returns a named list of raw_sims data frames (one per effect).
drawSimBatch <- function(
    getContribFun,
    effectSpecList,
    thetaHat,
    covTheta,
    type         = "changeProb",
    nsim         = 1000L,
    useCluster   = FALSE,
    nbrNodes     = 1L,
    clusterType  = c("PSOCK", "FORK"),
    cl      = NULL,
    batchSize    = 100L,
    batchDir     = "temp",
    prefix       = "simBatchB_b",
    keepBatch    = FALSE,
    verbose      = TRUE,
    gcEachBatch  = FALSE,
    gcEachSim    = FALSE,
    egoNormalize = TRUE
) {
    N <- length(effectSpecList)
    if (N == 0L) stop("'effectSpecList' must not be empty.")
    eff_names <- names(effectSpecList)
    if (is.null(eff_names) || any(eff_names == ""))
        stop("All elements of 'effectSpecList' must be named.")

    nbrNodes <- as.integer(nbrNodes)
    ct       <- match.arg(clusterType, c("PSOCK", "FORK"))
    useFork  <- isTRUE(useCluster) && (ct == "FORK") &&
                    (nbrNodes > 1L) && (.Platform$OS.type != "windows")
    usePSOCK <- isTRUE(useCluster) && (ct == "PSOCK") && (nbrNodes > 1L)

    nbatches <- ceiling(nsim / batchSize)
    # When there is only one batch and we are not asked to keep it, skip
    # disk I/O entirely and hold the result in memory.
    skip_disk <- (nbatches == 1L && !keepBatch)

    if (!skip_disk && !dir.exists(batchDir))
        dir.create(batchDir, recursive = TRUE)

    mode_label <- if (useFork) paste0("FORK (", nbrNodes, " cores)")
                  else if (usePSOCK) paste0("PSOCK (", nbrNodes, " workers)")
                  else "sequential"

    if (verbose) {
        message("Batch-shared simulation: nsim=", nsim, ", effects=", N,
                ", batchSize=", batchSize, ", batches=", nbatches,
                ", mode=", mode_label)
    }

    t_loop_start   <- proc.time()
    sims_computed  <- 0L  # sims actually run this session (not loaded from disk)
    batches_done   <- 0L  # completed batches (for ETA)

    # ---- Batch loop with resume support -----------------------------------
    # If a batch file already exists on disk (from a previous crashed run),
    # skip recomputing it.  This makes the server pipeline resumable.
    batch_pattern <- sprintf("^%s\\d{3}\\.rds$", prefix)

    in_memory_batch <- NULL  # used only when skip_disk == TRUE

    for (b in seq_len(nbatches)) {
        first_sim <- 1L + (b - 1L) * batchSize
        last_sim  <- min(b * batchSize, nsim)
        batch_idx <- first_sim:last_sim
        bsz <- length(batch_idx)

        fname <- file.path(batchDir, sprintf("%s%03d.rds", prefix, b))
        if (file.exists(fname)) {
            if (verbose) {
                message(sprintf("  Batch %d/%d: found on disk, skipping.",
                                b, nbatches))
            }
            next
        }

        if (verbose) {
            message(sprintf("  Batch %d/%d: sims %d-%d (%d draws, %s)",
                            b, nbatches, first_sim, last_sim, bsz, mode_label))
        }

        t_batch_start <- proc.time()

        # bpe[[j]][[k]]: per-draw result for effect j, k-th sim in this batch
        bpe <- vector("list", N)
        for (j in seq_len(N)) bpe[[j]] <- vector("list", bsz)

        # Function to run one draw and return per-effect results
        run_one_draw <- function(k) {
            i <- batch_idx[k]
            theta_sim <- MASS::mvrnorm(1L, mu = thetaHat, Sigma = covTheta)
            cc <- getContribFun(theta_sim)
            baseline <- predictProbability(cc, theta_sim, type,
                                           returnComponents = TRUE)
            theta_use <- baseline$theta_use

            per_effect <- vector("list", N)
            for (j in seq_len(N)) {
                spec      <- effectSpecList[[j]]
                unit_pred <- do.call(
                    spec$diffFun,
                    c(list(changeContributions = cc, theta_use = theta_use,
                           baseline = baseline),
                      spec$diffArgs)
                )
                out_nm <- spec$outcomeName
                lvl    <- if (!is.null(spec$level))     spec$level     else "period"
                cond   <- if (!is.null(spec$condition)) spec$condition else NULL

                row_j <- agg(out_nm, unit_pred,
                             level = lvl, condition = cond,
                             sum_fun = mean, na.rm = TRUE,
                             egoNormalize = egoNormalize)
                row_j[["sim"]] <- i

                massCols <- intersect(c("massCreation", "massDissolution"),
                                      names(unit_pred))
                for (mc in massCols) {
                    mc_row <- agg(mc, unit_pred,
                                  level = lvl, condition = NULL,
                                  sum_fun = mean, na.rm = TRUE,
                                  egoNormalize = egoNormalize)
                    mc_by  <- intersect(
                        getGroupVars(level = lvl, condition = NULL),
                        intersect(names(row_j), names(mc_row))
                    )
                    if (length(mc_by) > 0L) {
                        row_j <- merge(row_j, mc_row, by = mc_by,
                                       all.x = TRUE, sort = FALSE)
                    } else {
                        row_j[[mc]] <- mc_row[[mc]]
                    }
                }
                per_effect[[j]] <- row_j
            }
            per_effect
        }

        if (usePSOCK) {
            # Export needed objects to workers each batch
            parallel::clusterExport(
                cl,
                c("getContribFun", "effectSpecList", "thetaHat", "covTheta",
                  "type", "N", "bsz", "batch_idx", "egoNormalize"),
                envir = environment()
            )
            draw_results <- parallel::parLapply(cl, seq_len(bsz), run_one_draw)
            for (k in seq_len(bsz)) {
                for (j in seq_len(N)) bpe[[j]][[k]] <- draw_results[[k]][[j]]
            }
            sims_computed <- sims_computed + bsz
        } else if (useFork) {
            # Set thread limits in parent BEFORE forking — forked children
            # inherit the parent's BLAS/OpenMP thread pool at fork time.
            old_omp  <- Sys.getenv("OMP_NUM_THREADS",      unset = NA)
            old_blas <- Sys.getenv("OPENBLAS_NUM_THREADS",  unset = NA)
            old_mkl  <- Sys.getenv("MKL_NUM_THREADS",       unset = NA)
            Sys.setenv(OMP_NUM_THREADS      = "1",
                       OPENBLAS_NUM_THREADS  = "1",
                       MKL_NUM_THREADS       = "1")
            draw_results <- parallel::mclapply(
                seq_len(bsz), run_one_draw,
                mc.cores = nbrNodes, mc.set.seed = TRUE,
                mc.preschedule = FALSE
            )
            # Restore parent thread settings
            if (is.na(old_omp))  Sys.unsetenv("OMP_NUM_THREADS")      else Sys.setenv(OMP_NUM_THREADS      = old_omp)
            if (is.na(old_blas)) Sys.unsetenv("OPENBLAS_NUM_THREADS") else Sys.setenv(OPENBLAS_NUM_THREADS = old_blas)
            if (is.na(old_mkl))  Sys.unsetenv("MKL_NUM_THREADS")     else Sys.setenv(MKL_NUM_THREADS      = old_mkl)
            # Validate mclapply results — detect OOM-killed workers
            failed <- vapply(draw_results, function(x) {
                inherits(x, "try-error") || is.null(x)
            }, logical(1L))
            if (any(failed)) {
                n_fail <- sum(failed)
                stop(sprintf(
                    "mclapply: %d of %d workers failed (likely OOM or signal). ",
                    n_fail, bsz),
                    "Consider reducing nbrNodes or batchSize.")
            }
            for (k in seq_len(bsz)) {
                for (j in seq_len(N)) bpe[[j]][[k]] <- draw_results[[k]][[j]]
            }
            sims_computed <- sims_computed + bsz
        } else {
            for (k in seq_len(bsz)) {
                t_sim_start <- proc.time()
                per_effect  <- run_one_draw(k)
                sims_computed <- sims_computed + 1L
                if (verbose) {
                    i <- batch_idx[k]
                    t_draw    <- (proc.time() - t_sim_start)[["elapsed"]]
                    t_elapsed <- (proc.time() - t_loop_start)[["elapsed"]]
                    eta_sec   <- if (sims_computed < nsim)
                        (t_elapsed / sims_computed) * (nsim - sims_computed)
                        else 0
                    message(sprintf("    [%d/%d] %.1fs/draw | ETA ~%.0fs",
                                    i, nsim, t_draw, eta_sec))
                }
                for (j in seq_len(N)) bpe[[j]][[k]] <- per_effect[[j]]
                if (gcEachSim) gc(verbose = FALSE)
            }
        }

        # Save batch to disk (or hold in memory if single-batch run)
        if (skip_disk) {
            in_memory_batch <- bpe
        } else {
            saveRDS(bpe, file = fname)
        }
        rm(bpe)
        if (gcEachBatch) gc(verbose = FALSE)
        batches_done <- batches_done + 1L
        if (verbose) {
            t_batch   <- (proc.time() - t_batch_start)[["elapsed"]]
            t_elapsed <- (proc.time() - t_loop_start)[["elapsed"]]
            draw_rate <- if (sims_computed > 0) sims_computed / t_elapsed else NA_real_
            remaining <- nsim - last_sim
            eta_sec   <- if (remaining > 0 && !is.na(draw_rate) && draw_rate > 0)
                             remaining / draw_rate else 0
            eta_str <- if (eta_sec > 3600) sprintf("%.1f h", eta_sec / 3600)
                       else if (eta_sec > 60) sprintf("%.0f min", eta_sec / 60)
                       else sprintf("%.0f s", eta_sec)
            mem_mb <- if (gcEachBatch) {
                # gc() was just called above — reuse its summary for reporting
                tryCatch(
                    round(sum(gc(reset = FALSE, verbose = FALSE)[, 2]) / 1024, 0),
                    error = function(e) NA_integer_
                )
            } else {
                NA_integer_  # skip gc() when disabled — it's expensive on large heaps
            }
            message(sprintf(
                "  Done batch %d/%d (%d%%) | %.1fs (%.2f draws/s) | ETA %s | mem ~%s MB",
                b, nbatches, round(100 * last_sim / nsim),
                t_batch, if (!is.na(draw_rate)) draw_rate else 0,
                eta_str, if (!is.na(mem_mb)) as.character(mem_mb) else "?"))
            # After the first real batch, report per-core throughput
            # (useful to compare across runs / node counts)
            if (batches_done == 1L && (useFork || usePSOCK)) {
                rate_per_core <- draw_rate / nbrNodes
                message(sprintf(
                    "  Throughput (batch 1): %.2f draws/s total, %.2f draws/s/core (%d cores)",
                    draw_rate, rate_per_core, nbrNodes))
            }
        }
    } # batch loop

    t_total <- (proc.time() - t_loop_start)[["elapsed"]]
    if (verbose) {
        total_str <- if (t_total > 3600) sprintf("%.1f h", t_total / 3600)
                     else if (t_total > 60) sprintf("%.1f min", t_total / 60)
                     else sprintf("%.0f s", t_total)
        final_rate <- if (sims_computed > 0) sims_computed / t_total else 0
        message(sprintf(
            "Simulation complete: %d draws in %s (%.2f draws/s, %s)",
            sims_computed, total_str, final_rate, mode_label))
        message("Combining batches...")
    }

    # Read batch files and flatten per-effect lists
    per_effect_raw <- vector("list", N)
    for (j in seq_len(N)) per_effect_raw[[j]] <- vector("list", nsim)

    if (!is.null(in_memory_batch)) {
        # Single-batch path: no disk I/O
        for (j in seq_len(N)) per_effect_raw[[j]] <- in_memory_batch[[j]]
    } else {
        batch_files <- sort(list.files(
            batchDir,
            pattern    = batch_pattern,
            full.names = TRUE
        ))
        for (b in seq_along(batch_files)) {
            bp    <- readRDS(batch_files[[b]])
            first <- 1L + (b - 1L) * batchSize
            last  <- min(b * batchSize, nsim)
            for (j in seq_len(N)) {
                per_effect_raw[[j]][first:last] <- bp[[j]]
            }
            if (!keepBatch) file.remove(batch_files[[b]])
        }
    }

    # Collapse each effect's list of per-draw rows into one data frame
    # Column-wise unlist is faster than do.call(rbind, ...) on many small dfs
    out <- lapply(per_effect_raw, function(draws) {
        nms <- names(draws[[1L]])
        cols <- setNames(vector("list", length(nms)), nms)
        for (nm in nms) {
            cols[[nm]] <- unlist(lapply(draws, `[[`, nm), use.names = FALSE)
        }
        attr(cols, "row.names") <- .set_row_names(length(cols[[1L]]))
        class(cols) <- "data.frame"
        cols
    })
    names(out) <- eff_names
    out
}


openCluster <- function(nbrNodes, clusterType, cl, export_vars, export_envir, verbose) {
  clusterType <- if (nbrNodes > 1) match.arg(clusterType, c("PSOCK", "FORK")) else "FORK"
  usePSOCK <- (clusterType == "PSOCK" && nbrNodes > 1)
  useFORK  <- (clusterType == "FORK"  && nbrNodes > 1 && .Platform$OS.type != "windows")

  if (clusterType == "FORK" && capabilities("X11") && dev.cur() == 1)
    warning("An X11 graphics device is open. This can cause parallel processing to hang. ...")

  cluster_created <- FALSE
  cl <- cl
  if (usePSOCK && is.null(cl)) {
    cl <- parallel::makeCluster(nbrNodes, type = "PSOCK")
    cluster_created <- TRUE
    if (verbose) message("Created new parallel cluster with ", nbrNodes, " cores.")
  }
  if (usePSOCK) {
    parallel::clusterExport(cl, export_vars, envir = export_envir)
    # Limit inner threading on workers so BLAS/OpenMP don't oversubscribe
    # when multiple PSOCK workers run concurrently.
    parallel::clusterEvalQ(cl, {
      Sys.setenv(OMP_NUM_THREADS = "1",
                 OPENBLAS_NUM_THREADS = "1",
                 MKL_NUM_THREADS = "1")
    })
  }
  list(cl = cl, cluster_created = cluster_created, usePSOCK = usePSOCK, useFORK = useFORK)
}

closeCluster <- function(cli, verbose) {
  if (cli$cluster_created) {
    parallel::stopCluster(cli$cl)
    if (verbose) message("Stopped internal cluster.")
  }
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

getEffectNamesNoRate <- function(effects, depvar) {
    include <- effects[["include"]]
    includedEffects <- effects[include, ]
    keep <- includedEffects[["type"]] != "rate" & includedEffects[["name"]] == depvar
    inc <- includedEffects[keep, ]
    sn <- numberIntShortNames(inc[["shortName"]])
    n <- length(sn)
    i1 <- if (!is.null(inc[["interaction1"]])) inc[["interaction1"]] else rep("", n)
    i2 <- if (!is.null(inc[["interaction2"]])) inc[["interaction2"]] else rep("", n)
    cs <- effectCovarSuffix(i1, i2)
    snWithCovar <- ifelse(cs == "", sn, paste(sn, cs, sep = "_"))
    # changeStats names: depvar_shortName (no type suffix), deduplicated.
    changeStats <- unique(paste(inc[["name"]], snWithCovar, sep = "_"))
    changeStats
}

# Look up the centering mean for the covariate underlying a given effect.
# Returns the centering mean (numeric) if the effect is covariate-based and
# centered, otherwise 0.  Used to translate user-supplied contrast values
# from the raw (uncentered) scale to the centered scale used internally.
getCovCenteringMean <- function(effectName, effects, data, depvar) {
  inc <- effects[effects$include & effects$type != "rate" &
                   effects$name == depvar, ]
  # Build changeStats base names aligned with inc rows.
  sn <- numberIntShortNames(inc[["shortName"]])
  n <- length(sn)
  i1 <- if (!is.null(inc[["interaction1"]])) inc[["interaction1"]] else rep("", n)
  i2 <- if (!is.null(inc[["interaction2"]])) inc[["interaction2"]] else rep("", n)
  cs <- effectCovarSuffix(i1, i2)
  snWithCovar <- ifelse(cs == "", sn, paste(sn, cs, sep = "_"))
  rowchangeStatsNames <- paste(inc[["name"]], snWithCovar, sep = "_")
  idx <- match(effectName, rowchangeStatsNames)
  if (is.na(idx)) return(0)
  covarName <- inc$interaction1[idx]
  if (is.na(covarName) || covarName == "") return(0)
  # Look up in all covariate stores
  cov <- data[["cCovars"]][[covarName]]
  if (is.null(cov)) cov <- data[["vCovars"]][[covarName]]
  if (is.null(cov)) cov <- data[["dycCovars"]][[covarName]]
  if (is.null(cov)) cov <- data[["dyvCovars"]][[covarName]]
  if (is.null(cov)) return(0)
  # coCovar/varCovar always store the raw mean; check centered flag.
  # dycCovars/dyvCovars store 0 when not centered.
  if (isTRUE(attr(cov, "centered"))) {
    m <- attr(cov, "mean")
    if (!is.null(m) && is.finite(m)) return(m)
  }
  0
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

# Extract group columns as a named list for grouped_agg_from_cols.
# Integer and numeric (double) columns pass through directly.
# Character/factor columns are integer-encoded; their decode levels are
# returned so the caller can restore original values in the C++ output.
# Returns list(cols = named list of int/double vectors,
#              decode = named list: NULL for passthrough, character levels for encoded)
extractGroupCols <- function(data, group_vars) {
  cols <- vector("list", length(group_vars))
  names(cols) <- group_vars
  decode <- vector("list", length(group_vars))
  names(decode) <- group_vars
  needs_decode <- FALSE
  for (j in seq_along(group_vars)) {
    col <- data[[group_vars[j]]]
    if (is.integer(col) || is.double(col)) {
      cols[[j]] <- col
      decode[j] <- list(NULL)
    } else {
      # Character/factor: integer-encode
      vals_num <- suppressWarnings(as.numeric(as.character(col)))
      f <- if (all(!is.na(vals_num))) {
        factor(col, levels = as.character(sort(unique(vals_num))))
      } else {
        factor(col)
      }
      cols[[j]] <- as.integer(f)
      decode[[j]] <- levels(f)
      needs_decode <- TRUE
    }
  }
  list(cols = cols, decode = decode, needs_decode = needs_decode)
}

# Restore original values in columns that were factor-encoded by extractGroupCols.
decodeResultCols <- function(res, decode) {
  for (nm in names(decode)) {
    if (!is.null(decode[[nm]])) {
      res[[nm]] <- decode[[nm]][res[[nm]]]
    }
  }
  res
}

# Encode group columns as a contiguous integer matrix for grouped_agg_cpp.
# Non-integer/non-numeric columns are factor-encoded.
# Returns list(G = integer matrix, decode = list of level vectors per column).
encodeGroupKeys <- function(data, group_vars) {
  n <- if (is.data.frame(data)) nrow(data) else length(data[[1L]])
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
  nc <- ncol(key_mat)
  out <- vector("list", nc)
  names(out) <- group_vars
  for (j in seq_len(nc)) {
    col <- key_mat[, j]
    out[[j]] <- if (!is.null(decode[[j]])) decode[[j]][col] else col
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

  enc <- encodeGroupKeys(data, pre_agg_vars)
  ord <- do.call(order, lapply(seq_len(ncol(enc$G)), function(j) enc$G[, j]))
  res <- grouped_agg_cpp(data[[outcomeName]][ord], enc$G[ord, , drop = FALSE],
                         na_rm = na.rm, do_mean = TRUE)
  out <- decodeGroupKeys(res$key, pre_agg_vars, enc$decode)
  out[[outcomeName]] <- res$value
  out
}


agg <- function(outcomeName,
                data,
                level = "none",
                condition = NULL,
                sum_fun = mean,
                na.rm = TRUE,
                egoNormalize = TRUE) {
  if (length(outcomeName) != 1L) stop("'outcomeName' must be a single column name.")
  if (!is.null(condition)) condition <- resolveEffectName(condition, names(data))
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

  # ---- Rcpp fast path for simple mean/sum ----
  is_simple_mean <- identical(sum_fun, mean) || identical(sum_fun, base::mean)
  if (is_simple_mean && length(group_vars) > 0) {
    enc <- encodeGroupKeys(data, group_vars)
    ord <- do.call(order, lapply(seq_len(ncol(enc$G)), function(j) enc$G[, j]))
    res <- grouped_agg_cpp(data[[outcomeName]][ord], enc$G[ord, , drop = FALSE],
                           na_rm = na.rm, do_mean = TRUE)
    out <- decodeGroupKeys(res$key, group_vars, enc$decode)
    out[[outcomeName]] <- res$value
    attr(out, "row.names") <- .set_row_names(length(res$value))
    class(out) <- "data.frame"
    return(out)
  }

  # ---- Base R path ----
  if (length(group_vars) == 0) {
    output <- expand_summary(sum_fun(data[[outcomeName]], na.rm = na.rm))
    output <- as.data.frame(output)
    return(output)
  }

  # General base R path (custom sum_fun or multi-column output)
  # Coerce to data.frame for split() / interaction() which need [.data.frame
  if (!is.data.frame(data)) {
    attr(data, "row.names") <- .set_row_names(length(data[[1L]]))
    class(data) <- "data.frame"
  }
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

# Convert a raw contribution matrix to change-statistic space.
#
# The raw contribMat has one column per (effect x type) combination, with
# composite names like "depvar_shortName_eval", "depvar_shortName_creation".
# In change-statistic space, creation/endow are NOT independent statistics —
# they are the eval statistic masked by direction (density).  This function
# collapses type-variant columns into one column per underlying effect
# (eval-space values) and builds direction-dependent theta.
#
# When called without theta, returns only the collapsed matrix + density
# (backward-compatible sign-flip path for display).
#
# contribMat:   N x nRaw numeric matrix (raw contributions from C++).
# effectNames:  character vector of length nRaw (composite names).
# theta:        named numeric vector (optional). When supplied, also returns
#               thetaEff (direction-dependent effective theta).
#
# Returns a list:
#   csMat     — N x nchangeStats matrix (eval-space change statistics).
#   csNames   — character(nchangeStats): changeStats names without type suffix.
#   density   — integer(N): +1 (creation), -1 (dissolution), 0 (no-change).
#   thetaEff  — nchangeStats x 2 matrix with columns "creation", "dissolution"
#               (only when theta is supplied). nochange direction is always 0.
##@contribToChangeStats Contribution-to-CS conversion (changeStats)
contribToChangeStats <- function(contribMat, effectNames, theta = NULL) {
  nRaw <- length(effectNames)
  stopifnot(ncol(contribMat) == nRaw)

  # ---- Parse composite names into (base, type) ----
  # Composite format: "depvar_shortName[_covar]_type"
  # Type is always the last underscore-separated segment, but only when it is
  # one of the known types: eval, endow, creation.  Otherwise the whole name
  # is the base (for simplified names used in tests or by the user).
  knownTypes <- c("eval", "endow", "creation")
  lastSeg <- sub("^.*_", "", effectNames)           # last segment
  isKnownType <- lastSeg %in% knownTypes
  types <- ifelse(isKnownType, lastSeg, "eval")     # treat plain names as eval
  bases <- ifelse(isKnownType, sub("_[^_]+$", "", effectNames), effectNames)

  # ---- Extract density column ----
  densityIdx <- grep("density", bases, fixed = TRUE)[1L]
  if (is.na(densityIdx)) stop("No density column found in effectNames")
  density <- as.integer(contribMat[, densityIdx])
  density[is.na(density)] <- 0L

  # ---- Dissolution-row sign flip (signed → pure change statistics) ----
  # The C++ stores signed contributions: on dissolution rows (d=-1) both
  # eval and endow contributions equal -calculateContribution(alter).
  # Negate NON-DENSITY columns on dissolution rows to recover the pure
  # (unsigned) change statistics Δs.  The density column keeps its ±1
  # sign so it carries the direction d and can be used for conditioning.
  # Utility is computed as: d * (θ_density + Δs_rest %*% θ_rest).
  neg <- density == -1L
  if (any(neg)) {
    nonDens <- setdiff(seq_len(ncol(contribMat)), densityIdx)
    contribMat[neg, nonDens] <- -contribMat[neg, nonDens]
  }

  # ---- Group columns by changeStats base name ----
  # ALL columns (including density) become changeStats effects.
  # Density stays ±1 in csMat and serves as both the direction indicator
  # and a column available for conditioning (condition = "density").

  # Unique changeStats names, preserving order of first appearance.
  allBases <- bases
  allTypes <- types
  changeStatsOrder <- unique(allBases)
  nChangeStats <- length(changeStatsOrder)

  # Build changeStats matrix: for each changeStats effect, take the eval column
  # if it exists; otherwise combine creation + endow.
  csMat <- matrix(0, nrow = nrow(contribMat), ncol = nChangeStats)
  colnames(csMat) <- changeStatsOrder

  # Also build thetaEff if theta is supplied.
  if (!is.null(theta)) {
    thetaEff <- matrix(0, nrow = nChangeStats, ncol = 2L,
                       dimnames = list(changeStatsOrder, c("creation", "dissolution")))
  }

  for (k in seq_len(nChangeStats)) {
    base_k <- changeStatsOrder[k]
    members <- which(allBases == base_k)
    memberTypes <- allTypes[members]
    memberRawIdx <- members  # indices into effectNames/contribMat

    hasEval     <- "eval"     %in% memberTypes
    hasCreation <- "creation" %in% memberTypes
    hasEndow    <- "endow"    %in% memberTypes

    # changeStats column: pure (unsigned) change statistic (Δs),
    # except density which keeps its ±1 sign.
    if (hasEval) {
      csMat[, k] <- contribMat[, memberRawIdx[memberTypes == "eval"]]
    } else if (hasCreation && hasEndow) {
      # Creation and endow are structurally zero on the "other" direction
      # rows, so summing recovers the signed Δs for both directions.
      crIdx <- memberRawIdx[memberTypes == "creation"]
      enIdx <- memberRawIdx[memberTypes == "endow"]
      csMat[, k] <- contribMat[, crIdx] + contribMat[, enIdx]
    } else if (hasCreation) {
      csMat[, k] <- contribMat[, memberRawIdx[memberTypes == "creation"]]
    } else if (hasEndow) {
      csMat[, k] <- contribMat[, memberRawIdx[memberTypes == "endow"]]
    }

    # Direction-dependent theta.
    #
    # csMat stores the pure (unsigned) change statistics Δs.  Because the
    # C++ stores c_eval = c_endow (same sign) on dissolution rows, eval
    # and endow go in the SAME direction.  The utility formula is:
    #   u = d × Δs × θ_combined
    # where d is applied separately in calculateUtility, and:
    #   creation:    θ_combined = θ_eval + θ_creation
    #   dissolution: θ_combined = θ_eval + θ_endow
    # A positive θ_endow makes dissolution utility more negative (d = -1),
    # protecting existing ties.
    if (!is.null(theta)) {
      th_eval     <- if (hasEval)     theta[effectNames[memberRawIdx[memberTypes == "eval"]]]     else 0
      th_creation <- if (hasCreation) theta[effectNames[memberRawIdx[memberTypes == "creation"]]] else 0
      th_endow    <- if (hasEndow)    theta[effectNames[memberRawIdx[memberTypes == "endow"]]]    else 0
      thetaEff[k, "creation"]    <- th_eval + th_creation
      thetaEff[k, "dissolution"] <- th_eval + th_endow
    }
  }

  # Find density's index in the changeStats matrix.
  denschangeStatsIdx <- match(bases[densityIdx], changeStatsOrder)

  # changeStatsMap: lightweight mapping for buildThetaEff (theta-independent).
  changeStatsMap <- list(
    effectNames = effectNames,
    bases       = allBases,
    types       = allTypes,
    changeStatsOrder  = changeStatsOrder
  )

  out <- list(csMat = csMat, csNames = changeStatsOrder,
              densityIdx = denschangeStatsIdx, density = density,
              changeStatsMap = changeStatsMap)
  if (!is.null(theta)) out$thetaEff <- thetaEff
  out
}

# Build direction-dependent thetaEff matrix from theta and changeStatsMap.
# O(p) — no matrix work, just theta arithmetic.
# theta: named numeric vector aligned to changeStatsMap$effectNames.
# changeStatsMap: from contribToChangeStats()$changeStatsMap.
buildThetaEff <- function(theta, changeStatsMap) {
  changeStatsOrder  <- changeStatsMap$changeStatsOrder
  nchangeStats  <- length(changeStatsOrder)
  effectNames <- changeStatsMap$effectNames
  allBases    <- changeStatsMap$bases
  allTypes    <- changeStatsMap$types

  thetaEff <- matrix(0, nrow = nchangeStats, ncol = 2L,
                     dimnames = list(changeStatsOrder, c("creation", "dissolution")))
  for (k in seq_len(nchangeStats)) {
    base_k      <- changeStatsOrder[k]
    members     <- which(allBases == base_k)
    memberTypes <- allTypes[members]
    hasEval     <- "eval"     %in% memberTypes
    hasCreation <- "creation" %in% memberTypes
    hasEndow    <- "endow"    %in% memberTypes
    th_eval     <- if (hasEval)     theta[effectNames[members[memberTypes == "eval"]]]     else 0
    th_creation <- if (hasCreation) theta[effectNames[members[memberTypes == "creation"]]] else 0
    th_endow    <- if (hasEndow)    theta[effectNames[members[memberTypes == "endow"]]]    else 0
    thetaEff[k, "creation"]    <- th_eval + th_creation
    thetaEff[k, "dissolution"] <- th_eval + th_endow
  }
  thetaEff
}

# Attach effect contribution columns to a named list (or data.frame).
# Accepts either:
#   (a) changeStats output from contribToChangeStats(): csNames + csMat
#   (b) raw effectNames + contrib matrix (legacy path)
# flip is ignored when changeStats output is provided (already in eval-space).
attachContributions <- function(out, effectNames, contrib, flip = TRUE) {
  if (flip) {
    cs <- contribToChangeStats(contrib, effectNames)
    contrib <- cs$csMat
    effectNames <- cs$csNames
  }
  nc <- ncol(contrib)
  if (nc == 0L) return(out)
  extra <- vector("list", nc)
  colNms <- if (!is.null(colnames(contrib))) colnames(contrib) else effectNames
  names(extra) <- colNms
  for (j in seq_len(nc)) extra[[j]] <- contrib[, j]
  if (is.data.frame(out)) {
    attr(extra, "row.names") <- .set_row_names(nrow(contrib))
    class(extra) <- "data.frame"
    cbind(out, extra)
  } else {
    c(out, extra)
  }
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
  # With changeStats names (no type suffix), just pass through.
  # Strip _eval/_endow/_creation suffix if user explicitly supplied one.
  sub("_(eval|endow|creation)$", "", condition)
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

  massCols <- intersect(c("massCreation", "massDissolution"), names(x))
  mass_related <- character(0)
  if (length(massCols) > 0L) {
    mass_related <- grep(
      paste0("^(", paste(massCols, collapse = "|"), ")($|_)"),
      names(x), value = TRUE
    )
  }

  main_cols <- setdiff(names(x), mass_related)
  print.data.frame(x[, main_cols, drop = FALSE], row.names = FALSE, ...)

  for (mc in massCols) {
    mc_uc_pattern <- paste0("^", mc, "_")
    mc_uc_cols <- grep(mc_uc_pattern, names(x), value = TRUE)
    level_cols <- intersect(c("period", "group"), names(x))
    mc_all <- c(level_cols, mc, mc_uc_cols)
    mc_df <- unique(x[, mc_all, drop = FALSE])
    names(mc_df) <- sub(mc_uc_pattern, "", names(mc_df))
    cat("\n  ", mc, ":\n")
    print.data.frame(mc_df, row.names = FALSE, ...)
  }

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