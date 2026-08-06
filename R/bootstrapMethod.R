# ---- Batch-shared simulation loop ----------------------------------------
#
# Runs nsim parametric-bootstrap draws. Each draw calls estimator(theta_sim)
# which returns a named list of data.frames (one per effect). Results are
# accumulated in batches and optionally written to disk for resume support.
#
# estimator: closure built by makeEstimatorFun; must carry attr("eff_names").
#
# Returns a named list of raw_sims data frames (one per effect).
drawSimBatch <- function(
    estimator,
    thetaHat,
    covTheta,
    nsim         = 1000L,
    useCluster   = FALSE,
    nbrNodes     = 1L,
    clusterType  = c("PSOCK", "FORK"),
    cl           = NULL,
    batchSize    = 100L,
    batchDir     = "temp",
    prefix       = "simBatchB_b",
    keepBatch    = FALSE,
    drawSeed     = NULL,
    verbose      = TRUE,
    gcEachBatch  = FALSE,
    gcEachSim    = FALSE
) {
    eff_names <- attr(estimator, "eff_names")
    if (is.null(eff_names))
        stop("estimator must have an 'eff_names' attribute.")
    n_effects <- length(eff_names)

    nbrNodes <- as.integer(nbrNodes)
    clusterType       <- match.arg(clusterType, c("PSOCK", "FORK"))
    useFork  <- isTRUE(useCluster) && (clusterType == "FORK") &&
                    (nbrNodes > 1L) && (.Platform$OS.type != "windows")
    usePSOCK <- isTRUE(useCluster) && (clusterType == "PSOCK") && (nbrNodes > 1L)

    # ---- Ensure PSOCK cluster exists before the batch loop ---------------
    # openCluster creates a new cluster when cl=NULL and usePSOCK; it also
    # sets OMP/BLAS/MKL thread limits on workers.  closeCluster stops it
    # only if it was created here (not if the caller supplied their own).
    if (usePSOCK) {
        cli_psock <- openCluster(nbrNodes, clusterType, cl,
            export_vars  = character(0),
            export_envir = environment(),
            verbose      = verbose
        )
        cl <- cli_psock$cl
        on.exit(closeCluster(cli_psock, verbose = verbose), add = TRUE)
    }

    # Drawn before the first estimator() call: a dynamic estimator reaches
    # siena07, which calls set.seed() and would reset the stream between draws.
    drawThetas <- function()
        MASS::mvrnorm(nsim, mu = thetaHat, Sigma = covTheta)
    thetaDraws <- if (is.null(drawSeed)) drawThetas()
                  else .withSeed(drawSeed, drawThetas())
    if (nsim == 1L)
        thetaDraws <- matrix(thetaDraws, nrow = 1L,
                             dimnames = list(NULL, names(thetaHat)))

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
        message("Batch-shared simulation: nsim=", nsim, ", effects=", n_effects,
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
        bpe <- vector("list", n_effects)
        for (j in seq_len(n_effects)) bpe[[j]] <- vector("list", bsz)

        # Run one parametric-bootstrap draw: sample theta, call estimator.
        run_one_draw <- function(k) {
            i <- batch_idx[k]
            theta_sim <- thetaDraws[i, ]
            results_named <- estimator(theta_sim)
            per_effect <- vector("list", n_effects)
            for (j in seq_len(n_effects)) {
                row_j <- results_named[[eff_names[j]]]
                row_j[["sim"]] <- i
                per_effect[[j]] <- row_j
            }
            per_effect
        }

        if (usePSOCK) {
            # Export needed objects to workers each batch
            export_names <- c("estimator", "thetaDraws",
                              "n_effects", "bsz", "batch_idx", "eff_names")
            parallel::clusterExport(
                cl, export_names, envir = environment()
            )
            draw_results <- parallel::parLapply(cl, seq_len(bsz), run_one_draw)
            for (k in seq_len(bsz)) {
                for (j in seq_len(n_effects)) bpe[[j]][[k]] <- draw_results[[k]][[j]]
            }
            sims_computed <- sims_computed + bsz
        } else if (useFork) {
            # Pre-fork diagnostics: warn if any live chains survive in heap
            if (batches_done == 0L) {
                assertChainsFreed(verbose = verbose)
                memReport("pre-fork (uncertainty)", verbose = verbose)
            }
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
                for (j in seq_len(n_effects)) bpe[[j]][[k]] <- draw_results[[k]][[j]]
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
                for (j in seq_len(n_effects)) bpe[[j]][[k]] <- per_effect[[j]]
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
    per_effect_raw <- vector("list", n_effects)
    for (j in seq_len(n_effects)) per_effect_raw[[j]] <- vector("list", nsim)

    if (!is.null(in_memory_batch)) {
        # Single-batch path: no disk I/O
        for (j in seq_len(n_effects)) per_effect_raw[[j]] <- in_memory_batch[[j]]
    } else {
        batch_files <- sort(list.files(
            batchDir,
            pattern    = batch_pattern,
            full.names = TRUE
        ))
        # The resume feature keys on filename alone, so a stale batch from an
        # interrupted run -- or from a concurrent run sharing batchDir and
        # prefix -- looks indistinguishable from one of ours.  Assigning it
        # positionally used to recycle silently (a warning at most) and
        # quietly corrupt the draws.  Fail loudly instead.
        if (length(batch_files) != nbatches)
            stop(sprintf(paste0(
                "Found %d batch file(s) in '%s' but this run has %d batch(es). ",
                "Stale or foreign batches from another run would corrupt the ",
                "draws. Clear the directory, or give this run its own ",
                "batchDir/prefix."),
                length(batch_files), batchDir, nbatches), call. = FALSE)
        for (b in seq_along(batch_files)) {
            bp    <- readRDS(batch_files[[b]])
            first <- 1L + (b - 1L) * batchSize
            last  <- min(b * batchSize, nsim)
            expected <- last - first + 1L
            if (length(bp) != n_effects)
                stop(sprintf(paste0(
                    "Batch file '%s' holds %d effect(s); this run expects %d. ",
                    "It belongs to a different run -- clear '%s'."),
                    basename(batch_files[[b]]), length(bp), n_effects,
                    batchDir), call. = FALSE)
            for (j in seq_len(n_effects)) {
                if (length(bp[[j]]) != expected)
                    stop(sprintf(paste0(
                        "Batch file '%s' holds %d draw(s) for effect %d; ",
                        "this run expects %d. It belongs to a different run ",
                        "-- clear '%s'."),
                        basename(batch_files[[b]]), length(bp[[j]]), j,
                        expected, batchDir), call. = FALSE)
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


makeUncertaintySummarizer <- function(
    return_sd = TRUE,
    return_ci = TRUE,
    return_mean = FALSE,
    return_median = FALSE,
    ci_interval = c(0.025, 0.975)
) {
  ci_interval <- as.numeric(ci_interval)
  if (length(ci_interval) != 2L || ci_interval[1] >= ci_interval[2]) {
    stop("'ciInterval' must be a length-2 increasing vector, e.g. c(0.025, 0.975).")
  }
  if (any(!is.finite(ci_interval)) || any(ci_interval <= 0) || any(ci_interval >= 1)) {
    stop("'ciInterval' must be strictly between 0 and 1.")
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
      # median (0.5) is derived, not user-specified; CI bounds come from ci_interval
      qprobs <- if (isTRUE(return_median)) {
        c(ci_interval[1], 0.5, ci_interval[2])
      } else {
        ci_interval
      }
      if (n >= 2L) {
        qu <- unname(stats::quantile(x, probs = qprobs, na.rm = FALSE))
      } else {
        qu <- rep(NA_real_, length(qprobs))
      }
      if (isTRUE(return_median)) out[["Median"]] <- qu[2]
      if (isTRUE(return_ci)) {
        # column names derive from the interval bounds, e.g. q_025 / q_975
        out[[sprintf("q_%03d", round(ci_interval[1] * 1000))]] <- qu[1]
        out[[sprintf("q_%03d", round(ci_interval[2] * 1000))]] <- qu[length(qu)]
      }
    }

    out
  }
}



