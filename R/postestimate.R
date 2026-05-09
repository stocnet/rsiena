# --------------------------------------------------------------------------
# sienaPostestimate
#
# Callers supply:
#   contribFun    – hat contribFun (built by caller: static/preloaded/per_batch)
#   nChainBatches – number of chain-batches consumed by contribFun for hat
#   type, rateParams – forwarded to makeEstimatorFun
#   specs         – named list of per-effect spec entries
#   thetaHat, covTheta – point estimate + covariance
#
# Two regimes for hat vs uncertainty:
#   Static:  contribFun is theta-independent; same closure used for both.
#            No n3 involved; batching irrelevant.
#   Dynamic: hat uses caller-supplied contribFun (preloaded/per_batch)
#            with n3PointEst total chains batched in n3BatchSize chunks.
#            Uncertainty builds a fresh per_batch contribFun from dynArgs
#            with n3 total chains per theta draw (one siena07 per draw).
#            Chains depend on theta — hat chains cannot be reused here.
# --------------------------------------------------------------------------

sienaPostestimate <- function(
    contribFun,
    nChainBatches,
    type,
    rateParams = NULL,
    rateIdx    = NULL,
    specs,
    thetaHat,
    covTheta,
    dynamic  = FALSE,
    dynArgs  = NULL,
    n3       = NULL,
    n3BatchSize = NULL,
    useChangeContributions = FALSE,
    uncertainty = TRUE,
    uncertaintyMode = c("bootstrap", "delta", "deltaFull"),
    nsim        = 1000L,
    batchSize   = 100L,
    useCluster  = FALSE,
    nbrNodes    = 1L,
    clusterType = c("PSOCK", "FORK"),
    cl          = NULL,
    batchDir    = "temp",
    prefix      = "simBatch_b",
    keepBatch   = FALSE,
    verbose     = TRUE,
    na.rm        = TRUE,
    egoNormalize = TRUE,
    uncertaintySd         = TRUE,
    uncertaintyCi         = TRUE,
    uncertaintyMean       = FALSE,
    uncertaintyMedian     = FALSE,
    uncertaintyProbs      = c(0.025, 0.5, 0.975),
    uncertaintyMcse       = FALSE,
    uncertaintymcseBatches = NULL,
    decisionDetails = NULL,
    saveDir         = NULL,
    gcEachBatch     = FALSE,
    gcEachSim       = FALSE,
    preloadedChains = NULL
) {
    clusterType     <- match.arg(clusterType)
    uncertaintyMode <- match.arg(uncertaintyMode)

    isDeltaMode <- uncertaintyMode %in% c("delta", "deltaFull")
    isFullMode  <- uncertaintyMode == "deltaFull"

    # ---- Route ----
    # Three paths share point estimates (expects) but diverge on what
    # chains/estimators they need.  Each path is self-contained below.
    if (!uncertainty)
        return(.runHatOnly(
            contribFun             = contribFun,
            nChainBatches          = nChainBatches,
            specs                  = specs,
            thetaHat               = thetaHat,
            type                   = type,
            rateParams             = rateParams,
            rateIdx                = rateIdx,
            verbose                = verbose,
            nbrNodes               = nbrNodes,
            useChangeContributions = useChangeContributions,
            decisionDetails        = decisionDetails,
            saveDir                = saveDir
        ))

    if (isDeltaMode)
        return(.runDeltaPath(
            contribFun      = contribFun,
            nChainBatches   = nChainBatches,
            specs           = specs,
            thetaHat        = thetaHat,
            covTheta        = covTheta,
            type            = type,
            rateParams      = rateParams,
            rateIdx         = rateIdx,
            verbose         = verbose,
            nbrNodes        = nbrNodes,
            isFullMode      = isFullMode,
            dynamic         = dynamic,
            dynArgs         = dynArgs,
            preloadedChains = preloadedChains,
            n3              = n3,
            n3BatchSize     = n3BatchSize,
            useChangeContributions = useChangeContributions,
            decisionDetails = decisionDetails,
            saveDir         = saveDir
        ))

    .runBootstrapPath(
        contribFun             = contribFun,
        nChainBatches          = nChainBatches,
        specs                  = specs,
        thetaHat               = thetaHat,
        covTheta               = covTheta,
        type                   = type,
        rateParams             = rateParams,
        rateIdx                = rateIdx,
        verbose                = verbose,
        nbrNodes               = nbrNodes,
        dynamic                = dynamic,
        dynArgs                = dynArgs,
        n3                     = n3,
        n3BatchSize            = n3BatchSize,
        useChangeContributions = useChangeContributions,
        nsim                   = nsim,
        batchSize              = batchSize,
        useCluster             = useCluster,
        clusterType            = clusterType,
        cl                     = cl,
        batchDir               = batchDir,
        prefix                 = prefix,
        keepBatch              = keepBatch,
        na.rm                  = na.rm,
        egoNormalize           = egoNormalize,
        uncertaintySd          = uncertaintySd,
        uncertaintyCi          = uncertaintyCi,
        uncertaintyMean        = uncertaintyMean,
        uncertaintyMedian      = uncertaintyMedian,
        uncertaintyProbs       = uncertaintyProbs,
        uncertaintyMcse        = uncertaintyMcse,
        uncertaintymcseBatches = uncertaintymcseBatches,
        decisionDetails        = decisionDetails,
        saveDir                = saveDir,
        gcEachBatch            = gcEachBatch,
        gcEachSim              = gcEachSim
    )
}

# --------------------------------------------------------------------------
# .runHatOnly — point-estimate-only path (uncertainty = FALSE)
#
# Builds hatEstimatorFun, computes expects, attaches attrs, saves, returns.
# --------------------------------------------------------------------------
.runHatOnly <- function(contribFun, nChainBatches, specs, thetaHat,
                        type, rateParams, rateIdx, verbose, nbrNodes,
                        useChangeContributions, decisionDetails, saveDir) {
    hatEstimatorFun <- makeEstimatorFun(specs, contribFun, nChainBatches,
                                        type, rateParams, rateIdx, verbose,
                                        mc.cores = nbrNodes)
    expects <- hatEstimatorFun(thetaHat,
                               useChangeContributions = useChangeContributions)
    rm(hatEstimatorFun, contribFun)
    gc(verbose = FALSE)

    results <- lapply(names(specs), function(nm) {
        res <- attachPostestAttrs(expects[[nm]], specs[[nm]]$metadata)
        if (!is.null(decisionDetails) && !is.null(decisionDetails[[nm]]))
            attr(res, "decisionDetails") <- decisionDetails[[nm]]
        res
    })
    names(results) <- names(specs)
    if (!is.null(saveDir))
        for (nm in names(results))
            saveRDS(results[[nm]], file.path(saveDir, paste0(nm, ".rds")))
    results
}

# --------------------------------------------------------------------------
# .runDeltaPath — full delta / deltaFull path
#
# 1. Optionally simulates chains once (dynamic case) and wraps in a memory
#    store so that hatEstimatorFun (and the FD gradient) share the same chains.
# 2. Computes point estimates via hatEstimatorFun.
# 3. Calls deltaMethodUncertainty using the same estimator for the gradient.
# 4. Assembles per-effect result data.frames with delta_se / delta_q025/975.
# --------------------------------------------------------------------------
.runDeltaPath <- function(contribFun, nChainBatches, specs, thetaHat, covTheta,
                           type, rateParams, rateIdx, verbose, nbrNodes,
                           isFullMode, dynamic, dynArgs, preloadedChains,
                           n3, n3BatchSize, useChangeContributions,
                           decisionDetails, saveDir) {
    # Step 1 — replace contribFun with a memory-backed store (dynamic only).
    # Static delta: contribFun is already theta-independent; no chain sim needed.
    delta_wide <- NULL
    ssc_sum    <- NULL
    if (dynamic && !is.null(dynArgs)) {
        deltaState    <- .initDeltaMode(
            isFullMode      = isFullMode,
            dynamic         = dynamic,
            dynArgs         = dynArgs,
            preloadedChains = preloadedChains,
            contribFun      = contribFun,
            nChainBatches   = nChainBatches,
            thetaHat        = thetaHat,
            n3              = n3,
            n3BatchSize     = n3BatchSize
        )
        contribFun    <- deltaState$contribFun
        nChainBatches <- deltaState$nChainBatches
        delta_wide    <- deltaState$delta_wide
        ssc_sum       <- deltaState$ssc_sum
    }

    # Step 2 — point estimates.
    hatEstimatorFun <- makeEstimatorFun(specs, contribFun, nChainBatches,
                                        type, rateParams, rateIdx, verbose,
                                        mc.cores = nbrNodes)

    # Step 3+4 — hat + Jacobian in one sweep, then SE / CI assembly.
    # For accumulated or rate-weighted specs the analytical path is
    # unavailable/unsafe; fall back to FD inside deltaMethodUncertainty
    # (hat will be computed there too).
    has_accumulated <- any(vapply(specs, function(s) isTRUE(s$accumulated),
                                  logical(1L)))
    has_rateweight <- any(vapply(specs, function(s) isTRUE(s$rateWeight),
                   logical(1L)))
    if (isFullMode && (!dynamic || is.null(ssc_sum)))
        stop("uncertaintyMode = 'deltaFull' requires dynamic = TRUE and ",
             "score statistics (includeScores = TRUE in dynArgs).")

    if (!has_accumulated && !has_rateweight) {
        # Analytical path: one estimator call returns both hat and Jacobian.
        # Suppress warnings here — if it fails, deltaMethodUncertainty will
        # detect the failure itself and emit the appropriate warning.
        precomputed <- tryCatch(
            hatEstimatorFun(thetaHat, mode = "jacobian"),
            error = function(e) NULL)
    } else {
        precomputed <- NULL
    }
    expects <- if (!is.null(precomputed)) precomputed$hat
               else hatEstimatorFun(thetaHat,
                                    useChangeContributions = useChangeContributions)
    memReport("post-hat + Jacobian", verbose = verbose >= 2)

    delta_uncert <- deltaMethodUncertainty(
        wide        = delta_wide,
        estimator   = hatEstimatorFun,
        ssc_sum     = ssc_sum,
        thetaHat    = thetaHat,
        covTheta    = covTheta,
        specs       = specs,
        type        = type,
        fullMode    = isFullMode,
        precomputed = precomputed
    )
    rm(hatEstimatorFun, contribFun)
    gc(verbose = FALSE)

    results <- vector("list", length(specs))
    names(results) <- names(specs)
    for (nm in names(specs)) {
        spec  <- specs[[nm]]
        df    <- expects[[nm]]
        du    <- delta_uncert[[nm]]
        df[["delta_se"]] <- du$SE_delta
        if (isFullMode) df[["delta_full_se"]] <- as.numeric(du$SE_deltaFull)
        se_use <- as.numeric(if (isFullMode) du$SE_deltaFull else du$SE_delta)
        oc     <- spec$outcomeName
        if (!is.null(df[[oc]]) && length(se_use) %in% c(1L, nrow(df))) {
            df[["delta_q025"]] <- df[[oc]] - qnorm(0.975) * se_use
            df[["delta_q975"]] <- df[[oc]] + qnorm(0.975) * se_use
        }
        res <- attachPostestAttrs(df, spec$metadata)
        attr(res, "delta_jacobians") <- du[c("J_cond", "J_full",
                                              "baseline", "ssc_colMeans")]
        if (isFullMode && isTRUE(attr(du$SE_deltaFull, "fallback")))
            attr(res, "delta_full_fallback") <- TRUE
        if (!is.null(decisionDetails) && !is.null(decisionDetails[[nm]]))
            attr(res, "decisionDetails") <- decisionDetails[[nm]]
        results[[nm]] <- res
    }
    if (!is.null(saveDir))
        for (nm in names(results))
            saveRDS(results[[nm]], file.path(saveDir, paste0(nm, ".rds")))
    results
}

# --------------------------------------------------------------------------
# .runBootstrapPath — parametric-bootstrap path
#
# 1. Computes point estimates (hat path, using caller-supplied contribFun).
# 2. Builds a separate estimatorFun for draws (dynamic: fresh chains per draw;
#    static: same contribFun as hat).
# 3. Draws nsim bootstrap samples via drawSimBatch.
# 4. Aggregates via aggregatePostEstimation.
# --------------------------------------------------------------------------
.runBootstrapPath <- function(contribFun, nChainBatches, specs, thetaHat,
                               covTheta, type, rateParams, rateIdx, verbose,
                               nbrNodes, dynamic, dynArgs, n3, n3BatchSize,
                               useChangeContributions, nsim, batchSize,
                               useCluster, clusterType, cl, batchDir, prefix,
                               keepBatch, na.rm, egoNormalize,
                               uncertaintySd, uncertaintyCi, uncertaintyMean,
                               uncertaintyMedian, uncertaintyProbs,
                               uncertaintyMcse, uncertaintymcseBatches,
                               decisionDetails, saveDir, gcEachBatch, gcEachSim) {
    # Step 1 — point estimates from the caller-supplied (hat) contribFun.
    hatEstimatorFun <- makeEstimatorFun(specs, contribFun, nChainBatches,
                                        type, rateParams, rateIdx, verbose,
                                        mc.cores = nbrNodes)
    expects <- hatEstimatorFun(thetaHat,
                               useChangeContributions = useChangeContributions)
    rm(hatEstimatorFun)
    gc(verbose = FALSE)
    memReport("post-hat (chains freed)", verbose = verbose >= 2)

    # Step 2 — build draw estimator.
    # Dynamic: each bootstrap draw needs its own fresh chains (theta-dependent),
    #   so build a chainStore_simulate that re-runs siena07 per draw.
    # Static: contribFun is theta-independent; re-use the same closure as hat.
    if (dynamic && !is.null(dynArgs)) {
        uncertBatchN3 <- if (!is.null(n3BatchSize)) min(n3BatchSize, n3) else n3
        uncertDynArgs <- dynArgs
        uncertDynArgs$n3 <- uncertBatchN3
        uncertStore <- chainStore_simulate(uncertDynArgs, uncertBatchN3, n3)
        uncertContribFun <- makeContribFun(store   = uncertStore,
                                           effects = dynArgs$effects,
                                           depvar  = dynArgs$depvar,
                                           data    = dynArgs$data)
        estimatorFun <- makeEstimatorFun(specs, uncertContribFun,
                                         uncertStore$nBatches,
                                         type, rateParams, rateIdx,
                                         verbose = verbose)
    } else {
        # Static: contribFun is theta-independent; same closure is correct.
        estimatorFun <- makeEstimatorFun(specs, contribFun, nChainBatches,
                                         type, rateParams, rateIdx,
                                         verbose = verbose)
    }
    rm(contribFun)

    # Step 3 — bootstrap draws.
    uncertainty_summary_fun <- makeUncertaintySummarizer(
        return_sd     = uncertaintySd,
        return_ci     = uncertaintyCi,
        return_mean   = uncertaintyMean,
        return_median = uncertaintyMedian,
        probs         = uncertaintyProbs,
        return_mcse   = uncertaintyMcse,
        mcseBatches   = uncertaintymcseBatches
    )
    keepBatch_internal <- if (!is.null(saveDir)) TRUE else keepBatch
    raw_sims_list <- drawSimBatch(
        estimator   = estimatorFun,
        thetaHat    = thetaHat,
        covTheta    = covTheta,
        nsim        = nsim,
        batchSize   = batchSize,
        useCluster  = useCluster,
        nbrNodes    = nbrNodes,
        clusterType = clusterType,
        cl          = cl,
        batchDir    = batchDir,
        prefix      = prefix,
        keepBatch   = keepBatch_internal,
        verbose     = verbose,
        gcEachBatch = gcEachBatch,
        gcEachSim   = gcEachSim
    )

    # Step 4 — aggregate.
    t_agg_start <- proc.time()
    results <- aggregatePostEstimation(
        expects                 = expects,
        raw_sims_list           = raw_sims_list,
        specs                   = specs,
        uncertainty_summary_fun = uncertainty_summary_fun,
        na.rm                   = na.rm,
        egoNormalize            = egoNormalize,
        decisionDetails         = decisionDetails,
        saveDir                 = saveDir
    )
    if (verbose >= 2)
        message(sprintf("Aggregation complete: %.2fs",
                        (proc.time() - t_agg_start)[["elapsed"]]))
    if (!is.null(saveDir) && !keepBatch) {
        batch_pattern <- sprintf("^%s\\d{3}\\.rds$", prefix)
        bf <- list.files(batchDir, pattern = batch_pattern, full.names = TRUE)
        for (f in bf) file.remove(f)
    }
    results
}

# --------------------------------------------------------------------------
# chainStore — backend-agnostic batched access to change contributions.
#
# A chainStore is a list (class "chainStore") with:
#   $mode      – character label ("memory", "disk", "simulate")
#   $nBatches  – integer, number of batches
#   $nTotal    – integer, total chain iterations (NA for simulate)
#   $getBatch  – function(batchIdx, theta = NULL) → raw chain list (or
#                for simulate mode, a wide contrib struct)
#   $cleanup   – function() → frees temp resources
#
# The "simulate" backend is special: getBatch returns an already-flattened
# wide struct (from getDynamicChangeContributions), whereas "memory" and
# "disk" return raw chain lists that makeContribFun passes through
# flattenAndEnrichWide.
#
# Future backends (Arrow IPC, DuckDB, Rcpp XPtr) implement the same
# interface — see docs/design/chain_storage_backends.md.
# --------------------------------------------------------------------------

chainStore_memory <- function(chains, batchSize) {
    nTotal   <- length(chains)
    nBatches <- ceiling(nTotal / batchSize)
    structure(list(
        mode     = "memory",
        nBatches = nBatches,
        nTotal   = nTotal,
        getBatch = function(batchIdx, theta = NULL) {
            bsz      <- ceiling(nTotal / nBatches)
            startIdx <- (batchIdx - 1L) * bsz + 1L
            endIdx   <- min(batchIdx * bsz, nTotal)
            chains[startIdx:endIdx]
        },
        cleanup  = function() invisible(NULL)
    ), class = "chainStore")
}

chainStore_disk <- function(chains, batchSize, dir = tempdir(),
                            compress = FALSE, depvar = NULL, verbose = FALSE) {
    nTotal   <- length(chains)
    nBatches <- ceiling(nTotal / batchSize)
    files    <- character(nBatches)
    for (b in seq_len(nBatches)) {
        bsz      <- ceiling(nTotal / nBatches)
        startIdx <- (b - 1L) * bsz + 1L
        endIdx   <- min(b * bsz, nTotal)
        slice    <- chains[startIdx:endIdx]
        # Filter by depvar at write time so disk files contain only
        # the relevant network's ministeps.
        if (!is.null(depvar)) {
            slice <- lapply(slice, function(chain_nit)
                lapply(chain_nit, function(chain_group)
                    lapply(chain_group, function(chainPeriod)
                        Filter(function(ms) {
                            nm <- attr(ms, "networkName")
                            is.null(nm) || as.character(nm) %in% depvar
                        }, chainPeriod))))
        }
        files[b] <- tempfile(sprintf("chains_b%03d_", b), tmpdir = dir,
                             fileext = ".rds")
        saveRDS(slice, files[b], compress = compress)
    }
    if (verbose)
        message(sprintf("chainStore_disk: wrote %d batch files (%.1f MB total)",
                        nBatches, sum(file.size(files)) / 1e6))
    structure(list(
        mode     = "disk",
        nBatches = nBatches,
        nTotal   = nTotal,
        getBatch = function(batchIdx, theta = NULL) {
            readRDS(files[[batchIdx]])
        },
        cleanup  = function() {
            unlink(files[file.exists(files)])
            invisible(NULL)
        }
    ), class = "chainStore")
}

chainStore_simulate <- function(dynArgs, batchSize, n3Total) {
    nBatches <- ceiling(n3Total / batchSize)
    batchDynArgs <- dynArgs
    structure(list(
        mode     = "simulate",
        nBatches = nBatches,
        nTotal   = n3Total,
        getBatch = function(batchIdx, theta = NULL) {
            # Compute the actual n3 for this batch.  The last batch may be
            # smaller than batchSize when n3Total is not a multiple of batchSize.
            actual_n3 <- min(batchSize, n3Total - (batchIdx - 1L) * batchSize)
            this_args <- batchDynArgs
            this_args$n3 <- actual_n3
            this_args$useChangeContributions <- FALSE
            do.call(getDynamicChangeContributions,
                    c(list(theta = theta), this_args))
        },
        cleanup  = function() invisible(NULL)
    ), class = "chainStore")
}

# --------------------------------------------------------------------------
# .detectVmaxGB — detect R's effective memory ceiling in GB.
#
# Uses R_MAX_VSIZE if set; falls back to physical RAM via sysctl (macOS) or
# /proc/meminfo (Linux); final fallback is 16 GB conservative.
# --------------------------------------------------------------------------
.detectVmaxGB <- function() {
    tryCatch({
        vmax <- Sys.getenv("R_MAX_VSIZE", unset = "")
        if (nzchar(vmax)) {
            if (grepl("[gG][bB]?$", vmax)) {
                as.numeric(sub("[gG][bB]?$", "", vmax))
            } else {
                as.numeric(vmax) / 1024^3
            }
        } else {
            physGB <- tryCatch({
                if (Sys.info()[["sysname"]] == "Darwin") {
                    hw <- system("sysctl -n hw.memsize", intern = TRUE)
                    as.numeric(hw) / 1024^3
                } else if (file.exists("/proc/meminfo")) {
                    ln <- grep("^MemTotal", readLines("/proc/meminfo"),
                               value = TRUE)[1L]
                    as.numeric(sub("[^0-9]+", "", ln)) / 1024^2
                } else {
                    NA_real_
                }
            }, error = function(e) NA_real_)
            if (!is.na(physGB) && physGB > 0) physGB else 16
        }
    }, error = function(e) 16)
}

# --------------------------------------------------------------------------
# .checkDynMemory — memory safety guard for dynamic chain paths.
#
# Checks peak memory estimate against vmaxGB (stop/warning on breach), then
# applies the FORK-multiplier guard, potentially auto-reducing nbrNodes.
#
# Arguments:
#   data, depvar, effects — passed to estimateDynMemory
#   n3_per_batch  — chains per hat batch (point-estimate path)
#   n3_uncert     — chains per uncertainty-draw batch
#   useCluster, nbrNodes, clusterType — parallelism settings
#   uncertainty   — whether uncertainty draws are requested
#   verbose       — emit memory estimate messages
#
# Returns: list(nbrNodes = <adjusted>, vmaxGB = <detected>)
# --------------------------------------------------------------------------
.checkDynMemory <- function(data, depvar, effects,
                             n3_per_batch, n3_uncert,
                             useCluster = FALSE, nbrNodes = 1L,
                             clusterType = "PSOCK", uncertainty = TRUE,
                             verbose = FALSE) {
    vmaxGB <- .detectVmaxGB()
    memEst <- estimateDynMemory(data, depvar, effects, n3_per_batch)

    if (memEst$estGB_peak > vmaxGB) {
        stop(sprintf(
            paste0("Estimated peak memory (%.1f GB) for dynamic forward ",
            "simulation exceeds R_MAX_VSIZE (%.0f GB).\n",
            "  Network: %d actors, %d periods, %d effects, ",
            "mean rate %.1f, n3 = %d -> ~%.0fM rows.\n",
            "Reduce n3/n3BatchSize, increase R_MAX_VSIZE, or reduce chains."),
            memEst$estGB_peak, vmaxGB,
            memEst$nActor, memEst$nPer, memEst$nEff,
            memEst$meanRate, n3_per_batch, memEst$estRows / 1e6),
            call. = FALSE)
    } else if (memEst$estGB_peak > vmaxGB * 0.6) {
        warning(sprintf(
            paste0("Estimated peak memory %.1f GB (%.0f%% of R_MAX_VSIZE). ",
            "Network: %d actors, %d periods, mean rate %.1f, n3 = %d.\n",
            "Memory pressure likely; consider reducing n3/n3BatchSize."),
            memEst$estGB_peak,
            100 * memEst$estGB_peak / vmaxGB,
            memEst$nActor, memEst$nPer,
            memEst$meanRate, n3_per_batch),
            call. = FALSE, immediate. = TRUE)
    }

    ct      <- match.arg(clusterType, c("PSOCK", "FORK"))
    useFork <- isTRUE(useCluster) && ct == "FORK" && nbrNodes > 1L
    if (useFork && uncertainty) {
        uncertMem   <- estimateDynMemory(data, depvar, effects, n3_uncert)
        parentGB    <- memEst$estGB_peak + 2
        forkTotalGB <- uncertMem$estGB_peak * nbrNodes + parentGB
        forkSafeFrac <- 0.80
        if (forkTotalGB > vmaxGB * forkSafeFrac && nbrNodes > 1L) {
            safeGB    <- vmaxGB * forkSafeFrac - parentGB
            safeNodes <- max(1L, as.integer(floor(safeGB / uncertMem$estGB_peak)))
            warning(sprintf(
                paste0("FORK memory safety: estimated %.1f GB ",
                "(%d workers \u00d7 %.1f GB + %.1f GB parent) ",
                "exceeds %.0f%% of %.0f GB available.\n",
                "Auto-reducing nbrNodes from %d to %d to prevent OOM."),
                forkTotalGB, nbrNodes, uncertMem$estGB_peak,
                parentGB, 100 * forkSafeFrac, vmaxGB,
                nbrNodes, safeNodes),
                call. = FALSE, immediate. = TRUE)
            nbrNodes    <- safeNodes
            forkTotalGB <- uncertMem$estGB_peak * nbrNodes + parentGB
        }
        if (forkTotalGB > vmaxGB * 0.6) {
            warning(sprintf(
                paste0("FORK memory pressure: estimated %.1f GB ",
                "(%d workers \u00d7 %.1f GB + parent) is %.0f%% of %.0f GB.\n",
                "Consider reducing nbrNodes or n3."),
                forkTotalGB, nbrNodes, uncertMem$estGB_peak,
                100 * forkTotalGB / vmaxGB, vmaxGB),
                call. = FALSE, immediate. = TRUE)
        }
        if (verbose) {
            message(sprintf(
                "Fork memory estimate: %d workers \u00d7 %.1f GB + %.1f GB parent = %.1f GB total (of %.0f GB)",
                nbrNodes, uncertMem$estGB_peak, parentGB, forkTotalGB, vmaxGB))
        }
    }

    if (verbose) {
        message(sprintf(
            "Dynamic memory estimate: ~%.0fM rows, %.1f GB contrib, %.1f GB peak (of %.0f GB limit)",
            memEst$estRows / 1e6, memEst$estGB_contrib,
            memEst$estGB_peak, vmaxGB))
    }

    list(nbrNodes = nbrNodes, vmaxGB = vmaxGB)
}

# --------------------------------------------------------------------------
# .buildDynChainStore — shared dynamic chain store builder.
#
# Selects and constructs the appropriate chainStore backend, then wraps it
# in a contribFun via makeContribFun.
#
#   chains != NULL             → disk backend (caller has resolved chains)
#   chains == NULL, n3Batch < n3Total → batched simulate
#   chains == NULL, n3Batch == n3Total → all-at-once simulate (1 batch)
#
# Caller is responsible for rm(chains); gc() after this call to free RAM
# (chainStore_disk writes chains to disk; the original list is then redundant).
#
# Returns: list(store, contribFun, nChainBatches, mode)
#
# chainStoreMode: "auto" (default) picks memory when peak fits in RAM,
#   "disk" always serialises to temp files, "memory" always keeps in RAM.
#   Auto uses estimateDynMemory with a conservative 5× raw-chain multiplier.
# --------------------------------------------------------------------------
.buildDynChainStore <- function(chains, dynArgs, n3Total, n3Batch,
                                 depvar, effects, data, verbose,
                                 chainStoreMode = "auto",
                                 vmaxGB = NULL, nbrNodes = 1L) {
    chainStoreMode <- match.arg(chainStoreMode, c("auto", "disk", "memory"))
    nChainBatches <- ceiling(n3Total / n3Batch)
    if (!is.null(chains)) {
        useMemory <- if (chainStoreMode == "memory") {
            TRUE
        } else if (chainStoreMode == "disk") {
            FALSE
        } else { # "auto"
            if (!is.null(vmaxGB) && vmaxGB > 0) {
                memEst <- estimateDynMemory(data, depvar, effects, n3Batch)
                # raw chains in R are ~5x the flattened-peak estimate
                estTotalGB <- memEst$estGB_peak * 5 * max(1L, as.integer(nbrNodes))
                estTotalGB < vmaxGB * 0.7
            } else {
                FALSE  # no memory info → safe default is disk
            }
        }
        if (useMemory) {
            store      <- chainStore_memory(chains, n3Batch)
            contribFun <- makeContribFun(store = store, effects = effects,
                                          depvar = depvar, data = data)
            list(store = store, contribFun = contribFun,
                 nChainBatches = nChainBatches, mode = "memory")
        } else {
            store      <- chainStore_disk(chains, n3Batch,
                                          depvar = depvar, verbose = verbose)
            contribFun <- makeContribFun(store = store, effects = effects,
                                          depvar = depvar, data = data)
            list(store = store, contribFun = contribFun,
                 nChainBatches = nChainBatches, mode = "disk")
        }
    } else if (n3Batch < n3Total) {
        bArgs      <- dynArgs; bArgs$n3 <- n3Batch
        store      <- chainStore_simulate(bArgs, n3Batch, n3Total)
        contribFun <- makeContribFun(store = store, effects = effects,
                                      depvar = depvar, data = data)
        list(store = store, contribFun = contribFun,
             nChainBatches = nChainBatches, mode = "simulate")
    } else {
        aArgs      <- dynArgs; aArgs$n3 <- n3Total
        store      <- chainStore_simulate(aArgs, n3Total, n3Total)
        contribFun <- makeContribFun(store = store, effects = effects,
                                      depvar = depvar, data = data)
        list(store = store, contribFun = contribFun,
             nChainBatches = 1L, mode = "all_at_once")
    }
}

# --------------------------------------------------------------------------
# makeContribFun
#
# Wraps a chainStore into a closure compatible with makeEstimatorFun.
# --------------------------------------------------------------------------
# .resolveBuiltinJac — map a built-in predictFun to its paired Jacobian fn.
#
# Returns the corresponding *Jac function for the two built-in predictFuns,
# or NULL for any other function (signals FD fallback in evalBatchJacobian).
# Callers can override per-spec by setting spec$jacobianFun directly.
# --------------------------------------------------------------------------
.resolveBuiltinJac <- function(predictFun) {
  if (identical(predictFun, predictProbability)) return(predictProbabilityJac)
  if (identical(predictFun, predictFirstDiff))   return(predictFirstDiffJac)
  NULL
}

# --------------------------------------------------------------------------
# makeEstimatorFun
#
# Unified replacement for makePartialEstimator + makeBatchedMultiEstimator.
# Takes a named list of specs (one per effect) and a contribFun, returns a
# closure function(theta) → named list of data.frames.
#
# Each spec must have:
#   predictFun(changeContributions, theta, ...) → data.frame
#   predictArgs:   additional named args forwarded to predictFun
#   outcomeName:   column in predictFun output to aggregate
#   level:         aggregation level ("none", "period", "ego", ...)
#   condition:     optional conditioning column name
#   accumulated:   logical — use aggAccumulatedSumCount instead of aggSumCount
#   na.rm:         logical
#   egoNormalize:  logical
#   rateWeight:    logical — multiply by period rate (static path only)
#   dynamic:       logical — TRUE when chains come from siena07
# --------------------------------------------------------------------------
makeEstimatorFun <- function(specs, contribFun, nBatches,
    type, rateParams = NULL, rateIdx = NULL, verbose = FALSE, mc.cores = 1L) {
  force(specs); force(contribFun); force(nBatches)
  force(type); force(rateParams); force(rateIdx); force(verbose); force(mc.cores)
  N         <- length(specs)
  eff_names <- names(specs)
  sep       <- "\x1f"

  non_acc <- which(!vapply(specs, function(s) isTRUE(s$accumulated),
                           logical(1L)))
  pred_group_id <- integer(N)
  {
    next_gid <- 1L
    for (j in seq_len(N)) {
      gid <- 0L
      for (k in seq_len(j - 1L)) {
        if (identical(specs[[j]]$predictFun,    specs[[k]]$predictFun) &&
            identical(specs[[j]]$predictArgs,   specs[[k]]$predictArgs) &&
            identical(specs[[j]]$massContrasts, specs[[k]]$massContrasts)) {
          gid <- pred_group_id[k]
          break
        }
      }
      if (gid == 0L) { gid <- next_gid; next_gid <- next_gid + 1L }
      pred_group_id[j] <- gid
    }
    n_pred_groups <- next_gid - 1L
  }

  # Bundle spec-derived constants for the top-level helpers defined below.
  evalConfig <- list(
    specs         = specs,
    N             = N,
    eff_names     = eff_names,
    non_acc       = non_acc,
    pred_group_id = pred_group_id,
    n_pred_groups = n_pred_groups,
    sep           = sep,
    type          = type,
    rateParams    = rateParams,
    rateIdx       = rateIdx,
    verbose       = verbose,
    mc.cores      = mc.cores
  )

  # ---- one_batch: load cc, build context, call .evalBatch once or more ----
  #
  # When perturbations is non-NULL (delta path), .evalBatch is called once per
  # perturbation theta on the same cc and batchCtx — structural caches are
  # built once and reused for all (2*nParams + 1) theta values.
  # mode = "jacobian": also computes the analytical Jacobian in one sweep.
  one_batch <- function(b, theta, useChangeContributions,
                        perturbations = NULL, mode = "outcome") {
    t_ob_start <- if (verbose >= 2) proc.time() else NULL
    cc <- tryCatch(
      contribFun(theta, b, nBatches,
                 useChangeContributions = useChangeContributions),
      error = function(e) {
        if (grepl("cannot allocate|vector memory exhausted|cannot coerce",
                  e$message, ignore.case = TRUE))
          stop(sprintf(
            "Out of memory in chain batch %d/%d. Reduce n3BatchSize.",
            b, nBatches), call. = FALSE)
        stop(e)
      })
    if (verbose >= 2) {
      t_contrib <- (proc.time() - t_ob_start)[["elapsed"]]
      message(sprintf("    [contrib %d/%d] %.2fs", b, nBatches, t_contrib))
    }

    batchCtx <- .prepBatchContext(cc, evalConfig)

    if (mode == "jacobian") {
      br <- .evalBatch(cc, theta, batchCtx, evalConfig, jacobian = TRUE)
      return(list(hat = br$hat, jac = br$jac, jac_eff_names = br$eff_names))
    }

    hat <- .evalBatch(cc, theta, batchCtx, evalConfig)

    if (is.null(perturbations))
      return(list(hat = hat))

    # Perturbation evals are independent (same cc/batchCtx, different theta).
    # Strip the C++ probability cache (computed at theta_hat) so that
    # predictProbability recomputes utilities from contribMat at the perturbed
    # theta rather than silently returning stale cached values.
    cc_pert <- cc
    cc_pert$changeUtility    <- NULL
    cc_pert$changeProbability <- NULL

    # Use mclapply with FORK when mc.cores > 1 for ~nParams-fold speedup.
    if (mc.cores > 1L && length(perturbations) > 1L) {
      perts <- parallel::mclapply(perturbations,
                                  function(tp) .evalBatch(cc_pert, tp, batchCtx, evalConfig),
                                  mc.cores = min(mc.cores, length(perturbations)))
    } else {
      perts <- lapply(perturbations, function(tp) .evalBatch(cc_pert, tp, batchCtx, evalConfig))
    }
    list(hat = hat, perturbations = perts)
    # cc and batchCtx freed on return
  }

  # ---- estimator: batch loop + accumulator reduction ----
  estimator <- function(theta, perturbations = NULL,
                        useChangeContributions = FALSE, mode = "outcome") {

    use_parallel <- mc.cores > 1L && nBatches > 1L &&
                    is.null(perturbations) && mode == "outcome"

    if (use_parallel) {
      assertChainsFreed(verbose = FALSE)
      all_batch_results <- parallel::mclapply(
        seq_len(nBatches),
        function(b) one_batch(b, theta, useChangeContributions,
                              perturbations = NULL),
        mc.cores = mc.cores)
      for (b in seq_along(all_batch_results)) {
        if (inherits(all_batch_results[[b]], "error"))
          stop(sprintf("Point-estimate batch %d failed: %s",
                       b, conditionMessage(all_batch_results[[b]])))
      }
    } else {
      all_batch_results <- vector("list", nBatches)
      for (b in seq_len(nBatches)) {
        all_batch_results[[b]] <- one_batch(b, theta, useChangeContributions,
                                            perturbations = perturbations,
                                            mode = mode)
        gc(verbose = FALSE)
      }
    }

    # Reduce hat partials across batches.
    accums_hat <- vector("list", N)
    for (br in all_batch_results) {
      for (j in seq_len(N))
        accums_hat[[j]] <- mergePartialAccum(accums_hat[[j]], br$hat[[j]])
    }
    results_hat <- .finaliseAccums(accums_hat, evalConfig)

    # Jacobian mode: reduce K_eff gradient-column accumulators, zero-pad
    # to full nParams, return list(hat, jac).
    if (mode == "jacobian") {
      jac_eff_names <- all_batch_results[[1L]]$jac_eff_names
      theta_names_j <- names(theta)
      accums_jac <- vector("list", N)
      for (br in all_batch_results) {
        if (!is.null(br$jac)) {
          for (j in seq_len(N)) {
            if (!is.null(br$jac[[j]])) {
              K_j <- length(br$jac[[j]])
              if (is.null(accums_jac[[j]]))
                accums_jac[[j]] <- vector("list", K_j)
              for (k in seq_len(K_j))
                accums_jac[[j]][[k]] <- mergePartialAccum(
                  accums_jac[[j]][[k]], br$jac[[j]][[k]])
            }
          }
        }
      }
      results_jac <- .finaliseJac(accums_jac, theta_names_j, jac_eff_names, evalConfig)
      return(list(hat = results_hat, jac = results_jac))
    }

    if (is.null(perturbations))
      return(results_hat)

    # Reduce per-perturbation partials across batches.
    pert_names   <- names(perturbations)
    accums_perts <- setNames(vector("list", length(pert_names)), pert_names)
    for (pn in pert_names) accums_perts[[pn]] <- vector("list", N)
    for (br in all_batch_results) {
      for (pn in pert_names) {
        for (j in seq_len(N))
          accums_perts[[pn]][[j]] <- mergePartialAccum(
            accums_perts[[pn]][[j]], br$perturbations[[pn]][[j]])
      }
    }
    results_perts <- setNames(vector("list", length(pert_names)), pert_names)
    for (pn in pert_names)
      results_perts[[pn]] <- .finaliseAccums(accums_perts[[pn]], evalConfig)

    list(hat = results_hat, perturbations = results_perts)
  }

  attr(estimator, "eff_names") <- eff_names
  estimator
}

# ==========================================================================
# Top-level helpers for makeEstimatorFun
#
# Defined BELOW makeEstimatorFun per "caller before callee" convention.
# Each takes an `evalConfig` bundle (built inside makeEstimatorFun) instead of
# closing over the factory's local variables.
# ==========================================================================

# --------------------------------------------------------------------------
# .prepBatchContext — θ-independent per-batch setup
#
# Builds structural frame, keep mask, ego_id_cols, sort/mass caches, and
# cache keys from a freshly flattened cc.  All outputs depend only on the
# chain data (cc) and the fixed specs — NOT on theta.  Sharing this across
# all perturbation evals on the same cc avoids redundant work.
# --------------------------------------------------------------------------
.prepBatchContext <- function(cc, evalConfig) {
  specs    <- evalConfig$specs
  N        <- evalConfig$N
  non_acc  <- evalConfig$non_acc
  sep      <- evalConfig$sep

  cs      <- cc$changeStats
  density <- cs$density
  perm    <- cc$permitted
  keep     <- density != 0L & (if (is.null(perm)) TRUE else perm)
  all_kept <- all(keep)

  # Structural frame: group-column vectors from cc (CoW references).
  structural <- groupColsList(cc)

  # Add condition columns from csMat.
  for (j in seq_len(N)) {
    cond <- specs[[j]]$condition
    if (!is.null(cond)) {
      resolved <- resolveEffectName(cond, cs$csNames)
      for (cn in resolved) {
        if (is.null(structural[[cn]]))
          structural[[cn]] <- cs$csMat[, cn]
      }
    }
  }

  # Apply keep mask to structural frame (once).
  if (!all_kept) {
    for (nm in names(structural))
      structural[[nm]] <- structural[[nm]][keep]
  }

  ego_id_cols <- detectEgoUnit(structural)

  # Build spec_cache_key, mass_cache_key, spec_ego_cache_key.
  spec_cache_key     <- character(N)
  mass_cache_key     <- character(N)
  spec_ego_cache_key <- character(N)
  for (j in non_acc) {
    spec <- specs[[j]]
    cond <- spec$condition
    resolved_cond <- if (!is.null(cond))
      resolveEffectName(cond, cs$csNames) else NULL
    cond_str <- paste(sort(
      if (is.null(resolved_cond)) character(0L)
      else as.character(resolved_cond)), collapse = ",")
    spec_cache_key[j] <- paste(spec$level, cond_str,
                               spec$na.rm, spec$egoNormalize,
                               sep = sep)
    if (isTRUE(spec$massContrasts))
      mass_cache_key[j] <- paste(spec$level, "",
                                 spec$na.rm, spec$egoNormalize,
                                 sep = sep)
    if (isTRUE(spec$egoNormalize)) {
      pagg <- unique(c(ego_id_cols,
                       getGroupVars(level = spec$level),
                       as.character(if (!is.null(resolved_cond))
                         resolved_cond else character(0L))))
      spec_ego_cache_key[j] <- paste(pagg, collapse = sep)
    }
  }

  # Build sort_caches for all unique grouping configurations.
  sort_caches <- list()
  for (gk in unique(spec_cache_key[non_acc])) {
    idx1  <- which(spec_cache_key == gk)[1L]
    spec1 <- specs[[idx1]]
    resolved_cond1 <- if (!is.null(spec1$condition))
      resolveEffectName(spec1$condition, cs$csNames) else NULL
    gvars <- getGroupVars(level = spec1$level, condition = resolved_cond1)
    sort_caches[[gk]] <- buildAggCache(
      structural   = structural,
      group_vars   = gvars,
      ego_id_cols  = ego_id_cols,
      egoNormalize = isTRUE(spec1$egoNormalize),
      na.rm        = isTRUE(spec1$na.rm))
  }
  # Mass caches (level but condition=NULL).
  for (j in non_acc) {
    mk <- mass_cache_key[j]
    if (nchar(mk) > 0L && is.null(sort_caches[[mk]])) {
      spec <- specs[[j]]
      gvars_mass <- getGroupVars(level = spec$level, condition = NULL)
      sort_caches[[mk]] <- buildAggCache(
        structural   = structural,
        group_vars   = gvars_mass,
        ego_id_cols  = ego_id_cols,
        egoNormalize = isTRUE(spec$egoNormalize),
        na.rm        = isTRUE(spec$na.rm))
    }
  }

  list(
    keep               = keep,
    all_kept           = all_kept,
    structural         = structural,
    ego_id_cols        = ego_id_cols,
    sort_caches        = sort_caches,
    spec_cache_key     = spec_cache_key,
    mass_cache_key     = mass_cache_key,
    spec_ego_cache_key = spec_ego_cache_key,
    cs                 = cs
  )
}

# --------------------------------------------------------------------------
# .evalBatch — θ-dependent prediction + aggregation
#
# Given a flattened cc and a .prepBatchContext output (batchCtx), compute
# the baseline probabilities at theta and run the predict-group loop.
# All within-call sharing (predict-group dedup, ego-vals stage-1 cache)
# is preserved.  Structural caches come from batchCtx and are shared
# across multiple calls with different theta on the same cc.
# --------------------------------------------------------------------------
.evalBatch <- function(cc, theta, batchCtx, evalConfig, jacobian = FALSE) {
  specs         <- evalConfig$specs
  N             <- evalConfig$N
  n_pred_groups <- evalConfig$n_pred_groups
  pred_group_id <- evalConfig$pred_group_id
  sep           <- evalConfig$sep
  type          <- evalConfig$type
  rateParams    <- evalConfig$rateParams
  rateIdx       <- evalConfig$rateIdx
  verbose       <- evalConfig$verbose

  keep               <- batchCtx$keep
  all_kept           <- batchCtx$all_kept
  structural         <- batchCtx$structural
  sort_caches        <- batchCtx$sort_caches
  spec_cache_key     <- batchCtx$spec_cache_key
  mass_cache_key     <- batchCtx$mass_cache_key
  spec_ego_cache_key <- batchCtx$spec_ego_cache_key

  baseline <- predictProbability(cc, theta, type, returnComponents = TRUE)

  # ---- stream specs: compute outcome → aggregate → discard ----
  t_stream_start  <- if (verbose >= 2) proc.time() else NULL
  t_predict_specs <- 0
  t_agg_specs     <- 0
  batch_partials     <- vector("list", N)
  outcomes_by_group  <- vector("list", n_pred_groups)
  ego_vals_store     <- list()

  # Jacobian initialisation: validate contribMat, extract changeProb, allocate.
  # changeProb uses the C++ cache at theta-hat (free) or recomputes from contribMat.
  if (jacobian) {
    if (is.null(cc$contribMat))
      stop(".evalBatch(jacobian=TRUE): cc$contribMat is NULL; ",
           "keepContribMat must be TRUE in makeContribFun for the analytical path.")
    K_eff        <- ncol(cc$contribMat)
    cs           <- batchCtx$cs
    eff_names_cc <- cc$effectNames
    if (is.null(eff_names_cc) || length(eff_names_cc) == 0L)
      eff_names_cc <- colnames(cc$contribMat)
    if (is.null(eff_names_cc) || length(eff_names_cc) == 0L)
      eff_names_cc <- cs$csNames
    theta_use    <- theta[eff_names_cc]
    changeProb   <- if (!is.null(cc$changeProbability) &&
                         !all(is.na(cc$changeProbability))) {
      as.numeric(cc$changeProbability)
    } else {
      utility <- as.numeric(cc$contribMat %*% theta_use)
      calculateChangeProb(utility, cc$group_id)
    }
    density   <- as.numeric(cs$density)
    batch_jac <- vector("list", N)
  }

  for (j in seq_len(N)) {
    spec <- specs[[j]]

    t_sp <- if (verbose >= 2) proc.time() else NULL
    gid <- pred_group_id[j]
    if (is.null(outcomes_by_group[[gid]])) {
      extra <- list(changeContributions = cc,
                    theta              = theta,
                    baseline           = baseline,
                    outcomesOnly       = TRUE)
      if (!is.null(spec[["massContrasts"]]))
        extra$massContrasts <- isTRUE(spec$massContrasts)
      outcomes <- do.call(spec$predictFun, c(extra, spec$predictArgs))
      outcomes_by_group[[gid]] <- outcomes
    } else {
      outcomes <- outcomes_by_group[[gid]]
    }
    if (verbose >= 2)
      t_predict_specs <- t_predict_specs + (proc.time() - t_sp)[["elapsed"]]

    t_ag <- if (verbose >= 2) proc.time() else NULL
    if (is.data.frame(outcomes)) {
      # Fall back for predict functions that don't support outcomesOnly.
      unit_pred <- outcomes
      if (isTRUE(spec$accumulated)) {
        batch_partials[[j]] <- list(main =
          aggAccum(spec$outcomeName, unit_pred,
            level = spec$level, condition = spec$condition,
            na.rm = spec$na.rm, accumulated = TRUE))
      } else {
        batch_partials[[j]] <- list(main =
          aggAccum(spec$outcomeName, unit_pred,
            level = spec$level, condition = spec$condition,
            na.rm = spec$na.rm,
            egoNormalize = spec$egoNormalize))
      }
      if (verbose >= 2)
        t_agg_specs <- t_agg_specs + (proc.time() - t_ag)[["elapsed"]]
      next
    }

    # --- outcomesOnly path ---
    outcome_vec <- outcomes[[spec$outcomeName]]

    # rateWeight: static path multiplies by per-period rate.
    if (isTRUE(spec$rateWeight) && !isTRUE(spec$dynamic)) {
      rates_cur <- if (!is.null(rateIdx)) theta[rateIdx] else rateParams
      if (!is.null(rates_cur) && !is.null(cc$period))
        outcome_vec <- outcome_vec * rates_cur[cc$period]
    }

    if (isTRUE(spec$accumulated)) {
      acc_data <- structural
      acc_data[[spec$outcomeName]] <- if (all_kept) outcome_vec
                                      else outcome_vec[keep]
      attr(acc_data, "row.names") <- .set_row_names(length(structural[[1L]]))
      class(acc_data) <- "data.frame"
      batch_partials[[j]] <- list(main =
        aggAccum(spec$outcomeName, acc_data,
          level = spec$level, condition = spec$condition,
          na.rm = spec$na.rm, accumulated = TRUE))
    } else {
      kept_vec <- if (all_kept) outcome_vec else outcome_vec[keep]

      # Stage-1 ego-normalization caching across specs sharing the same
      # prediction outcome and canonical pre_agg_vars.
      cache_j      <- sort_caches[[spec_cache_key[j]]]
      pre_ego_vals <- NULL
      ek           <- spec_ego_cache_key[j]
      if (ek != "" && !is.null(cache_j$ego_scatter_idx)) {
        ev_key <- paste(gid, spec$outcomeName, ek, sep = sep)
        if (is.null(ego_vals_store[[ev_key]])) {
          s1 <- scatter_agg_1d(kept_vec, cache_j$ego_scatter_idx,
                               cache_j$ego_nGroups, cache_j$na.rm)
          ego_vals_store[[ev_key]] <- ifelse(
            s1$count > 0L, s1$sum / s1$count, NA_real_)
        }
        pre_ego_vals <- ego_vals_store[[ev_key]]
      }

      batch_partials[[j]] <- list(main =
        aggWithCache(spec$outcomeName, kept_vec, cache_j,
                     pre_ego_vals = pre_ego_vals))

      massCols <- intersect(c("massCreation", "massDissolution"),
                            names(outcomes))
      if (length(massCols) > 0L) {
        mc_cache    <- sort_caches[[mass_cache_key[j]]]
        partial_mass <- list()
        for (mc_col in massCols) {
          mc_vec <- if (all_kept) outcomes[[mc_col]]
                    else outcomes[[mc_col]][keep]
          partial_mass[[mc_col]] <- aggWithCache(mc_col, mc_vec, mc_cache)
        }
        batch_partials[[j]]$mass <- partial_mass
      }

      # Jacobian: same pass as prediction, changeProb already computed above.
      if (jacobian) {
        oc      <- spec$outcomeName
        cache_j <- sort_caches[[spec_cache_key[j]]]
        jFun    <- if (!is.null(spec$jacobianFun)) spec$jacobianFun
                   else .resolveBuiltinJac(spec$predictFun)
        if (!is.null(jFun)) {
          Jp_row <- jFun(cc = cc, theta = theta, changeProb = changeProb,
                         oc = oc, density = density,
                         pa = spec$predictArgs, cs = cs)
          if (!is.null(Jp_row)) {
            Jp_kept    <- if (all_kept) Jp_row else Jp_row[keep, , drop = FALSE]
            jac_k_list <- vector("list", K_eff)
            for (k in seq_len(K_eff))
              jac_k_list[[k]] <- list(main = aggWithCache(oc, Jp_kept[, k], cache_j))
            batch_jac[[j]] <- jac_k_list
          }
        }
      }
    }
    if (verbose >= 2)
      t_agg_specs <- t_agg_specs + (proc.time() - t_ag)[["elapsed"]]
  }

  if (verbose >= 2) {
    t_stream <- (proc.time() - t_stream_start)[["elapsed"]]
    message(sprintf(
      "    [stream] %.2fs (predict=%.2fs, agg=%.2fs, other=%.2fs)",
      t_stream, t_predict_specs, t_agg_specs,
      t_stream - t_predict_specs - t_agg_specs))
  }
  if (jacobian) list(hat = batch_partials, jac = batch_jac, eff_names = eff_names_cc)
  else batch_partials
}

# --------------------------------------------------------------------------
# .finaliseAccums — partial accumulators → named list of data.frames
# --------------------------------------------------------------------------
.finaliseAccums <- function(accums, evalConfig) {
  specs     <- evalConfig$specs
  N         <- evalConfig$N
  eff_names <- evalConfig$eff_names

  results <- setNames(vector("list", N), eff_names)
  for (j in seq_len(N)) {
    spec        <- specs[[j]]
    results[[j]] <- finalizeAccum(accums[[j]],
                                           spec$outcomeName,
                                           spec$level)
  }
  results
}

# --------------------------------------------------------------------------
# .finaliseJac — Jacobian accumulators → named list of [n_out x nParams] matrices
#
# accums_jac[[j]]: list of K_eff partial accums (one per effect column).
# theta_names:     character(nParams) from names(thetaHat).
# eff_names_j:     cc$effectNames — K_eff subset of theta_names.
# Returns: named list; non-effect columns are zero.
# --------------------------------------------------------------------------
.finaliseJac <- function(accums_jac, theta_names, eff_names_j, evalConfig) {
  specs     <- evalConfig$specs
  N         <- evalConfig$N
  eff_names <- evalConfig$eff_names

  nParams <- length(theta_names)
  results_jac <- setNames(vector("list", N), eff_names)
  eff_idx <- match(eff_names_j, theta_names)
  eff_idx <- eff_idx[!is.na(eff_idx)]
  K_eff   <- length(eff_idx)
  for (j in seq_len(N)) {
    if (is.null(accums_jac[[j]])) { results_jac[j] <- list(NULL); next }
    if (K_eff == 0L) { results_jac[j] <- list(NULL); next }
    spec  <- specs[[j]]
    oc    <- spec$outcomeName
    col1  <- finalizeAccum(accums_jac[[j]][[1L]], oc, spec$level)
    n_out <- nrow(col1)
    J     <- matrix(0, nrow = n_out, ncol = nParams,
                    dimnames = list(NULL, theta_names))
    J[, eff_idx[1L]] <- col1[[oc]]
    for (k in seq(2L, K_eff)) {
      col_k        <- finalizeAccum(accums_jac[[j]][[k]], oc, spec$level)
      J[, eff_idx[k]] <- col_k[[oc]]
    }
    results_jac[[j]] <- J
  }
  results_jac
}

# --------------------------------------------------------------------------
# aggregatePostEstimation
#
# Extracted per-effect aggregation loop used by the new sienaPostestimate
# and by marginalEffects when the stored-chains case needs separate
# expects vs draw estimators.
#
# expects:              named list of data.frames (point estimates, one per effect)
# raw_sims_list:        named list of stacked draw data.frames (from drawSimBatch)
# specs:                named list — same structure as passed to makeEstimatorFun
# uncertainty_summary_fun: function(x, na.rm) → named list of scalars
# na.rm, egoNormalize:  passed through to agg/aggMulti
# decisionDetails:      optional named list; entries attached as attr if present
# saveDir:              optional path — saves each result as <nm>.rds
# --------------------------------------------------------------------------
aggregatePostEstimation <- function(expects, raw_sims_list, specs,
    uncertainty_summary_fun, na.rm = TRUE, egoNormalize = TRUE,
    decisionDetails = NULL, saveDir = NULL) {

  results <- vector("list", length(specs))
  names(results) <- names(specs)

  for (nm in names(specs)) {
    spec     <- specs[[nm]]
    raw_j    <- raw_sims_list[[nm]]
    expect_j <- expects[[nm]]

    massCols <- intersect(c("massCreation", "massDissolution"), names(raw_j))

    if (is.null(spec$condition) && length(massCols) > 0L) {
      all_uncert <- aggMeanSim(
        c(spec$outcomeName, massCols), raw_j,
        level        = spec$level,
        condition    = NULL,
        sum_fun      = uncertainty_summary_fun,
        na.rm        = na.rm,
        egoNormalize = egoNormalize)
      uncert_j <- all_uncert[[spec$outcomeName]]
    } else {
      uncert_j <- aggUncertainty(
        spec$outcomeName, raw_j,
        level        = spec$level,
        condition    = spec$condition,
        sum_fun      = uncertainty_summary_fun,
        na.rm        = na.rm,
        egoNormalize = egoNormalize)
      if (length(massCols) > 0L)
        all_uncert <- aggMeanSim(massCols, raw_j,
          level        = spec$level,
          condition    = NULL,
          sum_fun      = uncertainty_summary_fun,
          na.rm        = na.rm,
          egoNormalize = egoNormalize)
    }

    result_j <- mergeEstimates(expect_j, uncert_j,
                               level     = spec$level,
                               condition = spec$condition)

    for (mc in massCols) {
      mc_uncert <- all_uncert[[mc]]
      level_by  <- intersect(
        getGroupVars(level = spec$level, condition = NULL),
        names(mc_uncert))
      uc_cols <- setdiff(names(mc_uncert), level_by)
      for (uc in uc_cols)
        names(mc_uncert)[names(mc_uncert) == uc] <- paste0(mc, "_", uc)
      if (length(level_by) > 0L) {
        result_j <- merge(result_j, mc_uncert,
          by = level_by, all.x = TRUE, sort = FALSE)
      } else {
        for (col in setdiff(names(mc_uncert), level_by))
          result_j[[col]] <- mc_uncert[[col]]
      }
    }

    result_j <- attachPostestAttrs(result_j, spec$metadata)

    if (!is.null(decisionDetails) && !is.null(decisionDetails[[nm]]))
      attr(result_j, "decisionDetails") <- decisionDetails[[nm]]

    if (!is.null(saveDir))
      saveRDS(result_j, file.path(saveDir, paste0(nm, ".rds")))

    results[[nm]] <- result_j
  }

  results
}

# --------------------------------------------------------------------------
# combinePostestResults — opt-in wide merge for specs sharing (level, condition)
#
# When multiple effects in `results` were computed at the same (level, condition)
# and the user opts in via combineSameLevel = TRUE, merge them into a single
# wider data.frame rather than returning N separate ones.
#
# `results`: named list of data.frames (output of aggregatePostEstimation or
#            sienaPostestimate). Each element may carry class-level attributes.
# `specs`:   named list of spec entries (same keys as results).
#
# Returns: named list of data.frames. Each group-of-same-level specs is replaced
# by a single wider data.frame keyed on group columns; elements that are alone at
# their (level, condition) stay as single-element lists (unchanged).
# --------------------------------------------------------------------------
combinePostestResults <- function(results, specs) {
  if (length(results) == 0L) return(results)

  # Group specs by (level, condition) key.
  group_key <- vapply(names(specs), function(nm) {
    sp  <- specs[[nm]]
    lv  <- if (!is.null(sp$level)) sp$level else "none"
    cnd <- if (!is.null(sp$condition)) paste(sp$condition, collapse = ",") else ""
    paste0(lv, "|", cnd)
  }, character(1L))

  out <- list()
  for (gk in unique(group_key)) {
    nms <- names(group_key)[group_key == gk]
    if (length(nms) == 1L) {
      out[[nms]] <- results[[nms]]
    } else {
      # Identify the structural (group) columns shared by all frames in group.
      # These are columns that appear in all frames with identical values.
      first <- results[[nms[1L]]]
      sp    <- specs[[nms[1L]]]
      lv    <- if (!is.null(sp$level)) sp$level else "none"
      cnd   <- if (!is.null(sp$condition))
                  resolveEffectName(sp$condition, names(first)) else NULL
      key_cols <- getGroupVars(level = lv, condition = cnd)

      # Merge all frames on key_cols.
      merged <- Reduce(function(a, b) {
        if (length(key_cols) > 0L)
          merge(a, b, by = key_cols, all = TRUE, sort = FALSE)
        else
          cbind(a, b[, setdiff(names(b), names(a)), drop = FALSE])
      }, lapply(nms, function(nm) results[[nm]]))

      # Name the combined element after the group's effects joined by "+".
      combined_nm <- paste(nms, collapse = "+")
      out[[combined_nm]] <- merged
    }
  }
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

# ======================================================================
# SECTION 3: Aggregation
# ======================================================================

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
  # For sim-aware uncertainty aggregation, use aggUncertainty() instead.
  if (egoNormalize) {
    ego_id_cols <- detectEgoUnit(data)
    extra <- setdiff(ego_id_cols, group_vars)
    if (length(extra) > 0) {
      data <- preAggEgo(data, outcomeName, group_vars, ego_id_cols, na.rm)
    }
  }

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

# Aggregate multiple outcome columns sharing the expensive preAggEgo step.
# Returns a named list of data.frames, one per outcomeName.
aggMulti <- function(outcomeNames, data, level = "none", condition = NULL,
                     sum_fun = mean, na.rm = TRUE, egoNormalize = TRUE) {
  if (length(outcomeNames) <= 1L) {
    return(setNames(
      list(agg(outcomeNames[1L], data, level, condition, sum_fun, na.rm,
               egoNormalize)),
      outcomeNames
    ))
  }
  group_vars <- getGroupVars(level = level, condition = condition)

  # Share ego-level pre-aggregation across all outcome columns:
  # encode group keys and order once, then run grouped_agg_cpp per column.
  # For sim-aware uncertainty aggregation, use aggMultiUncertainty() instead.
  if (egoNormalize && length(group_vars) > 0) {
    ego_id_cols <- detectEgoUnit(data)
    extra <- setdiff(ego_id_cols, group_vars)
    if (length(extra) > 0) {
      pre_agg_vars <- unique(c(group_vars, ego_id_cols))
      enc  <- encodeGroupKeys(data, pre_agg_vars)
      ord  <- do.call(order,
                 lapply(seq_len(ncol(enc$G)), function(j) enc$G[, j]))
      Gord <- enc$G[ord, , drop = FALSE]

      # First column: also extract key structure.
      res0   <- grouped_agg_cpp(data[[outcomeNames[1L]]][ord], Gord,
                                na_rm = na.rm, do_mean = TRUE)
      preagg <- decodeGroupKeys(res0$key, pre_agg_vars, enc$decode)
      preagg[[outcomeNames[1L]]] <- res0$value

      # Remaining columns: same key, different values.
      for (nm in outcomeNames[-1L]) {
        res_i <- grouped_agg_cpp(data[[nm]][ord], Gord,
                                 na_rm = na.rm, do_mean = TRUE)
        preagg[[nm]] <- res_i$value
      }

      # Outer aggregation on the small pre-aggregated data.
      results <- lapply(outcomeNames, function(nm)
        agg(nm, preagg, level = level, condition = condition,
            sum_fun = sum_fun, na.rm = na.rm, egoNormalize = FALSE))
      return(setNames(results, outcomeNames))
    }
  }

  # Fallback: individual calls.
  results <- lapply(outcomeNames, function(nm)
    agg(nm, data, level, condition, sum_fun, na.rm, egoNormalize))
  setNames(results, outcomeNames)
}

# aggMean: canonical multi-column mean aggregation (always takes a vector of names).
aggMean <- aggMulti

# ---- Sim-aware uncertainty aggregation ----------------------------------------
# collapseSimEgo: two-step ego normalization for simulation draws.
#   Step 1: within-ego means per {sim, ego, group_vars} (average across alters)
#   Step 2: average across egos within {sim, group_vars}
# Returns data with one row per sim × group combination.
# outcomeNames: one or more column names to collapse (vectorised for efficiency).
collapseSimEgo <- function(outcomeNames, data, group_vars, na.rm) {
  ego_id_cols <- detectEgoUnit(data)
  extra <- setdiff(ego_id_cols, group_vars)
  if (length(extra) == 0L) return(data)

  sim_group_vars <- unique(c("sim", group_vars))
  pre_agg_vars   <- unique(c(sim_group_vars, ego_id_cols))

  # Step 1: within-ego means per sim (across alters)
  enc  <- encodeGroupKeys(data, pre_agg_vars)
  ord  <- do.call(order, lapply(seq_len(ncol(enc$G)), function(j) enc$G[, j]))
  Gord <- enc$G[ord, , drop = FALSE]
  res0 <- grouped_agg_cpp(data[[outcomeNames[1L]]][ord], Gord,
                          na_rm = na.rm, do_mean = TRUE)
  preagg <- decodeGroupKeys(res0$key, pre_agg_vars, enc$decode)
  preagg[[outcomeNames[1L]]] <- res0$value
  for (nm in outcomeNames[-1L]) {
    preagg[[nm]] <- grouped_agg_cpp(data[[nm]][ord], Gord,
                                    na_rm = na.rm, do_mean = TRUE)$value
  }

  # Step 2: average across egos within {sim, group_vars}
  enc2  <- encodeGroupKeys(preagg, sim_group_vars)
  ord2  <- do.call(order, lapply(seq_len(ncol(enc2$G)), function(j) enc2$G[, j]))
  Gord2 <- enc2$G[ord2, , drop = FALSE]
  res2_0 <- grouped_agg_cpp(preagg[[outcomeNames[1L]]][ord2], Gord2,
                            na_rm = na.rm, do_mean = TRUE)
  result <- decodeGroupKeys(res2_0$key, sim_group_vars, enc2$decode)
  result[[outcomeNames[1L]]] <- res2_0$value
  for (nm in outcomeNames[-1L]) {
    result[[nm]] <- grouped_agg_cpp(preagg[[nm]][ord2], Gord2,
                                    na_rm = na.rm, do_mean = TRUE)$value
  }
  attr(result, "row.names") <- .set_row_names(length(res2_0$value))
  class(result) <- "data.frame"
  result
}

# aggUncertainty: sim-aware aggregation for uncertainty estimation.
# Collapses sim draws through ego normalization, then applies sum_fun across
# sims per group.  Falls back to agg() when data has no "sim" column.
aggUncertainty <- function(outcomeName, data, level = "none", condition = NULL,
                           sum_fun = mean, na.rm = TRUE, egoNormalize = TRUE) {
  if (!("sim" %in% names(data)))
    return(agg(outcomeName, data, level, condition, sum_fun, na.rm, egoNormalize))
  if (!is.null(condition)) condition <- resolveEffectName(condition, names(data))
  if (egoNormalize) {
    group_vars <- getGroupVars(level = level, condition = condition)
    data <- collapseSimEgo(outcomeName, data, group_vars, na.rm)
  }
  agg(outcomeName, data, level = level, condition = condition,
      sum_fun = sum_fun, na.rm = na.rm, egoNormalize = FALSE)
}

# aggMultiUncertainty: multi-column variant of aggUncertainty.
aggMultiUncertainty <- function(outcomeNames, data, level = "none",
                                condition = NULL, sum_fun = mean,
                                na.rm = TRUE, egoNormalize = TRUE) {
  if (length(outcomeNames) <= 1L)
    return(setNames(
      list(aggUncertainty(outcomeNames[1L], data, level, condition,
                          sum_fun, na.rm, egoNormalize)),
      outcomeNames))
  if (!("sim" %in% names(data)))
    return(aggMulti(outcomeNames, data, level, condition, sum_fun, na.rm,
                    egoNormalize))
  if (!is.null(condition)) condition <- resolveEffectName(condition, names(data))
  if (egoNormalize) {
    group_vars <- getGroupVars(level = level, condition = condition)
    data <- collapseSimEgo(outcomeNames, data, group_vars, na.rm)
  }
  aggMulti(outcomeNames, data, level = level, condition = condition,
           sum_fun = sum_fun, na.rm = na.rm, egoNormalize = FALSE)
}

# aggMeanSim: canonical sim-aware uncertainty aggregation (always multi-col).
aggMeanSim <- aggMultiUncertainty

# ---- Batched aggregation helpers ------------------------------------------------
# These support n3-batched accumulation: processing chains in smaller batches
# to avoid OOM from materializing all chains at once.

# --------------------------------------------------------------------------
# buildAggCache / aggWithCache — pre-sort once, aggregate many
#
# When multiple specs share the same (level, condition, na.rm, egoNormalize,
# structural frame), encodeGroupKeys + order() are identical across specs.
# buildAggCache does that work once; aggWithCache applies it per outcome vector.
#
# This turns O(N × n log n) into O(n log n + N × n).
# --------------------------------------------------------------------------

# Build a reusable aggregation cache for a given structural frame + grouping.
#
# structural: named list of integer/numeric vectors (the grouping columns from cc),
#             already subsetted by `keep` if applicable.  NOT a data.frame.
# group_vars: character vector of column names to group by (from getGroupVars).
# ego_id_cols: character vector of ego-unit columns (from detectEgoUnit).
# egoNormalize: whether to pre-aggregate within ego units.
# na.rm: passed to grouped_agg_cpp.
#
# Returns a list with class "aggCache" containing pre-computed sort permutations.
buildAggCache <- function(structural, group_vars, ego_id_cols,
                          egoNormalize = TRUE, na.rm = TRUE) {
  cache <- list(group_vars = group_vars, na.rm = na.rm,
                egoNormalize = egoNormalize)

  needs_ego <- egoNormalize && length(setdiff(ego_id_cols, group_vars)) > 0L

  if (needs_ego) {
    # Stage 1: ego pre-aggregation (choice-level → ego-level means).
    # Use scatter-accumulate instead of sort-permutation: precompute a
    # row→group label vector (ego_scatter_idx) so that aggWithCache
    # can do a sequential scan of vals rather than a random-access permutation.
    # Sequential reads are ~100× faster for 24M-row data (L3-miss vs sequential).
    pre_agg_vars <- unique(c(ego_id_cols, group_vars))
    ego_enc  <- encodeGroupKeys(structural, pre_agg_vars)
    ego_ord  <- do.call(order,
                  lapply(seq_len(ncol(ego_enc$G)), function(j) ego_enc$G[, j]))
    ego_G_sorted <- ego_enc$G[ego_ord, , drop = FALSE]

    # Build scatter index: ego_scatter_idx[i] = ego-group label for original row i
    cache$ego_scatter_idx <- build_scatter_idx(ego_G_sorted, ego_ord)

    # Run a dummy aggregation to determine ego-group count and key matrix.
    dummy <- rep(1.0, length(ego_ord))
    ego_res <- grouped_agg_cpp(dummy[ego_ord], ego_G_sorted,
                               na_rm = FALSE, do_mean = TRUE)
    cache$ego_nGroups <- nrow(ego_res$key)

    # Stage 2: scatter index on ego-level data (nEgoGroups rows).
    if (length(group_vars) > 0L) {
      ego_structural <- decodeGroupKeys(ego_res$key, pre_agg_vars, ego_enc$decode)
      main_enc  <- encodeGroupKeys(ego_structural, group_vars)
      main_ord  <- do.call(order,
                    lapply(seq_len(ncol(main_enc$G)), function(j) main_enc$G[, j]))
      main_G_sorted <- main_enc$G[main_ord, , drop = FALSE]

      # level_scatter_idx[g] = level-group label for ego-group g
      cache$level_scatter_idx <- build_scatter_idx(main_G_sorted, main_ord)

      # Precompute decoded group-key data frame for the output (one row per level group).
      dummy_ego <- rep(1.0, nrow(main_enc$G))
      level_res <- grouped_agg_cpp(dummy_ego[main_ord], main_G_sorted,
                                   na_rm = FALSE, do_mean = TRUE)
      cache$level_nGroups <- nrow(level_res$key)
      cache$level_key_df  <- decodeGroupKeys(level_res$key, group_vars,
                                             main_enc$decode)
    }
  } else if (length(group_vars) > 0L) {
    # No ego normalization — scatter directly on the raw structural frame.
    enc <- encodeGroupKeys(structural, group_vars)
    ord <- do.call(order,
              lapply(seq_len(ncol(enc$G)), function(j) enc$G[, j]))
    G_sorted <- enc$G[ord, , drop = FALSE]

    cache$main_scatter_idx <- build_scatter_idx(G_sorted, ord)

    dummy <- rep(1.0, length(ord))
    main_res <- grouped_agg_cpp(dummy[ord], G_sorted, na_rm = FALSE, do_mean = TRUE)
    cache$main_nGroups <- nrow(main_res$key)
    cache$main_key_df  <- decodeGroupKeys(main_res$key, group_vars, enc$decode)
  }
  # else: no group_vars and no ego normalization → scalar aggregation

  class(cache) <- "aggCache"
  cache
}

# Aggregate a single outcome vector using a pre-computed cache.
# Returns a data.frame with group columns + "{outcomeName}_sum" + "{outcomeName}_n".
aggWithCache <- function(outcomeName, vals, cache, pre_ego_vals = NULL) {
  na.rm <- cache$na.rm
  sumCol <- paste0(outcomeName, "_sum")
  cntCol <- paste0(outcomeName, "_n")

  if (!is.null(cache$ego_scatter_idx)) {
    # Stage 1: scatter-accumulate choice-level → ego-level means.
    # Sequential scan of vals (cache-friendly) + random writes to small
    # ego-bucket accumulator (fits in L2 cache).
    if (is.null(pre_ego_vals)) {
      s1 <- scatter_agg_1d(vals, cache$ego_scatter_idx,
                           cache$ego_nGroups, na.rm)
      ego_counts <- s1$count
      ego_means  <- ifelse(ego_counts > 0L, s1$sum / ego_counts, NA_real_)
    } else {
      ego_means <- pre_ego_vals
    }
    vals <- ego_means
  }

  if (!is.null(cache$level_scatter_idx)) {
    # Stage 2: scatter on ego-level data (tiny; always fast).
    s2  <- scatter_agg_1d(vals, cache$level_scatter_idx,
                          cache$level_nGroups, na.rm)
    out <- cache$level_key_df
    out[[sumCol]] <- s2$sum
    out[[cntCol]] <- as.double(s2$count)
  } else if (!is.null(cache$main_scatter_idx)) {
    # No ego normalization — scatter directly on raw data.
    s   <- scatter_agg_1d(vals, cache$main_scatter_idx,
                          cache$main_nGroups, na.rm)
    out <- cache$main_key_df
    out[[sumCol]] <- s$sum
    out[[cntCol]] <- as.double(s$count)
  } else {
    # Scalar aggregation (no grouping, e.g. level = "none").
    if (na.rm) vals <- vals[!is.na(vals)]
    out <- list()
    out[[sumCol]] <- sum(vals)
    out[[cntCol]] <- as.double(length(vals))
  }

  attr(out, "row.names") <- .set_row_names(length(out[[1L]]))
  class(out) <- "data.frame"
  out
}

# Returns sum + count per group rather than the mean, for batch accumulation.
# Output columns: group_vars + "{outcomeName}_sum" + "{outcomeName}_n".
aggSumCount <- function(outcomeName, data, level = "none", condition = NULL,
                        na.rm = TRUE, egoNormalize = TRUE) {
  if (!is.null(condition)) condition <- resolveEffectName(condition, names(data))
  group_vars <- getGroupVars(level = level, condition = condition)
  sumCol <- paste0(outcomeName, "_sum")
  cntCol <- paste0(outcomeName, "_n")

  if (egoNormalize) {
    ego_id_cols <- detectEgoUnit(data)
    extra <- setdiff(ego_id_cols, group_vars)
    if (length(extra) > 0)
      data <- preAggEgo(data, outcomeName, group_vars, ego_id_cols, na.rm)
  }

  if (length(group_vars) > 0) {
    enc  <- encodeGroupKeys(data, group_vars)
    ord  <- do.call(order,
               lapply(seq_len(ncol(enc$G)), function(j) enc$G[, j]))
    vals <- data[[outcomeName]][ord]
    Gord <- enc$G[ord, , drop = FALSE]

    res_sum <- grouped_agg_cpp(vals, Gord, na_rm = na.rm, do_mean = FALSE)
    not_na  <- as.double(!is.na(vals))
    res_cnt <- grouped_agg_cpp(not_na, Gord, na_rm = FALSE, do_mean = FALSE)

    out <- decodeGroupKeys(res_sum$key, group_vars, enc$decode)
    out[[sumCol]] <- res_sum$value
    out[[cntCol]] <- res_cnt$value
  } else {
    vals <- data[[outcomeName]]
    if (na.rm) vals <- vals[!is.na(vals)]
    out <- list()
    out[[sumCol]] <- sum(vals)
    out[[cntCol]] <- length(vals)
  }

  attr(out, "row.names") <- .set_row_names(length(out[[1L]]))
  class(out) <- "data.frame"
  out
}

# === FIX B: new function — delete to revert (see one_batch comment) =========
# Vectorized variant of aggSumCount: aggregates multiple outcome columns using
# a single sort + encodeGroupKeys pass per step (ego-normalize and group-level).
# outcomeVecs : named list of numeric vectors, all length == nrow(data).
# data        : shared structural data frame supplying grouping/condition columns.
# Returns     : named list of data.frames matching the aggSumCount output format.
batchAggSumCount <- function(outcomeVecs, data, level = "none",
                              condition = NULL, na.rm = TRUE,
                              egoNormalize = TRUE) {
  if (length(outcomeVecs) == 0L) return(list())
  outcomeNames <- names(outcomeVecs)
  for (nm in outcomeNames) data[[nm]] <- outcomeVecs[[nm]]
  if (!is.null(condition))
    condition <- resolveEffectName(condition, names(data))
  group_vars <- getGroupVars(level = level, condition = condition)

  if (egoNormalize) {
    ego_id_cols <- detectEgoUnit(data)
    extra <- setdiff(ego_id_cols, group_vars)
    if (length(extra) > 0L)
      data <- preAggEgoMulti(data, outcomeNames, group_vars, ego_id_cols, na.rm)
  }

  result_list <- setNames(vector("list", length(outcomeNames)), outcomeNames)

  if (length(group_vars) > 0L) {
    enc   <- encodeGroupKeys(data, group_vars)
    ord   <- do.call(order,
               lapply(seq_len(ncol(enc$G)), function(j) enc$G[, j]))
    G_ord  <- enc$G[ord, , drop = FALSE]
    key_df <- NULL
    for (nm in outcomeNames) {
      vals    <- data[[nm]][ord]
      res_sum <- grouped_agg_cpp(vals, G_ord, na_rm = na.rm, do_mean = FALSE)
      res_cnt <- grouped_agg_cpp(as.double(!is.na(vals)), G_ord,
                                  na_rm = FALSE, do_mean = FALSE)
      if (is.null(key_df))
        key_df <- as.data.frame(
          decodeGroupKeys(res_sum$key, group_vars, enc$decode),
          stringsAsFactors = FALSE)
      out <- key_df
      out[[paste0(nm, "_sum")]] <- res_sum$value
      out[[paste0(nm, "_n")]]   <- res_cnt$value
      attr(out, "row.names") <- .set_row_names(nrow(key_df))
      class(out) <- "data.frame"
      result_list[[nm]] <- out
    }
  } else {
    for (nm in outcomeNames) {
      vals <- data[[nm]]
      if (na.rm) vals <- vals[!is.na(vals)]
      out <- list()
      out[[paste0(nm, "_sum")]] <- sum(vals)
      out[[paste0(nm, "_n")]]   <- length(vals)
      attr(out, "row.names") <- .set_row_names(1L)
      class(out) <- "data.frame"
      result_list[[nm]] <- out
    }
  }
  result_list
}

# Like aggAccumulated() but returns sum + count at the final step.
# Output columns: final_group_vars + "{outcomeName}_sum" + "{outcomeName}_n".
aggAccumulatedSumCount <- function(outcomeName, data, level = "period",
                                   condition = NULL, na.rm = TRUE) {
  if (!"chain" %in% names(data))
    stop("Accumulated aggregation requires dynamic (chain-level) data.")
  if (!is.null(condition))
    warning("'condition' is ignored for accumulated aggregation.")

  sumCol <- paste0(outcomeName, "_sum")
  cntCol <- paste0(outcomeName, "_n")

  # Step 1: Ego-normalize (mean over alters within each ministep).
  ego_id_cols <- detectEgoUnit(data)
  step1 <- preAggEgo(data, outcomeName,
                     group_vars = character(0),
                     ego_id_cols = ego_id_cols,
                     na.rm = na.rm)

  # Step 2: Sum across ministeps within (period, chain, ego).
  acc_group <- intersect(c("group", "period", "chain", "ego"), names(step1))
  enc  <- encodeGroupKeys(step1, acc_group)
  ord  <- do.call(order,
             lapply(seq_len(ncol(enc$G)), function(j) enc$G[, j]))
  res  <- grouped_agg_cpp(step1[[outcomeName]][ord],
                           enc$G[ord, , drop = FALSE],
                           na_rm = na.rm, do_mean = FALSE)
  step2 <- decodeGroupKeys(res$key, acc_group, enc$decode)
  step2[[outcomeName]] <- res$value

  # Step 3: Return sum + count per final-level group (not mean).
  final_group <- switch(level,
    "ego"    = intersect(c("group", "period", "ego"), names(step2)),
    "period" = intersect(c("group", "period"), names(step2)),
    "none"   = character(0),
    intersect(c("group", "period"), names(step2))
  )

  if (length(final_group) > 0) {
    enc2 <- encodeGroupKeys(step2, final_group)
    ord2 <- do.call(order,
               lapply(seq_len(ncol(enc2$G)), function(j) enc2$G[, j]))
    vals2 <- step2[[outcomeName]][ord2]
    G2    <- enc2$G[ord2, , drop = FALSE]

    res_sum <- grouped_agg_cpp(vals2, G2, na_rm = na.rm, do_mean = FALSE)
    not_na  <- as.double(!is.na(vals2))
    res_cnt <- grouped_agg_cpp(not_na, G2, na_rm = FALSE, do_mean = FALSE)

    out <- decodeGroupKeys(res_sum$key, final_group, enc2$decode)
    out[[sumCol]] <- res_sum$value
    out[[cntCol]] <- res_cnt$value
  } else {
    vals <- step2[[outcomeName]]
    if (na.rm) vals <- vals[!is.na(vals)]
    out <- list()
    out[[sumCol]] <- sum(vals)
    out[[cntCol]] <- length(vals)
  }

  attr(out, "row.names") <- .set_row_names(length(out[[1L]]))
  class(out) <- "data.frame"
  out
}

# aggAccum: unified accumulation function (replaces aggSumCount + aggAccumulatedSumCount).
# accumulated=FALSE → sum+count per group (batched mean).
# accumulated=TRUE  → sum ministep-chains before sum+count (dynamic accumulated ME).
aggAccum <- function(outcomeName, data, level = "none", condition = NULL,
                     na.rm = TRUE, egoNormalize = TRUE, accumulated = FALSE) {
  if (accumulated)
    aggAccumulatedSumCount(outcomeName, data, level = level,
                           condition = condition, na.rm = na.rm)
  else
    aggSumCount(outcomeName, data, level = level, condition = condition,
                na.rm = na.rm, egoNormalize = egoNormalize)
}

# Merge two sum+count accumulators by group columns, adding sums and counts.
# If prev is NULL, returns curr unchanged.
mergeBatchAccum <- function(prev, curr) {
  if (is.null(prev)) return(curr)

  all_cols   <- names(curr)
  sc_cols    <- grep("_sum$|_n$", all_cols, value = TRUE)
  group_cols <- setdiff(all_cols, sc_cols)

  if (length(group_cols) > 0) {
    merged <- merge(prev, curr, by = group_cols, all = TRUE,
                    suffixes = c(".prev", ".curr"))
    for (sc in sc_cols) {
      col_prev <- paste0(sc, ".prev")
      col_curr <- paste0(sc, ".curr")
      merged[[sc]] <- ifelse(is.na(merged[[col_prev]]), 0,
                             merged[[col_prev]]) +
                      ifelse(is.na(merged[[col_curr]]), 0,
                             merged[[col_curr]])
      merged[[col_prev]] <- NULL
      merged[[col_curr]] <- NULL
    }
    merged
  } else {
    for (sc in sc_cols) {
      prev[[sc]] <- prev[[sc]] + curr[[sc]]
    }
    prev
  }
}

# Convert sum+count accumulator to mean: {name} = {name}_sum / {name}_n.
# Drops the _sum and _n columns.
finalizeBatchMean <- function(accum, outcomeName) {
  sumCol <- paste0(outcomeName, "_sum")
  cntCol <- paste0(outcomeName, "_n")
  accum[[outcomeName]] <- accum[[sumCol]] / accum[[cntCol]]
  accum[[sumCol]] <- NULL
  accum[[cntCol]] <- NULL
  accum
}

mergePartialAccum <- function(prev, curr) {
  if (is.null(prev)) return(curr)
  result <- list(main = mergeBatchAccum(prev$main, curr$main))
  if (!is.null(curr$mass)) {
    result$mass <- list()
    for (mc in names(curr$mass))
      result$mass[[mc]] <- mergeBatchAccum(prev$mass[[mc]], curr$mass[[mc]])
  }
  result
}

finalizePartialResult <- function(accum, outcomeName, level) {
  main_result <- finalizeBatchMean(accum$main, outcomeName)
  if (!is.null(accum$mass)) {
    for (mc in names(accum$mass)) {
      mc_df <- finalizeBatchMean(accum$mass[[mc]], mc)
      mc_by <- intersect(
        getGroupVars(level = level, condition = NULL),
        intersect(names(main_result), names(mc_df)))
      if (length(mc_by) > 0L) {
        main_result <- merge(main_result, mc_df,
          by = mc_by, all.x = TRUE, sort = FALSE)
      } else {
        main_result[[mc]] <- mc_df[[mc]]
      }
    }
  }
  main_result
}

# finalizeAccum: canonical name for finalizePartialResult.
finalizeAccum <- finalizePartialResult

# --------------------------------------------------------------------------
# makeContribFun
#
# Replaces makeChainSource. Returns function(theta, batchIdx, nBatches) → cc.
# Modes:
#   "static"     – theta-independent; calls getContribFun(theta) once.
#   "preloaded"  – stored chains, theta-independent; subsets rawChains per batch.
#   "per_batch"  – fresh siena07 simulation per (theta, batch).
# keepContribMat: when TRUE, cc$contribMat is NOT dropped (needed for RI).
# --------------------------------------------------------------------------
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

# === FIX B: new function — delete to revert (see one_batch comment) =========
# Like preAggEgo but computes within-ego means for multiple outcome columns in
# a single sort + encodeGroupKeys pass.  All outcomeNames must already exist in
# data as numeric columns with the same row ordering.
preAggEgoMulti <- function(data, outcomeNames, group_vars, ego_id_cols, na.rm) {
  pre_agg_vars <- unique(c(group_vars, ego_id_cols))
  enc     <- encodeGroupKeys(data, pre_agg_vars)
  ord     <- do.call(order, lapply(seq_len(ncol(enc$G)), function(j) enc$G[, j]))
  G_sorted <- enc$G[ord, , drop = FALSE]
  key_mat <- NULL
  avgs    <- setNames(vector("list", length(outcomeNames)), outcomeNames)
  for (nm in outcomeNames) {
    res        <- grouped_agg_cpp(data[[nm]][ord], G_sorted,
                                  na_rm = na.rm, do_mean = TRUE)
    avgs[[nm]] <- res$value
    if (is.null(key_mat)) key_mat <- res$key
  }
  out <- as.data.frame(decodeGroupKeys(key_mat, pre_agg_vars, enc$decode),
                       stringsAsFactors = FALSE)
  for (nm in outcomeNames) out[[nm]] <- avgs[[nm]]
  out
}


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
# ======================================================================
# SECTION 4: Utilities
# ======================================================================

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

getEffectMetaNoRate <- function(effects, depvar) {
  reg <- buildEffectNameRegistry(effects,
                   depvar = depvar,
                   includeOnly = TRUE,
                   append_parm = FALSE)
  reg <- reg[reg$effect_type != "rate", , drop = FALSE]

  if (nrow(reg) == 0L) {
    return(list(
      base_names = character(0),
      interaction_types = setNames(character(0), character(0))
    ))
  }

  base_names <- unique(reg$base_name)

  # Collapse creation/endow/eval variants to one interactionType per base effect.
  if (!all(is.na(reg$interaction_type))) {
    split_types <- split(reg$interaction_type, reg$base_name)
    int_types <- vapply(base_names, function(nm) {
      vals <- unique(split_types[[nm]])
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (length(vals) == 0L) return(NA_character_)
      if (length(vals) > 1L) {
        stop("Inconsistent interactionType across type variants for effect '",
             nm, "': ", paste(vals, collapse = ", "), call. = FALSE)
      }
      vals[[1L]]
    }, character(1L), USE.NAMES = FALSE)
    names(int_types) <- base_names
  } else {
    int_types <- setNames(rep(NA_character_, length(base_names)), base_names)
  }

  list(base_names = base_names, interaction_types = int_types)
}

# Backward-compatible accessor: keep a single implementation source in
# getEffectMetaNoRate(), and expose names where legacy callers expect it.
getEffectNamesNoRate <- function(effects, depvar) {
  getEffectMetaNoRate(effects, depvar)$base_names
}

# Effect-name resolver used by post-estimation helpers.
# Resolve against an explicit candidate vector (caller selects keyspace).
resolveEffectName <- function(effectName,
                candidateNames,
                              allowLegacySuffix = TRUE,
                              strict = TRUE) {
    if (is.null(effectName)) return(NULL)
  if (!is.character(candidateNames))
    stop("'candidateNames' must be a character vector.")
  vals <- unique(as.character(candidateNames))

    resolveOne <- function(nm) {
        if (nm %in% vals) return(nm)

        nmStripped <- if (allowLegacySuffix)
            sub("_(eval|endow|creation|rate)$", "", nm) else nm
        if (nmStripped %in% vals) return(nmStripped)

        if (!grepl("_", nmStripped, fixed = TRUE)) {
            intFamilies <- c("unspInt", "behUnspInt", "contUnspInt")
            numPart <- if (nmStripped %in% intFamilies) "(\\d+)?" else ""
            pattern <- paste0("(^|_)", nmStripped, numPart, "(_[^_]+)*$")
            m <- grep(pattern, vals, perl = TRUE, value = TRUE)
            m <- unique(m)
            if (length(m) == 1L) return(m)
            if (length(m) > 1L && strict) {
                stop("Effect '", nm, "' is ambiguous. Matches: ",
                     paste(m, collapse = ", "), call. = FALSE)
            }
            if (length(m) > 1L && !strict) return(m[1L])
        } else {
            m <- grep(paste0("(^|_)", nmStripped, "$"), vals,
                      perl = TRUE, value = TRUE)
            m <- unique(m)
            if (length(m) == 1L) return(m)
            if (length(m) > 1L && strict) {
                stop("Effect '", nm, "' is ambiguous. Matches: ",
                     paste(m, collapse = ", "), call. = FALSE)
            }
            if (length(m) > 1L && !strict) return(m[1L])
        }

        stop("Effect '", nm, "' not found in available names: ",
             paste(vals, collapse = ", "), call. = FALSE)
    }

    vapply(effectName, resolveOne, character(1L), USE.NAMES = FALSE)
}

# Look up the centering mean for the covariate underlying a given effect.
# Returns the centering mean (numeric) if the effect is covariate-based and
# centered, otherwise 0.  Used to translate user-supplied contrast values
# from the raw (uncentered) scale to the centered scale used internally.
getCovCenteringMean <- function(effectName, effects, data, depvar) {
  inc <- effects[effects$include & effects$type != "rate" &
                   effects$name == depvar, ]
  # Build changeStats base names aligned with inc rows.
  sn <- numberIntShortNames(inc[["shortName"]], key = intKey(inc))
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
contribToChangeStats <- function(contribMat, effectNames, theta = NULL,
                                 direction = NULL) {
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
  # density may be absent in upOnly/downOnly models where it is collinear
  # with the rate and excluded from the specification.
  hasDensity <- !is.na(densityIdx)

  # Unique changeStats names, preserving order of first appearance.
  allBases <- bases
  allTypes <- types
  changeStatsOrder <- unique(allBases)
  nChangeStats <- length(changeStatsOrder)

  # ---- Build csMat: sign-flip + merge eval/endow/creation columns ----
  if (hasDensity) {
    density <- as.integer(contribMat[, densityIdx])
    density[is.na(density)] <- 0L
  } else if (!is.null(direction)) {
    # Caller-supplied direction: upOnly/downOnly models where direction is
    # determined by the constraint, not a density column.
    density <- as.integer(direction)
  } else {
    # No density effect and no direction supplied — cannot determine
    # creation vs dissolution.  Use NA so downstream consumers fail
    # visibly rather than silently producing wrong results.
    density <- rep(NA_integer_, nrow(contribMat))
  }

  # Dissolution-row sign flip is applied per column below during csMat
  # construction (avoids bulk COW copy of contribMat).
  neg <- density == -1L
  anyNeg <- any(neg)

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

    # Extract column(s) from contribMat.  R extracts a fresh vector,
    # so the per-column sign flip below is always in-place (no COW).
    if (hasEval) {
      col <- contribMat[, memberRawIdx[memberTypes == "eval"]]
    } else if (hasCreation && hasEndow) {
      # Creation and endow are structurally zero on the "other" direction
      # rows, so summing recovers the signed Δs for both directions.
      crIdx <- memberRawIdx[memberTypes == "creation"]
      enIdx <- memberRawIdx[memberTypes == "endow"]
      col <- contribMat[, crIdx] + contribMat[, enIdx]
    } else if (hasCreation) {
      col <- contribMat[, memberRawIdx[memberTypes == "creation"]]
    } else if (hasEndow) {
      col <- contribMat[, memberRawIdx[memberTypes == "endow"]]
    }

    # Sign-flip non-density columns on dissolution rows.
    isDensity <- grepl("density", base_k, fixed = TRUE)
    if (!isDensity && anyNeg) {
      col[neg] <- -col[neg]
    }
    csMat[, k] <- col

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
  denschangeStatsIdx <- if (hasDensity) match(bases[densityIdx], changeStatsOrder) else NULL

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
# --------------------------------------------------------------------------
# attachPostestAttrs — attach metadata attributes to a postestimation result.
#
# Stores all metadata fields as attributes and infers conditionCols from the
# output column names. Does NOT assign an S3 class — each caller is
# responsible for setting class() on the returned data.frame.
# --------------------------------------------------------------------------
attachPostestAttrs <- function(result, metadata = NULL) {
  if (!is.null(metadata)) {
    # conditionCols should come from spec metadata, not column-name subtraction.
    if (!is.null(metadata$condition)) {
      cond_cols <- tryCatch(
        resolveEffectName(metadata$condition, names(result)),
        error = function(e) character(0)
      )
      cond_cols <- intersect(cond_cols, names(result))
      if (length(cond_cols) > 0)
        metadata$conditionCols <- cond_cols
    }
    for (nm in names(metadata)) {
      attr(result, nm) <- metadata[[nm]]
    }
  }
  result
}

##@print.sienaPrediction S3 print
summarizeValue <- function(x, na.rm = TRUE){
  if(na.rm){
    x <- x[! is.na(x)]
  }
  qu <- unname(quantile(x, probs = c(.025,0.5,.975)))
  mn <- mean(x)
  se <- sd(x)
  n <- length(x)
  as.list(c("Mean" = mn, "SE" = se, "Median" = qu[2], "q_025" = qu[1], "q_975" = qu[3]))
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


# --------------------------------------------------------------------------
# Memory monitoring utilities
# --------------------------------------------------------------------------

# getRSS: current process Resident Set Size in MB (cross-platform)
getRSS <- function() {
    tryCatch({
        os <- Sys.info()[["sysname"]]
        if (os == "Darwin" || grepl("BSD", os)) {
            # macOS / BSD: ps reports RSS in KB
            txt <- system2("ps", c("-o", "rss=", "-p", Sys.getpid()),
                           stdout = TRUE, stderr = FALSE)
            as.numeric(trimws(txt)) / 1024
        } else if (file.exists("/proc/self/status")) {
            # Linux: VmRSS line in KB
            lines <- readLines("/proc/self/status", warn = FALSE)
            rss_ln <- grep("^VmRSS:", lines, value = TRUE)[1L]
            as.numeric(sub("[^0-9]+", "", rss_ln)) / 1024
        } else {
            NA_real_
        }
    }, error = function(e) NA_real_)
}

# memReport: snapshot of RSS + R heap (via gc), for verbose logging
memReport <- function(label = "", verbose = TRUE) {
    if (!verbose) return(invisible(NULL))
    rss_mb <- getRSS()
    gc_info <- gc(verbose = FALSE, reset = FALSE)
    heap_mb <- round(sum(gc_info[, 2]), 0)  # Mb used (column 2: current)
    rss_str <- if (!is.na(rss_mb)) sprintf("%.0f", rss_mb) else "?"
    message(sprintf("[mem] %s — RSS: %s MB, R heap: %d MB",
                    label, rss_str, heap_mb))
    invisible(list(rss_mb = rss_mb, heap_mb = heap_mb))
}

# assertChainsFreed: verify no live changeContributions in the calling
# frame chain before fork.  Warns (not errors) so production runs proceed.
assertChainsFreed <- function(verbose = FALSE) {
  # Best-effort diagnostics only: never fail the computation.
  get_cc <- function(x) {
    tryCatch({
      if (!is.list(x)) return(NULL)
      x[["changeContributions"]]
    }, error = function(e) NULL)
  }

    for (fi in rev(seq_len(sys.nframe() - 1L))) {
        env <- sys.frame(fi)
        for (nm in ls(env)) {
      tryCatch({
        obj <- get(nm, envir = env, inherits = FALSE)
        cc <- get_cc(obj)
        if (is.list(cc) && length(cc) > 0L) {
          sz_mb <- tryCatch(
            as.numeric(object.size(cc)) / 1e6,
            error = function(e) NA_real_)
          msg <- sprintf(
            paste0("Pre-fork warning: '%s' in frame %d still holds ",
                 "live changeContributions (%d chains, ~%.0f MB). ",
                 "These will be copy-on-write duplicated by fork()."),
            nm, fi, length(cc),
            if (!is.na(sz_mb)) sz_mb else 0)
          if (verbose) message(msg)
          warning(msg, call. = FALSE, immediate. = TRUE)
        }
      }, error = function(e) NULL)
        }
    }
    # Also check .GlobalEnv
    for (nm in ls(.GlobalEnv)) {
    tryCatch({
      obj <- get(nm, envir = .GlobalEnv, inherits = FALSE)
      cc <- get_cc(obj)
      if (is.list(cc) && length(cc) > 0L) {
        sz_mb <- tryCatch(
          as.numeric(object.size(cc)) / 1e6,
          error = function(e) NA_real_)
        msg <- sprintf(
          paste0("Pre-fork warning: '%s' in .GlobalEnv still holds ",
               "live changeContributions (%d chains, ~%.0f MB). ",
               "These will be copy-on-write duplicated by fork()."),
          nm, length(cc),
          if (!is.na(sz_mb)) sz_mb else 0)
        if (verbose) message(msg)
        warning(msg, call. = FALSE, immediate. = TRUE)
      }
    }, error = function(e) NULL)
    }
    invisible(NULL)
}

planBatch <- function(
  data,
  depvar,
  nsim,
  nbrNodes = 1L,
  useCluster = FALSE,
  dynamic = FALSE,
  n3 = NULL,
  unitBudget = 2.5e8,
  dynamicMinistepFactor = 10,
  memoryScale = NULL,
  clusterType = NULL,
  useChangeContributions = FALSE
) {
  dv    <- data$depvars[[depvar]]
  dvDim <- dim(dv)
  nEgo   <- if (!is.null(dvDim) && length(dvDim) >= 1L) as.integer(dvDim[1]) else as.integer(length(dv))
  nChoice <- if (!is.null(dvDim) && length(dvDim) >= 2L) as.integer(dvDim[2]) else 1L
  nPer   <- if (!is.null(dvDim) && length(dvDim) >= 3L) max(1L, as.integer(dvDim[3] - 1L)) else 1L
  units  <- as.numeric(nEgo) * as.numeric(nChoice) * as.numeric(nPer)

  effective_workers <- if (isTRUE(useCluster)) max(1L, as.integer(nbrNodes)) else 1L

  if (dynamic) {
    n3Val        <- if (is.null(n3)) 1L else max(1L, as.integer(n3))
    unitsPerCall <- units * as.numeric(dynamicMinistepFactor) * as.numeric(n3Val)
  } else {
    n3Val        <- 1L
    unitsPerCall <- units
  }

  unitsPerAgg <- max(1.0, units * as.numeric(n3Val) *
                       if (dynamic) as.numeric(dynamicMinistepFactor) else 1.0)

  budgetForAgg <- as.numeric(unitBudget) - effective_workers * unitsPerCall

  if (budgetForAgg <= 0) {
    if (effective_workers > 1L) {
      warning(sprintf(
        "Memory budget (%.0f units) may be insufficient for %d parallel worker(s) at %.0f units each. Consider reducing nbrNodes or increasing batchUnitBudget.",
        unitBudget, effective_workers, unitsPerCall
      ))
    }
    maxBatch <- max(1L, as.integer(floor(as.numeric(unitBudget) / unitsPerAgg)))
  } else {
    maxBatch <- max(1L, as.integer(floor(budgetForAgg / unitsPerAgg)))
  }

  batch_raw <- min(as.integer(nsim), maxBatch)

  if (!is.null(memoryScale) && as.integer(memoryScale) > 1L)
    batch_raw <- max(1L, as.integer(floor(batch_raw / as.integer(memoryScale))))

  k <- if (isTRUE(useCluster) && as.integer(nbrNodes) > 1L) as.integer(nbrNodes) else 1L
  if (k > 1L && as.integer(nsim) >= k) {
    b2        <- (batch_raw %/% k) * k
    batch_raw <- if (b2 < k) k else b2
  }
  min(max(1L, batch_raw), as.integer(nsim))
}