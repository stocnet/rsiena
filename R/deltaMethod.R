## deltaMethod.R — delta-method standard errors for postestimation QoIs
##
## Implements two uncertainty modes (cf. docs/design/postestimation/deltaMethodDesign.md):
##
##   "delta"     — frozen-chain SE: sqrt(grad_cond' Sigma grad_cond)
##                 Channels (A) parameter uncertainty + (B) chain MC noise.
##                 Works for both static and dynamic contributions.
##
##   "deltaFull" — adds Dolby-corrected REINFORCE term (channel C, path-
##                 distribution sensitivity). Dynamic contributions only;
##                 requires ssc_sum from getDynamicChangeContributions(...,
##                 includeScores = TRUE).
##
## Internal functions (not exported):
##   .dolbyRegrCoef          — optimal Dolby baseline coefficients (cf. phase1.r:150-189)
##   .evalSpecScalar         — evaluate spec QoI at theta on frozen wide struct
##   gradCondDelta           — central-FD gradient of Q w.r.t. theta
##   gradReinforceDolby      — Dolby-corrected REINFORCE gradient
##   seDelta                 — delta-method SE from gradient + covariance
##   deltaMethodUncertainty  — orchestrator; returns SE list for all specs
##
## ─────────────────────────────────────────────────────────────────────────────

## ── Orchestrator ─────────────────────────────────────────────────────────────
##
## Compute delta-method uncertainty for all specs.
##
## Per-row Jacobian: for each spec, the conditional FD Jacobian J_cond has
## one row per output row (R rows when the spec is aggregated to level=
## "period"/"ego" or has a `condition`). SE is computed per row as
## sqrt(J_i Σ J_i').
##
## REINFORCE channel (deltaFull): bucketed per-(chain, cell), covering both
## multi-row aggregated specs and accumulated specs.  Where the bucketed path
## is unavailable we issue a one-shot warning and fall back to the conditional
## SE.  That fallback is well-defined but omits the path-distribution channel;
## it is NOT simply a lower bound — the corresponding variance term is a
## covariance, 2*Cov(.,.), and can take either sign.
##
## @param wide         frozen wide struct (returnWide=TRUE from
##                     getDynamicChangeContributions); must contain
##                     `changeStats`, `chain`, `permitted` (can be NULL),
##                     `changeUtility` / `changeProbability` (will be cleared).
## @param ssc_sum      [n3 x nParams] per-chain wave-summed score matrix, or NULL.
##                     Obtained via includeScores=TRUE.  NULL => grad_re = 0.
## @param thetaHat     named numeric vector, length nParams.
## @param covTheta     [nParams x nParams] covariance matrix of thetaHat.
## @param specs        named list of spec entries (same as passed to
##                     makeEstimatorFun; each needs predictFun, predictArgs,
##                     outcomeName).
## @param type         model type string (passed to predictProbability).
## @param fullMode     logical; if TRUE compute Dolby REINFORCE correction
##                     where applicable (requires ssc_sum != NULL).
## @param eps          FD step size.
##
## @return Named list (one entry per spec name), each a list with:
##   J_cond, J_full      — per-row Jacobians [n_out x nParams]
##   SE_delta            — per-row SE (length n_out); conditional channel only
##   SE_deltaFull        — per-row SE (length n_out); cond + REINFORCE where
##                         supported, else equal to SE_delta with attr "fallback"
##   baseline            — Dolby baseline coefficients (scalar specs only; NULL otherwise)
##   ssc_colMeans        — colMeans(ssc_sum) diagnostic; NULL if no ssc_sum
##   Q_hat               — point estimate (scalar for R=1; vector for R>1)
##
deltaMethodUncertainty <- function(wide, estimator, ssc_sum, thetaHat, covTheta,
                                   specs, type,
                                   fullMode = FALSE,
                                   eps = 1e-5,
                                   precomputed = NULL,
                                   verbose = FALSE) {

  nParams     <- length(thetaHat)
  theta_names <- names(thetaHat)

  ## Choose Jacobian method.
  ##   Analytical: one batch sweep, O(n * K) — fast (~7 s for Glasgow dynamic).
  ##               Works for all non-accumulated specs; requires cc$contribMat.
  ##   FD fallback: 2*K batch sweeps, O(n * K * n_sort) — slow (~1200 s).
  ##               Always correct; used when any spec is accumulated,
  ##               rate-weighted, or when the analytical path is explicitly
  ##               disabled.
  has_accumulated <- any(vapply(specs, function(s) isTRUE(s$accumulated),
                                logical(1L)))
  has_rateweight <- any(vapply(specs, function(s) isTRUE(s$rateWeight),
                               logical(1L)))
  if (has_rateweight) {
    if (verbose >= 1) message(
      "rateWeight detected: forcing finite-difference Jacobian because ",
      "analytical rateWeight Jacobian correction is not yet implemented."
    )
  }
  use_analytical  <- !(has_accumulated || has_rateweight)

  # If the caller already ran mode="jacobian" (analytical path), reuse the
  # result directly — avoids a redundant full estimation sweep.
  if (!is.null(precomputed) && use_analytical) {
    jac_result     <- precomputed
    use_analytical <- TRUE
  } else if (use_analytical) {
    jac_result <- tryCatch(
      jacobianCondAnalytical(estimator, thetaHat, specs),
      error = function(e) {
        if (verbose >= 1) message("Analytical Jacobian failed (", conditionMessage(e),
                "); falling back to finite-difference Jacobian.")
        NULL
      })
    if (is.null(jac_result))
      use_analytical <- FALSE
  }
  if (!use_analytical) {
    jac_result <- jacobianCondDelta(estimator, thetaHat, specs, eps)
  } else {
    # evalBatchJacobian returns NULL for specs it cannot handle analytically
    # (e.g. direction-split effects, interactions, predictSecondDiff).
    # Fill those in via FD on just the NULL subset.
    null_specs <- Filter(function(sn) is.null(jac_result$jac[[sn]]),
                         names(jac_result$jac))
    if (length(null_specs) > 0L) {
      if (verbose >= 1) message(
        "Analytical Jacobian unavailable for ", length(null_specs),
        " spec(s); using finite-difference fallback for: ",
        paste(null_specs, collapse = ", ")
      )
      fd_result <- jacobianCondDelta(estimator, thetaHat,
                                     specs[null_specs], eps)
      for (sn in null_specs)
        jac_result$jac[[sn]] <- fd_result$jac[[sn]]
    }
  }

  jac <- jac_result$jac
  hat <- jac_result$hat

  ## ssc_sum availability for REINFORCE.
  have_scores <- !is.null(ssc_sum) && nrow(ssc_sum) > 1L
  reinforce_requested <- isTRUE(fullMode) && have_scores

  ## Names of specs that fell back to the conditional SE; warned once, after
  ## the loop, so the message can name every affected spec rather than only
  ## the first.
  fallback_specs <- character(0)

  results <- setNames(vector("list", length(specs)), names(specs))

  for (specName in names(specs)) {
    spec        <- specs[[specName]]
    outcomeName <- spec$outcomeName

    Q_vec <- hat[[specName]][[outcomeName]]
    n_out <- length(Q_vec)
    J_cond <- jac[[specName]]                      # [n_out x nParams]

    SE_delta <- seDeltaRows(J_cond, covTheta)

    ## Default: full == conditional (with attribute marking fallback path).
    J_full       <- J_cond
    SE_deltaFull <- SE_delta
    baseline     <- NULL
    fallback     <- TRUE

    if (reinforce_requested) {
      ## Unified REINFORCE path for all specs (scalar and multi-row).
      ## .evalSpecChainBuckets handles the n_out==1 case with key_df = one-row
      ## placeholder; gradReinforceRowwise handles [1×R] just as well as [C×R].
      bucket_res <- NULL
      if (!is.null(wide)) {
        cond_resolved_b <- if (!is.null(spec$condition))
          resolveEffectName(spec$condition, wide$changeStats$csNames) else NULL
        gvars_b  <- getGroupVars(level = spec$level, condition = cond_resolved_b)
        key_cols <- intersect(gvars_b, names(hat[[specName]]))
        key_df_b <- if (length(key_cols) > 0L)
          hat[[specName]][, key_cols, drop = FALSE]
        else
          data.frame(.row = seq_len(n_out))
        bucket_res <- tryCatch(
          .evalSpecChainBuckets(spec, wide, thetaHat, type, key_df_b),
          error = function(e) NULL)
      }

      if (!is.null(bucket_res)) {
        ssc_aligned  <- ssc_sum[bucket_res$chain_ids, , drop = FALSE]
        g_score      <- gradReinforceRowwise(bucket_res$Y, bucket_res$N,
                                             ssc_aligned)
        J_full       <- J_cond + g_score
        dimnames(J_full) <- dimnames(J_cond)
        SE_deltaFull <- seDeltaRows(J_full, covTheta)
        ## Dolby baseline diagnostic (scalar specs only; NULL for multi-row).
        if (n_out == 1L) {
          N_c      <- pmax(bucket_res$N[, 1L], 1)
          Q_c      <- bucket_res$Y[, 1L] / N_c
          Qsk_mat  <- outer(Q_c, rep(1, nParams)) * ssc_aligned
          baseline <- .dolbyRegrCoef(Qsk_mat, ssc_aligned)
          names(baseline) <- theta_names
        }
        fallback <- FALSE
      } else {
        fallback_specs <- c(fallback_specs, specName)
      }
    }

    results[[specName]] <- list(
      J_cond       = J_cond,
      J_full       = J_full,
      SE_delta     = SE_delta,
      SE_deltaFull = structure(SE_deltaFull, fallback = fallback),
      baseline     = baseline,
      ssc_colMeans = if (have_scores) colMeans(ssc_sum) else NULL,
      Q_hat        = Q_vec
    )
  }

  ## Deliberately a warning(), not a verbose-gated message(): unlike the
  ## FD-Jacobian fallbacks this returns a *different quantity* — an SE that
  ## omits the path-distribution channel — rather than a slower route to the
  ## same number.  Note it is not simply an under-estimate: the corresponding
  ## variance term is 2*Cov(.,.) and can take either sign.
  if (length(fallback_specs) > 0L)
    warning(
      "uncertaintyMode = 'deltaFull': bucketed REINFORCE unavailable for ",
      "spec(s) ", paste0("'", fallback_specs, "'", collapse = ", "),
      ". Reported SE for these is the conditional (frozen-chain) delta SE, ",
      "which omits the path-distribution channel. Common causes: a custom ",
      "predictFun without outcomesOnly support, or a 'condition' that cannot ",
      "be resolved against the chain data -- note that a conditioned ",
      "ACCUMULATED spec does stratify its point estimate, but the REINFORCE ",
      "term cannot yet be bucketed per stratum, so only its SE falls back ",
      "here. For an SE covering all channels use ",
      "uncertaintyMode = 'bootstrap'.",
      call. = FALSE)

  results
}



# --------------------------------------------------------------------------
# .initDeltaMode — set up chain store for delta / deltaFull modes
#
# Simulates chains once at thetaHat, wraps them in a chainStore_memory, and
# returns an updated contribFun (with keepContribMat = TRUE).  Also builds
# delta_wide when isFullMode = TRUE.
#
# Only called when isDeltaMode = TRUE and dynamic = TRUE; caller guards both.
# --------------------------------------------------------------------------
.initDeltaMode <- function(isFullMode, dynamic, dynArgs,
                           preloadedChains, contribFun, nChainBatches,
                           thetaHat, n3, n3BatchSize) {
    if (!dynamic || is.null(dynArgs))
        return(list(contribFun    = contribFun,
                    nChainBatches = nChainBatches,
                    delta_wide    = NULL,
                    ssc_sum       = NULL))

    # Caller-supplied preloadedChains takes priority over dynArgs$preloadedChains
    if (!is.null(preloadedChains) && is.null(dynArgs$preloadedChains))
        dynArgs$preloadedChains <- preloadedChains

    raw_chains <- dynArgs$preloadedChains
    ssc_sum    <- NULL
    if (is.null(raw_chains)) {
        sim_args                        <- dynArgs
        sim_args$theta                  <- thetaHat
        sim_args$n3                     <- n3
        sim_args$useChangeContributions <- FALSE
        sim_args$returnWide             <- FALSE  # default; explicit for clarity
        sim_args$includeScores          <- isFullMode
        raw_chains <- do.call(getDynamicChangeContributions, sim_args)
        ssc_sum    <- attr(raw_chains, "ssc_sum")
        attr(raw_chains, "ssc_sum") <- NULL
    } else if (isFullMode) {
        stop("uncertaintyMode = 'deltaFull' requires simulated chains with ",
             "includeScores = TRUE; preloaded chains lack ssc_sum.")
    }

    bsz   <- if (!is.null(n3BatchSize)) min(n3BatchSize, length(raw_chains))
             else length(raw_chains)

    # deltaFull must pair Q_chain and ssc_sum from the same simulated chains.
    # Build the wide object from raw_chains directly, rather than re-simulating.
    delta_wide <- NULL
    if (isFullMode) {
      if (is.null(ssc_sum))
        stop("uncertaintyMode = 'deltaFull' requires score statistics.")
      delta_wide <- flattenAndEnrichWide(raw_chains,
                         effects = dynArgs$effects,
                         depvar  = dynArgs$depvar,
                         data    = dynArgs$data)
      delta_wide$ssc_sum <- ssc_sum
    }

    store <- chainStore_memory(raw_chains, bsz)
    rm(raw_chains)

    contribFun <- makeContribFun(store          = store,
                                 effects        = dynArgs$effects,
                                 depvar         = dynArgs$depvar,
                                 data           = dynArgs$data,
                                 keepContribMat = TRUE)
    nChainBatches <- store$nBatches

    list(contribFun    = contribFun,
         nChainBatches = nChainBatches,
         delta_wide    = delta_wide,
         ssc_sum       = ssc_sum)
}

## ── accumulatedByChain — per-(chain, cell) sum/count for accumulated specs ───
##
## Seam between the accumulated aggregation and the bucketed REINFORCE
## consumer.  Returns the same shape aggWithCache() returns for the
## non-accumulated path:
##
##     chain, <final group cols...>, <outcome>_sum, <outcome>_n
##
## Accumulated aggregation already groups by chain when summing over ministeps
## (`acc_group` in aggAccumulatedSumCount); passing extraGroup = "chain" simply
## stops the final stage from collapsing that away, so the per-chain
## accumulated values survive to the output.
##
## NOTE: this is the seam the planned aggregation unification replaces (see
## dev-notes/accumulated-deltafull-patch-plan.md Sec. 4) — the signature is
## intended to outlive the body.
##
## @param spec       spec entry (outcomeName, level, na.rm).
## @param structural structural columns for kept rows (chain, group, period,
##                   ego, ministep).
## @param vals       outcome values for kept rows (same length/order).
## @return data.frame, or NULL if the accumulated path is unavailable.
accumulatedByChain <- function(spec, structural, vals) {
  acc_data <- structural
  acc_data[[spec$outcomeName]] <- vals
  attr(acc_data, "row.names") <- .set_row_names(length(structural[[1L]]))
  class(acc_data) <- "data.frame"
  tryCatch(
    aggAccumulatedSumCount(spec$outcomeName, acc_data,
                           level      = spec$level,
                           condition  = NULL,   # unsupported; main path warns
                           na.rm      = isTRUE(spec$na.rm),
                           extraGroup = "chain"),
    error = function(e) NULL)
}

## ── Chain-bucket evaluation for bucketed REINFORCE ───────────────────────────
##
## .evalSpecChainBuckets — evaluate spec on delta_wide at thetaHat and
## aggregate by (chain, cell) to obtain [n_chains × R] sum/count matrices.
##
## Non-accumulated specs use the same buildAggCache / aggWithCache machinery as
## the main estimation path, with "chain" prepended to group_vars.  This means
## egoNormalize is handled identically: Stage 1 pre-aggregates by
##   unique(c(ego_id_cols, c("chain", gvars)))   (= chain-scoped ego means)
## and Stage 2 collapses to (chain, cell).
##
## Accumulated specs route through accumulatedByChain() instead (below), which
## sums over ministeps within (chain, ego) before the final sum/count.  Both
## producers return the same shape, so everything downstream is shared.
##
## Returns NULL if predictFun returns a data.frame (no outcomesOnly support),
## or required columns are absent.
##
## @param spec     spec entry (predictFun, predictArgs, outcomeName, level,
##                 condition, egoNormalize, accumulated).
## @param wide     delta_wide from flattenAndEnrichWide; must have chain,
##                 period, ego, permitted, changeStats.
## @param theta    named numeric (thetaHat).
## @param type     model type string passed to predictProbability.
## @param key_df   data.frame [R × k] of cell key columns in the row order
##                 matching Q_hat_r (from hat[[specName]][, gvars, drop=FALSE]).
##
## @return list(Y, N, chain_ids) or NULL.
##   Y, N        — [n_chains × R] double matrices (row = chain, col = cell).
##   chain_ids   — integer vector of chain labels mapping Y rows to ssc_sum rows.

.evalSpecChainBuckets <- function(spec, wide, theta, type, key_df) {
  if (is.null(wide$chain))      return(NULL)

  R <- nrow(key_df)
  if (R == 0L) return(NULL)

  ## Evaluate predictFun with outcomesOnly on delta_wide at theta.
  wide_eval                   <- wide
  wide_eval$changeUtility     <- NULL
  wide_eval$changeProbability <- NULL
  baseline <- predictProbability(wide_eval, theta, type, returnComponents = TRUE)
  extra    <- c(list(changeContributions = wide_eval,
                     theta               = theta,
                     baseline            = baseline,
                     outcomesOnly        = TRUE),
                spec$predictArgs)
  outcomes <- tryCatch(do.call(spec$predictFun, extra), error = function(e) NULL)
  if (is.null(outcomes) || is.data.frame(outcomes)) return(NULL)
  ov <- outcomes[[spec$outcomeName]]
  if (is.null(ov) || length(ov) == 0L) return(NULL)

  ## Keep mask — same logic as .evalSpecScalar.
  density <- wide$changeStats$density
  perm    <- wide$permitted
  keep    <- density != 0L & (if (is.null(perm)) TRUE else perm)

  ## Structural frame (kept rows) — includes chain, group, period, ego, ministep.
  structural_k <- groupColsList(wide, keep = which(keep))

  ## Add condition column from changeStats$csMat if needed.
  cond_resolved <- if (!is.null(spec$condition))
    resolveEffectName(spec$condition, wide$changeStats$csNames) else NULL
  if (!is.null(cond_resolved)) {
    cm <- wide$changeStats$csMat
    if (is.null(cm) || !all(cond_resolved %in% colnames(cm))) return(NULL)
    for (cn in cond_resolved)
      structural_k[[cn]] <- cm[, cn][keep]
  }

  ## Group vars for output cells; prepend "chain" to get per-(chain, cell) agg.
  gvars       <- getGroupVars(level = spec$level, condition = cond_resolved)
  gvars_chain <- c("chain", gvars)
  ego_id_cols <- detectEgoUnit(structural_k)
  egoNorm     <- isTRUE(spec$egoNormalize)
  na_rm       <- isTRUE(spec$na.rm)

  ## Two producers, one consumer: everything below only needs bc_agg to carry
  ## chain, <gvars>, <outcome>_sum, <outcome>_n.
  ##
  ## Accumulated: sum over ministeps within (chain, ego) before the final
  ## sum/count, with chain retained (see accumulatedByChain above).
  ## Otherwise: one buildAggCache call, which handles egoNormalize correctly —
  ## pre_agg_vars = unique(c(ego_id_cols, gvars_chain)) groups by all ego-id
  ## columns first (chain-scoped ego means), then collapses to (chain, cell).
  bc_agg <- if (isTRUE(spec$accumulated)) {
    accumulatedByChain(spec, structural_k, ov[keep])
  } else {
    cache_bc <- tryCatch(
      buildAggCache(structural_k, group_vars = gvars_chain,
                    ego_id_cols = ego_id_cols,
                    egoNormalize = egoNorm, na.rm = na_rm),
      error = function(e) NULL)
    if (is.null(cache_bc)) NULL
    else aggWithCache(spec$outcomeName, ov[keep], cache_bc)
  }
  if (is.null(bc_agg)) return(NULL)

  sumCol <- paste0(spec$outcomeName, "_sum")
  cntCol <- paste0(spec$outcomeName, "_n")
  if (!all(c(sumCol, cntCol, "chain") %in% names(bc_agg))) return(NULL)

  chains_all <- sort(unique(bc_agg$chain))
  n_chains   <- length(chains_all)
  if (n_chains == 0L) return(NULL)

  ## Match bc_agg rows to (chain_idx, cell_idx) for matrix filling.
  chain_idx_m <- match(bc_agg$chain, chains_all)

  if (length(gvars) > 0L) {
    ref_parts <- lapply(gvars, function(v) key_df[[v]])
    if (any(vapply(ref_parts, is.null, logical(1L)))) return(NULL)
    ref_keys <- if (length(ref_parts) == 1L) as.character(ref_parts[[1L]])
                else do.call(paste, c(ref_parts, list(sep = "\r")))
    row_parts <- lapply(gvars, function(v) bc_agg[[v]])
    if (any(vapply(row_parts, is.null, logical(1L)))) return(NULL)
    row_keys <- if (length(row_parts) == 1L) as.character(row_parts[[1L]])
                else do.call(paste, c(row_parts, list(sep = "\r")))
    cell_idx <- match(row_keys, ref_keys)
  } else {
    cell_idx <- rep(1L, nrow(bc_agg))
  }

  valid_m <- !is.na(cell_idx)
  if (!any(valid_m)) return(NULL)

  ## Accumulate rather than assign: bc_agg may carry finer grouping columns
  ## than gvars (the accumulated path retains "group" from its level switch),
  ## in which case several bc_agg rows map to the same (chain, cell).  Summing
  ## sums and counts is the correct collapse; assigning would silently keep
  ## only the last row.  For the non-accumulated path there is exactly one row
  ## per (chain, cell), so this is equivalent to the previous behaviour.
  Y <- matrix(0, nrow = n_chains, ncol = R)
  N <- matrix(0, nrow = n_chains, ncol = R)
  for (i in which(valid_m)) {
    Y[chain_idx_m[i], cell_idx[i]] <- Y[chain_idx_m[i], cell_idx[i]] +
                                        bc_agg[[sumCol]][i]
    N[chain_idx_m[i], cell_idx[i]] <- N[chain_idx_m[i], cell_idx[i]] +
                                        bc_agg[[cntCol]][i]
  }

  list(Y = Y, N = N, chain_ids = as.integer(chains_all))
}


## ── Per-row Dolby REINFORCE gradient (bucketed multi-cell specs) ─────────────
##
## For each output cell r, computes the Dolby-corrected REINFORCE gradient:
##   A_c_r    = (Y_c_r − Q_hat_r * N_c_r) / mean_c(N_c_r)
##   g_score_r = gradReinforceDolby(A_c_r, ssc_aligned)
##
## @param Y           [n_chains × R] per-chain cell sum matrix.
## @param N           [n_chains × R] per-chain cell count matrix.
## @param ssc_aligned [n_chains × nParams] score matrix (rows match Y rows).
##
## @return [R × nParams] score gradient matrix (unnamed columns).
gradReinforceRowwise <- function(Y, N, ssc_aligned) {
  R         <- ncol(Y)
  n_chains  <- nrow(Y)
  nParams   <- ncol(ssc_aligned)

  Q_hat_r  <- colSums(Y, na.rm = TRUE) / pmax(colSums(N, na.rm = TRUE), 1)
  mean_N_r <- colMeans(N, na.rm = TRUE)
  ## Avoid division by zero for empty cells.
  mean_N_r[mean_N_r == 0] <- 1

  ## A_c_r = (Y_c_r − Q_hat_r * N_c_r) / mean_c(N_c_r)   [n_chains × R]
  A <- (Y - outer(rep(1, n_chains), Q_hat_r) * N) /
         outer(rep(1, n_chains), mean_N_r)

  g_score <- matrix(0, nrow = R, ncol = nParams)
  for (r in seq_len(R))
    g_score[r, ] <- gradReinforceDolby(A[, r], ssc_aligned)
  g_score
}


## ── Dolby regression coefficient ─────────────────────────────────────────────
##
## Optimal scalar baseline for the k-th parameter:
##   baseline[k] = Cov(stat[, k], scores[, k]) / Var(scores[, k])
##
## `stat`   — [n3 x nParams] matrix (e.g. Q_c * s_k^{(c)}, one col per parameter)
## `scores` — [n3 x nParams] matrix of per-chain score contributions
##
## Lifted from phase1.r:150-189.  NA / zero-variance guards match the original.
.dolbyRegrCoef <- function(stat, scores) {
  nParams  <- ncol(scores)
  baseline <- numeric(nParams)
  for (i in seq_len(nParams)) {
    v <- var(scores[, i])
    if (is.finite(v) && v > 0)
      baseline[i] <- cov(stat[, i], scores[, i]) / v
  }
  baseline[!is.finite(baseline)] <- 0
  baseline
}


## ── Evaluate a spec's QoI at `theta` on frozen wide struct ───────────────────
##
## Returns list(Q_scalar, Q_chain):
##   Q_scalar — scalar mean(fd[keep], na.rm=TRUE)
##   Q_chain  — per-chain mean; named numeric vector; names are chain indices
##
## IMPORTANT: `wide$changeUtility` / `wide$changeProbability` are cleared so
## that `predictProbability` recomputes R-side from `changeStats$csMat` at the
## supplied theta.  This is required for FD gradient correctness: the C++
## pre-computed values were stored at theta-hat and must not be reused when
## theta is perturbed.
.evalSpecScalar <- function(spec, wide, theta, type) {
  ## Clear C++ pre-computed probabilities — force R-side recomputation at theta.
  wide_eval <- wide
  wide_eval$changeUtility     <- NULL
  wide_eval$changeProbability <- NULL

  baseline <- predictProbability(wide_eval, theta, type,
                                  returnComponents = TRUE)
  extra    <- c(list(changeContributions = wide_eval,
                     theta               = theta,
                     baseline            = baseline,
                     outcomesOnly        = TRUE),
                spec$predictArgs)
  outcomes <- do.call(spec$predictFun, extra)

  ## Extract outcome vector (outcomesOnly = TRUE returns a named list).
  ov <- if (is.list(outcomes) && !is.data.frame(outcomes))
    outcomes[[spec$outcomeName]]
  else
    outcomes[[spec$outcomeName]]

  ## Keep mask: exclude density=0 (no-change) rows and non-permitted rows.
  density <- wide$changeStats$density
  perm    <- wide$permitted
  keep    <- density != 0L & (if (is.null(perm)) TRUE else perm)
  ov_k    <- ov[keep]
  chain_k <- wide$chain[keep]

  Q_scalar <- mean(ov_k, na.rm = TRUE)
  Q_chain  <- tapply(ov_k, chain_k, mean, na.rm = TRUE)
  list(Q_scalar = Q_scalar, Q_chain = as.numeric(Q_chain),
       chain_ids = as.integer(names(Q_chain)))
}


## ── Central finite-difference Jacobian (all specs at once) ──────────────────
##
## For each spec, returns the [n_out x nParams] Jacobian
##
##   J[i, k] = (outcome[i](theta+eps_k) - outcome[i](theta-eps_k)) / (2*eps)
##
## Builds the full 2*nParams FD grid as a `perturbations` named list and
## calls `estimator(thetaHat, perturbations = ...)` ONCE.  Inside the
## estimator each batch b is loaded once; `prepBatchContext(cc_b)` is called
## once; `evalBatch` is called (2*nParams + 1) times on the same cc_b with
## different theta values.  This means θ-independent structural caches are
## built once per batch regardless of the number of parameters.
##
## `estimator` — the closure returned by makeEstimatorFun.  Must accept the
##               `perturbations` argument added in the Commit-1 refactor.
## `thetaHat`  — named numeric vector of length nParams.
## `specs`     — named list of specs (each needs outcomeName field).
## `eps`       — finite-difference step.
##
## Returns: list(jac, hat)
##   jac — named list (one per spec) of [n_out x nParams] matrices.
##   hat — named list of data.frames at thetaHat (reused by caller).
##
jacobianCondDelta <- function(estimator, thetaHat, specs, eps = 1e-5) {
  nParams     <- length(thetaHat)
  spec_names  <- names(specs)
  theta_names <- names(thetaHat)

  ## Build the full FD grid as a named perturbations list.
  ## Names: "p_1".."p_nParams" for +eps, "m_1".."m_nParams" for -eps.
  perts <- vector("list", 2L * nParams)
  pnames <- character(2L * nParams)
  for (k in seq_len(nParams)) {
    tp <- thetaHat; tp[k] <- tp[k] + eps
    tm <- thetaHat; tm[k] <- tm[k] - eps
    perts[[2L * k - 1L]] <- tp
    perts[[2L * k]]      <- tm
    pnames[[2L * k - 1L]] <- paste0("p_", k)
    pnames[[2L * k]]      <- paste0("m_", k)
  }
  names(perts) <- pnames

  ## Single estimator call: one batch sweep, all FD evals on shared cc_b.
  ev <- estimator(thetaHat, perturbations = perts)

  hat <- ev$hat

  ## Assemble Jacobian matrices from the returned perturbation results.
  jac <- setNames(vector("list", length(spec_names)), spec_names)
  for (specName in spec_names) {
    oc    <- specs[[specName]]$outcomeName
    n_out <- length(hat[[specName]][[oc]])
    J     <- matrix(NA_real_, nrow = n_out, ncol = nParams,
                    dimnames = list(NULL, theta_names))
    for (k in seq_len(nParams)) {
      ep_oc <- ev$perturbations[[paste0("p_", k)]][[specName]][[oc]]
      em_oc <- ev$perturbations[[paste0("m_", k)]][[specName]][[oc]]
      J[, k] <- (ep_oc - em_oc) / (2 * eps)
    }
    jac[[specName]] <- J
  }

  list(jac = jac, hat = hat)
}


## ── Analytical per-row Jacobian (non-accumulated specs only) ─────────────────
##
## Calls estimator(thetaHat, mode = "jacobian") which executes one batch sweep
## and returns J[i, k] = ∂Q_i/∂θ_k analytically via the softmax derivative.
##
## Requires keepContribMat = TRUE in the contribFun closure (set automatically
## by sienaPostestimate when isDeltaMode = TRUE).
##
## Returns: list(jac, hat) — same shape as jacobianCondDelta.
##   jac — named list (one per spec) of [n_out × nParams] matrices.
##              Accumulated specs → NULL (caller should fall back to FD).
##   hat — named list of data.frames at thetaHat.
##
jacobianCondAnalytical <- function(estimator, thetaHat, specs) {
  ev       <- estimator(thetaHat, mode = "jacobian")
  hat <- ev$hat
  jac_raw  <- ev$jac   # named list, each [n_out × nParams] or NULL

  jac <- setNames(vector("list", length(specs)), names(specs))
  for (specName in names(specs))
    jac[specName] <- list(jac_raw[[specName]])

  list(jac = jac, hat = hat)
}


## ── Dolby-corrected REINFORCE gradient ───────────────────────────────────────
##
## grad_re[k] = mean(Q_c * s_k^{(c)}) - baseline_k * mean(s_k^{(c)})
##            = E[Q S_k] - baseline_k * E[S_k]
##
## where baseline_k = Cov(Q_c * s_k, s_k) / Var(s_k)   (Dolby optimal baseline)
##
## `Q_chain`     — per-chain Q values, length n3 (indexed by chain_ids)
## `ssc_aligned` — [n3_eff x nParams] aligned score matrix (rows match Q_chain)
##
gradReinforceDolby <- function(Q_chain, ssc_aligned) {
  nParams <- ncol(ssc_aligned)
  Qv      <- as.numeric(Q_chain)
  vapply(seq_len(nParams), function(k) {
    sk         <- ssc_aligned[, k]
    Qsk        <- Qv * sk
    v          <- var(sk)
    baseline_k <- if (is.finite(v) && v > 0) cov(Qsk, sk) / v else 0
    mean(Qsk) - baseline_k * mean(sk)
  }, numeric(1L))
}


## ── Delta-method SE (per-row, vectorised) ──────────────────────────────────
##
## Input J: [n_out x nParams] Jacobian.
## Output: numeric vector of length n_out with SE_i = sqrt(J[i,] %*% Σ %*% J[i,]).
##
## Computed as rowSums((J %*% Σ) * J) — equivalent to diag(J %*% Σ %*% t(J))
## but O(n_out * nParams^2) instead of O(n_out^2 * nParams).
##
seDeltaRows <- function(J, covTheta) {
  M <- J %*% covTheta            # [n_out x nParams]
  v <- rowSums(M * J)
  sqrt(pmax(v, 0))               # guard against tiny negative from FP
}

## Scalar-form SE (kept for REINFORCE Dolby path which still works on scalars).
seDelta <- function(g, covTheta) {
  as.numeric(sqrt(t(g) %*% covTheta %*% g))
}