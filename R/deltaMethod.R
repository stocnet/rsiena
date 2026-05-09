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


## ── Orchestrator ─────────────────────────────────────────────────────────────
##
## Compute delta-method uncertainty for all specs.
##
## Per-row Jacobian: for each spec, the conditional FD Jacobian J_cond has
## one row per output row (R rows when the spec is aggregated to level=
## "period"/"ego" or has a `condition`). SE is computed per row as
## sqrt(J_i Σ J_i').
##
## REINFORCE channel (deltaFull): currently implemented only for scalar
## specs (R = 1, i.e. level="none" without condition). For multi-row
## aggregated specs we issue a one-shot warning and fall back to the
## conditional SE — this is well-defined and conservative (lower bound on
## SE for that spec). Per-bucket REINFORCE is a planned extension.
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
                                   precomputed = NULL) {

  nParams     <- length(thetaHat)
  theta_names <- names(thetaHat)

  ## Choose Jacobian method.
  ##   Analytical: one batch sweep, O(n * K) — fast (~7 s for Glasgow dynamic).
  ##               Works for all non-accumulated specs; requires cc$contribMat.
  ##   FD fallback: 2*K batch sweeps, O(n * K * n_sort) — slow (~1200 s).
  ##               Always correct; used when any spec is accumulated or when
  ##               the analytical path is explicitly disabled.
  has_accumulated <- any(vapply(specs, function(s) isTRUE(s$accumulated),
                                logical(1L)))
  use_analytical  <- !has_accumulated

  # If the caller already ran mode="jacobian" (analytical path), reuse the
  # result directly — avoids a redundant full estimation sweep.
  if (!is.null(precomputed)) {
    jac_result     <- precomputed
    use_analytical <- TRUE
  } else if (use_analytical) {
    jac_result <- tryCatch(
      jacobianCondAnalytical(estimator, thetaHat, specs),
      error = function(e) {
        warning("Analytical Jacobian failed (", conditionMessage(e),
                "); falling back to finite-difference Jacobian.", call. = FALSE)
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
      warning(
        "Analytical Jacobian unavailable for ", length(null_specs),
        " spec(s); using finite-difference fallback for: ",
        paste(null_specs, collapse = ", "),
        call. = FALSE
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

  ## One-shot warning state for fallback on multi-row specs.
  fallback_warned <- FALSE

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
      if (n_out == 1L) {
        ## Scalar spec: full Dolby REINFORCE on per-chain Q.
        hat_chain <- .evalSpecScalar(spec, wide, thetaHat, type)
        Q_chain   <- hat_chain$Q_chain
        chain_ids <- hat_chain$chain_ids
        if (length(Q_chain) > 1L && length(chain_ids) > 0L) {
          ssc_aligned <- ssc_sum[chain_ids, , drop = FALSE]
          grad_cond   <- as.numeric(J_cond[1L, ])
          names(grad_cond) <- theta_names

          grad_re      <- gradReinforceDolby(Q_chain, ssc_aligned)
          names(grad_re) <- theta_names
          grad_full    <- grad_cond + grad_re
          SE_deltaFull <- seDelta(grad_full, covTheta)

          J_full       <- matrix(grad_full, nrow = 1L,
                                 dimnames = list(NULL, theta_names))

          Qsk_mat  <- outer(Q_chain, rep(1, nParams)) * ssc_aligned
          baseline <- .dolbyRegrCoef(Qsk_mat, ssc_aligned)
          names(baseline) <- theta_names
          fallback <- FALSE
        }
      } else {
        ## Multi-row aggregated spec: per-bucket REINFORCE not yet implemented.
        ## Fall back to conditional SE and warn (one-shot).
        ## TODO(deltaFull-bucketed): implement per-bucket Dolby REINFORCE for
        ## level={period,ego} and condition=* specs by computing a
        ## [n3 x R] per-chain × per-bucket Y matrix at thetaHat, with
        ## row-count weighting matching the Q aggregation. Fall back for
        ## accumulated/egoNormalize specs (different inner reduction).
        if (!fallback_warned) {
          warning(
            "uncertaintyMode = 'deltaFull' is not yet implemented for ",
            "aggregated specs (level != 'none' or with `condition`). ",
            "Falling back to conditional delta SE for those specs. ",
            "Note: conditional SE only captures parameter uncertainty + ",
            "frozen-chain noise; it does NOT include path-distribution ",
            "sensitivity. For publication-quality SE on aggregated outputs, ",
            "use uncertaintyMode = 'batch' (parametric bootstrap).",
            call. = FALSE)
          fallback_warned <- TRUE
        }
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

  results
}

