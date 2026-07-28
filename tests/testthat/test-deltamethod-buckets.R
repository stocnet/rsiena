## test-deltamethod-buckets.R
##
## Tests for the unified bucketed REINFORCE path in deltaMethodUncertainty.
##
## The key invariant: for a scalar spec (n_out==1) with constant per-chain
## counts (equal N_c across all chains), the SE returned by
## deltaMethodUncertainty (which uses .evalSpecChainBuckets +
## gradReinforceRowwise) must be identical (to machine precision) to the SE
## that the old .evalSpecScalar + gradReinforceDolby path would have returned.
##
## This is the "before/after" check: the reference ("before") is computed
## explicitly in the test using .evalSpecScalar + gradReinforceDolby; the
## "after" is the SE returned by deltaMethodUncertainty via the unified path.
## Algebraic equality for equal N_c is exact, so tolerance is 1e-12.

context("deltaFull bucketed REINFORCE — orchestrator-level parity")

make_wide_dyn_bucket <- function(n_chains = 3L, n_egos = 4L, seed = 42L) {
  set.seed(seed)
  n_alts  <- 3L
  n_total <- n_chains * n_egos * n_alts

  chain_vec <- rep(seq_len(n_chains), each = n_egos * n_alts)
  ego_vec   <- rep(rep(seq_len(n_egos), each = n_alts), times = n_chains)
  ms_vec    <- ego_vec

  density_col <- rep(c(1L, -1L, 0L), times = n_chains * n_egos)
  recip_col   <- sample(c(0L, 1L), n_total, replace = TRUE)
  recip_col[density_col == 0L] <- 0L

  contribMat <- matrix(c(density_col, recip_col), ncol = 2L,
                       dimnames = list(NULL, c("density", "recip")))

  group_id <- as.integer(factor(paste(chain_vec, ego_vec, sep = "_")))

  wide <- list(
    chain             = chain_vec,
    group             = rep(1L, n_total),
    period            = rep(1L, n_total),
    ego               = ego_vec,
    ministep          = ms_vec,
    choice            = density_col,
    group_id          = group_id,
    permitted         = NULL,
    contribMat        = contribMat,
    effectNames       = c("density", "recip"),
    changeUtility     = NULL,
    changeProbability = NULL,
    changeStats       = NULL
  )
  wide$changeStats <- contribToChangeStats(wide$contribMat, wide$effectNames)
  wide
}

## Deterministic mock predictFun — outcome depends only on change stats.
mock_pred_bucket <- function(changeContributions, theta, baseline = NULL,
                             outcomesOnly = FALSE, ...) {
  cs <- changeContributions$changeStats
  ov <- as.numeric(cs$density * 0.3 + cs$csMat[, "recip"] * 0.1)
  if (outcomesOnly) return(list(me_out = ov))
  data.frame(me_out = ov)
}

make_scalar_spec <- function() {
  list(
    predictFun   = mock_pred_bucket,
    predictArgs  = list(),
    outcomeName  = "me_out",
    level        = "none",
    condition    = NULL,
    egoNormalize = FALSE,
    na.rm        = FALSE,
    accumulated  = FALSE
  )
}

make_fd_estimator_scalar <- function(spec_name = "spec", outcome_name = "me_out") {
  function(theta, perturbations = NULL,
           useChangeContributions = FALSE, mode = "outcome") {
    eval_at <- function(th) {
      val <- th[["density"]] + 2 * th[["recip"]]
      setNames(list(data.frame(me_out = val)), spec_name)
    }
    hat <- eval_at(theta)
    if (is.null(perturbations)) return(hat)
    perts_out <- lapply(perturbations, eval_at)
    list(hat = hat, perturbations = perts_out)
  }
}

test_that("deltaMethodUncertainty scalar deltaFull SE matches old .evalSpecScalar path", {
  skip_on_cran()

  wide     <- make_wide_dyn_bucket()
  spec     <- make_scalar_spec()
  theta    <- c(density = -2, recip = 1)
  type     <- "changeProb"
  nParams  <- length(theta)
  n_chains <- length(unique(wide$chain))

  set.seed(77L)
  ssc_sum <- matrix(rnorm(n_chains * nParams), nrow = n_chains, ncol = nParams,
                    dimnames = list(NULL, names(theta)))

  ## Fixed J_cond and covariance.
  J_cond   <- matrix(c(0.1, -0.05), nrow = 1L,
                     dimnames = list(NULL, names(theta)))
  covTheta <- diag(c(0.01, 0.02))
  Q_hat    <- 0.25

  ## Inject known Jacobian and hat via precomputed — bypasses estimator call.
  precomp <- list(
    jac = list(myspec = J_cond),
    hat = list(myspec = data.frame(me_out = Q_hat))
  )

  res <- deltaMethodUncertainty(
    wide        = wide,
    estimator   = NULL,
    ssc_sum     = ssc_sum,
    thetaHat    = theta,
    covTheta    = covTheta,
    specs       = list(myspec = spec),
    type        = type,
    fullMode    = TRUE,
    precomputed = precomp
  )

  expect_false(attr(res$myspec$SE_deltaFull, "fallback"),
    info = "bucketed path must succeed for this valid scalar spec")

  old_path    <- .evalSpecScalar(spec, wide, theta, type)
  ssc_aligned <- ssc_sum[old_path$chain_ids, , drop = FALSE]
  grad_re_old <- gradReinforceDolby(old_path$Q_chain, ssc_aligned)
  J_full_old  <- matrix(as.numeric(J_cond[1L, ]) + grad_re_old, nrow = 1L,
                        dimnames = dimnames(J_cond))
  SE_old      <- seDeltaRows(J_full_old, covTheta)

  expect_equal(as.numeric(res$myspec$SE_deltaFull), SE_old, tolerance = 1e-12,
    info = "unified bucketed path must give same SE as old evalSpecScalar path")
})

test_that("deltaMethodUncertainty: accumulated spec uses bucketed REINFORCE", {
  skip_on_cran()

  wide  <- make_wide_dyn_bucket()
  spec  <- make_scalar_spec()
  spec$accumulated <- TRUE
  theta    <- c(density = -2, recip = 1)
  n_chains <- length(unique(wide$chain))
  est <- make_fd_estimator_scalar(spec_name = "acc")
  set.seed(77L)
  ssc_sum  <- matrix(rnorm(n_chains * 2L), nrow = n_chains,
                     dimnames = list(NULL, names(theta)))

  expect_no_warning(
    res <- deltaMethodUncertainty(
      wide = wide, estimator = est, ssc_sum = ssc_sum,
      thetaHat = theta, covTheta = diag(c(0.01, 0.02)),
      specs = list(acc = spec), type = "changeProb",
      fullMode = TRUE
    )
  )

  expect_false(attr(res$acc$SE_deltaFull, "fallback"),
    info = "accumulated spec must now take the bucketed REINFORCE path")
  expect_false(isTRUE(all.equal(as.numeric(res$acc$SE_deltaFull),
                                as.numeric(res$acc$SE_delta))),
    info = "deltaFull SE must differ from the conditional SE once REINFORCE runs")
})

test_that("accumulatedByChain: Y/N reconstruct the accumulated point estimate", {
  skip_on_cran()

  wide <- make_wide_dyn_bucket()
  spec <- make_scalar_spec()
  spec$accumulated <- TRUE
  theta <- c(density = -2, recip = 1)

  ## Point estimate via the ordinary accumulated aggregation.
  keep <- wide$changeStats$density != 0L &
          (if (is.null(wide$permitted)) TRUE else wide$permitted)
  wide_eval                   <- wide
  wide_eval$changeUtility     <- NULL
  wide_eval$changeProbability <- NULL
  baseline <- RSiena:::predictProbability(wide_eval, theta, "changeProb",
                                          returnComponents = TRUE)
  outcomes <- do.call(spec$predictFun,
    c(list(changeContributions = wide_eval, theta = theta,
           baseline = baseline, outcomesOnly = TRUE), spec$predictArgs))
  ov <- outcomes[[spec$outcomeName]]

  structural_k <- RSiena:::groupColsList(wide, keep = which(keep))
  acc_data <- structural_k
  acc_data[[spec$outcomeName]] <- ov[keep]
  attr(acc_data, "row.names") <- .set_row_names(length(structural_k[[1L]]))
  class(acc_data) <- "data.frame"

  ref <- RSiena:::aggAccumulatedSumCount(spec$outcomeName, acc_data,
                                         level = spec$level, na.rm = TRUE)
  Q_ref <- sum(ref[[paste0(spec$outcomeName, "_sum")]]) /
           sum(ref[[paste0(spec$outcomeName, "_n")]])

  ## Same quantity via the per-chain bucket path.
  bc <- RSiena:::accumulatedByChain(spec, structural_k, ov[keep])
  expect_false(is.null(bc), info = "accumulatedByChain must not return NULL")
  expect_true("chain" %in% names(bc),
    info = "chain must survive to the bucket output")

  Q_bucket <- sum(bc[[paste0(spec$outcomeName, "_sum")]]) /
              sum(bc[[paste0(spec$outcomeName, "_n")]])

  expect_equal(Q_bucket, Q_ref,
    info = "per-chain buckets must reconstruct the accumulated point estimate")

  ## Retaining chain must refine, never coarsen: at least as many rows.
  expect_gte(nrow(bc), nrow(ref))
})

## ── 3. No chain column → falls back to conditional SE ────────────────────────

test_that("deltaMethodUncertainty: no chain column forces fallback", {
  skip_on_cran()

  wide       <- make_wide_dyn_bucket()
  wide$chain <- NULL          # simulate static wide
  spec  <- make_scalar_spec()
  theta <- c(density = -2, recip = 1)
  n_chains <- 3L
  est <- make_fd_estimator_scalar(spec_name = "s")
  set.seed(77L)
  ssc_sum <- matrix(rnorm(n_chains * 2L), nrow = n_chains,
                    dimnames = list(NULL, names(theta)))

  expect_warning(
    res <- deltaMethodUncertainty(
      wide = wide, estimator = est, ssc_sum = ssc_sum,
      thetaHat = theta, covTheta = diag(c(0.01, 0.02)),
      specs = list(s = spec), type = "changeProb",
      fullMode = TRUE
    ),
    regexp = "bucketed REINFORCE unavailable for spec\\(s\\) 's'"
  )

  expect_true(attr(res$s$SE_deltaFull, "fallback"))
})

test_that("deltaFull fallback warning names every affected spec", {
  skip_on_cran()

  wide       <- make_wide_dyn_bucket()
  wide$chain <- NULL                     # forces fallback for all specs
  theta      <- c(density = -2, recip = 1)
  set.seed(5L)
  ssc_sum <- matrix(rnorm(3L * 2L), nrow = 3L,
                    dimnames = list(NULL, names(theta)))

  spec_a <- make_scalar_spec(); spec_b <- make_scalar_spec()
  est <- make_fd_estimator_scalar(spec_name = "alpha")

  w <- tryCatch(
    deltaMethodUncertainty(
      wide = wide, estimator = est, ssc_sum = ssc_sum,
      thetaHat = theta, covTheta = diag(c(0.01, 0.02)),
      specs = list(alpha = spec_a, beta = spec_b),
      type = "changeProb", fullMode = TRUE),
    warning = function(w) conditionMessage(w))

  expect_true(is.character(w), info = "a warning must be raised")
  expect_match(w, "'alpha'")
  expect_match(w, "'beta'")
  expect_match(w, "uncertaintyMode = 'bootstrap'")
  expect_false(grepl("'batch'", w, fixed = TRUE))
})

