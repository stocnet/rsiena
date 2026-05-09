## test-deltamethod-buckets.R
##
## Unit tests for .evalSpecChainBuckets / gradReinforceRowwise.
##
## Focus: scalar parity — for R=1 (no grouping vars), .evalSpecChainBuckets
## returns Y/N consistent with .evalSpecScalar Q_chain; gradReinforceRowwise
## agrees with gradReinforceDolby (equal-N_c fixture, so equality is exact).
##
## Algebraic guarantee (constant N_c = K per chain):
##   Q_hat_r = sum(Y)/sum(N) = mean(Q_chain)
##   A[c]    = (Y[c] - Q_hat_r * K) / K = Q_chain[c] - mean(Q_chain)
##   grad_new[k] = mean(A*sk) - cov(A*sk, sk)/var(sk) * mean(sk)
##              = grad_old[k]    (add/subtract mean(Q_chain)*mean(sk) terms cancel)
##
## These tests are the "before" check: they pass against the existing code.
## After scalar unification they remain the correctness anchor.

context("deltaFull bucketed REINFORCE — scalar parity")

## ── Synthetic dynamic wide fixture ───────────────────────────────────────────
##
## n_chains chains × n_egos egos × 3 alternatives (add/drop/no-change).
## The no-change rows are excluded by the keep mask (density == 0), leaving
## exactly 2 * n_egos rows per chain → constant N_c across all chains.
make_wide_dyn_bucket <- function(n_chains = 3L, n_egos = 4L, seed = 42L) {
  set.seed(seed)
  n_alts  <- 3L                           # add / drop / no-change
  n_total <- n_chains * n_egos * n_alts

  chain_vec <- rep(seq_len(n_chains), each = n_egos * n_alts)
  ego_vec   <- rep(rep(seq_len(n_egos), each = n_alts), times = n_chains)
  ms_vec    <- ego_vec                    # 1 ministep per ego (simple)

  ## Density column: +1 (add), -1 (drop), 0 (no-change) — repeating per ego.
  density_col <- rep(c(1L, -1L, 0L), times = n_chains * n_egos)
  ## Recip column: random 0/1; zero on no-change rows.
  recip_col <- sample(c(0L, 1L), n_total, replace = TRUE)
  recip_col[density_col == 0L] <- 0L

  contribMat <- matrix(c(density_col, recip_col), ncol = 2L,
                       dimnames = list(NULL, c("density", "recip")))

  ## group_id: each (chain, ego) combination forms its own softmax group.
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

## Mock predictFun: deterministic function of change stats; ignores theta.
## outcome = density * 0.3 + recip * 0.1
mock_pred_bucket <- function(changeContributions, theta, baseline = NULL,
                             outcomesOnly = FALSE, ...) {
  cs <- changeContributions$changeStats
  ov <- as.numeric(cs$density * 0.3 + cs$csMat[, "recip"] * 0.1)
  if (outcomesOnly) return(list(me_out = ov))
  data.frame(me_out = ov)
}

make_scalar_spec_bucket <- function() {
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

## ── 1. Y/N from bucketed path matches Q_chain from scalar path ───────────────

test_that("scalar spec: .evalSpecChainBuckets Y/N matches .evalSpecScalar Q_chain", {
  skip_on_cran()

  wide  <- make_wide_dyn_bucket()
  spec  <- make_scalar_spec_bucket()
  theta <- c(density = -2, recip = 1)
  type  <- "changeProb"

  res_s  <- .evalSpecScalar(spec, wide, theta, type)
  key_df <- data.frame(.row = 1L)
  res_b  <- .evalSpecChainBuckets(spec, wide, theta, type, key_df)

  expect_false(is.null(res_b),
    info = ".evalSpecChainBuckets must return non-NULL for a valid scalar spec")
  expect_equal(nrow(res_b$Y), 3L, info = "3 chains expected")
  expect_equal(ncol(res_b$Y), 1L, info = "1 output cell (scalar spec)")

  ## Chain IDs must match between the two paths (both sort unique chains).
  expect_equal(res_b$chain_ids, res_s$chain_ids)

  ## Per-chain ratio Y/N must equal per-chain mean from .evalSpecScalar.
  Q_bucket <- as.numeric(res_b$Y[, 1L] / pmax(res_b$N[, 1L], 1L))
  expect_equal(Q_bucket, res_s$Q_chain, tolerance = 1e-10)
})

## ── 2. Gradient parity: gradReinforceRowwise == gradReinforceDolby ────────────
##
## When N_c is constant across chains, the ratio-influence formula in
## gradReinforceRowwise collapses algebraically to gradReinforceDolby.
## Tolerance is machine-precision tight.

test_that("scalar parity: gradReinforceRowwise[1,] equals gradReinforceDolby (equal N_c)", {
  skip_on_cran()

  wide  <- make_wide_dyn_bucket()
  spec  <- make_scalar_spec_bucket()
  theta <- c(density = -2, recip = 1)
  type  <- "changeProb"

  res_s <- .evalSpecScalar(spec, wide, theta, type)
  key_df <- data.frame(.row = 1L)
  res_b  <- .evalSpecChainBuckets(spec, wide, theta, type, key_df)

  ## Synthetic score matrix aligned to chain_ids.
  n_chains <- length(res_b$chain_ids)
  nParams  <- 2L
  set.seed(77L)
  ssc_mat <- matrix(rnorm(n_chains * nParams), nrow = n_chains, ncol = nParams)

  grad_old <- gradReinforceDolby(res_s$Q_chain, ssc_mat)
  grad_new <- as.numeric(gradReinforceRowwise(res_b$Y, res_b$N, ssc_mat)[1L, ])

  expect_equal(grad_new, grad_old, tolerance = 1e-10,
    info = "ratio-influence and direct-Q formulas must agree for equal N_c")
})

## ── 3. Bucketed path returns NULL for accumulated spec ───────────────────────

test_that(".evalSpecChainBuckets returns NULL for accumulated spec", {
  skip_on_cran()

  wide   <- make_wide_dyn_bucket()
  spec   <- make_scalar_spec_bucket()
  spec$accumulated <- TRUE
  theta  <- c(density = -2, recip = 1)
  type   <- "changeProb"
  key_df <- data.frame(.row = 1L)

  res <- .evalSpecChainBuckets(spec, wide, theta, type, key_df)
  expect_null(res, info = "accumulated spec must return NULL")
})

## ── 4. Bucketed path returns NULL when wide has no chain column ───────────────

test_that(".evalSpecChainBuckets returns NULL when wide has no chain column", {
  skip_on_cran()

  wide  <- make_wide_dyn_bucket()
  wide$chain <- NULL                  # simulate static wide (no chain col)
  spec  <- make_scalar_spec_bucket()
  theta <- c(density = -2, recip = 1)
  type  <- "changeProb"
  key_df <- data.frame(.row = 1L)

  res <- .evalSpecChainBuckets(spec, wide, theta, type, key_df)
  expect_null(res, info = "no chain column must return NULL")
})
