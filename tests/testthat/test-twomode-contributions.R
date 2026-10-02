# Regression test: static change contributions for two-mode (bipartite) networks.
#
# The choice set open to an ego differs by network type:
#   one-mode  n alters, the diagonal (a == e) serving as "no change"
#   two-mode  m receivers plus a trailing "no change" slot, i.e. m + 1
#
# Before the fix the static path sized and iterated the per-ego array by the
# number of EGOS in both cases. For bipartite data that is wrong whenever
# m != n: with m > n the trailing receivers and the no-change slot were
# silently dropped; with m < n the loop ran past the receivers.
#
# Three pieces had to change together, so this test guards all of them:
#   - StatisticCalculator: nChoices = oneMode ? egos : alters + 1, used both for
#     the contributions array and for the permitted cache
#   - siena07setup: size by the dependent variable the effect belongs to rather
#     than by the first one
#   - contributions.R: sienaSetupDataForCpp(includeBipartite = TRUE). With
#     FALSE the bipartite depvar never reaches C++ and R *aborts* inside
#     C_setupModelOptions -- an abort, not a catchable error, so a regression
#     here takes the whole test session down rather than failing one expectation.
#
# The 12 x 10 fixture is the two-mode setup from checkEffects.R (same seed and
# construction), so it matches the data the effects are normally checked with.

testthat::skip_on_cran()

make_twomode <- function(nS, nR) {
  set.seed(12321)
  wave1 <- matrix(0, nS, nR)
  for (i in seq_len(nS * nR)) {
    wave1[i] <- sample(c(0, 1), 1, prob = c(0.7, 0.3))
  }
  wave2 <- wave1
  for (i in seq_len(nS * nR)) {
    wave2[i] <- abs(wave1[i] - 0.5 + sample(c(-0.5, 0.5), 1, prob = c(0.3, 0.7)))
  }
  senders    <- sienaNodeSet(nS, nodeSetName = "senders")
  recipients <- sienaNodeSet(nR, nodeSetName = "recipients")
  network <- sienaDependent(array(c(wave1, wave2), dim = c(nS, nR, 2)),
                            type = "bipartite",
                            nodeSet = c("senders", "recipients"),
                            allowOnly = FALSE)
  sienaDataCreate(network, nodeSets = list(senders, recipients))
}

fit_quick <- function(dat, eff) {
  algo <- sienaAlgorithmCreate(projname = NULL, nsub = 1L, n3 = 20L, seed = 42L)
  siena07(algo, data = dat, effects = eff,
          batch = TRUE, silent = TRUE, verbose = FALSE)
}

# m < n is the checkEffects fixture; m > n is the orientation that silently
# truncated, and newparallel.R exercises both orientations for the same reason.
for (dims in list(c(nS = 12L, nR = 10L), c(nS = 10L, nR = 12L))) {

  nS <- dims[["nS"]]; nR <- dims[["nR"]]

  test_that(sprintf("two-mode %dx%d: choice set is m + 1", nS, nR), {
    dat <- make_twomode(nS, nR)
    eff <- getEffects(dat)
    fit <- fit_quick(dat, eff)

    cc <- getStaticChangeContributions(ans = fit, data = dat, effects = eff)

    # structure: [[depvar]][[effect]][[ego]] -> numeric of length nChoices
    expect_true(is.list(cc))
    perEffect <- cc[[1L]]
    expect_gt(length(perEffect), 0L)

    perEgo <- perEffect[[1L]]
    expect_length(perEgo, nS)                        # one entry per sender
    expect_length(perEgo[[1L]], nR + 1L)             # m receivers + no change
    expect_true(all(vapply(perEgo, length, integer(1)) == nR + 1L))

    # the fix removed a branch that wrote NaN into unreachable slots
    expect_false(any(vapply(perEgo, function(x) any(is.nan(x)), logical(1))))
  })
}

test_that("one-mode choice set is unchanged by the two-mode fix", {
  n <- 8L
  set.seed(2)
  arr <- array(stats::rbinom(n * n * 2L, 1L, 0.3), dim = c(n, n, 2L))
  arr[, , 1L][diag(n) == 1] <- 0
  arr[, , 2L][diag(n) == 1] <- 0
  dat <- sienaDataCreate(net = sienaDependent(arr))
  eff <- getEffects(dat)
  fit <- fit_quick(dat, eff)

  cc <- getStaticChangeContributions(ans = fit, data = dat, effects = eff)
  perEgo <- cc[[1L]][[1L]]

  expect_length(perEgo, n)
  expect_length(perEgo[[1L]], n)   # diagonal serves as "no change": n, not n + 1
})
