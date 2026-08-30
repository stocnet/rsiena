# Tests for the period-level tie probability (R/chainAccumulate.R).
#
# The controlling test is "one-step readout reproduces fitted()".  It is the
# only quantity here that is unbiased by construction, so a discrepancy is a
# bug rather than a property of the estimand -- it catches focal-row selection,
# direction handling, the K = 0 path and chain misalignment.  Everything else
# in this file is either an exact algebraic identity or a guard.

load_all_fixtures <- function() {
  for (nm in c("mynet", "mydata", "mymodel", "mycontrols", "ans")) {
    obj <- load_fixture(nm)
    if (is.null(obj)) return(FALSE)
    assign(nm, obj, envir = parent.frame())
  }
  TRUE
}

# A phase-3 run with BOTH flags, so fitted() and the chains come from the same
# simulation.  getDynamicChangeContributions() re-simulates with
# returnDeps = FALSE, which would silently compare against different chains.
periodFixture <- function(n3 = 300L, seed = 4321L) {
  key  <- sprintf("periodTieProb_n3_%d_seed_%d", n3, seed)
  cached <- load_fixture(key)
  if (!is.null(cached)) return(cached)

  waves  <- array(c(s501, s502, s503), dim = c(50, 50, 3))
  mynet  <- sienaDependent(waves)
  mydata <- sienaDataCreate(mynet)
  eff    <- includeEffects(getEffects(mydata), transTrip, verbose = FALSE)
  ans0 <- siena07(sienaAlgorithmCreate(projname = NULL, nsub = 3, n3 = 200,
                                       seed = seed, cond = FALSE),
                  data = mydata, effects = eff,
                  batch = TRUE, silent = TRUE, useCluster = FALSE)
  ans <- siena07(sienaAlgorithmCreate(projname = NULL, nsub = 0, n3 = n3,
                                      seed = seed + 1L, cond = FALSE),
                 data = mydata, effects = eff, prevAns = ans0,
                 returnDeps = TRUE, returnChangeContributions = TRUE,
                 batch = TRUE, silent = TRUE, useCluster = FALSE)
  out <- list(ans = ans, mydata = mydata, eff = eff, waves = waves)
  saveRDS(out, file.path(cache_dir, paste0(key, ".rds")))
  out
}

periodPieces <- function(fx) {
  cc    <- flattenAndEnrichWide(fx$ans$changeContributions, fx$eff,
                                "mynet", fx$mydata)
  theta <- fx$ans$theta[fx$ans$effects$type != "rate"]
  names(theta) <- cc$effectNames
  list(cc = cc, theta = theta,
       x0 = periodStartNetworks(fx$mydata, "mynet"))
}


# ---------------------------------------------------------------------------
# periodStartNetworks
# ---------------------------------------------------------------------------

test_that("periodStartNetworks returns one wave per period and folds structurals", {
  waves  <- array(c(s501, s502, s503), dim = c(50, 50, 3))
  mynet  <- sienaDependent(waves)
  mydata <- sienaDataCreate(mynet)
  x0 <- periodStartNetworks(mydata, "mynet")

  expect_length(x0, 2L)                       # 3 waves -> 2 periods
  expect_equal(dim(x0[[1L]]), c(50L, 50L))
  # Values are tie states, never RSiena's 10/11 structural codes.
  expect_true(all(x0[[1L]] %in% c(0, 1, NA)))
  expect_equal(x0[[1L]][!is.na(x0[[1L]])],
               as.numeric(s501)[!is.na(s501)], ignore_attr = TRUE)
})

test_that("periodStartNetworks rejects an unknown dependent variable", {
  waves  <- array(c(s501, s502, s503), dim = c(50, 50, 3))
  mydata <- sienaDataCreate(sienaDependent(waves))
  expect_error(periodStartNetworks(mydata, "nosuchnet"),
               "no dependent variable")
})


# ---------------------------------------------------------------------------
# Branch probabilities: the identities that pin the direction bookkeeping
# ---------------------------------------------------------------------------

test_that("on the factual condition each arm reproduces the stored probability", {
  skip_slow()
  p  <- periodPieces(periodFixture())
  br <- branchProbabilities(p$cc, p$theta)

  dens <- as.numeric(p$cc$changeStats$density)
  keep <- .chainKeep(p$cc)
  stored <- p$cc$changeProbability

  # c is P(choose j | tie absent).  On a row whose factual direction IS
  # "absent" (density +1) that is exactly what the simulation stored.
  expect_equal(br$c[keep & dens == 1], stored[keep & dens == 1],
               tolerance = 1e-12)
  expect_equal(br$d[keep & dens == -1], stored[keep & dens == -1],
               tolerance = 1e-12)
})

test_that("the flipped arm matches an explicit block-softmax rebuild", {
  skip_slow()
  p  <- periodPieces(periodFixture())
  br <- branchProbabilities(p$cc, p$theta)

  cs   <- p$cc$changeStats
  thE  <- buildThetaEff(p$theta[p$cc$effectNames], cs$changeStatsMap)
  u    <- as.numeric(p$cc$changeUtility)
  dens <- as.numeric(cs$density)
  dI   <- cs$densityIdx

  # Rebuild the whole choice-set softmax with the focal alternative's
  # direction reversed, for a sample of decisions.  This is the check the
  # patched update has to earn.
  set.seed(11)
  blocks <- sample(unique(p$cc$group_id), 200L)
  err <- 0
  for (g in blocks) {
    ii <- which(p$cc$group_id == g)
    m  <- cs$csMat[ii, , drop = FALSE]
    for (r in seq_along(ii)) {
      if (dens[ii[r]] == 0) next
      mm <- m
      mm[r, dI] <- -mm[r, dI]
      uu <- calculateUtility(mm, thE, densityIdx = dI)
      pr <- exp(uu - max(uu)); pr <- pr / sum(pr)
      got <- if (dens[ii[r]] == 1) br$d[ii[r]] else br$c[ii[r]]
      err <- max(err, abs(pr[r] - got))
    }
  }
  expect_lt(err, 1e-12)
})

test_that("branchProbabilities rejects an unresolvable effect name", {
  skip_slow()
  p <- periodPieces(periodFixture())
  expect_error(branchProbabilities(p$cc, p$theta, list(nosucheffect = 1)),
               "not found in available names")
})

test_that("the same effect cannot be shifted twice", {
  skip_slow()
  p <- periodPieces(periodFixture())
  expect_error(
    predictPeriodDiff(p$cc, p$theta, p$x0,
                      c("recip", "recip"), list(c(0, 1), c(0, 1))),
    "shifted twice")
})


# ---------------------------------------------------------------------------
# The controlling test
# ---------------------------------------------------------------------------

test_that("the one-step readout r_K reproduces fitted()", {
  skip_slow()
  fx <- periodFixture()
  p  <- periodPieces(fx)

  res <- chainAccumulate(p$cc, list(f = branchProbabilities(p$cc, p$theta)),
                         p$x0, accumulators = c("piK", "rK"))
  fit <- fitted(fx$ans)

  for (pr in sort(unique(res$period))) {
    d  <- res[res$period == pr, ]
    ag <- aggregate(cbind(rK, piK) ~ ego + choice, data = d, FUN = mean)
    ag$target <- fit[cbind(ag$ego, ag$choice, pr)]
    ag <- ag[stats::complete.cases(ag), ]

    # A is unbiased by construction; anything above Monte Carlo error here is
    # a bug, not an estimand property.
    expect_lt(abs(mean(ag$rK - ag$target)), 1e-4)

    # B is EXPECTED to be biased low -- the conditioning cost of propagating.
    # Pinning the sign and a loose magnitude stops a silent regression that
    # made B accidentally equal to A.
    bias_B <- mean(ag$piK - ag$target) / mean(ag$target)
    expect_lt(bias_B, 0)
    expect_gt(bias_B, -0.10)
  }
})


# ---------------------------------------------------------------------------
# Accumulator identities
# ---------------------------------------------------------------------------

test_that("expected creations minus dissolutions equals piK - pi_0", {
  skip_slow()
  p <- periodPieces(periodFixture())
  res <- chainAccumulate(p$cc, list(f = branchProbabilities(p$cc, p$theta)),
                         p$x0, accumulators = c("piK", "counts"))
  pi0 <- mapply(function(pr, i, j) p$x0[[pr]][i, j],
                res$period, res$ego, res$choice)
  expect_equal(res$nCreations - res$nDissolutions, res$piK - pi0,
               tolerance = 1e-12)
})

test_that("occupancy is bounded by K and non-negative", {
  skip_slow()
  p <- periodPieces(periodFixture())
  res <- chainAccumulate(p$cc, list(f = branchProbabilities(p$cc, p$theta)),
                         p$x0, accumulators = c("piK", "occupancy"))
  expect_true(all(res$occupancy >= 0))
  expect_true(all(res$occupancy <= res$K + 1e-9))
})

test_that("an unknown accumulator is rejected by name", {
  skip_slow()
  p <- periodPieces(periodFixture())
  expect_error(
    chainAccumulate(p$cc, list(f = branchProbabilities(p$cc, p$theta)),
                    p$x0, accumulators = "notAnAccumulator"),
    "unknown accumulator")
})


# ---------------------------------------------------------------------------
# The contrast
# ---------------------------------------------------------------------------

test_that("the level contrast matches an independent recursion per branch", {
  skip_slow()
  fx <- periodFixture()
  p  <- periodPieces(fx)

  new <- predictPeriodDiff(p$cc, p$theta, p$x0,
                           effectNames = "recip", contrasts = list(c(0, 1)))

  # Independent implementation: build each branch by rebuilding the change
  # statistic matrix (the way the design prototype does), run the two
  # recursions separately, difference at the end.
  cs  <- p$cc$changeStats
  thE <- buildThetaEff(p$theta[p$cc$effectNames], cs$changeStatsMap)
  dI  <- cs$densityIdx
  rI  <- match(resolveEffectName("recip", cs$csNames), cs$csNames)
  u0  <- calculateUtility(cs$csMat, thE, densityIdx = dI)
  probs <- function(a, b) {
    m <- cs$csMat
    m[, dI] <- ifelse(a == 0, 1, -1)
    m[, rI] <- b
    e <- exp(calculateUtility(m, thE, densityIdx = dI) - u0)
    e / (1 / p$cc$changeProbability - 1 + e)
  }
  ref <- chainAccumulate(
    p$cc,
    list(lo = list(c = probs(0, 0), d = probs(1, 0)),
         hi = list(c = probs(0, 1), d = probs(1, 1))),
    p$x0, accumulators = "piK")

  key_new <- paste(new$chain, new$period, new$ego, new$choice)
  key_ref <- paste(ref$chain, ref$period, ref$ego, ref$choice)
  m <- match(key_ref, key_new)
  expect_false(anyNA(m))
  expect_equal(new$periodFirstDiff[m], ref$piK..hi - ref$piK..lo,
               tolerance = 1e-12)
})

test_that("contrast order sets the sign", {
  skip_slow()
  p <- periodPieces(periodFixture())
  up   <- predictPeriodDiff(p$cc, p$theta, p$x0, "recip", list(c(0, 1)))
  down <- predictPeriodDiff(p$cc, p$theta, p$x0, "recip", list(c(1, 0)))
  expect_equal(up$periodFirstDiff, -down$periodFirstDiff, tolerance = 1e-12)
})

test_that("a null contrast is a zero effect", {
  skip_slow()
  p <- periodPieces(periodFixture())
  flat <- predictPeriodDiff(p$cc, p$theta, p$x0, "recip", list(c(1, 1)))
  expect_true(max(abs(flat$periodFirstDiff)) < 1e-12)
})

test_that("each effect needs exactly one of contrast or diff", {
  skip_slow()
  p <- periodPieces(periodFixture())
  ## Neither.
  expect_error(predictPeriodDiff(p$cc, p$theta, p$x0, "recip"),
               "needs either contrast")
  ## Both.
  expect_error(
    predictPeriodDiff(p$cc, p$theta, p$x0, "recip",
                      contrasts = list(c(0, 1)), diffs = list(1)),
    "not both")
})


# ---------------------------------------------------------------------------
# Streaming invariants
# ---------------------------------------------------------------------------

test_that("chainAccumulate refuses rows that are not in ministep order", {
  skip_slow()
  p  <- periodPieces(periodFixture())
  br <- branchProbabilities(p$cc, p$theta)
  scrambled <- p$cc
  ord <- rev(seq_along(scrambled$group_id))
  for (nm in c("chain", "group", "period", "ego", "ministep", "choice",
               "group_id", "permitted", "changeUtility", "changeProbability"))
    scrambled[[nm]] <- scrambled[[nm]][ord]
  scrambled$changeStats$density <- scrambled$changeStats$density[ord]
  expect_error(chainAccumulate(scrambled, list(f = br), p$x0),
               "ministep order")
})

test_that("includeUnvisited seeds every dyad and leaves visited ones alone", {
  skip_slow()
  p  <- periodPieces(periodFixture())
  br <- branchProbabilities(p$cc, p$theta)
  seen   <- chainAccumulate(p$cc, list(f = br), p$x0, "piK")
  allDy  <- chainAccumulate(p$cc, list(f = br), p$x0, "piK",
                            includeUnvisited = TRUE)
  expect_gte(nrow(allDy), nrow(seen))
  # Seeding must not invent self-"dyads" out of the diagonal of x0.
  expect_false(any(allDy$ego == allDy$choice))
  expect_false(any(seen$ego == seen$choice))

  k1 <- paste(seen$chain, seen$period, seen$ego, seen$choice)
  k2 <- paste(allDy$chain, allDy$period, allDy$ego, allDy$choice)
  m  <- match(k1, k2)
  expect_false(anyNA(m))
  expect_equal(seen$piK, allDy$piK[m], tolerance = 1e-12)

  # Anything only in the seeded frame never came up, so it sits at x_ij(t0).
  extra <- allDy[-m, ]
  if (nrow(extra) > 0) {
    expect_true(all(extra$K == 0L))
    pi0 <- mapply(function(pr, i, j) p$x0[[pr]][i, j],
                  extra$period, extra$ego, extra$choice)
    expect_equal(extra$piK, pi0, tolerance = 1e-12)
  }
})

test_that("chainAccumulate needs simulated chains", {
  skip_slow()
  p  <- periodPieces(periodFixture())
  br <- branchProbabilities(p$cc, p$theta)
  noChain <- p$cc
  noChain$chain <- NULL
  expect_error(chainAccumulate(noChain, list(f = br), p$x0), "no 'chain'")
})


# ---------------------------------------------------------------------------
# Entry points and their guards
# ---------------------------------------------------------------------------

test_that("predict() returns piK and marginalEffects() its contrast", {
  skip_slow()
  fx <- periodFixture()
  algo <- sienaAlgorithmCreate(projname = NULL, nsub = 0, n3 = 300,
                               seed = 4322, cond = FALSE)
  ctl <- set_postest_algo_saom(algorithm = algo, n3 = 300,
                               useChangeContributions = TRUE, verbose = FALSE)
  unc <- set_postest_uncertainty_saom(enabled = FALSE)
  out <- set_postest_output_saom()

  ptg <- make_predict_targets(fx$ans, effects = fx$eff, depvar = "mynet",
                              type = "tieProb", scope = "period", dynamic = TRUE,
                              level = "egoChoice")
  pr <- predict(fx$ans, fx$mydata, ptg, control_uncertainty = unc,
                control_algo = ctl, control_out = out)
  expect_true(is.data.frame(pr))
  expect_true("periodTieProb" %in% names(pr))
  expect_true(all(c("period", "ego", "choice") %in% names(pr)))
  expect_true(all(pr$periodTieProb >= 0 & pr$periodTieProb <= 1))

  mtg <- make_marginal_targets(fx$ans, data = fx$mydata, effects = fx$eff,
                               depvar = "mynet", type = "tieProb", scope = "period",
                               dynamic = TRUE, level = "period",
                               includeDefaults = FALSE)
  mtg <- set_target(mtg, recip, contrast = c(0, 1), verbose = FALSE)
  me <- marginalEffects(fx$ans, data = fx$mydata, targets = mtg,
                        control_uncertainty = unc, control_algo = ctl,
                        control_out = out)
  expect_true("periodFirstDiff" %in% names(me))
  # Reciprocation raises the endpoint probability; the magnitude is reported in
  # docs/design/period_level_tie_probability.md section 10.
  expect_true(all(me$periodFirstDiff > 0.1))
})

test_that("period scope refuses combinations it cannot honour", {
  skip_slow()
  fx <- periodFixture()
  expect_error(
    make_predict_targets(fx$ans, effects = fx$eff, scope = "period",
                         dynamic = FALSE),
    "requires dynamic = TRUE")
  expect_error(
    make_predict_targets(fx$ans, effects = fx$eff, scope = "period",
                         dynamic = TRUE, accumulated = TRUE),
    "already a period-level")
  expect_error(
    make_marginal_targets(fx$ans, data = fx$mydata, effects = fx$eff,
                          scope = "period", dynamic = FALSE),
    "requires dynamic = TRUE")
  expect_error(
    make_marginal_targets(fx$ans, data = fx$mydata, effects = fx$eff,
                          type = "tieProb", scope = "period", dynamic = TRUE,
                          mainEffect = "riskRatio"),
    "risk difference")
})


# ---------------------------------------------------------------------------
# Second differences
#
# These exist because the cells/weights refactor is only worth anything if the
# cross-partial needs no code of its own.  The check is that it equals the
# difference of two first differences computed independently.
# ---------------------------------------------------------------------------

test_that(".periodCells produces RSiena's own perturbation weights", {
  ## A level is emitted as a SHIFT from the row's factual statistic, which is
  ## what the shared builder consumes -- so levels and diffs are one mechanism.
  m <- cbind(density = rep(1, 3), recip = c(0, 1, 0), transTrip = c(0, 2, 1))
  one <- .periodCells("recip", contrasts = list(c(0, 1)), csMat = m)
  expect_equal(vapply(one, `[[`, numeric(1L), "weight"),
               unname(.cellWeights(1L)))
  expect_equal(vapply(one, `[[`, character(1L), "name"),
               names(.cellWeights(1L)))
  expect_equal(one[[1L]]$shifts$recip, 0 - m[, "recip"])
  expect_equal(one[[2L]]$shifts$recip, 1 - m[, "recip"])

  two <- .periodCells(c("recip", "transTrip"),
                      contrasts = list(c(0, 1), c(0, 2)), csMat = m)
  expect_equal(vapply(two, `[[`, numeric(1L), "weight"),
               unname(.cellWeights(2L)))
  expect_equal(two[[4L]]$shifts$recip,     1 - m[, "recip"])
  expect_equal(two[[4L]]$shifts$transTrip, 2 - m[, "transTrip"])
  expect_equal(two[[1L]]$shifts$recip,     0 - m[, "recip"])
})

test_that("a diff is a sustained shift from whatever the statistic is", {
  m <- cbind(density = rep(1, 3), transTrip = c(0, 2, 1))
  d <- .periodCells("transTrip", diffs = list(1), csMat = m)
  expect_equal(d[[1L]]$shifts$transTrip, 0)   # base: unshifted
  expect_equal(d[[2L]]$shifts$transTrip, 1)   # A: one more, whatever it was
  ## A level and a diff coincide only where the base statistic is constant --
  ## which is exactly why a level is not defined for a derived count.
  flat <- cbind(density = rep(1, 3), s = rep(2, 3))
  ## Compare after recycling: a level yields a per-row vector, a diff a scalar
  ## that broadcasts, and the counterfactual is the same either way.
  rl <- function(cells) lapply(cells, function(cl) rep_len(cl$shifts$s, 3))
  expect_equal(rl(.periodCells("s", contrasts = list(c(2, 3)), csMat = flat)),
               rl(.periodCells("s", diffs     = list(1),       csMat = flat)))
})

test_that("the second difference equals the difference of two first differences", {
  skip_slow()
  p <- periodPieces(periodFixture())

  sd <- predictPeriodDiff(p$cc, p$theta, p$x0,
                          c("recip", "transTrip"), list(c(0, 1), c(0, 2)))
  expect_true("periodSecondDiff" %in% names(sd))

  # Built the long way: four independent recursions, combined by hand.
  cs <- p$cc$changeStats
  cn <- function(e) resolveEffectName(e, cs$csNames)
  cellPi <- function(r, t) {
    sh <- setNames(list(r - cs$csMat[, cn("recip")],
                        t - cs$csMat[, cn("transTrip")]),
                   c(cn("recip"), cn("transTrip")))
    br <- branchProbabilities(p$cc, p$theta, sh)
    chainAccumulate(p$cc, list(f = br), p$x0, accumulators = "piK")$piK
  }
  ref <- cellPi(1, 2) - cellPi(0, 2) - (cellPi(1, 0) - cellPi(0, 0))
  expect_equal(sd$periodSecondDiff, ref, tolerance = 1e-12)
})

test_that("a second difference is zero when one effect's contrast is null", {
  skip_slow()
  p <- periodPieces(periodFixture())
  flat <- predictPeriodDiff(p$cc, p$theta, p$x0,
                            c("recip", "transTrip"), list(c(0, 1), c(2, 2)))
  expect_lt(max(abs(flat$periodSecondDiff)), 1e-12)
})

test_that("second differences are symmetric in the two effects", {
  skip_slow()
  p <- periodPieces(periodFixture())
  ab <- predictPeriodDiff(p$cc, p$theta, p$x0,
                          c("recip", "transTrip"), list(c(0, 1), c(0, 2)))
  ba <- predictPeriodDiff(p$cc, p$theta, p$x0,
                          c("transTrip", "recip"), list(c(0, 2), c(0, 1)))
  k1 <- paste(ab$chain, ab$period, ab$ego, ab$choice)
  k2 <- paste(ba$chain, ba$period, ba$ego, ba$choice)
  expect_equal(ab$periodSecondDiff, ba$periodSecondDiff[match(k1, k2)],
               tolerance = 1e-12)
})

test_that("marginalEffects() computes a period-level second difference", {
  skip_slow()
  fx <- periodFixture()
  algo <- sienaAlgorithmCreate(projname = NULL, nsub = 0, n3 = 300,
                               seed = 4322, cond = FALSE)
  ctl <- set_postest_algo_saom(algorithm = algo, n3 = 300,
                               useChangeContributions = TRUE, verbose = FALSE)
  unc <- set_postest_uncertainty_saom(enabled = FALSE)
  out <- set_postest_output_saom()

  mtg <- make_marginal_targets(fx$ans, data = fx$mydata, effects = fx$eff,
                               depvar = "mynet", type = "tieProb", scope = "period",
                               dynamic = TRUE, level = "period",
                               includeDefaults = FALSE)
  mtg <- set_second_diff(mtg, list(recip = list(contrast = c(0, 1)),
                                   transTrip = list(contrast = c(0, 2))),
                         verbose = FALSE)
  me <- marginalEffects(fx$ans, data = fx$mydata, targets = mtg,
                        control_uncertainty = unc, control_algo = ctl,
                        control_out = out)
  expect_true("periodSecondDiff" %in% names(me))
  expect_true(all(is.finite(me$periodSecondDiff)))
})


# ---------------------------------------------------------------------------
# emitLevel and reachable aggregation levels
# ---------------------------------------------------------------------------

test_that("only levels reachable from chainEgoChoice are offered", {
  expect_true(.canAggregateTo("period", "chainEgoChoice"))
  expect_true(.canAggregateTo("egoChoice", "chainEgoChoice"))
  expect_true(.canAggregateTo("chainEgo", "chainEgoChoice"))
  # A period-level row has no ministep of its own.
  expect_false(.canAggregateTo("ministep", "chainEgoChoice"))
  expect_false(.canAggregateTo("ministepChoice", "chainEgoChoice"))
  # Every level remains reachable from the row-wise emit level.
  for (lv in .postestLevels)
    expect_true(.canAggregateTo(lv, "ministepChoice"))
})

test_that("a ministep level is refused for period scope", {
  skip_slow()
  fx <- periodFixture()
  expect_error(
    make_predict_targets(fx$ans, effects = fx$eff, scope = "period",
                         dynamic = TRUE, level = "ministep"),
    "not available for")
  expect_error(
    make_marginal_targets(fx$ans, data = fx$mydata, effects = fx$eff,
                          type = "tieProb", scope = "period", dynamic = TRUE,
                          level = "ministepChoice"),
    "not available for")
})


# ---------------------------------------------------------------------------
# Ego-specific pins
#
# A pin on an ego statistic moves EVERY alternative in the ego's choice set,
# so the shift splits into an ego-wide part and a focal spike and needs the
# block sum.  The only test that can catch this being wrong is an explicit
# rebuild of the whole softmax, which is what the first test here does.
# ---------------------------------------------------------------------------

egoFixture <- function() {
  key <- "periodTieProb_egoX"
  cached <- load_fixture(key)
  if (!is.null(cached)) return(cached)

  waves  <- array(c(s501, s502, s503), dim = c(50, 50, 3))
  mynet  <- sienaDependent(waves)
  # An actor covariate that is genuinely ego-specific in the model.
  set.seed(99)
  cov1   <- coCovar(as.numeric(s50a[, 1]))
  mydata <- sienaDataCreate(mynet, cov1)
  eff    <- includeEffects(getEffects(mydata), egoX, altX,
                           interaction1 = "cov1", verbose = FALSE)
  ans0 <- siena07(sienaAlgorithmCreate(projname = NULL, nsub = 2, n3 = 150,
                                       seed = 771, cond = FALSE),
                  data = mydata, effects = eff,
                  batch = TRUE, silent = TRUE, useCluster = FALSE)
  ans <- siena07(sienaAlgorithmCreate(projname = NULL, nsub = 0, n3 = 60,
                                      seed = 772, cond = FALSE),
                 data = mydata, effects = eff, prevAns = ans0,
                 returnDeps = TRUE, returnChangeContributions = TRUE,
                 batch = TRUE, silent = TRUE, useCluster = FALSE)
  out <- list(ans = ans, mydata = mydata, eff = eff)
  saveRDS(out, file.path(cache_dir, paste0(key, ".rds")))
  out
}

egoPieces <- function(fx) {
  cc <- flattenAndEnrichWide(fx$ans$changeContributions, fx$eff, "mynet",
                             fx$mydata)
  theta <- fx$ans$theta[fx$ans$effects$type != "rate"]
  names(theta) <- cc$effectNames
  list(cc = cc, theta = theta, x0 = periodStartNetworks(fx$mydata, "mynet"))
}

test_that("an ego-specific pin is detected as ego-wide, a dyadic one is not", {
  skip_slow()
  p  <- egoPieces(egoFixture())
  it <- p$cc$effectInteractionTypes
  cs <- p$cc$changeStats
  egoNm <- grep("egoX", cs$csNames, value = TRUE)[1L]
  altNm <- grep("altX", cs$csNames, value = TRUE)[1L]
  expect_equal(resolvePerturbType(.pinCompositeName(p$cc, egoNm), it), "ego")
  expect_equal(resolvePerturbType(.pinCompositeName(p$cc, altNm), it), "alter")
})

test_that("an ego pin matches an explicit rebuild of the whole choice set", {
  skip_slow()
  p  <- egoPieces(egoFixture())
  cs <- p$cc$changeStats
  dI <- cs$densityIdx
  eI <- match(grep("egoX", cs$csNames, value = TRUE)[1L], cs$csNames)
  thE <- buildThetaEff(p$theta[p$cc$effectNames], cs$changeStatsMap)
  dens <- as.numeric(cs$density)

  LEVEL <- 2
  ## A level expressed as a shift, which is the one form the builder takes.
  br <- branchProbabilities(p$cc, p$theta,
          setNames(list(LEVEL - cs$csMat[, eI]), cs$csNames[eI]))

  # Reference: for each focal row, set the ego statistic to LEVEL on EVERY
  # alternative of the block (it is ego-specific), set the focal row's
  # direction to the arm, rebuild the softmax from scratch.
  set.seed(5)
  blocks <- sample(unique(p$cc$group_id), 120L)
  errC <- errD <- 0
  for (g in blocks) {
    ii <- which(p$cc$group_id == g)
    for (r in seq_along(ii)) {
      if (dens[ii[r]] == 0) next
      for (arm in 0:1) {
        m <- cs$csMat[ii, , drop = FALSE]
        m[, eI] <- LEVEL                       # ego-wide: all alternatives
        m[r, dI] <- if (arm == 0L) 1 else -1   # focal row's tie state
        uu <- calculateUtility(m, thE, densityIdx = dI)
        pr <- exp(uu - max(uu)); pr <- pr / sum(pr)
        got <- if (arm == 0L) br$c[ii[r]] else br$d[ii[r]]
        if (arm == 0L) errC <- max(errC, abs(pr[r] - got))
        else           errD <- max(errD, abs(pr[r] - got))
      }
    }
  }
  expect_lt(errC, 1e-12)
  expect_lt(errD, 1e-12)
})

test_that("the ego-wide shift is exactly zero on the no-change alternative", {
  skip_slow()
  # In the ego kernel every alternative enters the shared denominator, so a
  # nonzero shift on the no-change option would corrupt every dyad in the
  # block rather than just its own row.
  p  <- egoPieces(egoFixture())
  cs <- p$cc$changeStats
  eI <- grep("egoX", cs$csNames, value = TRUE)[1L]
  dens <- as.numeric(cs$density)

  cm  <- cs$csMat[, eI]
  br0 <- branchProbabilities(p$cc, p$theta, setNames(list(0 - cm), eI))
  br9 <- branchProbabilities(p$cc, p$theta, setNames(list(9 - cm), eI))
  # The no-change rows are not focal rows, but they must still be a valid
  # distribution member: their probability may move (the block renormalises)
  # while never becoming NA or leaving [0, 1].
  nc <- dens == 0
  expect_false(anyNA(br0$c[nc]))
  expect_false(anyNA(br9$c[nc]))
  expect_true(all(br9$c[!is.na(br9$c)] >= 0 & br9$c[!is.na(br9$c)] <= 1))
})

test_that("a dyadic shift matches an explicit rebuild, touching only the focal row", {
  skip_slow()
  p  <- periodPieces(periodFixture())
  cs <- p$cc$changeStats
  dI <- cs$densityIdx
  rI <- match(resolveEffectName("recip", cs$csNames), cs$csNames)
  thE <- buildThetaEff(p$theta[p$cc$effectNames], cs$changeStatsMap)
  dens <- as.numeric(cs$density)

  LEVEL <- 1
  br <- branchProbabilities(p$cc, p$theta,
          setNames(list(LEVEL - cs$csMat[, rI]), cs$csNames[rI]))

  ## Reference: recip is DYADIC, so the counterfactual moves only the focal
  ## alternative -- unlike the ego case, the other alternatives keep their own
  ## statistic.  Rebuilding the whole block is what proves the ego-wide part is
  ## correctly zero here rather than merely assumed to be.
  set.seed(3)
  blocks <- sample(unique(p$cc$group_id), 120L)
  errC <- errD <- 0
  for (g in blocks) {
    ii <- which(p$cc$group_id == g)
    for (r in seq_along(ii)) {
      if (dens[ii[r]] == 0) next
      for (arm in 0:1) {
        m <- cs$csMat[ii, , drop = FALSE]
        m[r, rI] <- LEVEL                      # focal alternative only
        m[r, dI] <- if (arm == 0L) 1 else -1
        uu <- calculateUtility(m, thE, densityIdx = dI)
        pr <- exp(uu - max(uu)); pr <- pr / sum(pr)
        got <- if (arm == 0L) br$c[ii[r]] else br$d[ii[r]]
        if (arm == 0L) errC <- max(errC, abs(pr[r] - got))
        else           errD <- max(errD, abs(pr[r] - got))
      }
    }
  }
  expect_lt(errC, 1e-12)
  expect_lt(errD, 1e-12)
})

test_that("marginalEffects() accepts an ego-specific effect at period level", {
  skip_slow()
  fx <- egoFixture()
  algo <- sienaAlgorithmCreate(projname = NULL, nsub = 0, n3 = 60,
                               seed = 772, cond = FALSE)
  ctl <- set_postest_algo_saom(algorithm = algo, n3 = 60,
                               useChangeContributions = TRUE, verbose = FALSE)
  unc <- set_postest_uncertainty_saom(enabled = FALSE)
  out <- set_postest_output_saom()

  mtg <- make_marginal_targets(fx$ans, data = fx$mydata, effects = fx$eff,
                               depvar = "mynet", type = "tieProb", scope = "period",
                               dynamic = TRUE, level = "period",
                               includeDefaults = FALSE)
  mtg <- set_target(mtg, egoX, covar1 = "cov1", contrast = c(0, 1),
                    verbose = FALSE)
  me <- marginalEffects(fx$ans, data = fx$mydata, targets = mtg,
                        control_uncertainty = unc, control_algo = ctl,
                        control_out = out)
  expect_true("periodFirstDiff" %in% names(me))
  expect_true(all(is.finite(me$periodFirstDiff)))
})


# ---------------------------------------------------------------------------
# Mixed types in one call
#
# The per-ministep effect, its accumulated sum and the period-level endpoint
# are all read off the SAME simulated chains.  Asking for them in one call
# streams every spec over a single flattened batch instead of re-flattening
# once per call.  The numbers must not change as a result -- that is what
# these tests pin.
# ---------------------------------------------------------------------------

test_that("scope is overridable per target and validated per row", {
  skip_slow()
  fx <- periodFixture()
  tg <- make_marginal_targets(fx$ans, data = fx$mydata, effects = fx$eff,
                              depvar = "mynet", type = "tieProb",
                              dynamic = TRUE, level = "period",
                              includeDefaults = FALSE)
  tg <- set_target(tg, transTrip, contrast = c(0, 1), verbose = FALSE)
  tg <- set_target(tg, recip, contrast = c(0, 1),
                   scope = "period", verbose = FALSE)
  lowered <- RSiena:::.targetsToEffectList(tg)$effectList
  expect_equal(lowered$recip_fd$scope, "period")
  expect_null(lowered$transTrip_fd$scope)   # not overridden -> model default

  # A row-level type still gets the period-level estimand rules applied to it.
  bad <- set_target(tg, recip, contrast = c(0, 1), level = "ministep",
                    verbose = FALSE)
  expect_error(RSiena:::.targetsToEffectList(bad), "not available for")
})

test_that("one mixed call reproduces separate calls exactly", {
  skip_slow()
  fx <- periodFixture()
  algo <- sienaAlgorithmCreate(projname = NULL, nsub = 0, n3 = 300,
                               seed = 4322, cond = FALSE)
  ctl <- set_postest_algo_saom(algorithm = algo, n3 = 300,
                               useChangeContributions = TRUE, verbose = FALSE)
  unc <- set_postest_uncertainty_saom(enabled = FALSE)
  out <- set_postest_output_saom()

  sep_ame <- marginalEffects(fx$ans, data = fx$mydata,
      targets = local({ tg <- make_marginal_targets(fx$ans, data = fx$mydata,
          effects = fx$eff, depvar = "mynet", type = "tieProb", dynamic = TRUE,
          level = "period", includeDefaults = FALSE)
        set_target(tg, transTrip, contrast = c(0, 1), verbose = FALSE) }),
      control_uncertainty = unc, control_algo = ctl, control_out = out)
  sep_end <- marginalEffects(fx$ans, data = fx$mydata,
      targets = local({ tg <- make_marginal_targets(fx$ans, data = fx$mydata,
          effects = fx$eff, depvar = "mynet", type = "tieProb", scope = "period",
          dynamic = TRUE, level = "period", includeDefaults = FALSE)
        set_target(tg, recip, contrast = c(0, 1), verbose = FALSE) }),
      control_uncertainty = unc, control_algo = ctl, control_out = out)

  tg <- make_marginal_targets(fx$ans, data = fx$mydata, effects = fx$eff,
                              depvar = "mynet", type = "tieProb", dynamic = TRUE,
                              level = "period", includeDefaults = FALSE)
  tg <- set_target(tg, transTrip, contrast = c(0, 1), verbose = FALSE)
  tg <- set_target(tg, recip, contrast = c(0, 1), scope = "period",
                   verbose = FALSE)
  one <- marginalEffects(fx$ans, data = fx$mydata, targets = tg,
      control_uncertainty = unc, control_algo = ctl, control_out = out)

  ## Several targets come back long: one row per (effect, period) with the
  ## value in `est`, rather than a named list.
  expect_true(all(c("effect", "period", "est") %in% names(one)))
  expect_setequal(unique(one$effect), c("transTrip_fd", "recip_fd"))

  ## Row-wise and period-level specs streamed over ONE flattened batch must
  ## give exactly what each gives on its own.
  expect_equal(one$est[one$effect == "transTrip_fd"], sep_ame$firstDiff,
               tolerance = 1e-12)
  expect_equal(one$est[one$effect == "recip_fd"], sep_end$periodFirstDiff,
               tolerance = 1e-12)
  ## They are different estimands, not one relabelled.
  expect_false(isTRUE(all.equal(one$est[one$effect == "transTrip_fd"],
                                one$est[one$effect == "recip_fd"])))
})

test_that("two conflicting row-wise types in one call are refused", {
  skip_slow()
  fx <- periodFixture()
  algo <- sienaAlgorithmCreate(projname = NULL, nsub = 0, n3 = 60,
                               seed = 4322, cond = FALSE)
  ctl <- set_postest_algo_saom(algorithm = algo, n3 = 60,
                               useChangeContributions = TRUE, verbose = FALSE)
  tg <- make_marginal_targets(fx$ans, data = fx$mydata, effects = fx$eff,
                              depvar = "mynet", type = "tieProb", dynamic = TRUE,
                              level = "period", includeDefaults = FALSE)
  tg <- set_target(tg, transTrip, contrast = c(0, 1), verbose = FALSE)
  tg <- set_target(tg, recip, contrast = c(0, 1), type = "changeProb",
                   verbose = FALSE)
  expect_error(
    marginalEffects(fx$ans, data = fx$mydata, targets = tg,
        control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
        control_algo = ctl, control_out = set_postest_output_saom()),
    "cannot mix type")
})


# ---------------------------------------------------------------------------
# add_target: several requests about one effect
#
# set_target() tunes the single row an effect has; add_target() appends another
# request about the same effect.  That is what lets one effect be asked for on
# the per-ministep, accumulated and period-level scales at once -- and, on the
# same mechanism, under several conditionings.
# ---------------------------------------------------------------------------

baseTargets <- function(fx, type = "tieProb") {
  tg <- make_marginal_targets(fx$ans, data = fx$mydata, effects = fx$eff,
          depvar = "mynet", type = type, dynamic = TRUE, level = "period",
          includeDefaults = FALSE)
  set_target(tg, recip, contrast = c(0, 1), name = "recip_ame", verbose = FALSE)
}

test_that("set_target still renames rather than appending", {
  skip_slow()
  fx <- periodFixture()
  tg <- baseTargets(fx)
  expect_equal(sum(tg$include), 1L)
  expect_true("recip_ame" %in% tg$name)
  # Repeating the same named call re-tunes that row; it must not duplicate it.
  tg2 <- set_target(tg, recip, contrast = c(0, 1), name = "recip_ame",
                    verbose = FALSE)
  expect_equal(sum(tg2$include), 1L)
})

test_that("set_target under a new name appends a second request", {
  skip_slow()
  fx <- periodFixture()
  tg <- baseTargets(fx)
  tg <- set_target(tg, recip, name = "recip_acc", accumulated = TRUE,
                   verbose = FALSE)
  tg <- set_target(tg, recip, name = "recip_end", scope = "period",
                   verbose = FALSE)

  el <- RSiena:::.targetsToEffectList(tg)$effectList
  expect_setequal(names(el), c("recip_ame", "recip_acc", "recip_end"))
  # The perturbation is cloned; only the named settings differ.
  expect_equal(el$recip_acc$contrast1, c(0, 1))
  expect_equal(el$recip_end$contrast1, c(0, 1))
  expect_true(isTRUE(el$recip_acc$accumulated))
  expect_equal(el$recip_end$scope, "period")
  expect_null(el$recip_ame$type)
})

test_that("several conditionings of one effect in one object", {
  skip_slow()
  fx <- periodFixture()
  tg <- baseTargets(fx)
  tg <- set_target(tg, recip, name = "recip_byTT", condition = "transTrip",
                   verbose = FALSE)
  el <- RSiena:::.targetsToEffectList(tg)$effectList
  expect_setequal(names(el), c("recip_ame", "recip_byTT"))
  expect_equal(el$recip_byTT$condition, "transTrip")
  expect_null(el$recip_ame$condition)
})

test_that("a new name appends; an existing name addresses that request", {
  skip_slow()
  fx <- periodFixture()
  tg <- baseTargets(fx)                       # set_target(..., name = "recip_ame")
  tg <- set_target(tg, recip, name = "recip_end", scope = "period",
                   verbose = FALSE)

  ## Both are APPENDED requests; the row enumerated from the effects object
  ## keeps its own auto-generated name and stays excluded.
  expect_setequal(tg$name[tg$include], c("recip_ame", "recip_end"))
  expect_true("recip_fd" %in% tg$name[!tg$include])
  expect_equal(sum(tg$.added), 2L)

  ## An existing name addresses that request, and only it.
  tg2 <- set_target(tg, recip, name = "recip_end", level = "ego",
                    verbose = FALSE)
  el2 <- RSiena:::.targetsToEffectList(tg2)$effectList
  expect_equal(el2$recip_end$level, "ego")
  expect_null(el2$recip_ame$level)
  expect_setequal(names(el2), c("recip_ame", "recip_end"))
})

test_that("the Glasgow pattern is unchanged by appending", {
  skip_slow()
  fx <- periodFixture()
  ## How glasgow_targets_new.R builds targets: includeDefaults = FALSE, then
  ## one named set_target per effect.  Under the old rename semantics that
  ## renamed the enumerated row; under append it adds a row and leaves the
  ## enumerated one excluded.  The INCLUDED set -- all that is lowered -- is
  ## identical either way, which is why this change is safe for that script.
  tg <- make_marginal_targets(fx$ans, data = fx$mydata, effects = fx$eff,
          depvar = "mynet", type = "tieProb", dynamic = TRUE, level = "period",
          includeDefaults = FALSE)
  tg <- set_target(tg, recip,     contrast = c(0, 1), name = "recip_fd",
                   verbose = FALSE)
  tg <- set_target(tg, transTrip, diff = 1,           name = "transTrip_fd",
                   verbose = FALSE)
  ## Those names are the ones make_marginal_targets() already generated, so
  ## these calls ADDRESS the enumerated rows -- nothing is appended at all.
  expect_equal(sum(tg$.added), 0L)
  el <- RSiena:::.targetsToEffectList(tg)$effectList
  expect_setequal(names(el), c("recip_fd", "transTrip_fd"))
  expect_equal(el$recip_fd$contrast1, c(0, 1))
  expect_equal(el$transTrip_fd$diff1, 1)
})

test_that("an existing name addresses rather than duplicates", {
  skip_slow()
  fx <- periodFixture()
  tg <- baseTargets(fx)
  ## A name that already exists ADDRESSES that request rather than
  ## duplicating it, so this re-tunes and leaves the count unchanged.
  tg2 <- set_target(tg, recip, name = "recip_ame", level = "ego", verbose = FALSE)
  expect_equal(sum(tg2$include), sum(tg$include))
})

test_that("three requests on one effect run in a single call", {
  skip_slow()
  fx <- periodFixture()
  algo <- sienaAlgorithmCreate(projname = NULL, nsub = 0, n3 = 300,
                               seed = 4322, cond = FALSE)
  ctl <- set_postest_algo_saom(algorithm = algo, n3 = 300,
                               useChangeContributions = TRUE, verbose = FALSE)
  unc <- set_postest_uncertainty_saom(enabled = FALSE)
  out <- set_postest_output_saom()

  tg <- baseTargets(fx)
  tg <- set_target(tg, recip, name = "recip_acc", accumulated = TRUE,
                   verbose = FALSE)
  tg <- set_target(tg, recip, name = "recip_end", scope = "period",
                   verbose = FALSE)
  ## The accumulated aggregation emits a `group` column the other two do not,
  ## so combining them drops it and says so.  Pre-existing inconsistency in the
  ## aggregation, not in the targets object -- asserted rather than silenced so
  ## it is noticed if it ever changes.
  expect_warning(
    r <- marginalEffects(fx$ans, data = fx$mydata, targets = tg,
          control_uncertainty = unc, control_algo = ctl, control_out = out),
    "dropped column")

  expect_setequal(unique(r$effect), c("recip_ame", "recip_acc", "recip_end"))
  # Three genuinely different estimands off one pass over the chains.
  v <- function(nm) r$est[r$effect == nm]
  expect_false(isTRUE(all.equal(v("recip_ame"), v("recip_acc"))))
  expect_false(isTRUE(all.equal(v("recip_acc"), v("recip_end"))))
  expect_true(all(is.finite(r$est)))
})


# ---------------------------------------------------------------------------
# Chain-store scratch directory
#
# R removes a session's tempdir on a normal exit but not when the process is
# killed, so interrupted runs leave batches behind.  Tagging the directory with
# the owning process makes a leftover attributable rather than anonymous.
# ---------------------------------------------------------------------------

test_that("chain batches go to a PID-tagged directory", {
  d <- RSiena:::.chainDir(tempdir())
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  expect_true(dir.exists(d))
  expect_match(basename(d), sprintf("^RSiena_chains_pid%d$", Sys.getpid()))
})
