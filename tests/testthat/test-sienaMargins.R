# testthat::skip_on_cran()

# -- Load fixtures from cache (built by helper-models.R) --
ans        <- load_fixture("ans")
mydata     <- load_fixture("mydata")
mymodel    <- load_fixture("mymodel")
mycontrols <- load_fixture("mycontrols")
ans2       <- load_fixture("ans2")
mydata2    <- load_fixture("mydata2")
mymodel2   <- load_fixture("mymodel2")
mycontrols2 <- load_fixture("mycontrols2")
ans3       <- load_fixture("ans3")
mydata3    <- load_fixture("mydata3")
mymodel3   <- load_fixture("mymodel3")
mycontrols3 <- load_fixture("mycontrols3")
# Full-mode fixtures (NULL in quick mode; those tests are guarded by skip_slow()/skip_if())
ans_ego        <- load_fixture("ans_ego")
mydata_ego     <- load_fixture("mydata_ego")
mymodel_ego    <- load_fixture("mymodel_ego")
mycontrols_ego <- load_fixture("mycontrols_ego")
ans_2int       <- load_fixture("ans_2int")
mydata_2int    <- load_fixture("mydata_2int")
mymodel_2int   <- load_fixture("mymodel_2int")
mycontrols_2int <- load_fixture("mycontrols_2int")

# ── Static tests ──────────────────────────────────────────────────────────────

test_that("marginalEffects static: firstDiff structure", {
  tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                             type = "tieProb", level = "egoChoice",
                             condition = c("recip", "density", "transTrip"),
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, transTrip, diff = 1))
  out <- marginalEffects(ans, mydata, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})



test_that("marginalEffects static: secondDiff structure", {
  tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                             type = "tieProb", level = "egoChoice",
                             condition = c("recip", "density", "transTrip"),
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_second_diff(tg,
                                         list(transTrip = list(diff = 1),
                                              recip     = list(contrast = c(0, 1)))))
  out <- marginalEffects(ans, mydata, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
  expect_true(is.data.frame(out))
  expect_true("secondDiff" %in% names(out))
})

test_that("marginalEffects static interaction, no uncertainty", {
  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             type = "tieProb", level = "period",
                             condition = "recip", includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, transTrip, diff = 1,
                                    interaction = TRUE,
                                    intEffectNames = "transRecTrip",
                                    modEffectNames = "recip"))
  out <- marginalEffects(ans2, mydata2, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
      control_algo = set_postest_algo_saom(verbose = FALSE))
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

test_that("marginalEffects static interaction with uncertainty", {
  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             type = "tieProb", level = "period",
                             condition = "recip", includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, transTrip, diff = 1,
                                    interaction = TRUE,
                                    intEffectNames = "transRecTrip",
                                    modEffectNames = "recip"))
  out <- marginalEffects(ans2, mydata2, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(nsim = 5),
      control_algo = set_postest_algo_saom(verbose = FALSE))
  expect_true(is.data.frame(out))
  expect_true("SE" %in% names(out))
})

test_that("marginalEffects static secondDiff with interactions, riskRatio", {
  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             type = "tieProb", level = "period",
                             mainEffect = "riskRatio", includeDefaults = FALSE)
  tg <- suppressMessages(set_second_diff(tg, list(
      recip     = list(contrast = c(0, 1), interaction = TRUE,
                       intEffectNames = "transRecTrip",
                       modEffectNames = "transTrip"),
      transTrip = list(diff = 1, interaction = TRUE,
                       intEffectNames = "transRecTrip",
                       modEffectNames = "recip"))))
  out <- marginalEffects(ans2, mydata2, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(nsim = 5),
      control_algo = set_postest_algo_saom(verbose = FALSE))
  expect_true(is.data.frame(out))
  expect_true("secondRiskRatio" %in% names(out))
})

test_that("marginalEffects static with custom interaction", {
  tg <- make_postest_targets(ans3, effects = mymodel3, depvar = "mynet3",
                             type = "tieProb", level = "period",
                             condition = "inPop", includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, recip, contrast = c(0, 1),
                                    interaction = TRUE,
                                    intEffectNames = "unspInt",
                                    modEffectNames = "inPop"))
  out <- marginalEffects(ans3, mydata3, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(nsim = 5),
      control_algo = set_postest_algo_saom(verbose = FALSE))
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

test_that("marginalEffects dynamic with custom interaction", {
  skip_slow()
  tg <- make_postest_targets(ans3, effects = mymodel3, depvar = "mynet3",
                             dynamic = TRUE,
                             type = "tieProb", level = "period",
                             condition = "inPop", includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, recip, contrast = c(0, 1),
                                    interaction = TRUE,
                                    intEffectNames = "unspInt",
                                    modEffectNames = "inPop"))
  out <- marginalEffects(ans3, mydata3, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
      control_algo = set_postest_algo_saom(
                                           algorithm = mycontrols3, n3 = 60,
                                           verbose = FALSE))
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

# ── Dynamic tests (marginalEffects dynamic = TRUE) ───────────────────────────────────

test_that("marginalEffects dynamic: firstDiff structure, uncertainty=FALSE", {
  skip_slow()
  tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                             dynamic = TRUE,
                             type = "tieProb", condition = "density",
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, transTrip, diff = 1))
  out <- marginalEffects(ans, mydata, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
      control_algo = set_postest_algo_saom(
                                           algorithm = mycontrols, n3 = 60))
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

test_that("marginalEffects dynamic: CI structure", {
  skip_slow()
  tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                             dynamic = TRUE,
                             type = "tieProb", condition = "density",
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, transTrip, diff = 1))
  out <- marginalEffects(ans, mydata, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(
          nsim = 20, sd = FALSE, ci = TRUE, simMean = TRUE, simMedian = TRUE),
      control_algo = set_postest_algo_saom(
                                           algorithm = mycontrols, n3 = 60,
                                           verbose = FALSE))
  expect_true(is.data.frame(out))
  expect_true(all(c("Mean", "q_025", "q_975", "Median") %in% names(out)))
  expect_false("SE" %in% names(out))
})



test_that("marginalEffects dynamic: interaction", {
  skip_slow()
  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             dynamic = TRUE,
                             type = "tieProb", level = "period",
                             condition = c("density", "recip"),
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, transTrip, diff = 1,
                                    interaction = TRUE,
                                    intEffectNames = "transRecTrip",
                                    modEffectNames = "recip"))
  out <- marginalEffects(ans2, mydata2, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
      control_algo = set_postest_algo_saom(
                                           algorithm = mycontrols2, n3 = 60))
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

# ── Ego-perturbation & behavior-DV integration tests ─────────────────────────

test_that("marginalEffects static: egoX auto-detects ego perturbation", {
  skip_slow()
  skip_if(is.null(ans_ego), "ans_ego not fitted (RSENA_FULL_TESTS=1)")
  tg <- make_postest_targets(ans_ego, effects = mymodel_ego,
                             depvar = "mynet_ego", type = "tieProb",
                             level = "period", condition = "density",
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, egoX, diff = 1))
  out <- marginalEffects(ans_ego, mydata_ego, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
  # Mass contrasts should be present for ego perturbation
  expect_true("massCreation" %in% names(out))
  expect_true("massDissolution" %in% names(out))
})

test_that("marginalEffects static: egoX mass contrasts sum to dyad-level diffs", {
  skip_slow()
  skip_if(is.null(ans_ego), "ans_ego not fitted (RSENA_FULL_TESTS=1)")
  tg <- make_postest_targets(ans_ego, effects = mymodel_ego,
                             depvar = "mynet_ego", type = "changeProb",
                             level = "none", includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, egoX, diff = 1))
  out <- marginalEffects(ans_ego, mydata_ego, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
  # For each ego×period, massCreation should equal sum of firstDiff
  # for creation rows (density_eval == 1) in that group.
  densName <- grep("density", names(out), value = TRUE, fixed = TRUE)[1]
  for (p in unique(out$period)) {
    for (e in unique(out$ego)) {
      rows <- out$period == p & out$ego == e
      cre <- rows & out[[densName]] == 1
      dis <- rows & out[[densName]] == -1
      if (any(cre)) {
        expect_equal(
          out$massCreation[which(rows)[1]],
          sum(out$firstDiff[cre]),
          tolerance = 1e-10, info = paste("ego", e, "period", p, "creation")
        )
      }
      if (any(dis)) {
        expect_equal(
          out$massDissolution[which(rows)[1]],
          sum(out$firstDiff[dis]),
          tolerance = 1e-10, info = paste("ego", e, "period", p, "dissolution")
        )
      }
    }
  }
})

test_that("marginalEffects static: egoX tieProb mass contrasts use changeProbDiff", {
  skip_slow()
  skip_if(is.null(ans_ego), "ans_ego not fitted (RSENA_FULL_TESTS=1)")
  # Compute both changeProb and tieProb versions
  tg_cp <- make_postest_targets(ans_ego, effects = mymodel_ego,
                                depvar = "mynet_ego", type = "changeProb",
                                level = "none", includeDefaults = FALSE)
  tg_cp <- suppressMessages(set_target(tg_cp, egoX, diff = 1))
  out_cp <- marginalEffects(ans_ego, mydata_ego, targets = tg_cp,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
  tg_tp <- make_postest_targets(ans_ego, effects = mymodel_ego,
                                depvar = "mynet_ego", type = "tieProb",
                                level = "none", includeDefaults = FALSE)
  tg_tp <- suppressMessages(set_target(tg_tp, egoX, diff = 1))
  out_tp <- marginalEffects(ans_ego, mydata_ego, targets = tg_tp,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
  # Mass contrasts should be identical regardless of type, because
  # they are always on the changeProbDiff scale
  expect_equal(out_tp$massCreation, out_cp$massCreation, tolerance = 1e-10)
  expect_equal(out_tp$massDissolution, out_cp$massDissolution, tolerance = 1e-10)
})

test_that("marginalEffects static: alter perturbType does NOT have mass columns", {
  skip_slow()
  skip_if(is.null(ans_ego), "ans_ego not fitted (RSENA_FULL_TESTS=1)")
  tg <- make_postest_targets(ans_ego, effects = mymodel_ego,
                             depvar = "mynet_ego", type = "tieProb",
                             level = "period", condition = "density",
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, egoX, diff = 1,
                                    perturbType = "alter"))
  out <- marginalEffects(ans_ego, mydata_ego, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
  # Mass contrasts should NOT be present for alter perturbation
  expect_false("massCreation" %in% names(out))
  expect_false("massDissolution" %in% names(out))
})

test_that("marginalEffects static: egoX interaction firstDiff", {
  skip_slow()
  skip_if(is.null(ans_ego), "ans_ego not fitted (RSENA_FULL_TESTS=1)")
  tg <- make_postest_targets(ans_ego, effects = mymodel_ego,
                             depvar = "mynet_ego", type = "tieProb",
                             level = "period",
                             condition = c("density", "transTrip"),
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, egoX, diff = 1,
                                    interaction = TRUE,
                                    intEffectNames = "unspInt",
                                    modEffectNames = "transTrip"))
  out <- marginalEffects(ans_ego, mydata_ego, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

test_that("marginalEffects dynamic: egoX auto-detects ego perturbation", {
  skip_slow()
  skip_if(is.null(ans_ego), "ans_ego not fitted (RSENA_FULL_TESTS=1)")
  tg <- make_postest_targets(ans_ego, effects = mymodel_ego,
                             depvar = "mynet_ego", type = "tieProb",
                             condition = "density", includeDefaults = FALSE, dynamic = TRUE)
  tg <- suppressMessages(set_target(tg, egoX, diff = 1))
  out <- marginalEffects(ans_ego, mydata_ego, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
      control_algo = set_postest_algo_saom(
                                           algorithm = mycontrols_ego,
                                           n3 = 60))
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

# NOTE: left in flat form deliberately -- see migration report. This test
# checks the early depvar-type guard in marginalEffects.sienaFit, which fires
# from `data`/`depvar` alone before `effects` is ever consulted. It relies on
# `depvar` disagreeing with what `effects` (mymodel_ego, built for
# "mynet_ego") has registered -- data here carries "mybeh_b", a behavior
# variable that does not exist in mymodel_ego at all. The object form cannot
# express this: make_postest_targets() requires depvar's effects to be
# present in `effects` and errors immediately ("No included effects for
# dependent variable 'mybeh_b'"), never reaching the behavior-DV guard the
# test is checking. Attempting to route around it by giving the targets
# object the model's real depvar ("mynet_ego") also fails the test's intent:
# marginalEffects() then reads attr(targets, "depvar") = "mynet_ego", looks up
# data[["depvars"]][["mynet_ego"]] (absent from mydata_b, which only has
# "mynet_b"/"mybeh_b"), gets NULL, and the behavior guard silently never
# fires. Either route changes what the test verifies, which the task
# forbids -- so this one call stays on the flat form.
test_that("marginalEffects: behavior DV stops with informative error", {
  skip_slow()
  skip_if(is.null(ans_ego), "ans_ego not fitted (RSENA_FULL_TESTS=1)")
  # Build a data object that has a behavior dependent variable,
  # so the depvar type guard is actually triggered.
  mynet_b <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
  mybeh_b <- sienaDependent(s50a, type = "behavior")
  mydata_b <- sienaDataCreate(mynet_b, mybeh_b)
  ## depvar names a behaviour variable, which the targets object cannot be
  ## built for either -- so the guard is checked wherever it now fires.
  expect_error(
    suppressMessages(marginalEffects(ans_ego, mydata_b,
      targets = make_postest_targets(ans_ego, effects = mymodel_ego,
                                     depvar = "mybeh_b", type = "tieProb",
                                     level = "period", condition = "density"),
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))),
    "[Bb]ehavio(u)?r|mybeh_b"
  )
})



test_that("two unspInt: firstDiff via unspInt1 (main effect = recip)", {
  skip_if(is.null(ans_2int), "ans_2int not fitted (RSENA_FULL_TESTS not set)")
  tg <- make_postest_targets(ans_2int, effects = mymodel_2int,
                             depvar = "mynet_2int", type = "tieProb",
                             level = "period",
                             condition = c("inPop", "density"),
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, recip, contrast = c(0, 1),
                                    interaction = TRUE,
                                    intEffectNames = "unspInt1",
                                    modEffectNames = "inPop"))
  out <- marginalEffects(ans_2int, mydata_2int, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
  expect_true(nrow(out) > 0L)
})

test_that("two unspInt: firstDiff via unspInt2 (main effect = recip)", {
  skip_if(is.null(ans_2int), "ans_2int not fitted (RSENA_FULL_TESTS not set)")
  tg <- make_postest_targets(ans_2int, effects = mymodel_2int,
                             depvar = "mynet_2int", type = "tieProb",
                             level = "period",
                             condition = c("outPop", "density"),
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, recip, contrast = c(0, 1),
                                    interaction = TRUE,
                                    intEffectNames = "unspInt2",
                                    modEffectNames = "outPop"))
  out <- marginalEffects(ans_2int, mydata_2int, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
  expect_true(nrow(out) > 0L)
})

test_that("two unspInt: secondDiff across unspInt1 and unspInt2", {
  skip_if(is.null(ans_2int), "ans_2int not fitted (RSENA_FULL_TESTS not set)")
  tg <- make_postest_targets(ans_2int, effects = mymodel_2int,
                             depvar = "mynet_2int", type = "tieProb",
                             level = "period",
                             condition = c("inPop", "outPop", "density"),
                             includeDefaults = FALSE)
  ## Same base effect crossed with itself, entering through two different
  ## unspInt terms -- so the list form's names repeat, which is allowed.
  tg <- suppressMessages(set_second_diff(tg, list(
      recip = list(contrast = c(0, 1), interaction = TRUE,
                   intEffectNames = "unspInt1", modEffectNames = "inPop"),
      recip = list(contrast = c(0, 1), interaction = TRUE,
                   intEffectNames = "unspInt2", modEffectNames = "outPop"))))
  out <- marginalEffects(ans_2int, mydata_2int, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
  expect_true(is.data.frame(out))
  expect_true("secondDiff" %in% names(out))
  expect_true(nrow(out) > 0L)
})

test_that("two unspInt: conditional prediction (tieProb) with uncertainty", {
  skip_if(is.null(ans_2int), "ans_2int not fitted (RSENA_FULL_TESTS not set)")
  tg <- make_postest_targets(ans_2int, effects = mymodel_2int,
                             depvar = "mynet_2int", type = "tieProb",
                             level = "period",
                             condition = c("inPop", "density"),
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, recip, contrast = c(0, 1),
                                    interaction = TRUE,
                                    intEffectNames = "unspInt1",
                                    modEffectNames = "inPop"))
  out <- marginalEffects(ans_2int, mydata_2int, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(nsim = 5),
      control_algo = set_postest_algo_saom(verbose = FALSE))
  expect_true(is.data.frame(out))
  expect_true("SE" %in% names(out))
})

# ── effectList: multi-effect batch path ──────────────────────────────────────

test_that("effectList batch: targets sharing (level, condition) merge into one data frame with an effect column", {
  set.seed(1)
  tg_recip <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                                   type = "tieProb", level = "period",
                                   condition = "density",
                                   includeDefaults = FALSE)
  tg_recip <- suppressMessages(set_target(tg_recip, recip,
                                          contrast = c(0, 1)))
  scalar_recip <- marginalEffects(ans, mydata, targets = tg_recip,
      control_uncertainty = set_postest_uncertainty_saom(nsim = 20),
      control_algo = set_postest_algo_saom(verbose = FALSE))

  tg_trans <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                                   type = "tieProb", level = "period",
                                   condition = "density",
                                   includeDefaults = FALSE)
  tg_trans <- suppressMessages(set_target(tg_trans, transTrip, diff = 1))
  scalar_trans <- marginalEffects(ans, mydata, targets = tg_trans,
      control_uncertainty = set_postest_uncertainty_saom(nsim = 20),
      control_algo = set_postest_algo_saom(verbose = FALSE))

  set.seed(1)
  tg_batch <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                                   type = "tieProb", level = "period",
                                   condition = "density",
                                   includeDefaults = FALSE)
  tg_batch <- suppressMessages(set_target(tg_batch, recip,
                                          contrast = c(0, 1),
                                          name = "recip_fd"))
  tg_batch <- suppressMessages(set_target(tg_batch, transTrip, diff = 1,
                                          name = "trans_fd"))
  batch <- marginalEffects(ans, mydata, targets = tg_batch,
      control_uncertainty = set_postest_uncertainty_saom(nsim = 20),
      control_algo = set_postest_algo_saom(verbose = FALSE))

  # recip_fd and trans_fd share (level, condition), so combineSameLevel =
  # TRUE (the default) merges them into a single data frame with an
  # `effect` column identifying which target each row belongs to, rather
  # than a two-element named list. Both are firstDiff-type quantities, so
  # their (differently-named) quantity columns are unified into `est`;
  # SE/q_025/q_975 are shared column names already and pass through as-is.
  expect_true(is.data.frame(batch))
  expect_setequal(unique(batch$effect), c("recip_fd", "trans_fd"))

  # Both effects carry the expected uncertainty columns (checked per
  # effect via the `effect` column, not merely present somewhere in the
  # frame).
  expect_true("est" %in% names(batch))
  expect_true("SE"  %in% names(batch))
  expect_true(all(is.finite(batch$est[batch$effect == "recip_fd"])))
  expect_true(all(is.finite(batch$SE[batch$effect == "recip_fd"])))
  expect_true(all(is.finite(batch$est[batch$effect == "trans_fd"])))
  expect_true(all(is.finite(batch$SE[batch$effect == "trans_fd"])))

  # Row counts match the scalar equivalents
  expect_equal(nrow(batch[batch$effect == "recip_fd", ]), nrow(scalar_recip))
  expect_equal(nrow(batch[batch$effect == "trans_fd", ]), nrow(scalar_trans))
})

# NOTE: left in flat form deliberately. A single set_target() on a `targets`
# object ALSO collapses to a bare data.frame now (marginalEffects.sienaFit
# unwraps any single-group `results` list via `if (length(results) == 1L)
# return(results[[1L]])`.  The `.single_effect` flag that used to drive this
# belonged to the flat effectName1 path and went with it in step 5c, so the
# unwrapping is now a property of the one remaining path.
test_that("a single target returns a data frame, not a one-element list", {
  tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                             type = "tieProb", level = "period",
                             condition = "density", includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, recip, contrast = c(0, 1)))
  out <- marginalEffects(ans, mydata, targets = tg,
    control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
    control_algo = set_postest_algo_saom(verbose = FALSE))
  expect_true(is.data.frame(out))
  expect_false(is.list(out) && !is.data.frame(out))
})

# "fd" (firstDiff) and "sd" (secondDiff) share the same (level, condition),
# so they merge under the default combineSameLevel = TRUE: each target's
# quantity column (firstDiff / secondDiff) is renamed to a common `est`
# column before rbind, and which kind of diff a row holds is carried by the
# `effect` column instead. One merged table, both kinds distinguished by
# `effect` -- not a two-element named list.
test_that("effectList: secondDiff mixed with firstDiff in one batch", {
  tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                             type = "tieProb", level = "period",
                             condition = "density", includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, recip, contrast = c(0, 1),
                                    name = "fd"))
  tg <- suppressMessages(set_second_diff(tg,
      list(transTrip = list(diff = 1), recip = list(contrast = c(0, 1))),
      name = "sd"))
  batch <- marginalEffects(ans, mydata, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE,
                                                          nsim = 10),
      control_algo = set_postest_algo_saom(verbose = FALSE))
  expect_true(is.data.frame(batch))
  expect_setequal(unique(batch$effect), c("fd", "sd"))
  expect_true("est" %in% names(batch))
  expect_true(nrow(batch[batch$effect == "fd", ]) > 0L)
  expect_true(nrow(batch[batch$effect == "sd", ]) > 0L)
  expect_true(all(is.finite(batch$est[batch$effect == "fd"])))
  expect_true(all(is.finite(batch$est[batch$effect == "sd"])))
})

# recip_fd and trans_fd share (level, condition) here too, so this batch
# merges into one data frame (effect column: recip_fd / trans_fd) exactly
# like the static case, rather than a named list -- the dynamic-specific
# thing being tested (forward sims shared across effects) is orthogonal to
# the shape and still exercised by running both targets in one call.
test_that("effectList dynamic: shared forward sims merge into one data frame", {
  skip_slow()
  tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                             dynamic = TRUE,
                             type = "tieProb", level = "period",
                             condition = "density", includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, recip, contrast = c(0, 1),
                                    name = "recip_fd"))
  tg <- suppressMessages(set_target(tg, transTrip, diff = 1,
                                    name = "trans_fd"))
  batch <- marginalEffects(ans, mydata, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE,
                                                          nsim = 10),
      control_algo = set_postest_algo_saom(
                                           algorithm = mycontrols, n3 = 50,
                                           verbose = FALSE))
  expect_true(is.data.frame(batch))
  expect_setequal(unique(batch$effect), c("recip_fd", "trans_fd"))
  expect_true("est" %in% names(batch))
  expect_true(all(is.finite(batch$est[batch$effect == "recip_fd"])))
  expect_true(all(is.finite(batch$est[batch$effect == "trans_fd"])))
})

