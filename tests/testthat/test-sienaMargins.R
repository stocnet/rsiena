# testthat::skip_on_cran()
# Models: ans, mydata, mymodel, mycontrols  — from helper-models.R (base)
#         ans2, mydata2, mymodel2            — from helper-models.R (transRecTrip)
#         ans3, mydata3, mymodel3            — from helper-models.R (unspInt)

# ── Static tests ──────────────────────────────────────────────────────────────

test_that("marginalEffects static: firstDiff structure", {
  out <- marginalEffects(
    object = ans, data = mydata,
    effectName1 = "transTrip", diff1 = 1,
    type = "tieProb", depvar = "mynet",
    level = "egoChoice",
    condition = c("recip", "density", "transTrip"),
    uncertainty = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})



test_that("marginalEffects static with uncertainty: MCSE columns, no CI", {
  out <- marginalEffects(
    object = ans, data = mydata,
    effectName1 = "recip", contrast1 = c(0, 1),
    type = "tieProb", depvar = "mynet",
    level = "period", condition = "density",
    nsim = 20, uncertainty = TRUE,
    uncertaintyMean = TRUE,
    uncertaintyMcse = TRUE, uncertaintymcseBatches = 4,
    uncertaintySd = TRUE, uncertaintyCi = FALSE,
    verbose = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true(all(c("Mean", "SE", "cases", "mcse_Mean", "mcse_SE") %in% names(out)))
  expect_false("q_025" %in% names(out))
})

test_that("marginalEffects static: secondDiff structure", {
  out <- marginalEffects(
    object = ans, data = mydata,
    effectName1 = "transTrip", diff1 = 1,
    effectName2 = "recip", contrast2 = c(0, 1),
    second = TRUE, type = "tieProb", depvar = "mynet",
    level = "egoChoice",
    condition = c("recip", "density", "transTrip"),
    uncertainty = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("secondDiff" %in% names(out))
})

test_that("marginalEffects static interaction, no uncertainty", {
  out <- marginalEffects(
    object = ans2, data = mydata2,
    effectName1 = "transTrip", diff1 = 1,
    interaction1 = TRUE, intEffectNames1 = "transRecTrip",
    modEffectNames1 = "recip",
    type = "tieProb", depvar = "mynet2",
    level = "period", condition = "recip",
    uncertainty = FALSE, verbose = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

test_that("marginalEffects static interaction with uncertainty", {
  out <- marginalEffects(
    object = ans2, data = mydata2,
    effectName1 = "transTrip", diff1 = 1,
    interaction1 = TRUE, intEffectNames1 = "transRecTrip",
    modEffectNames1 = "recip",
    type = "tieProb", depvar = "mynet2",
    level = "period", condition = "recip",
    nsim = 5, uncertainty = TRUE, verbose = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("SE" %in% names(out))
})

test_that("marginalEffects static secondDiff with interactions, riskRatio", {
  out <- marginalEffects(
    object = ans2, data = mydata2,
    effectName1 = "recip", contrast1 = c(0, 1),
    interaction1 = TRUE, intEffectNames1 = "transRecTrip",
    modEffectNames1 = "transTrip",
    second = TRUE,
    effectName2 = "transTrip", diff2 = 1,
    interaction2 = TRUE, intEffectNames2 = "transRecTrip",
    modEffectNames2 = "recip",
    type = "tieProb", depvar = "mynet2",
    level = "period", nsim = 5, uncertainty = TRUE,
    verbose = FALSE, mainEffect = "riskRatio"
  )
  expect_true(is.data.frame(out))
  expect_true("secondRiskRatio" %in% names(out))
})

test_that("marginalEffects static with custom interaction", {
  out <- marginalEffects(
    object = ans3, data = mydata3, effects = mymodel3,
    effectName1 = "recip", contrast1 = c(0, 1),
    interaction1 = TRUE, intEffectNames1 = "unspInt",
    modEffectNames1 = "inPop",
    type = "tieProb", depvar = "mynet3",
    level = "period", condition = "inPop",
    nsim = 5, uncertainty = TRUE, verbose = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

test_that("marginalEffects dynamic with custom interaction", {
  skip_slow()
  out <- marginalEffects(
    object = ans3, data = mydata3, effects = mymodel3,
    effectName1 = "recip", contrast1 = c(0, 1),
    interaction1 = TRUE, intEffectNames1 = "unspInt",
    modEffectNames1 = "inPop",
    type = "tieProb", depvar = "mynet3",
    level = "period", condition = "inPop",
    dynamic = TRUE, algorithm = mycontrols3, n3 = 60,
    uncertainty = FALSE, verbose = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

# ── Dynamic tests (marginalEffects dynamic = TRUE) ───────────────────────────────────

test_that("marginalEffects dynamic: firstDiff structure, uncertainty=FALSE", {
  skip_slow()
  out <- marginalEffects(
    object = ans, data = mydata,
    effectName1 = "transTrip", diff1 = 1,
    effects = mymodel, algorithm = mycontrols,
    dynamic = TRUE, n3 = 60,
    type = "tieProb", condition = "density",
    uncertainty = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

test_that("marginalEffects dynamic: MCSE + CI structure", {
  skip_slow()
  out <- marginalEffects(
    object = ans, data = mydata,
    effectName1 = "transTrip", diff1 = 1,
    effects = mymodel, algorithm = mycontrols,
    dynamic = TRUE, n3 = 60, nsim = 20,
    type = "tieProb", condition = "density",
    uncertainty = TRUE,
    uncertaintyMcse = TRUE, uncertaintymcseBatches = 4,
    uncertaintySd = FALSE, uncertaintyCi = TRUE,
    verbose = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true(all(c("Mean", "cases", "mcse_Mean", "q_025", "q_975", "Median") %in% names(out)))
  expect_false("SE" %in% names(out))
  expect_false("mcse_SE" %in% names(out))
})



test_that("marginalEffects dynamic: interaction", {
  skip_slow()
  out <- marginalEffects(
    object = ans2, data = mydata2,
    effectName1 = "transTrip", diff1 = 1,
    interaction1 = TRUE, intEffectNames1 = "transRecTrip",
    modEffectNames1 = "recip",
    effects = mymodel2, algorithm = mycontrols2,
    dynamic = TRUE, n3 = 60,
    type = "tieProb", level = "period",
    condition = c("density", "recip"),
    uncertainty = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

# ── Ego-perturbation & behavior-DV integration tests ─────────────────────────

test_that("marginalEffects static: egoX auto-detects ego perturbation", {
  skip_slow()
  skip_if(is.null(ans_ego), "ans_ego not fitted (RSENA_FULL_TESTS=1)")
  out <- marginalEffects(
    object = ans_ego, data = mydata_ego,
    effectName1 = "egoX", diff1 = 1,
    type = "tieProb", depvar = "mynet_ego",
    level = "period", condition = "density",
    uncertainty = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
  # Mass contrasts should be present for ego perturbation
  expect_true("massCreation" %in% names(out))
  expect_true("massDissolution" %in% names(out))
})

test_that("marginalEffects static: egoX mass contrasts sum to dyad-level diffs", {
  skip_slow()
  skip_if(is.null(ans_ego), "ans_ego not fitted (RSENA_FULL_TESTS=1)")
  out <- marginalEffects(
    object = ans_ego, data = mydata_ego,
    effectName1 = "egoX", diff1 = 1,
    type = "changeProb", depvar = "mynet_ego",
    level = "none",
    uncertainty = FALSE
  )
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
  out_cp <- marginalEffects(
    object = ans_ego, data = mydata_ego,
    effectName1 = "egoX", diff1 = 1,
    type = "changeProb", depvar = "mynet_ego",
    level = "none", uncertainty = FALSE
  )
  out_tp <- marginalEffects(
    object = ans_ego, data = mydata_ego,
    effectName1 = "egoX", diff1 = 1,
    type = "tieProb", depvar = "mynet_ego",
    level = "none", uncertainty = FALSE
  )
  # Mass contrasts should be identical regardless of type, because
  # they are always on the changeProbDiff scale
  expect_equal(out_tp$massCreation, out_cp$massCreation, tolerance = 1e-10)
  expect_equal(out_tp$massDissolution, out_cp$massDissolution, tolerance = 1e-10)
})

test_that("marginalEffects static: alter perturbType does NOT have mass columns", {
  skip_slow()
  skip_if(is.null(ans_ego), "ans_ego not fitted (RSENA_FULL_TESTS=1)")
  out <- marginalEffects(
    object = ans_ego, data = mydata_ego,
    effectName1 = "egoX", diff1 = 1,
    perturbType1 = "alter",
    type = "tieProb", depvar = "mynet_ego",
    level = "period", condition = "density",
    uncertainty = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
  # Mass contrasts should NOT be present for alter perturbation
  expect_false("massCreation" %in% names(out))
  expect_false("massDissolution" %in% names(out))
})

test_that("marginalEffects static: egoX interaction firstDiff", {
  skip_slow()
  skip_if(is.null(ans_ego), "ans_ego not fitted (RSENA_FULL_TESTS=1)")
  out <- marginalEffects(
    object = ans_ego, data = mydata_ego,
    effectName1 = "egoX", diff1 = 1,
    interaction1 = TRUE, intEffectNames1 = "unspInt",
    modEffectNames1 = "transTrip",
    type = "tieProb", depvar = "mynet_ego",
    level = "period", condition = c("density", "transTrip"),
    uncertainty = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

test_that("marginalEffects dynamic: egoX auto-detects ego perturbation", {
  skip_slow()
  skip_if(is.null(ans_ego), "ans_ego not fitted (RSENA_FULL_TESTS=1)")
  out <- marginalEffects(
    object = ans_ego, data = mydata_ego,
    effectName1 = "egoX", diff1 = 1,
    effects = mymodel_ego, algorithm = mycontrols_ego,
    dynamic = TRUE, n3 = 60,
    type = "tieProb", condition = "density",
    uncertainty = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

test_that("marginalEffects: behavior DV stops with informative error", {
  skip_slow()
  skip_if(is.null(ans_ego), "ans_ego not fitted (RSENA_FULL_TESTS=1)")
  # Build a data object that has a behavior dependent variable,
  # so the depvar type guard is actually triggered.
  mynet_b <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
  mybeh_b <- sienaDependent(s50a, type = "behavior")
  mydata_b <- sienaDataCreate(mynet_b, mybeh_b)
  expect_error(
    marginalEffects(
      object = ans_ego, data = mydata_b,
      effectName1 = "egoX", diff1 = 1,
      type = "tieProb", depvar = "mybeh_b",
      level = "period", condition = "density",
      uncertainty = FALSE
    ),
    "behavior"
  )
})



test_that("two unspInt: firstDiff via unspInt1 (main effect = recip)", {
  skip_if(is.null(ans_2int), "ans_2int not fitted (RSENA_FULL_TESTS not set)")
  out <- marginalEffects(
    object = ans_2int, data = mydata_2int, effects = mymodel_2int,
    effectName1 = "recip", contrast1 = c(0, 1),
    interaction1 = TRUE, intEffectNames1 = "unspInt1",
    modEffectNames1 = "inPop",
    type = "tieProb", depvar = "mynet_2int",
    level = "period", condition = c("inPop", "density"),
    uncertainty = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
  expect_true(nrow(out) > 0L)
})

test_that("two unspInt: firstDiff via unspInt2 (main effect = recip)", {
  skip_if(is.null(ans_2int), "ans_2int not fitted (RSENA_FULL_TESTS not set)")
  out <- marginalEffects(
    object = ans_2int, data = mydata_2int, effects = mymodel_2int,
    effectName1 = "recip", contrast1 = c(0, 1),
    interaction1 = TRUE, intEffectNames1 = "unspInt2",
    modEffectNames1 = "outPop",
    type = "tieProb", depvar = "mynet_2int",
    level = "period", condition = c("outPop", "density"),
    uncertainty = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
  expect_true(nrow(out) > 0L)
})

test_that("two unspInt: secondDiff across unspInt1 and unspInt2", {
  skip_if(is.null(ans_2int), "ans_2int not fitted (RSENA_FULL_TESTS not set)")
  out <- marginalEffects(
    object = ans_2int, data = mydata_2int, effects = mymodel_2int,
    effectName1 = "recip", contrast1 = c(0, 1),
    interaction1 = TRUE, intEffectNames1 = "unspInt1",
    modEffectNames1 = "inPop",
    second = TRUE,
    effectName2 = "recip", contrast2 = c(0, 1),
    interaction2 = TRUE, intEffectNames2 = "unspInt2",
    modEffectNames2 = "outPop",
    type = "tieProb", depvar = "mynet_2int",
    level = "period", condition = c("inPop", "outPop", "density"),
    uncertainty = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("secondDiff" %in% names(out))
  expect_true(nrow(out) > 0L)
})

test_that("two unspInt: conditional prediction (tieProb) with uncertainty", {
  skip_if(is.null(ans_2int), "ans_2int not fitted (RSENA_FULL_TESTS not set)")
  out <- marginalEffects(
    object = ans_2int, data = mydata_2int, effects = mymodel_2int,
    effectName1 = "recip", contrast1 = c(0, 1),
    interaction1 = TRUE, intEffectNames1 = "unspInt1",
    modEffectNames1 = "inPop",
    type = "tieProb", depvar = "mynet_2int",
    level = "period", condition = c("inPop", "density"),
    nsim = 5, uncertainty = TRUE, verbose = FALSE
  )
  expect_true(is.data.frame(out))
  expect_true("SE" %in% names(out))
})

# ── effectList: multi-effect batch path ──────────────────────────────────────

test_that("effectList returns named list with same structure as scalar calls", {
  set.seed(1)
  scalar_recip <- marginalEffects(
    object = ans, data = mydata,
    effectName1 = "recip", contrast1 = c(0, 1),
    type = "tieProb", level = "period", condition = "density",
    nsim = 20, uncertainty = TRUE, verbose = FALSE
  )
  scalar_trans <- marginalEffects(
    object = ans, data = mydata,
    effectName1 = "transTrip", diff1 = 1,
    type = "tieProb", level = "period", condition = "density",
    nsim = 20, uncertainty = TRUE, verbose = FALSE
  )

  set.seed(1)
  batch <- marginalEffects(
    object = ans, data = mydata,
    type = "tieProb", level = "period", condition = "density",
    nsim = 20, uncertainty = TRUE, verbose = FALSE,
    effectList = list(
      recip_fd  = list(effectName1 = "recip",     contrast1 = c(0, 1)),
      trans_fd  = list(effectName1 = "transTrip", diff1     = 1)
    )
  )

  # Returns a named list
  expect_true(is.list(batch))
  expect_setequal(names(batch), c("recip_fd", "trans_fd"))

  # Each element is a data frame with the expected uncertainty columns
  expect_true(is.data.frame(batch$recip_fd))
  expect_true(is.data.frame(batch$trans_fd))
  expect_true("firstDiff" %in% names(batch$recip_fd))
  expect_true("SE"        %in% names(batch$recip_fd))
  expect_true("firstDiff" %in% names(batch$trans_fd))
  expect_true("SE"        %in% names(batch$trans_fd))

  # Row counts match the scalar equivalents
  expect_equal(nrow(batch$recip_fd), nrow(scalar_recip))
  expect_equal(nrow(batch$trans_fd), nrow(scalar_trans))
})

test_that("effectList scalar call (single element) returns data frame, not list", {
  out <- marginalEffects(
    object = ans, data = mydata,
    effectName1 = "recip", contrast1 = c(0, 1),
    type = "tieProb", level = "period", condition = "density",
    nsim = 10, uncertainty = FALSE, verbose = FALSE
  )
  expect_true(is.data.frame(out))
  expect_false(is.list(out) && !is.data.frame(out))
})

test_that("effectList: secondDiff mixed with firstDiff in one batch", {
  batch <- marginalEffects(
    object = ans, data = mydata,
    type = "tieProb", level = "period", condition = "density",
    nsim = 10, uncertainty = FALSE, verbose = FALSE,
    effectList = list(
      fd   = list(effectName1 = "recip",    contrast1 = c(0, 1)),
      sd   = list(effectName1 = "transTrip", diff1 = 1,
                  effectName2 = "recip",    contrast2 = c(0, 1),
                  second = TRUE)
    )
  )
  expect_true(is.list(batch))
  expect_true("firstDiff"  %in% names(batch$fd))
  expect_true("secondDiff" %in% names(batch$sd))
})

test_that("effectList dynamic: shared forward sims, returns named list", {
  skip_slow()
  batch <- marginalEffects(
    object = ans, data = mydata,
    effects = mymodel, algorithm = mycontrols,
    type = "tieProb", level = "period", condition = "density",
    dynamic = TRUE, n3 = 50, nsim = 10,
    uncertainty = FALSE, verbose = FALSE,
    effectList = list(
      recip_fd = list(effectName1 = "recip",    contrast1 = c(0, 1)),
      trans_fd = list(effectName1 = "transTrip", diff1    = 1)
    )
  )
  expect_true(is.list(batch))
  expect_true("firstDiff" %in% names(batch$recip_fd))
  expect_true("firstDiff" %in% names(batch$trans_fd))
})

