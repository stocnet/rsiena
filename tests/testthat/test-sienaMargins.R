# testthat::skip_on_cran()

# Basic model: density + recip + transTrip
mynet <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata <- sienaDataCreate(mynet)
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, transTrip, name = "mynet")
mycontrols <- sienaAlgorithmCreate(projname = NULL, n3 = 50, cond = FALSE)
ans <- siena07(
  mycontrols,
  data = mydata,
  effects = mymodel,
  returnChangeContributions = TRUE,
  returnDataFrame = TRUE
)

# Model with two structural effects (for interaction/moderator tests)
mynet2 <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata2 <- sienaDataCreate(mynet2)
mymodel2 <- getEffects(mydata2)
mymodel2 <- includeEffects(mymodel2, transTrip, transRecTrip, name = "mynet2")
mycontrols2 <- sienaAlgorithmCreate(projname = NULL, n3 = 60, cond = FALSE)
ans2 <- siena07(
  mycontrols2,
  data = mydata2,
  effects = mymodel2,
  returnChangeContributions = TRUE,
  returnDataFrame = TRUE
)

# Model with a custom interaction (inPop * recip via unspInt)
mynet3 <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata3 <- sienaDataCreate(mynet3)
mymodel3 <- getEffects(mydata3)
mymodel3 <- includeEffects(mymodel3, inPop, name = "mynet3")
mymodel3 <- includeInteraction(mymodel3, recip, inPop, name = "mynet3")
mycontrols3 <- sienaAlgorithmCreate(projname = NULL, n3 = 60, cond = FALSE)
ans3 <- siena07(
  mycontrols3,
  data = mydata3,
  effects = mymodel3,
  returnChangeContributions = TRUE,
  returnDataFrame = TRUE
)

# ── Static tests ──────────────────────────────────────────────────────────────

test_that("sienaAME static: firstDiff structure (base R)", {
  with_mocked_bindings(
    {
      out <- sienaAME(
        ans = ans, data = mydata,
        effectName1 = "transTrip", diff1 = 1,
        type = "tieProb", depvar = "mynet",
        level = "egoChoice",
        condition = c("recip", "density", "transTrip"),
        uncertainty = FALSE
      )
      expect_true(is.data.frame(out))
      expect_false("data.table" %in% class(out))
      expect_true("firstDiff" %in% names(out))
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
})

test_that("sienaAME static: data.table output", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  out <- sienaAME(
    ans = ans, data = mydata,
    effectName1 = "transTrip", diff1 = 1,
    type = "tieProb", depvar = "mynet",
    level = "period", condition = "density",
    uncertainty = FALSE
  )
  expect_true("data.table" %in% class(out))
  expect_true("firstDiff" %in% names(out))
})

test_that("sienaAME static with uncertainty: MCSE columns, no CI (base R)", {
  with_mocked_bindings(
    {
      out <- sienaAME(
        ans = ans, data = mydata,
        effectName1 = "recip", contrast1 = c(0, 1),
        type = "tieProb", depvar = "mynet",
        level = "period", condition = "density",
        nsim = 20, uncertainty = TRUE,
        uncertainty_mcse = TRUE, uncertainty_mcse_batches = 4,
        uncertainty_sd = TRUE, uncertainty_ci = FALSE,
        verbose = FALSE
      )
      expect_true(is.data.frame(out))
      expect_true(all(c("Mean", "SE", "cases", "mcse_Mean", "mcse_SE") %in% names(out)))
      expect_false("q_025" %in% names(out))
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
})

test_that("sienaAME static: secondDiff structure (base R)", {
  with_mocked_bindings(
    {
      out <- sienaAME(
        ans = ans, data = mydata,
        effectName1 = "transTrip", diff1 = 1,
        effectName2 = "recip", contrast2 = c(0, 1),
        second = TRUE, type = "tieProb", depvar = "mynet",
        level = "egoChoice",
        condition = c("recip", "density", "transTrip"),
        uncertainty = FALSE
      )
      expect_true(is.data.frame(out))
      expect_true("secondDiff" %in% names(out))
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
})

test_that("sienaAME static interaction, no uncertainty (base R)", {
  with_mocked_bindings(
    {
      out <- sienaAME(
        ans = ans2, data = mydata2,
        effectName1 = "transTrip", diff1 = 1,
        interaction1 = TRUE, int_effectNames1 = "transRecTrip",
        mod_effectNames1 = "recip",
        type = "tieProb", depvar = "mynet2",
        level = "period", condition = "recip",
        uncertainty = FALSE, verbose = FALSE
      )
      expect_true(is.data.frame(out))
      expect_true("firstDiff" %in% names(out))
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
})

test_that("sienaAME static interaction with uncertainty (base R)", {
  with_mocked_bindings(
    {
      out <- sienaAME(
        ans = ans2, data = mydata2,
        effectName1 = "transTrip", diff1 = 1,
        interaction1 = TRUE, int_effectNames1 = "transRecTrip",
        mod_effectNames1 = "recip",
        type = "tieProb", depvar = "mynet2",
        level = "period", condition = "recip",
        nsim = 5, uncertainty = TRUE, verbose = FALSE
      )
      expect_true(is.data.frame(out))
      expect_true("Mean" %in% names(out))
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
})

test_that("sienaAME static secondDiff with interactions, riskRatio (base R)", {
  with_mocked_bindings(
    {
      out <- sienaAME(
        ans = ans2, data = mydata2,
        effectName1 = "recip", contrast1 = c(0, 1),
        interaction1 = TRUE, int_effectNames1 = "transRecTrip",
        mod_effectNames1 = "transTrip",
        second = TRUE,
        effectName2 = "transTrip", diff2 = 1,
        interaction2 = TRUE, int_effectNames2 = "transRecTrip",
        mod_effectNames2 = "recip",
        type = "tieProb", depvar = "mynet2",
        level = "period", nsim = 5, uncertainty = TRUE,
        verbose = FALSE, mainEffect = "riskRatio"
      )
      expect_true(is.data.frame(out))
      expect_true("secondRiskRatio" %in% names(out))
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
})

test_that("sienaAME static with custom interaction (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  out <- sienaAME(
    ans = ans3, data = mydata3, effects = mymodel3,
    effectName1 = "recip", contrast1 = c(0, 1),
    interaction1 = TRUE, int_effectNames1 = "unspInt",
    mod_effectNames1 = "inPop",
    type = "tieProb", depvar = "mynet3",
    level = "period", condition = "inPop",
    nsim = 5, uncertainty = FALSE, verbose = FALSE
  )
  expect_true("data.table" %in% class(out) || is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

# ── Dynamic tests (sienaAME dynamic = TRUE) ───────────────────────────────────

test_that("sienaAME dynamic: firstDiff structure, uncertainty=FALSE (base R)", {
  with_mocked_bindings(
    {
      out <- sienaAME(
        ans = ans, data = mydata,
        effectName1 = "transTrip", diff1 = 1,
        effects = mymodel, algorithm = mycontrols,
        dynamic = TRUE, n3 = 60, nsim = 5,
        type = "tieProb", condition = "density",
        uncertainty = FALSE
      )
      expect_true(is.data.frame(out))
      expect_false("data.table" %in% class(out))
      expect_true("firstDiff" %in% names(out))
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
})

test_that("sienaAME dynamic: MCSE + CI structure (base R)", {
  with_mocked_bindings(
    {
      out <- sienaAME(
        ans = ans, data = mydata,
        effectName1 = "transTrip", diff1 = 1,
        effects = mymodel, algorithm = mycontrols,
        dynamic = TRUE, n3 = 60, nsim = 20,
        type = "tieProb", condition = "density",
        uncertainty = TRUE,
        uncertainty_mcse = TRUE, uncertainty_mcse_batches = 4,
        uncertainty_sd = FALSE, uncertainty_ci = TRUE,
        verbose = FALSE
      )
      expect_true(is.data.frame(out))
      expect_true(all(c("Mean", "cases", "mcse_Mean", "q_025", "q_975", "Median") %in% names(out)))
      expect_false("SE" %in% names(out))
      expect_false("mcse_SE" %in% names(out))
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
})

test_that("sienaAME dynamic: data.table output", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  out <- sienaAME(
    ans = ans, data = mydata,
    effectName1 = "transTrip", diff1 = 1,
    effects = mymodel, algorithm = mycontrols,
    dynamic = TRUE, n3 = 60, nsim = 5,
    type = "tieProb", condition = "density"
  )
  expect_true("data.table" %in% class(out) || is.data.frame(out))
  expect_true("firstDiff" %in% names(out))
})

test_that("sienaAME dynamic: interaction (base R)", {
  with_mocked_bindings(
    {
      out <- sienaAME(
        ans = ans2, data = mydata2,
        effectName1 = "transTrip", diff1 = 1,
        interaction1 = TRUE, int_effectNames1 = "transRecTrip",
        mod_effectNames1 = "recip",
        effects = mymodel2, algorithm = mycontrols2,
        dynamic = TRUE, n3 = 60, nsim = 5,
        type = "tieProb", level = "period",
        condition = c("density", "recip"),
        uncertainty = FALSE
      )
      expect_true(is.data.frame(out))
      expect_false("data.table" %in% class(out))
      expect_true("firstDiff" %in% names(out))
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
})



