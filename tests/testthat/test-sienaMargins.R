# testthat::skip_on_cran()

# library(testthat)
# library(RSiena)

# Minimal RSiena setup for reproducible test
mynet <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata <- sienaDataCreate(mynet)
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, transTrip, name = "mynet")
mycontrols <- sienaAlgorithmCreate(projname=NULL, n3 = 50, cond = FALSE)
ans <- siena07(
  mycontrols,
  data = mydata,
  effects = mymodel,
  returnChangeContributions = TRUE,
  returnDataFrame = TRUE
)

test_that("sienaAME static (base R fallback)", {
  with_mocked_bindings(
    {
      ame_static <- sienaAME(
        ans = ans,
        data = mydata,
        effectName1 = "transTrip",
        diff1 = 1,
        useTieProb = TRUE,
        depvar = "mynet",
        level = "egoChoice",
        condition = c("recip","density","transTrip"),
        uncertainty = FALSE
      )
      expect_true(is.data.frame(ame_static))
      expect_false("data.table" %in% class(ame_static))
      expect_true("firstDiff" %in% names(ame_static))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )
})


test_that("sienaAME static (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  ame_static <- sienaAME(
    ans = ans,
    data = mydata,
    effectName1 = "transTrip",
    diff1 = 1,
    useTieProb = TRUE,
    depvar = "mynet",
    level = "egoChoice",
    condition = c("recip","density","transTrip"),
    uncertainty = FALSE
  )
  expect_true("data.table" %in% class(ame_static) || is.data.frame(ame_static))
  expect_true("firstDiff" %in% names(ame_static))
})

test_that("sienaAME static (base R fallback, uncertainty)", {
  with_mocked_bindings(
    {
      ame_static_uncert <- sienaAME(
        ans = ans,
        data = mydata,
        effectName1 = "recip",
        contrast1 = c(0, 1),
        useTieProb = FALSE,
        depvar = "mynet",
        level = "period",
        condition = c("density"),
        nsim = 10,  # keep small for speed
        uncertainty = TRUE
      )
      expect_true(is.data.frame(ame_static_uncert))
      expect_false("data.table" %in% class(ame_static_uncert))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )
})

test_that("sienaAME static (data.table, uncertainty)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  ame_static_uncert <- sienaAME(
    ans = ans,
    data = mydata,
    effectName1 = "recip",
    contrast1 = c(0, 1),
    useTieProb = FALSE,
    depvar = "mynet",
    level = "period",
    condition = c("density"),
    nsim = 10,  # keep small for speed
    uncertainty = TRUE
  )
  expect_true("data.table" %in% class(ame_static_uncert) || is.data.frame(ame_static_uncert))
})

test_that("sienaAME static (base R fallback, second difference)", {
  with_mocked_bindings(
    {
      ame_second <- sienaAME(
        ans = ans,
        data = mydata,
        effectName1 = "transTrip",
        effectName2 = "recip",
        diff1 = 1,
        contrast2 = c(0, 1),
        second = TRUE,
        useTieProb = TRUE,
        depvar = "mynet",
        level = "egoChoice",
        condition = c("recip","density","transTrip"),
        uncertainty = FALSE
      )
      expect_true(is.data.frame(ame_second))
      expect_false("data.table" %in% class(ame_second))
      expect_true("secondDiff" %in% names(ame_second))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )
})

test_that("sienaAME static (data.table, second difference)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  ame_second <- sienaAME(
    ans = ans,
    data = mydata,
    effectName1 = "transTrip",
    effectName2 = "recip",
    diff1 = 1,
    contrast2 = c(0, 1),
    second = TRUE,
    useTieProb = TRUE,
    depvar = "mynet",
    level = "egoChoice",
    condition = c("recip","density","transTrip"),
    uncertainty = FALSE
  )
  expect_true("data.table" %in% class(ame_second) || is.data.frame(ame_second))
  expect_true("secondDiff" %in% names(ame_second))
})


mynet2 <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata2 <- sienaDataCreate(mynet2)
mymodel2 <- getEffects(mydata2)
mymodel2 <- includeEffects(mymodel2, transTrip, transRecTrip, name = "mynet2")
mycontrols2 <- sienaAlgorithmCreate(projname=NULL, n3 = 60, cond = FALSE)
ans2 <- siena07(
  mycontrols2,
  data = mydata2,
  effects = mymodel2,
  returnChangeContributions = TRUE,
  returnDataFrame = TRUE
)

test_that("sienaAME static interaction (base R fallback)", {
  with_mocked_bindings(
    {
      ame_static_interaction <- sienaAME(
        ans = ans2,
        data = mydata2,
        effectName1 = "transTrip",
        diff1 = 1,
        interaction1 = TRUE,
        int_effectNames1 = "transRecTrip",
        mod_effectNames1 = "recip",
        useTieProb = TRUE,
        depvar = "mynet2",
        level = "period",
        condition = "recip",
        nsim = 10,
        uncertainty = TRUE,
        verbose = FALSE,
        mainEffect = "riskDifference"
      )
      expect_true(is.data.frame(ame_static_interaction))
      expect_false("data.table" %in% class(ame_static_interaction))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )
})

test_that("sienaAME static interaction (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  ame_static_interaction <- sienaAME(
    ans = ans2,
    data = mydata2,
    effectName1 = "transTrip",
    diff1 = 1,
    interaction1 = TRUE,
    int_effectNames1 = "transRecTrip",
    mod_effectNames1 = "recip",
    useTieProb = TRUE,
    depvar = "mynet2",
    level = "period",
    condition = "density",
    nsim = 10,
    uncertainty = TRUE,
    verbose = FALSE,
    mainEffect = "riskRatio"
  )
  expect_true("data.table" %in% class(ame_static_interaction) || is.data.frame(ame_static_interaction))
})

test_that("sienaAME static moderator (base R fallback)", {
  with_mocked_bindings(
    {
      ame_static_moderator <- sienaAME(
        ans = ans2,
        data = mydata2,
        effectName1 = "transTrip",
        diff1 = 1,
        interaction1 = TRUE,
        int_effectNames1 = "transRecTrip",
        mod_effectNames1 = "recip",
        second = TRUE,
        effectName2 = "recip",
        contrast2 = c(0,1),
        interaction2 = TRUE,
        int_effectNames2 = "transRecTrip",
        mod_effectNames2 = "transTrip",
        useTieProb = TRUE,
        depvar = "mynet2",
        level = "period",
        uncertainty = FALSE,
        verbose = FALSE,
        mainEffect = "riskDifference"
      )
      expect_true(is.data.frame(ame_static_moderator))
      expect_false("data.table" %in% class(ame_static_moderator))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )
})

test_that("sienaAME static moderator (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  ame_static_moderator <- sienaAME(
    ans = ans2,
    data = mydata2,
    effectName1 = "transTrip",
    diff1 = 1,
    interaction1 = TRUE,
    int_effectNames1 = "transRecTrip",
    mod_effectNames1 = "recip",
    second = TRUE,
    effectName2 = "recip",
    contrast2 = c(0,1),
    interaction2 = TRUE,
    int_effectNames2 = "transRecTrip",
    mod_effectNames2 = "transTrip",
    useTieProb = TRUE,
    depvar = "mynet2",
    level = "period",
    nsim = 10,
    uncertainty = FALSE,
    verbose = FALSE,
    mainEffect = "riskRatio"
  )
  expect_true("data.table" %in% class(ame_static_moderator) || is.data.frame(ame_static_moderator))
})

test_that("sienaAME static moderator2 (base R fallback, uncertainty)", {
  with_mocked_bindings(
    {
      ame_static_moderator2 <- sienaAME(
        ans = ans2,
        data = mydata2,
        effectName1 = "recip",
        contrast1 = c(0,1),
        interaction1 = TRUE,
        int_effectNames1 = "transRecTrip",
        mod_effectNames1 = "transTrip",
        second = TRUE,
        effectName2 = "transTrip",
        diff2 = 1,
        interaction2 = TRUE,
        int_effectNames2 = "transRecTrip",
        mod_effectNames2 = "recip",
        useTieProb = TRUE,
        depvar = "mynet2",
        level = "period",
        nsim = 5,
        uncertainty = TRUE,
        verbose = FALSE,
        mainEffect = "riskRatio"
      )
      expect_true(is.data.frame(ame_static_moderator2))
      expect_false("data.table" %in% class(ame_static_moderator2))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )
})

test_that("sienaAME static moderator2 (data.table, uncertainty)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  ame_static_moderator2 <- sienaAME(
    ans = ans2,
    data = mydata2,
    effectName1 = "recip",
    contrast1 = c(0,1),
    interaction1 = TRUE,
    int_effectNames1 = "transRecTrip",
    mod_effectNames1 = "transTrip",
    second = TRUE,
    effectName2 = "transTrip",
    diff2 = 1,
    interaction2 = TRUE,
    int_effectNames2 = "transRecTrip",
    mod_effectNames2 = "recip",
    useTieProb = TRUE,
    depvar = "mynet2",
    level = "period",
    nsim = 5,
    uncertainty = TRUE,
    verbose = FALSE,
    mainEffect = "riskRatio"
  )
  expect_true("data.table" %in% class(ame_static_moderator2) || is.data.frame(ame_static_moderator2))
})

test_that("sienaAME static interaction (base R fallback, uncertainty)", {
  with_mocked_bindings(
    {
      ame_static_interaction <- sienaAME(
        ans = ans2,
        data = mydata2,
        effectName1 = "transTrip",
        diff1 = 1,
        interaction1 = TRUE,
        int_effectNames1 = "transRecTrip",
        mod_effectNames1 = "recip",
        useTieProb = TRUE,
        depvar = "mynet2",
        level = "period",
        condition = "recip",
        nsim = 5,
        uncertainty = TRUE,
        verbose = FALSE
      )
      expect_true(is.data.frame(ame_static_interaction))
      expect_false("data.table" %in% class(ame_static_interaction))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )
})

test_that("sienaAME static interaction (data.table, uncertainty)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  ame_static_interaction <- sienaAME(
    ans = ans2,
    data = mydata2,
    effectName1 = "transTrip",
    diff1 = 1,
    interaction1 = TRUE,
    int_effectNames1 = "transRecTrip",
    mod_effectNames1 = "recip",
    useTieProb = TRUE,
    depvar = "mynet2",
    level = "period",
    condition = "recip",
    nsim = 5,
    uncertainty = TRUE,
    verbose = FALSE
  )
  expect_true("data.table" %in% class(ame_static_interaction) || is.data.frame(ame_static_interaction))
})


mynet3 <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata3 <- sienaDataCreate(mynet3)
mymodel3 <- getEffects(mydata3)
## inPop recip model
mymodel3 <- includeEffects(mymodel3, inPop, name = "mynet3")
mymodel3 <- includeInteraction(mymodel3, recip, inPop, name = "mynet3")

mycontrols3 <- sienaAlgorithmCreate(projname=NULL, n3 = 60, cond = FALSE)
ans3 <- siena07(
  mycontrols3,
  data = mydata3,
  effects = mymodel3,
  returnChangeContributions = TRUE,
  returnDataFrame = TRUE
)

test_that("predict.sienaFit with custom interactions (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)

  ame_static_interaction <- sienaAME(
    ans = ans3,
    data = mydata3,
    effects = mymodel3,
    effectName1 = "recip",
    contrast1 = c(0,1),
    interaction1 = TRUE,
    int_effectNames1 = "unspInt",
    mod_effectNames1 = "inPop",
    useTieProb = TRUE,
    depvar = "mynet3",
    level = "period",
    condition = "inPop",
    nsim = 5,
    uncertainty = TRUE,
    verbose = TRUE
  )
  expect_true("data.table" %in% class(ame_static_interaction) || 
    is.data.frame(ame_static_interaction))
})

test_that("sienaAME static moderator (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), 
    "data.table not available")
  library(data.table)
  ame_static_moderator <- sienaAME(
    ans = ans3,
    data = mydata3,
    effectName1 = "recip",
    contrast1 = c(0,1),
    interaction1 = TRUE,
    int_effectNames1 = "unspInt",
    mod_effectNames1 = "inPop",
    second = TRUE,
    effectName2 = "inPop",
    diff2 = 1,
    interaction2 = TRUE,
    int_effectNames2 = "unspInt",
    mod_effectNames2 = "recip",
    useTieProb = TRUE,
    depvar = "mynet3",
    level = "period",
    uncertainty = FALSE,
    verbose = FALSE,
    mainEffect = "riskDifference"
  )
  expect_true("data.table" %in% class(ame_static_moderator) || 
    is.data.frame(ame_static_moderator))
})


# ame_dynamic_moderator

# test1 <- sienaAME(
#     ans = ans,
#     data = mydata,
#     effectName1 = "recip",
#     contrast1 = c(0,1), #be more explicit about lower and hihger?
#     second = TRUE,
#     effectName2 = "transTrip",
#     diff2 = 1,
#     effects = mymodel,
#     useTieProb = TRUE,
#     depvar = "mynet1",
#     level = "period",
#     nsim = 100,  # keep small for speed
#     uncertainty = FALSE,
#     verbose = TRUE
# )

# test2 <- sienaAME(
#     ans = ans,
#     data = mydata,
#     effectName1 = "transTrip",
#     diff1 = 1,
#     second = TRUE,
#     effectName2 = "recip",
#     contrast2 = c(0,1), #be more explicit about lower and hihger?
#     effects = mymodel,
#     useTieProb = TRUE,
#     depvar = "mynet1",
#     level = "period",
#     nsim = 100,  # keep small for speed
#     uncertainty = TRUE,
#     verbose = TRUE
# )
