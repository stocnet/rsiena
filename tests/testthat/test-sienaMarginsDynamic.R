# testthat::skip_on_cran()

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

test_that("sienaAMEDynamic (base R fallback, uncertainty = FALSE)", {
  with_mocked_bindings(
    {
      margins_dynamic <- sienaAMEDynamic(
        ans = ans,
        data = mydata,
        effectName1 = "transTrip",
        diff1 = 1,
        effects = mymodel,
        algorithm = mycontrols,
        n3 = 60,
        nsim = 5,
        useTieProb = TRUE,
        condition = "density",
        uncertainty = FALSE
      )
      expect_true(is.data.frame(margins_dynamic))
      expect_false("data.table" %in% class(margins_dynamic))
      expect_true("firstDiff" %in% names(margins_dynamic))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )
})

test_that("sienaAMEDynamic (data.table, uncertainty = FALSE)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  margins_dynamic <- sienaAMEDynamic(
    ans = ans,
    data = mydata,
    effectName1 = "transTrip",
    diff1 = 1,
    effects = mymodel,
    algorithm = mycontrols,
    n3 = 60,
    nsim = 5,
    useTieProb = TRUE,
    condition = "density",
    uncertainty = FALSE
  )
  expect_true("data.table" %in% class(margins_dynamic) || is.data.frame(margins_dynamic))
  expect_true("firstDiff" %in% names(margins_dynamic))
})

test_that("sienaAMEDynamic (base R fallback)", {
  with_mocked_bindings(
    {
      margins_dynamic <- sienaAMEDynamic(
        ans = ans,
        data = mydata,
        effectName1 = "transTrip",
        diff1 = 1,
        effects = mymodel,
        algorithm = mycontrols,
        n3 = 60, # keep small for test speed
        nsim = 5, # keep small for test speed
        useTieProb = TRUE,
        condition = "density"
      )
      expect_true(is.data.frame(margins_dynamic))
      expect_false("data.table" %in% class(margins_dynamic))
      expect_true("firstDiff" %in% names(margins_dynamic))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )
})

test_that("sienaAMEDynamic (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  margins_dynamic <- sienaAMEDynamic(
    ans = ans,
    data = mydata,
    effectName1 = "transTrip",
    diff1 = 1,
    effects = mymodel,
    algorithm = mycontrols,
    n3 = 60, # keep small for test speed
    nsim = 5, # keep small for test speed
    useTieProb = TRUE,
    condition = "density"
  )
  expect_true("data.table" %in% class(margins_dynamic) || is.data.frame(margins_dynamic))
  expect_true("firstDiff" %in% names(margins_dynamic))
})

  mynet2 <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
  mydata2 <- sienaDataCreate(mynet2)
  mymodel2 <- getEffects(mydata2)
  mymodel2 <- includeEffects(mymodel2, transTrip, transRecTrip, name = "mynet2")
  mycontrols2 <- sienaAlgorithmCreate(projname=NULL, n3 = 50, cond = FALSE)
  ans2 <- siena07(
    mycontrols2,
    data = mydata2,
    effects = mymodel2,
    returnChangeContributions = TRUE,
    returnDataFrame = TRUE
  )

test_that("sienaAMEDynamic (base R fallback, interaction)", {
  with_mocked_bindings(
    {
      margins_dynamic_interaction <- sienaAMEDynamic(
        ans = ans2,
        data = mydata2,
        effectName1 = "transTrip",
        diff1 = 1,
        interaction1 = TRUE,
        int_effectNames1 = "transRecTrip",
        mod_effectNames1 = "recip",
        effects = mymodel2,
        algorithm = mycontrols2,
        n3 = 60,
        nsim = 5,
        useTieProb = TRUE,
        level = "period",
        condition = "density"
      )
      expect_true(is.data.frame(margins_dynamic_interaction))
      expect_false("data.table" %in% class(margins_dynamic_interaction))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )
})

test_that("sienaAMEDynamic (data.table, interaction)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  margins_dynamic_interaction <- sienaAMEDynamic(
    ans = ans2,
    data = mydata2,
    effectName1 = "transTrip",
    diff1 = 1,
    interaction1 = TRUE,
    int_effectNames1 = "transRecTrip",
    mod_effectNames1 = "recip",
    effects = mymodel2,
    algorithm = mycontrols2,
    n3 = 60,
    nsim = 5,
    useTieProb = TRUE,
    level = "period",
    condition = "density"
  )
  expect_true("data.table" %in% class(margins_dynamic_interaction) || is.data.frame(margins_dynamic_interaction))
})

mynet3 <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata3 <- sienaDataCreate(mynet3)
mymodel3 <- getEffects(mydata3)
mymodel3 <- includeEffects(mymodel3, transTrip, name = "mynet3")
mymodel3 <- includeEffects(mymodel3, outPop, name = "mynet3")
mymodel3 <- includeInteraction(mymodel3, recip, outPop, name = "mynet3")
mycontrols3 <- sienaAlgorithmCreate(projname=NULL, n3 = 50, cond = FALSE)
ans3 <- siena07(
  mycontrols3,
  data = mydata3,
  effects = mymodel3,
  returnChangeContributions = TRUE,
  returnDataFrame = TRUE
)

test_that("sienaAMEDynamic (custom interaction)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
    ame_dynamic_custom <- sienaAMEDynamic(
      ans = ans3,
      data = mydata3,
      effectName1 = "outPop",
      diff1 = 1,
      interaction1 = TRUE,
      int_effectNames1 = "unspInt",
      mod_effectNames1 = "recip",
      effects = mymodel3,
      algorithm = mycontrols3,
      n3 = 60,
      nsim = 5,
      useTieProb = TRUE,
      level = "period",
      condition = "density"
    )
  expect_true("data.table" %in% class(ame_dynamic_custom) || is.data.frame(ame_dynamic_custom))
  expect_true("firstDiff" %in% names(ame_dynamic_custom))
})

test_that("sienaAMEDynamic (base R fallback, custom interaction)", {
  with_mocked_bindings(
    {
      ame_dynamic_custom <- sienaAMEDynamic(
        ans = ans3,
        data = mydata3,
        effectName1 = "outPop",
        diff1 = 1,
        interaction1 = TRUE,
        int_effectNames1 = "unspInt",
        mod_effectNames1 = "recip",
        effects = mymodel3,
        algorithm = mycontrols3,
        n3 = 60,
        nsim = 5,
        useTieProb = TRUE,
        level = "period",
        condition = "density"
      )
      expect_true(is.data.frame(ame_dynamic_custom))
      expect_false("data.table" %in% class(ame_dynamic_custom))
      expect_true("firstDiff" %in% names(ame_dynamic_custom))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )
})
# test_sienaAMEDynamic_all_cases <- function(ans, mydata, mymodel, mycontrols) {
#   # Setup
#   mynet <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
#   mydata <- sienaDataCreate(mynet)
#   mymodel <- getEffects(mydata)
#   mymodel <- includeEffects(mymodel, transTrip, name = "mynet")
#   mycontrols <- sienaAlgorithmCreate(projname=NULL, n3 = 500, cond = FALSE)
#   ans <- siena07(
#     mycontrols,
#     data = mydata,
#     effects = mymodel,
#     returnChangeContributions = TRUE,
#     returnDataFrame = TRUE
#   )
#   results <- list()
#   cat("=== Standard First Difference ===\n")
#   results[["first_diff"]] <- sienaAMEDynamic(
#     ans = ans,
#     data = mydata,
#     effectName1 = "transTrip",
#     diff1 = 1,
#     effects = mymodel,
#     algorithm = mycontrols,
#     n3 = 50, # keep small for test speed
#     nsim = 10, # Keep small for test speed
#     useTieProb = TRUE,
#     condition = "density"
#   )
#   print(head(results[["first_diff"]]))

#   cat("=== First Difference with Interaction ===\n")
#   mynet2 <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
#   mydata2 <- sienaDataCreate(mynet2)
#   mymodel2 <- getEffects(mydata2)
#   mymodel2 <- includeEffects(mymodel2, transTrip, transRecTrip, name = "mynet2")
#   mycontrols2 <- sienaAlgorithmCreate(projname=NULL, n3 = 500, cond = FALSE)
#   ans2 <- siena07(
#     mycontrols2,
#     data = mydata2,
#     effects = mymodel2,
#     returnChangeContributions = TRUE,
#     returnDataFrame = TRUE
#   )
#   results[["first_diff_interaction"]] <- sienaAMEDynamic(
#     ans = ans2,
#     data = mydata2,
#     effectName1 = "transTrip",
#     diff1 = 1,
#     interaction1 = TRUE,
#     int_effectNames1 = "transRecTrip",
#     mod_effectNames1 = "recip",
#     effects = mymodel2,
#     algorithm = mycontrols2,
#     second = FALSE,
#     n3 = 50, # keep small for test speed
#     nsim = 10, # Keep small for test speed
#     useTieProb = TRUE,
#     level = "period",
#     condition = "density"
#   )
#   print(head(results[["first_diff_interaction"]]))

#     cat("=== Second Difference, mixed diff/contrast ===\n")
#   results[["second_diff_mixed"]] <- sienaAMEDynamic(
#     ans = ans2,
#     data = mydata2,
#     effectName1 = "transTrip", diff1 = 1,
#     interaction1 = TRUE,
#     int_effectNames1 = "transRecTrip",
#     mod_effectNames1 = "recip",
#     effectName2 = "recip", contrast2 = c(0,1),
#     interaction2 = TRUE,
#     int_effectNames2 = "transRecTrip",
#     mod_effectNames2 = "transTrip",
#     effects = mymodel2,
#     algorithm = mycontrols2,
#     second = TRUE,
#     n3 = 50, # keep small for test speed
#     nsim = 10, # Keep small for test speed
#     useTieProb = FALSE,
#     level = "period",
#     condition = "density"
#   )
#   print(head(results[["second_diff_mixed"]]))

#     cat("=== Test Summary ===\n")
#   for (nm in names(results)) {
#     cat(nm, ": Rows=", nrow(results[[nm]]), if (is.data.frame(results[[nm]])) paste("Cols=", ncol(results[[nm]])), "\n")
#   }
#   invisible(results)

#   ame_dynamic_interaction <- sienaAMEDynamic(
#     ans = ans2,
#     data = mydata2,
#     effectName1 = "transTrip",
#     diff1 = 1,
#     interaction1 = TRUE,
#     int_effectNames1 = "transRecTrip",
#     mod_effectNames1 = "recip",
#     effects = mymodel2,
#     algorithm = mycontrols2, # checks need to be added, that algorithm is algorithm object etc.
#     useTieProb = FALSE,
#     depvar = "mynet2",
#     level = "period",
#     condition = "density",
#     n3 = 50, # keep small for test speed
#     nsim = 10, # Keep small for test speed
#     uncertainty = TRUE,
#     verbose = TRUE
#   )
#   ame_dynamic_interaction

#   ame_dynamic_moderator <- sienaAMEDynamic(
#       ans = ans2,
#       data = mydata2,
#       effectName1 = "transTrip",
#       diff1 = 1,
#       interaction1 = TRUE,
#       int_effectNames1 = "transRecTrip",
#       mod_effectNames1 = "recip",
#       second = TRUE,
#       effectName2 = "recip",
#       contrast2 = c(0,1),
#       interaction2 = TRUE,
#       int_effectNames2 = "transRecTrip",
#       mod_effectNames2 = "transTrip",
#       effects = mymodel2,
#       algorithm = mycontrols2,
#       useTieProb = TRUE,
#       depvar = "mynet2", # add check that this is ANY of the depvars/netnames!
#       level = "period",
#       n3 = 50, # keep small for test speed
#       nsim = 10, # Keep small for test speed
#       batch_size = 5,
#       uncertainty = TRUE,
#       verbose = TRUE
#   )

#   ame_dynamic_moderator

#   ame_dynamic_moderator2 <- sienaAMEDynamic(
#       ans = ans2,
#       data = mydata2,
#       effectName1 = "recip",
#       contrast1 = c(0,1), #be more explicit about lower and hihger?
#       interaction1 = TRUE,
#       int_effectNames1 = "transRecTrip",
#       mod_effectNames1 = "transTrip",
#       second = TRUE,
#       effectName2 = "transTrip",
#       diff2 = 1,
#       interaction2 = TRUE,
#       int_effectNames2 = "transRecTrip",
#       mod_effectNames2 = "recip",
#       effects = mymodel2,
#       algorithm = mycontrols2,
#       useTieProb = TRUE,
#       depvar = "mynet2",
#       level = "period",
#       condition = "density",
#       n3 = 50,
#       nsim = 10,  # keep small for speed
#       batch_size = 5,
#       uncertainty = TRUE,
#       verbose = TRUE
#   )

#   ame_dynamic_moderator2

#   mynet3 <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
#   mydata3 <- sienaDataCreate(mynet3)
#   mymodel3 <- getEffects(mydata3)
#   ## outdegree recip model
#   mymodel3 <- includeEffects(mymodel3, transTrip, name = "mynet3")
#   mymodel3 <- includeEffects(mymodel3, outPop, name = "mynet3")
#   mymodel3 <- includeInteraction(mymodel3, recip, outPop, name = "mynet3")

#   mycontrols3 <- sienaAlgorithmCreate(projname=NULL, n3 = 500, cond = FALSE)
#   ans3 <- siena07(
#     mycontrols3,
#     data = mydata3,
#     effects = mymodel3,
#     returnChangeContributions = TRUE,
#     returnDataFrame = TRUE,
#     useCluster = TRUE, nbrNodes = 10, clusterType = "FORK"
#   )

#   ame_dynamic_moderator3 <- sienaAMEDynamic(
#       ans = ans3,
#       data = mydata3,
#       effectName1 = "transTrip",
#       diff1 = 1,
#       second = TRUE,
#       effectName2 = "outPop",
#       diff2 = 1,
#       effects = mymodel3,
#       algorithm = mycontrols3,
#       useTieProb = TRUE,
#       depvar = "mynet3",
#       level = "period",
#       n3 = 50, # keep small for test speed
#       nsim = 10, # Keep small for test speed
#       useCluster = TRUE, nbrNodes = 10, clusterType = "FORK",
#       batch_size = 10,
#       verbose = TRUE
#   )

#   ame_dynamic_moderator3

# }

# test_sienaAMEDynamic_all_cases()