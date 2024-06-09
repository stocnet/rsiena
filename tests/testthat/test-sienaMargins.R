skip_on_cran()

# library(RSiena) after installation or better devtools::load_all() for interactive use

# Run simple reciprocity GWESPFF model (but return chains)
# Should probably be extracted to a helper file.
# Alternatively the results could be stored and post calculations compared

mynet <- sienaDependent(array(c(s502, s503), dim = c(50, 50, 2)))
mydata <- sienaDataCreate(mynet)
mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, gwespFF, name = "mynet")
mycontrols <- sienaAlgorithmCreate(projname = NULL, seed = 42)
ans <- siena07(
  mycontrols,
  batch = TRUE, 
  silent = TRUE,
  data = mydata,
  effects = mymodel,
  returnChains = TRUE
)

# Test expectedChangeProbabilities
test_that("expectedChangeProbabilities calculates  the same expected probabilities as expectedRelativeImportance()", {
  # Read contributions
  contributions <- RSiena:::getChangeContributions(
    mycontrols,
    mydata,
    mymodel
  )

  # Run expectedRelativeImportance to compare results
  # Should be made independent of expectedRelativeImportance later
  RI <- RSiena:::expectedRelativeImportance(
    conts = contributions,
    effects = ans$effects,
    theta = ans$theta,
    thedata = mydata
  )

  expect_equal(
    RSiena:::expectedChangeProbabilities(
      conts = contributions,
      effects = ans$effects,
      theta = ans$theta,
      thedata = mydata
    ),
    RI$toggleProbabilities,
    ignore_attr = TRUE
  )
}
)

# Test expectedChangeDynamics
test_that("Test if expectedChangeDynamics creates a list of matrices of correct dimensions", {
  expectedProbChain <- RSiena:::expectedChangeDynamics(
    data = mydata,
    theta = c(ans$rate, ans$theta),
    algorithm = ans$x,
    effects = mymodel,
    depvar = "mynet"
  )

  # Not a real unit test yet. How to define true result?
  expect_type(
    expectedProbChain,
    "list"
  )

  expect_in(
    lapply(expectedProbChain[[1]], function(x) class(x)),
    list(c("matrix", "array"))
  )

  expect_in(
    lapply(expectedProbChain[[1]], function(x) dim(x)),
    list(c(dim(mynet)[1], 1))
  )
}
)