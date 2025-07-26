# skip_on_cran()

# # library(RSiena) after installation or better devtools::load_all() for interactive use

# # Run simple reciprocity GWESPFF model (but return chains)
# # Should probably be extracted to a helper file.
# # Alternatively the results could be stored and post calculations compared

# mynet <- sienaDependent(array(c(s502, s503), dim = c(50, 50, 2)))
# mydata <- sienaDataCreate(mynet)
# mymodel <- getEffects(mydata)
# mymodel <- setEffect(mymodel, transTrip, name = "mynet")
# mymodel <- includeEffects(mymodel, transRecTrip, name = "mynet")
# mycontrols <- sienaAlgorithmCreate(projname = NULL, seed = 42)
# ans <- siena07(
#   mycontrols,
#   batch = TRUE, 
#   silent = TRUE,
#   data = mydata,
#   effects = mymodel,
#   returnChains = TRUE
# )

# # Read contributions
# contributions <- RSiena:::getChangeContributions(
#   mycontrols,
#   mydata,
#   mymodel
# )

# # Run expectedRelativeImportance to compare results
# # Should be made independent of expectedRelativeImportance later
# RI <- RSiena:::expectedRelativeImportance(
#   conts = contributions,
#   effects = ans$effects,
#   theta = ans$theta,
#   thedata = mydata
# )

# probs <- RSiena:::expectedChangeProbabilities(
#       conts = contributions,
#       effects = ans$effects,
#       theta = ans$theta,
#       thedata = mydata
#     )

# # Test expectedChangeProbabilities
# test_that("expectedChangeProbabilities calculates the same expected probabilities as expectedRelativeImportance()", {
#   expect_equal(
#     probs,
#     RI$toggleProbabilities,
#     ignore_attr = TRUE
#   )
# }
# )

# # Test directMaringalEffect

# # Test sienaAME for observed networks for each actor
# test_that("Test if sienaAME creates a matrix for each actor and that the average of the actor AMEs is a double ", {
#       # Not a real unit test yet. How to define true result?
#   expect_in(
#     class(sienaAME(probs[1,,1], ans, ans$effects$shortName)),
#     c("matrix", "array")
#   )

#   actorAME <- t(apply(probs[,,1], 1, function(p) sienaAME(p, ans, "recip")))
#   averageAME <- mean(actorAME[,1])

#   expect_type(
#     averageAME,
#     "double")
# }
# )


# # Test expectedChangeDynamics

#   expectedProbChain <- RSiena:::expectedChangeDynamics(
#     data = mydata,
#     theta = c(ans$rate, ans$theta),
#     algorithm = ans$x,
#     effects = mymodel,
#     depvar = "mynet"
#   )

# test_that("Test if expectedChangeDynamics creates a list of matrices of correct dimensions", {
#   # Not a real unit test yet. How to define true result?
#   expect_type(
#     expectedProbChain,
#     "list"
#   )

#   expect_in(
#     lapply(expectedProbChain[[1]], function(x) class(x)),
#     list(c("matrix", "array"))
#   )

#   expect_in(
#     lapply(expectedProbChain[[1]], function(x) dim(x)),
#     list(c(dim(mynet)[1], 1))
#   )
# }
# )

# # Test sienaAME for chain
# test_that("Test if sienaAME also creates a matrix for each actor for chains and that the average of the actorAME is a double", {
#       # Not a real unit test yet. How to define true result?  
#   chainprobs <- expectedProbChain[[1]]
#   expect_in(
#     class(sienaAME(chainprobs[[1]], ans, ans$effects$shortName)),
#     c("matrix", "array")
#   )

#   actorAME <- sapply(chainprobs, function(p) sienaAME(p, ans, "recip"))
#   averageAME <- mean(actorAME[1,])

#   expect_type(
#     averageAME,
#     "double")
# }
# )
