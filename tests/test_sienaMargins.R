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
effectNames <- mymodel$shortName[mymodel$include & (mymodel$type != "rate")]

# Basic one-diff AME, point estimate only
cat("\nTesting sienaAME (static) ...\n")
try({
  ame_static <- sienaAME(
    ans = ans,
    data = mydata,
    effectName1 = "transTrip",
    diff1 = 1,
    effectNames = effectNames,
    effects = mymodel,
    useTieProb = TRUE,
    depvar = "mynet",
    uncertainty = FALSE,  # no sims for quick check
    verbose = TRUE
  )
  print(head(ame_static))
  stopifnot(is.data.frame(ame_static) || is.data.table(ame_static))
  stopifnot("firstDiff" %in% names(ame_static))
  stopifnot(sum(abs(ame_static$firstDiff)) > 0) # Should be nontrivial
  cat("sienaAME: SUCCESS\n")
})

# Optional: check with uncertainty/parallel (should be similar, but gives quantiles etc.)
cat("\nTesting sienaAME (static, uncertainty) ...\n")
try({
  ame_static_uncert <- sienaAME(
    ans = ans,
    data = mydata,
    effectName1 = "recip",
    contrast1 = c(0, 1),
    effectNames = effectNames,
    effects = mymodel,
    useTieProb = TRUE,
    depvar = "mynet",
    level = "period",
    condition = c("transTrip"),
    nsim = 100,  # keep small for speed
    uncertainty = TRUE,
    verbose = TRUE
  )
  print(ame_static_uncert)
   ame_static_uncert2 <- sienaAME(
    ans = ans,
    data = mydata,
    effectName1 = "transTrip",
    diff1 = 1,
    effectNames = effectNames,
    effects = mymodel,
    useTieProb = TRUE,
    depvar = "mynet",
    level = "period",
    condition = c("recip"),
    nsim = 100,  # keep small for speed
    uncertainty = TRUE,
    verbose = TRUE
  )
  # This should return something with quantiles/SD/mean
  stopifnot(is.data.frame(ame_static_uncert) || is.data.table(ame_static_uncert))
  cat("sienaAME (with uncertainty): SUCCESS\n")
})

cat("\nAll sienaAME static tests completed successfully.\n")


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
effectNames2 <- mymodel2$shortName[mymodel2$include & (mymodel2$type != "rate")]

ame_static_interaction <- sienaAME(
    ans = ans2,
    data = mydata2,
    effectName1 = "transTrip",
    diff1 = 1,
    interaction1 = TRUE,
    int_effectNames1 = "transRecTrip",
    mod_effectNames1 = "recip",
    effectNames = effectNames2,
    effects = mymodel2,
    useTieProb = FALSE,
    depvar = "mynet2",
    level = "period",
    condition = "density",
    nsim = 100,  # keep small for speed
    uncertainty = TRUE,
    verbose = TRUE
)
ame_static_interaction

ame_static_interaction2 <- sienaAME(
    ans = ans2,
    data = mydata2,
    effectName1 = "recip",
    contrast1 = c(0,1),
    interaction1 = TRUE,
    int_effectNames1 = "transRecTrip",
    mod_effectNames1 = "transTrip",
    effectNames = effectNames2,
    effects = mymodel2,
    useTieProb = TRUE,
    depvar = "mynet2",
    level = "period",
    condition = "transTrip",
    nsim = 100,  # keep small for speed
    uncertainty = TRUE,
    verbose = TRUE
)

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
    contrast2 = c(0,1), #be more explicit about lower and hihger?
    interaction2 = TRUE,
    int_effectNames2 = "transRecTrip",
    mod_effectNames2 = "transTrip",
    effectNames = effectNames2,
    effects = mymodel2,
    useTieProb = TRUE,
    depvar = "mynet2",
    level = "period",
    nsim = 100,  # keep small for speed
    uncertainty = TRUE,
    verbose = TRUE
)

ame_static_moderator

ame_static_moderator2 <- sienaAME(
    ans = ans2,
    data = mydata2,
    effectName1 = "recip",
    contrast1 = c(0,1), #be more explicit about lower and hihger?
    interaction1 = TRUE,
    int_effectNames1 = "transRecTrip",
    mod_effectNames1 = "transTrip",
    second = TRUE,
    effectName2 = "transTrip",
    diff2 = 1,
    interaction2 = TRUE,
    int_effectNames2 = "transRecTrip",
    mod_effectNames2 = "recip",
    effectNames = effectNames2,
    effects = mymodel2,
    useTieProb = TRUE,
    depvar = "mynet2",
    level = "period",
    condition = "density",
    nsim = 100,  # keep small for speed
    uncertainty = TRUE,
    verbose = TRUE
)

ame_static_moderator2

mynet3 <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata3 <- sienaDataCreate(mynet3)
mymodel3 <- getEffects(mydata3)
## outdegree recip model
mymodel3 <- includeEffects(mymodel3, transTrip, name = "mynet3")
mymodel3 <- includeEffects(mymodel3, outPop, name = "mynet3")
mymodel3 <- includeInteraction(mymodel3, recip, outPop, name = "mynet3")
effectNames3 <- mymodel3$shortName[mymodel3$include & (mymodel3$type != "rate")]

mycontrols3 <- sienaAlgorithmCreate(projname=NULL, n3 = 500, cond = FALSE)
ans3 <- siena07(
  mycontrols3,
  data = mydata3,
  effects = mymodel3,
  returnChangeContributions = TRUE,
  returnDataFrame = TRUE
)

ame_static_moderator <- sienaAME(
    ans = ans3,
    data = mydata3,
    effectName1 = "transTrip",
    diff1 = 1,
    second = TRUE,
    effectName2 = "outPop",
    diff2 = 1,
    effectNames = effectNames3,
    effects = mymodel3,
    algorithm = mycontrols3,
    useTieProb = TRUE,
    depvar = "mynet3",
    level = "period",
    n3 = 500,
    nsim = 100,  # keep small for speed
    uncertainty = TRUE,
    verbose = TRUE
)

ame_dynamic_moderator

test1 <- sienaAME(
    ans = ans,
    data = mydata,
    effectName1 = "recip",
    contrast1 = c(0,1), #be more explicit about lower and hihger?
    second = TRUE,
    effectName2 = "transTrip",
    diff2 = 1,
    effectNames = effectNames,
    effects = mymodel,
    useTieProb = TRUE,
    depvar = "mynet1",
    level = "period",
    nsim = 100,  # keep small for speed
    uncertainty = FALSE,
    verbose = TRUE
)

test2 <- sienaAME(
    ans = ans,
    data = mydata,
    effectName1 = "transTrip",
    diff1 = 1,
    second = TRUE,
    effectName2 = "recip",
    contrast2 = c(0,1), #be more explicit about lower and hihger?
    effectNames = effectNames,
    effects = mymodel,
    useTieProb = TRUE,
    depvar = "mynet1",
    level = "period",
    nsim = 100,  # keep small for speed
    uncertainty = TRUE,
    verbose = TRUE
)
