# Helper: shared model fits — sourced once per test session.
# All variables defined here are available to every test_that() block.
#
# Quick mode (default): only the three base models (ans, ans2, ans3) are fitted.
# Full mode (RSENA_FULL_TESTS=1): additionally fits predict-specific models
# (interaction-without-main, co-evolution, creation/maintenance).
#
# n3 = 50: minimum for 50-actor networks.

# Guard: run models only when called from a test runner, not standalone load_all().
# Detects devtools::test(), testthat::test_dir(), testthat::test_check() via call stack.
# Override with RSENA_LOAD_MODELS=1 for any other case.
.in_test_run <- function() {
  any(vapply(sys.calls(), function(x) {
    fn <- tryCatch(deparse(x[[1]], nlines = 1), error = function(e) "")
    # Match both namespace-qualified (devtools::test) and bare forms (R CMD check).
    grepl("(^|::)(test|test_dir|test_check)$", fn)
  }, logical(1))) || identical(Sys.getenv("RSENA_LOAD_MODELS"), "1")
}

if (.in_test_run()) {
print("helper-models.R is being sourced")
# ── Base model: density + recip + transTrip ──────────────────────────────────
mynet    <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata   <- sienaDataCreate(mynet)
mymodel  <- getEffects(mydata)
mymodel  <- includeEffects(mymodel, transTrip, name = "mynet")
mycontrols <- sienaAlgorithmCreate(projname = NULL, n3 = 50, cond = FALSE, seed = 42)
ans <- siena07(
  mycontrols,
  data    = mydata,
  effects = mymodel,
  returnDeps               = TRUE,
  returnChangeContributions = TRUE,
  returnDataFrame          = TRUE,
  silent                   = TRUE
)

# ── Model 2: add transRecTrip (interaction/moderator tests in sienaMargins) ──
mynet2    <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata2   <- sienaDataCreate(mynet2)
mymodel2  <- getEffects(mydata2)
mymodel2  <- includeEffects(mymodel2, transTrip, transRecTrip, name = "mynet2")
mycontrols2 <- sienaAlgorithmCreate(projname = NULL, n3 = 50, cond = FALSE, seed = 7)
ans2 <- siena07(
  mycontrols2,
  data    = mydata2,
  effects = mymodel2,
  returnChangeContributions = TRUE,
  returnDataFrame           = TRUE,
  silent                    = TRUE
)

# ── Model 3: inPop * recip via unspInt (custom interaction in sienaMargins) ──
mynet3    <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata3   <- sienaDataCreate(mynet3)
mymodel3  <- getEffects(mydata3)
mymodel3  <- includeEffects(mymodel3, inPop, name = "mynet3")
mymodel3  <- includeInteraction(mymodel3, recip, inPop, name = "mynet3")
mycontrols3 <- sienaAlgorithmCreate(projname = NULL, n3 = 50, cond = FALSE, seed = 13)
ans3 <- siena07(
  mycontrols3,
  data    = mydata3,
  effects = mymodel3,
  returnChangeContributions = TRUE,
  returnDataFrame           = TRUE,
  silent                    = TRUE
)

# ── Predict-specific models — only in full test mode ─────────────────────────
# Initialised to NULL so tests can call skip_slow() without "object not found" errors.
ans_ego           <- NULL
mydata_ego        <- NULL
mymodel_ego       <- NULL
mycontrols_ego    <- NULL
ans_int           <- NULL
mydata_int        <- NULL
mymodel_int       <- NULL
mycontrols_int    <- NULL
ans_int_uncond    <- NULL
ans_co            <- NULL
mydata_co         <- NULL
mymodel_co        <- NULL
mycontrols_co     <- NULL
ans_cm            <- NULL
mydata_cm         <- NULL
mymodel_cm        <- NULL
mycontrols_cm     <- NULL

if (identical(Sys.getenv("RSENA_FULL_TESTS"), "1")) {

  # ego-effect model: density + recip + egoX(mybeh) + egoX*transTrip interaction
  # (used by test-sienaMargins.R for ego-wide perturbation integration tests)
  mynet_ego   <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
  mybeh_ego   <- coCovar(s50a[, 1])
  mydata_ego  <- sienaDataCreate(mynet_ego, mybeh_ego)
  mymodel_ego <- getEffects(mydata_ego)
  mymodel_ego <- includeEffects(mymodel_ego, egoX, interaction1 = "mybeh_ego",
                                name = "mynet_ego")
  mymodel_ego <- includeEffects(mymodel_ego, transTrip, name = "mynet_ego")
  mymodel_ego <- includeInteraction(mymodel_ego, egoX, transTrip,
                                    interaction1 = c("mybeh_ego", ""),
                                    name = "mynet_ego")
  mycontrols_ego <- sienaAlgorithmCreate(projname = NULL, n3 = 50,
                                         cond = FALSE, seed = 99)
  ans_ego <- siena07(
    mycontrols_ego,
    data    = mydata_ego,
    effects = mymodel_ego,
    returnDeps               = TRUE,
    returnChangeContributions = TRUE,
    returnDataFrame           = TRUE,
    silent                    = TRUE
  )

  # interaction-without-main: recip * outPop unspInt WITHOUT outPop main effect
  # (used by test-contributions.R and test-predict.R)
  mynet_int   <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
  mydata_int  <- sienaDataCreate(mynet_int)
  mymodel_int <- getEffects(mydata_int)
  mymodel_int <- includeInteraction(mymodel_int, recip, outPop, name = "mynet_int")
  mycontrols_int <- sienaAlgorithmCreate(projname = NULL, seed = 42, n3 = 50)
  ans_int <- siena07(
    mycontrols_int,
    data    = mydata_int,
    effects = mymodel_int,
    returnChangeContributions = TRUE,
    returnDataFrame           = TRUE,
    silent                    = TRUE
  )

  # cond=FALSE version of the same model
  mycontrols_int_uncond <- sienaAlgorithmCreate(
    projname = NULL, seed = 42, n3 = 50, cond = FALSE
  )
  ans_int_uncond <- siena07(
    mycontrols_int_uncond,
    data    = mydata_int,
    effects = mymodel_int,
    silent  = TRUE
  )

  # two-network co-evolution (50 actors, 2 periods per network)
  mynet_a    <- sienaDependent(array(c(s501, s502), dim = c(50, 50, 2)))
  mynet_b    <- sienaDependent(array(c(s502, s503), dim = c(50, 50, 2)))
  mydata_co  <- sienaDataCreate(mynet_a, mynet_b)
  mymodel_co <- getEffects(mydata_co)
  mymodel_co <- includeEffects(mymodel_co, transTrip, name = "mynet_a")
  mymodel_co <- includeEffects(mymodel_co, transTrip, name = "mynet_b")
  mycontrols_co <- sienaAlgorithmCreate(projname = NULL, seed = 1, n3 = 50, cond = FALSE)
  ans_co <- siena07(
    mycontrols_co,
    data    = mydata_co,
    effects = mymodel_co,
    returnChangeContributions = TRUE,
    returnDataFrame           = FALSE,
    silent                    = TRUE
  )

  # creation / maintenance effects (recip eval+creation, transTrip eval+endow)
  mynet_cm   <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
  mydata_cm  <- sienaDataCreate(mynet_cm)
  mymodel_cm <- getEffects(mydata_cm)
  mymodel_cm <- includeEffects(mymodel_cm, recip,     type = "creation")
  mymodel_cm <- includeEffects(mymodel_cm, transTrip, type = "eval")
  mymodel_cm <- includeEffects(mymodel_cm, transTrip, type = "endow")
  mycontrols_cm <- sienaAlgorithmCreate(projname = NULL, seed = 3, n3 = 50, cond = FALSE)
  ans_cm <- siena07(
    mycontrols_cm,
    data    = mydata_cm,
    effects = mymodel_cm,
    returnChangeContributions = TRUE,
    returnDataFrame           = FALSE,
    silent                    = TRUE
  )
}

}
