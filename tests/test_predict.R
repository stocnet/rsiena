# testthat::skip_on_cran()

library(RSiena)


test_predict <- function() {
  cat("\n== Testing Prediction Pipeline ==\n")
  # Setup
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

  # --- Test sienaPredict and downstream simProb and softmax/agg (base R fallback)
#   suppressPackageStartupMessages(detach("package:data.table", unload = TRUE))
#   pred <- sienaPredict(
#     ans = ans,
#     data = mydata,
#      useuseTieProb = TRUE,
#      condition = "transTrip",
#      level = "period"
#   )
#   print(head(pred))
#   stopifnot(is.data.frame(pred) || "data.table" %in% class(pred))
#   stopifnot("changeProb" %in% colnames(pred))
#   summ <- agg("changeProb", pred, level = "period")
#   print(summ)
#   stopifnot(is.data.frame(summ) || "data.table" %in% class(summ))

  # --- If data.table is available, check with data.table
  if (requireNamespace("data.table", quietly = TRUE)) {
    library(data.table)
    cat("* Prediction pipeline (data.table):\n")
    pred_dt <- sienaPredict(
      ans = ans,
      data = mydata,
      useTieProb = TRUE,
      nsim = 50,
      condition = "transTrip",
      level = "period",
      uncertainty = TRUE
    )
    print(head(pred_dt))
    stopifnot("data.table" %in% class(pred_dt) || is.data.frame(pred_dt))
    # detach("package:data.table", unload = TRUE)
  }

  # --- Test sienaPredictDynamic and downstream simProbDynamic and softmax/agg (base R fallback)
#   suppressPackageStartupMessages(detach("package:data.table", unload = TRUE))
#   pred <- sienaPredictDynamic(
#     ans = ans,
#     data = mydata,
#     useuseTieProb = TRUE,
#      n3 = 50,
#      nsim = 50,
#      condition = "density"
#   )
#   print(head(pred))
#   stopifnot(is.data.frame(pred) || "data.table" %in% class(pred))
#   stopifnot("changeProb" %in% colnames(pred))
#   summ <- agg("changeProb", pred, level = "period")
#   print(summ)
#   stopifnot(is.data.frame(summ) || "data.table" %in% class(summ))

  # --- If data.table is available, check with data.table
  if (requireNamespace("data.table", quietly = TRUE)) {
    library(data.table)
    cat("* Prediction pipeline (data.table):\n")
    pred_dt <- sienaPredictDynamic(
      ans = ans,
      data = mydata,
      effects = mymodel,
      algorithm = mycontrols,
      useTieProb = TRUE,
      n3 = 50,
      nsim = 20,
      condition = "density",
      uncertainty = FALSE
    )
    print(head(pred_dt))
    stopifnot("data.table" %in% class(pred_dt) || is.data.frame(pred_dt))
    # detach("package:data.table", unload = TRUE)
  }

  cat("Prediction pipeline tests passed.\n")
}

test_predict()
