test_contributions <- function() {
  cat("\n== Testing Contribution Extraction ==\n")
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

#   # --- Test static contribution extraction
#   cat("* Static contributions (base R fallback):\n")
#   suppressPackageStartupMessages(detach("package:data.table", unload = TRUE))
#   stat_df <- calculateContribution(ans, mydata)
#   print(head(stat_df))
#   stopifnot(is.data.frame(stat_df) || inherits(stat_df, "data.table"))

  # --- If data.table is available, check with data.table
  if (requireNamespace("data.table", quietly = TRUE)) {
    library(data.table)
    cat("* Static contributions (data.table):\n")
    stat_dt <- calculateContribution(ans, mydata)
    print(head(stat_dt))
    stopifnot("data.table" %in% class(stat_dt))

    ## add comparison with "direct" data.frame output / c++ conversion
  }

  # --- Test dynamic contribution extraction
#   cat("* Dynamic contributions (base R fallback):\n")
#   dyn_df <- getChangeContributionsDynamic(
#     ans = ans,
#     data = mydata,
#     theta = c(ans$rate, ans$theta),
#     algorithm = mycontrols,
#     effects = mymodel,
#     depvar = "mynet",
#     returnDataFrame = TRUE
#   )
#   print(head(dyn_df))
#   stopifnot(is.matrix(dyn_df) || is.data.frame(dyn_df) || "data.table" %in% class(dyn_df))

  if (requireNamespace("data.table", quietly = TRUE)) {
    library(data.table)
    cat("* Dynamic contributions (data.table):\n")
    dyn_dt <- getChangeContributionsDynamic(
      ans = ans,
      data = mydata,
      theta = c(ans$rate,ans$theta),
      algorithm = mycontrols,
      effects = mymodel,
      depvar = "mynet",
      returnDataFrame = TRUE
    )
    print(head(dyn_dt))
    stopifnot("data.table" %in% class(dyn_dt))
    cat("* Dynamic contributions (data.table) with useChangeContributions = FALSE:\n")
    dyn_dt <- getChangeContributionsDynamic(
      ans = ans,
      data = mydata,
      theta = c(ans$rate, ans$theta),
      algorithm = mycontrols,
      effects = mymodel,
      depvar = "mynet",
      useChangeContributions = FALSE,
      n3 = 1000,
      returnDataFrame = TRUE
    )
    print(head(dyn_dt))
    stopifnot("data.table" %in% class(dyn_dt))
    cat("* Dynamic contributions (data.table) with new theta values:\n")
    new_theta <- MASS::mvrnorm(n=1,
                           mu = ans$theta,
                           Sigma = ans$covtheta)
    dyn_dt <- getChangeContributionsDynamic(
      ans = NULL,
      data = mydata,
      theta = c(ans$rate, new_theta),
      algorithm = mycontrols,
      effects = mymodel,
      depvar = "mynet",
      useChangeContributions = FALSE,
      n3 = 100,
      returnDataFrame = TRUE
    )
    print(head(dyn_dt))
    stopifnot("data.table" %in% class(dyn_dt))
  }

  cat("All contribution extractions passed.\n")
}

test_contributions()
