# Install/load RSiena if needed
# Uncomment if not installed:
# install.packages("RSiena")
# Or for development:
devtools::load_all(".")

library(RSiena)
library(data.table)

# Dummy data (replace with your actual s501, s502, s503)
set.seed(42)

mynet <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata <- sienaDataCreate(mynet)
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, transTrip, name = "mynet")
mycontrols <- sienaAlgorithmCreate(projname = NULL, n3 = 100, cond = FALSE)
ans <- siena07(
  mycontrols,
  data = mydata,
  effects = mymodel,
  returnChangeContributions = TRUE,
  returnDataFrame = TRUE,
  silent = TRUE
)

# Multicore setup
nbrNodes <- 6
useCluster <- TRUE
clusterType <- "PSOCK"

cat("Starting multicore predict.sienaFit...\n")
start_time <- Sys.time()
      pred_df <- predictDynamic(
        ans = ans,
        data = mydata,
        effects = mymodel,
        algorithm = mycontrols,
        useTieProb = TRUE,
        n3 = 100,
        nsim = 24,
        condition = "transTrip",
        level = "period",
        uncertainty = TRUE,
        useCluster = TRUE,
        clusterType = "FORK",
        nbrNodes = 6,
        silent = FALSE
      )

end_time <- Sys.time()
cat("Finished. Duration: ", end_time - start_time, "\n")

# Save output for inspection
saveRDS(pred_dt, "pred_dt.rds")
cat("Results saved to pred_dt.rds\n")

# Keep R session alive for monitoring
cat("Sleeping for 60 seconds so you can check Activity Monitor...\n")
Sys.sleep(60)