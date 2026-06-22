skip_on_cran()

################################################################################
################################################################################
### file test-totGwInAltDist2_nc.R
################################################################################
################################################################################

################################################################################
### Check totGwInAltDist2_nc effect ###
################################################################################

# Use internal data set to test totGwInAltDist2_nc effect
mynet <- sienaDependent(array(c(s502, s503), dim = c(50, 50, 2)))
mybeh <- sienaDependent(s50a[, 2:3], type = "behavior")
mydata <- sienaDataCreate(mynet, mybeh)
mymodel <- getEffects(mydata)
p <- 69
mymodel <- setEffect(mymodel, totGwInAltDist2_nc,
  name = "mybeh", interaction1 = "mynet", parameter = p)
mycontrols <- sienaAlgorithmCreate(projname = NULL, seed = 42)
ans <- siena07(
  mycontrols,
  batch = TRUE,
  silent = TRUE,
  data = mydata,
  effects = mymodel,
  returnChains = FALSE)


# Check the target statistics

test_that("Target statistics are correct", {
  alpha <- 0.69
  adj <- mynet[, , 1]
  beh <- mybeh[, , 2]
  weighted_alter_sum <- numeric(nrow(adj))

  for (i in seq_len(nrow(adj))) {
    for (h in which(adj[i, ] == 1)) {
      total <- sum(adj[, h] * beh) - beh[i]
      weighted_alter_sum[i] <- weighted_alter_sum[i] +
        exp(alpha) * (1 - (1 - exp(-alpha))^total)
    }
  }

  totGwInAltDist2_nc_target <- sum(beh * weighted_alter_sum)

  expect_equal(
    totGwInAltDist2_nc_target,
    ans$targets[7]
  )
}
)

mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, totGwAltDist2_nc,
  name = "mybeh", interaction1 = "mynet", parameter = p)
mycontrols <- sienaAlgorithmCreate(projname = NULL, seed = 42)
ans <- siena07(
  mycontrols,
  batch = TRUE,
  silent = TRUE,
  data = mydata,
  effects = mymodel,
  returnChains = FALSE)


# Check the target statistics

test_that("Target statistics are correct", {
  alpha <- 0.69
  adj <- mynet[, , 1]
  beh <- mybeh[, , 2]
  weighted_alter_sum <- numeric(nrow(adj))

for (i in seq_len(nrow(adj))) {
  for (j in which(adj[i,] == 1)) {
    total <- sum(adj[j, ] * beh * (seq_len(nrow(adj)) != i))
    weighted_alter_sum[i] <- weighted_alter_sum[i] +
      exp(alpha) * (1 - (1 - exp(-alpha))^total)
  }
}

  totGwAltDist2_nc_target <- sum(beh * weighted_alter_sum)

  expect_equal(
    totGwAltDist2_nc_target,
    ans$targets[7]
  )
}
)

# MISSING IN CHECK EFFECTS!