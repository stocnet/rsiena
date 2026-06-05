skip_on_cran()

################################################################################
################################################################################
### file test-actAlt.R
################################################################################
################################################################################

################################################################################
### Check actAlt and totActAlt effects ###
################################################################################

# Use internal data set to test actAlt and totActAlt effects
mynet <- sienaDependent(array(c(s502, s503), dim = c(50, 50, 2)))
mybeh <- sienaDependent(s50a[, 2:3], type = "behavior")
mydata <- sienaDataCreate(mynet, mybeh)
mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, actAlt,
  name = "mybeh", interaction1 = "mynet")
mymodel <- setEffect(mymodel, totActAlt,
  name = "mybeh", interaction1 = "mynet")
mycontrols <- sienaAlgorithmCreate(projname = NULL, seed = 42)
ans <- siena07(
  mycontrols,
  batch = TRUE,
  silent = TRUE,
  data = mydata,
  effects = mymodel,
  returnChains = FALSE)


# Check the target statistics

test_that("Target statistics are correct for actAlt", {
  adj <- mynet[, , 1]
  beh <- mybeh[, , 2]
  beh_cent <- beh - mean(beh)

  # actAlt: z_i(centered) * sum_j x_ij * outdeg(j) / outdeg(i)
  actAlt_target <- numeric(nrow(adj))
  for (i in seq_len(nrow(adj))) {
    outdeg_i <- sum(adj[i, ])
    if (outdeg_i > 0) {
      actAlt_target[i] <- beh_cent[i] * sum(adj[i, ] * colSums(adj)) / outdeg_i
    }
  }

  expect_equal(
    sum(actAlt_target),
    ans$targets[7],
    tolerance = 1e-10
  )
})

test_that("Target statistics are correct for totActAlt", {
  adj <- mynet[, , 1]
  beh <- mybeh[, , 2]
  beh_cent <- beh - mean(beh)

  # totActAlt: z_i(centered) * sum_j x_ij * outdeg(j)
  totActAlt_target <- numeric(nrow(adj))
  for (i in seq_len(nrow(adj))) {
    totActAlt_target[i] <- beh_cent[i] * sum(adj[i, ] * colSums(adj))
  }

  expect_equal(
    sum(totActAlt_target),
    ans$targets[8],
    tolerance = 1e-10
  )
})
