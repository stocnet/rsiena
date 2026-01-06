skip_on_cran()

################################################################################
################################################################################
### file test-totAInAltDist2_nc.R
################################################################################
################################################################################

################################################################################
### Check totAInAltDist2_nc effect ###
################################################################################

# Use internal data set to test totAInAltDist2_nc effect
mynet <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[, 2:3], type = "behavior")
mydata <- sienaDataCreate(mynet, mybeh)
mymodel <- getEffects(mydata)
p <- 1
mymodel <- setEffect(mymodel, totAInAltDist2_nc,
  name = "mybeh", interaction1 = "mynet")
mycontrols <- sienaAlgorithmCreate(projname = NULL, seed = 42)
ans <- siena07(
  mycontrols,
  batch = TRUE,
  silent = TRUE,
  data = mydata,
  effects = mymodel,
  returnChains = FALSE,
  thetaBound = 100
)


# Check the target statistics

test_that("Target statistics are correct", {
  # Calculate the mean behavior to test any behavior effect statistic
  mbh <- mean(mybeh)
  # Calculate c_p for testing the average group parameter
  c_p <- ifelse(p <= 0.5, 0, p - mbh)

  linearshape_target <- sum(mybeh[, , 2] - mbh)

  expect_equal(
    linearshape_target,
    ans$targets[5]
  )

  quadshape_target <- sum((mybeh[, , 2] - mbh)^2)

  expect_equal(
    quadshape_target,
    ans$targets[6]
  )

  divi <- function(x, y) {
    ifelse(y == 0, 0, x / y)
  }

  adj <- mynet[, , 1]
  beh <- mybeh[, , 2]
  n <- nrow(adj)
  # Precompute in-degrees for all j
  indeg <- colSums(adj)

  # For each i, we want to compute:
  # stat[i] = sum_{h != i} z_h * sum_j (adj[i, j] * adj[h, j] / (indeg[j] - adj[i, j]))

  stat <- numeric(n)
  for (i in 1:n) {
    # vector of adj[i, j] for all j
    xi <- adj[i, ]
    # Denominator for each j (vector)
    denom <- indeg - xi
    denom[denom == 0] <- NA  # Avoid division by zero

    # For each h != i, compute the sum over j
    for (h in setdiff(1:n, i)) {
      xhj <- adj[h, ]
      # Elementwise: adj[i, j] * adj[h, j] / denom[j]
      contrib <- xi * xhj / denom
      contrib[is.na(contrib)] <- 0
      stat[i] <- stat[i] + beh[h] * sum(contrib)
    }
  }


  totAInAltDist2_nc_target <- sum(beh * stat)

  expect_equal(
    totAInAltDist2_nc_target,
    ans$targets[7]
  )

}
)
