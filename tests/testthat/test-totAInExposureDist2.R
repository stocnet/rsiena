skip_on_cran()

################################################################################
### Use internal data set to test avInExposure effect ###
################################################################################

mynet <- sienaDependent(array(c(s501, s502), dim = c(50, 50, 2)))
sm1 <- 1 * (s50s[, 2] >= 2)
sm2 <- 1 * (s50s[, 3] >= 2)
sm2 <- pmax(sm1, sm2)
sm2[c(33, 28, 29, 44)] <- 1

mybeh <- sienaDependent(cbind(sm1, sm2), type = "behavior")
mydata <- sienaDataCreate(mynet, mybeh)

mymodel <- getEffects(mydata)

mycontrols <- sienaAlgorithmCreate(projname = NULL, seed = 1234, firstg = 0.1)
mycontrols2 <- sienaAlgorithmCreate(projname = NULL, seed = 1234, firstg = 0.001)

ans <- siena07(mycontrols, data = mydata, effects = mymodel,
               batch = TRUE, silent = TRUE)

mymodel0 <- updateTheta(mymodel, ans)

mymodel <- includeEffects(mymodel0, totAInExposureDist2, 
                          type = "rate", name = "mybeh", interaction1 = "mynet")

ans1 <- siena07(mycontrols2, data = mydata, effects = mymodel,
                batch = TRUE, silent = TRUE)

# Check the target statistics

test_that("Target statistics are correct", {

  divi <- function(x, y) {
          ifelse(y == 0, 0, x / y)
  }

adj <- mynet[, , 1]
beh <- mybeh[, , 1]
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

# To get the target statistic:
totAInExposureDist2_target <- sum((mybeh[, , 2] - mybeh[, , 1]) * stat)

  expect_equal(
    totAInExposureDist2_target,
    ans1$targets[5]
  )
})