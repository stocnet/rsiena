################################################################################
################################################################################
### file Test_anyInExposureDist2.R
################################################################################
################################################################################

################################################################################
### Check if RSiena can be loaded and the main estimation algorithm ###
################################################################################

# Use internal data set to test totTwoInStarAlt effect

mynet <- sienaDependent(array(c(s501, s502), dim = c(50, 50, 2)))
sm1 <- 1 * (s50s[, 2] >= 2)
sm2 <- 1 * (s50s[, 3] >= 2)
sm2 <- pmax(sm1, sm2)
sm2[c(33, 28, 29, 44)] <- 1

mybeh <- sienaDependent(cbind(sm1, sm2),
                        type = "behavior")
mydata <- sienaDataCreate(mynet, mybeh)

mymodel <- getEffects(mydata)
p <- 1
mymodel <- setEffect(mymodel, anyInExposureDist2, type= "rate", name = "mybeh", interaction1 = "mynet")

mycontrols <- sienaAlgorithmCreate(projname = NULL, seed = 1234, firstg = 0.01)

ans <- (siena07(
  mycontrols,
  batch = TRUE,
  silent = TRUE,
  data = mydata,
  effects = mymodel,
  returnChains = FALSE))

test_that("anyInExposureDist2 target statistic is correct", {

dist2behaviorsum <- t(mynet[, , 1]) %*% mybeh[, , 1]
n <- nrow(mynet[, , 1])

mat <- matrix(0, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    if (mynet[, , 1][i, j] == 1) {
      h_count <- dist2behaviorsum[j]
      # Exclude self if i also has behavior and i->j
      if (mynet[, , 1][i, j] == 1 && mybeh[, , 1][i] == 1) {
        h_count <- h_count - 1
      }
      mat[i, j] <- as.integer(h_count > 0)
    }
  }
}

mindist2beh <- rowSums(mat)

anyInExposureDist2_target <- sum((mybeh[, , 2] - mybeh[, , 1]) * mindist2beh)

  expect_equal(
    anyInExposureDist2_target,
    ans$targets[5]
  )

}
)
