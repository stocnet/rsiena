skip_on_cran()

################################################################################
################################################################################
### file Test_totInExposureDist2.R
################################################################################
################################################################################

################################################################################
### Use internal data set to test totInExposureDist2 effect ###
################################################################################


mynet <- sienaDependent(array(c(s501, s502), dim = c(50, 50, 2)))
sm1 <- 1 * (s50s[, 2] >= 2)
sm2 <- 1 * (s50s[, 3] >= 2)
sm2 <- pmax(sm1, sm2)
sm2[c(33, 28, 29, 44)] <- 1

mybeh <- sienaDependent(cbind(sm1, sm2),
                        type = "behavior")
mydata <- sienaDataCreate(mynet, mybeh)

mymodel <- getEffects(mydata)

mycontrols <- sienaAlgorithmCreate(projname = NULL, seed = 1234, firstg = 0.1)
mycontrols2 <- sienaAlgorithmCreate(projname = NULL, seed = 1234, firstg = 0.001)

ans <- siena07(mycontrols, data = mydata, effects = mymodel,
               batch = TRUE, silent = TRUE)

mymodel0 <- updateTheta(mymodel, ans)

mymodel <- includeEffects(mymodel0, totInExposureDist2, 
                          type = "rate", name = "mybeh", interaction1 = "mynet")

ans1 <- siena07(mycontrols2, data = mydata, effects = mymodel,
                batch = TRUE, silent = TRUE) # terrible conversion


# Check the target statistics

test_that("Target statistics are correct", {

  instars <- function(mat){
    # make an empty version of the matrix where we will store instar values
    matrix_vals <- mat
    matrix_vals[] <- 0
    # loop over the actors in the network, comparing pair-wise their values 
    for (i in 1:nrow(mat)){
      for (j in 1:nrow(mat)){
        a_row <- mat[i, c(-i, -j)]
        b_row <- mat[j, c(-i, -j)]
        matrix_vals[i, j] <- sum(a_row * b_row) # take sum of instars
      }
    }
    diag(matrix_vals) <- 0
    # return results
    return(matrix_vals)
  }

  instarmat <- instars(mynet[, , 1])

  weighted <- (mybeh[, , 1]) %*% instarmat


  totInExposureDist2_target <- sum((mybeh[, , 2] - mybeh[, , 1]) * weighted)

  expect_equal(
    totInExposureDist2_target,
    ans1$targets[5]
  )

}
)
