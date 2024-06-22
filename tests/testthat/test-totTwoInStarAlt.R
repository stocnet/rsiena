################################################################################
################################################################################
### file Test-neweffect.R
################################################################################
################################################################################

################################################################################
### Check totTwoInStarAlt effect ###
################################################################################

# Use internal data set to test totTwoInStarAlt effect
mynet <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[, 2:3], type = "behavior")
mydata <- sienaDataCreate(mynet, mybeh)
mymodel <- getEffects(mydata)
p <- 1
mymodel <- setEffect(mymodel, totTwoInStarAlt, name = "mybeh", interaction1 = "mynet")
mycontrols <- sienaAlgorithmCreate(projname = NULL, seed = 42)
ans <- siena07(
  mycontrols,
  batch = TRUE,
  silent = TRUE,
  data = mydata,
  effects = mymodel,
  returnChains = FALSE
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

  divi <- function(x, y) {
    ifelse(y == 0, 0, x / y)
  }

  instarmat <- instars(mynet[, , 1])
  weighted <- divi(
    (mybeh[, , 2]) %*% instarmat,
    (mybeh[, , 2])  %*%  (instarmat > 0)
  )
  totTwoInStarAlt_target <- sum((mybeh[, , 2]) * weighted)

  expect_equal(
    totTwoInStarAlt_target,
    ans$targets[7]
  )

}
)
