skip_on_cran()

################################################################################
################################################################################
### file test-totGwdspFFAlt_nc.R
################################################################################
################################################################################

################################################################################
### Check totGwdspFFAlt_nc effect ###
################################################################################

# Use internal data set to test totGwdspFFAlt_nc effect
mynet <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[, 2:3], type = "behavior")
mydata <- sienaDataCreate(mynet, mybeh)
mymodel <- getEffects(mydata)
p <- 69
mymodel <- setEffect(mymodel, totGwdspFFAlt_nc,
  name = "mybeh", interaction1 = "mynet", parameter = p)
mymodel <- includeEffects(mymodel, linear_nc, quad_nc,
    name = "mybeh")
mymodel <- includeEffects(mymodel, linear, quad,
    name = "mybeh", include = FALSE)
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
  # Calculate the mean behavior to test any behavior effect statistic

  linearshape_target <- sum(mybeh[, , 2])

  expect_equal(
    linearshape_target,
    ans$targets[5]
  )

  quadshape_target <- sum((mybeh[, , 2])^2)

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

  # just use instars2 <- adj %*% t(adj) and set the diagonal to 0 ?

  adj <- mynet[, , 1]
  alpha <- 0.69
  twopath_mat <- forward_twopaths(adj)
  weight_mat <- exp(alpha) * (1 - (1 - exp(-alpha))^twopath_mat)
  weighted_beh <- (mybeh[, , 2]) %*% weight_mat
  totGwdspFFAlttarget <- sum((mybeh[, , 2]) * weighted_beh)

  expect_equal(
    totGwdspFFAlttarget,
    ans$targets[7]
  )
}
)
