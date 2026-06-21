skip_on_cran()

################################################################################
################################################################################
### file test-indegTotGroup.R
################################################################################
################################################################################

################################################################################
### Use internal data set to test indegTotGroup effect
################################################################################

## Test for three networks without shape effects

mynet <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- sienaDependent(s50a[,1:3], type="behavior")
mydata <- sienaDataCreate(mynet, mybeh)
mymodel <- getEffects(mydata)
p <- 1
mymodel <- includeEffects (mymodel, linear, quad,  name='mybeh', include=FALSE)

mymodel <- setEffect(mymodel,indegTotGroup, 
                     name='mybeh', interaction1='mynet', parameter=p)
mycontrols <- sienaAlgorithmCreate(projname=NULL)

ans <- siena07(
  mycontrols,
  batch = TRUE,
  silent = TRUE,
  data = mydata,
  effects = mymodel,
  returnChains = FALSE
)


# Check the target statistics
mbh <- mean(mybeh)
c_p <- ifelse(p <= 0.5, 0, p - mbh)

test_that("Centered target statistics are correct", {
  # Wave 2
  z_w2 <- mybeh[, 1, 2] - mbh
  deg_w1 <- colSums(mynet[,,1])
  w_sum_w2 <- sum(z_w2 * deg_w1)
  net_weight_w1 <- sum(deg_w1)
  
  # Wave 3
  z_w3 <- mybeh[, 1, 3] - mbh
  deg_w2 <- colSums(mynet[,,2])
  w_sum_w3 <- sum(z_w3 * deg_w2)
  net_weight_w2 <- sum(deg_w2)

  # Substantively: Outer behavior * (Inner weighted sum - Total network weight * c_p)
  indegedTotale_target <- sum( z_w2 * (w_sum_w2 - net_weight_w1 * c_p) ) +
                          sum( z_w3 * (w_sum_w3 - net_weight_w2 * c_p) )

  expect_equal(
    indegedTotale_target,
    ans$targets[7]
  )
})

mymodel <- getEffects(mydata)
p <- 1
mymodel <- includeEffects (mymodel, linear, quad,  name='mybeh', include=FALSE)

mymodel <- setEffect(mymodel,indegTotGroup_nc, 
                     name='mybeh', interaction1='mynet', parameter=p)
mycontrols <- sienaAlgorithmCreate(projname=NULL)

ans <- siena07(
  mycontrols,
  batch = TRUE,
  silent = TRUE,
  data = mydata,
  effects = mymodel,
  returnChains = FALSE
)


test_that("Non-centered target statistics are correct", {
  # Wave 2
  z_w2 <- mybeh[, 1, 2]
  deg_w1 <- colSums(mynet[,,1])
  w_sum_w2 <- sum(z_w2 * deg_w1)
  net_weight_w1 <- sum(deg_w1)
  
  # Wave 3
  z_w3 <- mybeh[, 1, 3]
  deg_w2 <- colSums(mynet[,,2])
  w_sum_w3 <- sum(z_w3 * deg_w2)
  net_weight_w2 <- sum(deg_w2)

  # Substantively: Outer behavior * (Inner weighted sum - Total network weight * c_p)
  indegedTotale_target <- sum( z_w2 * (w_sum_w2) ) +
                          sum( z_w3 * (w_sum_w3) )

  expect_equal(
    indegedTotale_target,
    ans$targets[7]
  )
})

mymodel <- getEffects(mydata)
p <- 1
mymodel <- includeEffects (mymodel, linear, quad,  name='mybeh', include=FALSE)

mymodel <- setEffect(mymodel,indegTotGroup_noEgo, 
                     name='mybeh', interaction1='mynet', parameter=p)
mycontrols <- sienaAlgorithmCreate(projname=NULL)

ans <- siena07(
  mycontrols,
  batch = TRUE,
  silent = TRUE,
  data = mydata,
  effects = mymodel,
  returnChains = FALSE
)

test_that("Centered Target statistics excluding Ego are correct", {
  z_w2 <- mybeh_centered[, 1, 2]
  z_w3 <- mybeh_centered[, 1, 3]
  
  deg_w1 <- colSums(mynet[,,1])
  deg_w2 <- colSums(mynet[,,2])

  global_sum_w2 <- sum(z_w2 * deg_w1)
  global_sum_w3 <- sum(z_w3 * deg_w2)

  indegedTotale_noEgo_target <- sum(
    z_w2 * (
      (global_sum_w2 - (z_w2 * deg_w1)) - ((sum(deg_w1) - deg_w1) * c_p)
    )
  ) +
  sum(
    z_w3 * (
      (global_sum_w3 - (z_w3 * deg_w2)) - ((sum(deg_w2) - deg_w2) * c_p)
    )
  )

  expect_equal(
    indegedTotale_noEgo_target,
    ans$targets[7]
  )
})

mymodel <- getEffects(mydata)
p <- 1
mymodel <- includeEffects (mymodel, linear, quad,  name='mybeh', include=FALSE)

mymodel <- setEffect(mymodel,indegTotGroup_nc_noEgo, 
                     name='mybeh', interaction1='mynet', parameter=p)
mycontrols <- sienaAlgorithmCreate(projname=NULL)

ans <- siena07(
  mycontrols,
  batch = TRUE,
  silent = TRUE,
  data = mydata,
  effects = mymodel,
  returnChains = FALSE
)

test_that("Non-Centered Target statistics excluding Ego are correct", {
  z_w2 <- mybeh[, 1, 2]
  z_w3 <- mybeh[, 1, 3]
  
  deg_w1 <- colSums(mynet[,,1])
  deg_w2 <- colSums(mynet[,,2])

  global_sum_w2 <- sum(z_w2 * deg_w1)
  global_sum_w3 <- sum(z_w3 * deg_w2)

  indegedTotale_noEgo_target <- sum(
    z_w2 * (
      (global_sum_w2 - (z_w2 * deg_w1))
    )
  ) +
  sum(
    z_w3 * (
      (global_sum_w3 - (z_w3 * deg_w2))
    )
  )

  expect_equal(
    indegedTotale_noEgo_target,
    ans$targets[7]
  )
})