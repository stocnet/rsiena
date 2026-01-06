skip_on_cran()

################################################################################
################################################################################
### file test-indegAvGroup.R
################################################################################
################################################################################

################################################################################
### Use internal data set to test indegAvGroup effect
################################################################################

## Test for three networks without shape effects

mynet <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- sienaDependent(s50a[,1:3], type="behavior")
mydata <- sienaDataCreate(mynet, mybeh)
mymodel <- getEffects(mydata)
p <- 1
mymodel <- includeEffects (mymodel, linear, quad,  name='mybeh', include=FALSE)

mymodel <- setEffect(mymodel,indegAvGroup, 
                     name='mybeh', interaction1='mynet', parameter=p)
mymodel
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

test_that("Target statistics are correct", {
  # Calculate the mean behavior to test any behavior effect statistic
  mbh <- mean(mybeh)
  # Calculate c_p for testing the average group parameter
  c_p <- ifelse(p <= 0.5, 0, p - mbh)

  indegedAverage_target <- sum((mybeh[,,2]-mbh) * (( sum((mybeh[,,2]-mbh) * colSums(mynet[,,1])) / sum(colSums(mynet[,,1])) ) - c_p )) + 
  sum((mybeh[,,3]-mbh) * (( sum((mybeh[,,3]-mbh) * colSums(mynet[,,2])) / sum(colSums(mynet[,,2])) ) - c_p ))

  expect_equal(
    indegedAverage_target,
    ans$targets[7]
  )

}
)
