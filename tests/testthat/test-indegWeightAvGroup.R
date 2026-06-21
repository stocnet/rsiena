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

# Calculate the mean behavior to test any behavior effect statistic
mbh <- mean(mybeh)
# Calculate c_p for testing the average group parameter
c_p <- ifelse(p <= 0.5, 0, p - mbh)

mybeh_centered <- mybeh - mbh

test_that("Centered Target statistics are correct", {

  indegedAverage_target <- sum(
    (mybeh_centered[,,2]) * 
      (
        (sum((mybeh_centered[,,2]) * colSums(mynet[,,1])) / sum(colSums(mynet[,,1]))) 
          - c_p 
      )
  ) +
  sum(
    (mybeh_centered[,,3]) * 
      (
        (sum((mybeh_centered[,,3]) * colSums(mynet[,,2])) / sum(colSums(mynet[,,2])) )
        - c_p )
  ) 

  expect_equal(
    indegedAverage_target,
    ans$targets[7]
  )

}
)

mymodel <- getEffects(mydata)
mymodel <- includeEffects (mymodel, linear, quad,  name='mybeh', include=FALSE)
mymodel <- setEffect(mymodel,indegAvGroup_nc, 
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

test_that("Non-centered garget statistics are correct", {

  indegedAverage_nc_target <- sum(
    (mybeh[,,2]) * 
        (sum((mybeh[,,2]) * colSums(mynet[,,1])) / sum(colSums(mynet[,,1]))) 
  ) +
  sum(
    (mybeh[,,3]) * 
        (sum((mybeh[,,3]) * colSums(mynet[,,2])) / sum(colSums(mynet[,,2])) )
  ) 

  expect_equal(
    indegedAverage_nc_target,
    ans$targets[7]
  )

}
)

mymodel <- getEffects(mydata)
p <- 1
mymodel <- includeEffects (mymodel, linear, quad,  name='mybeh', include=FALSE)
mymodel <- setEffect(mymodel,indegAvGroup_noEgo, 
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

test_that("Non-centered garget statistics excluding ego are correct", {

indegedAverage_noEgo_target <- sum(
    (mybeh_centered[,,2]) * 
      (
        ((sum((mybeh_centered[,,2]) * colSums(mynet[,,1])) - (mybeh_centered[,,2]) * colSums(mynet[,,1]))  / (sum(colSums(mynet[,,1])) - colSums(mynet[,,1]))) 
          - c_p 
      )
  ) +
  sum(
    (mybeh_centered[,,3]) * 
      (
        ((sum((mybeh_centered[,,3]) * colSums(mynet[,,2])) - (mybeh_centered[,,3]) * colSums(mynet[,,2])) / (sum(colSums(mynet[,,2])) - colSums(mynet[,,2])))
        - c_p )
  )


  expect_equal(
    indegedAverage_noEgo_target,
    ans$targets[7]
  )

}
)