test_that("basic effects have correct target statistics", {
  mynet1 <- sienaDependent(array(c(tmp3, tmp4),dim=c(32, 32, 2)))
  mydata <- sienaDataCreate(mynet1)
  myeff<- getEffects(mydata)
  expect_equal(c(RSiena:::getTargets(mydata, myeff)), c(51,129,58))
})
