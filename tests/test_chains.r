# probably should devtools::skip_on_cran()
# for interactive testing use
# devtools::load_all()
library(RSiena)

# Run simple test model ----

mynet <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata <- sienaDataCreate(mynet)
mymodel <- getEffects(mydata)
## TransitiveTriplets model
mymodel <- includeEffects(mymodel, transTrip, name = "mynet")

# Test returnChangeContributions when running siena07 directly ----

mycontrols <- sienaAlgorithmCreate(projname="test", n3 = 500, cond = FALSE)
ans <- siena07(
  mycontrols,
  data = mydata,
  effects = mymodel,
  returnChangeContributions = TRUE
)

length(ans$changeContributions)
head(ans$changeContributions[[1]][[1]][[1]])

# Test returnChangeContributions when setting msub = 0 and prevAns = ans, using batch mode

mycontrols$nsub <- 0
ans2 <- siena07(
  mycontrols,
  data=mydata,
  effects=mymodel,
  batch=TRUE,
  prevAns = ans,
  returnChangeContributions = TRUE)

length(ans2$changeContributions)
head(ans2$changeContributions[[1]][[1]][[1]])

# Test sienaRIDynamics ----

## Unconditional Estimation

### Use changeContributions from ans
RIDynamics1 <- sienaRIDynamics(data=mydata, ans=ans, useChangeContributions=TRUE, intervalsPerPeriod=10)
RIDynamics1
plot(RIDynamics1)
# if not exported to namespace use plot.sienaRIDynamics(RIDynamics1)

### Unconditional Estimation with model rerun - reduce n3 to speed up
mycontrols2 <- sienaAlgorithmCreate(nsub=2, n3=50, cond=FALSE)
mynet2 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata2 <- sienaDataCreate(mynet2)
mymoldel2 <- getEffects(mydata2)
mymoldel2 <- includeEffects(mymoldel2, density, recip, transTies, nbrDist2)

ans2 <- siena07(mycontrols2, data=mydata2, effects=mymoldel2, batch=TRUE)

RIDynamics2 <- sienaRIDynamics(mydata2, ans=ans2)
RIDynamics2
plot(RIDynamics2)
# if not exported to namespace use plot.sienaRIDynamics(RIDynamics2)


### Don't use ans but previously estimated coefficients

RIDynamics3 <- sienaRIDynamics(data=mydata2, theta=c(ans2$rate,ans2$theta),
             algorithm=mycontrols2, effects=mymoldel2, intervalsPerPeriod=10)
RIDynamics3
plot(RIDynamics3)
# if not exported to namespace use plot.sienaRIDynamics(RIDynamics3)

## Conditional Estimation
mycontrols3 <- sienaAlgorithmCreate(nsub=2, n3=50, cond=TRUE)
ans3 <- siena07(mycontrols3, data=mydata2, effects=mymoldel2, batch=TRUE)

RIDynamics4 <- sienaRIDynamics(mydata2, ans=ans3, effects = mymoldel2)
RIDynamics4
plot(RIDynamics4)
# if not exported to namespace use plot.sienaRIDynamics(RIDynamics4)
