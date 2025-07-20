devtools::load_all()
mynet <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata <- sienaDataCreate(mynet)
mymodel <- getEffects(mydata)
## TransitiveTriplets model
mymodel <- includeEffects(mymodel, transTrip, name = "mynet")

mycontrols <- sienaAlgorithmCreate(projname="test", n3 = 1000, cond = FALSE)
ans <- siena07(
  mycontrols,
  data = mydata,
  effects = mymodel, 
  returnChangeContributions = TRUE
  )

length(ans$changeContributions)


mycontrols$nsub <- 0
ans2 <- siena07(mycontrols, data=mydata, effects=mymodel, batch=TRUE, prevAns = ans, 
returnChangeContributions = TRUE)
length(ans2$changeContributions)

devtools::load_all()

myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=500, cond=TRUE)
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, density, recip, transTies, nbrDist2)
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)

## Not run: 
RIDynamics1 <- sienaRIDynamics(mydata, ans=ans, effects = myeff)
RIDynamics1



plot(RIDynamics1)

myalgorithm2 <- sienaAlgorithmCreate(nsub=2, n3=50, cond=TRUE)
ans2 <- siena07(myalgorithm2, data=mydata, effects=myeff, batch=TRUE)

RIDynamics2 <- sienaRIDynamics(mydata, ans=ans2, effects=myeff)
RIDynamics2

RIDynamics3 <- sienaRIDynamics(data=mydata, theta=c(ans2$rate,ans2$theta),
             algorithm=myalgorithm2, effects=myeff, intervalsPerPeriod=4)
RIDynamics3


mynet <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata <- sienaDataCreate(mynet)
mymodel <- getEffects(mydata)
## TransitiveTriplets model
mymodel <- includeEffects(mymodel, transTrip, name = "mynet")

alg <- sienaAlgorithmCreate(projname="test", n3 = 60)
alg$returnChangeContributions <- TRUE
ans <- siena07(alg, data=mydata, effects=mymodel)
str(ans$changeContributions)
length(ans$changeContributions[[1]][[1]])
head(ans$changeContributions[[1]][[1]][[1]])
