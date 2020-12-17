################################################################################
################################################################################
### file checkEffects.R
################################################################################
################################################################################

################################################################################
### Procedure used for checking:
### use an internally available data set (two waves) without missings.
### See if estimation with default algorithm will converge,
### and check the target statistic.
################################################################################

library(RSiena)

############################################
### check from.w.ind
############################################

mynet1 <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaNet(array(c(s502, s503), dim=c(50, 50, 2)))
mynet3 <- sienaNet(array(c(s503, s501), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet1, mynet2, mynet3)
mymodel <- sienaModelCreate(projname=NULL, seed=1234)
myeff <- getEffects(mydata)

myeff <- includeEffects(myeff,from.w.ind, name='mynet2',
				interaction1='mynet1', interaction2='mynet1')
ans2 <- siena07(mymodel, data=mydata, effects=myeff)
ans2
ans2$targets

tp <- s501 %*% t(s501)
diag(s501) <- 0
outd <- diag(rowSums(s501))
ind <- diag(colSums(s501))
tpw <- s501 %*% outd %*% t(s501)
tpw2 <- s501 %*% ind %*% t(s501)
diag(tpw) <- 0
diag(tpw2) <- 0
sum(tpw*s503)
sum(tpw2*s503) # OK


(myeff <- setEffect(myeff,from.w.ind, name='mynet2',
				interaction1='mynet1', interaction2='mynet1', parameter=2))
ans2 <- siena07(mymodel, data=mydata, effects=myeff)
ans2
ans2$targets
ind3 <- diag(sqrt(colSums(s501)))
tpw3 <- s501 %*% ind3 %*% t(s501)
diag(tpw3) <- 0
sum(tpw3*s503) # OK


myeff <- getEffects(mydata)
(myeff <- setEffect(myeff,from.w.ind, name='mynet2',
				interaction1='mynet1', interaction2='mynet1', parameter=-1))
(ans2 <- siena07(mymodel, data=mydata, effects=myeff))
ans2$targets
ind3 <- diag(1/(colSums(s501) + 1))
tpw3 <- s501 %*% ind3 %*% t(s501)
diag(tpw3) <- 0
sum(tpw3*s503) # OK


myeff <- getEffects(mydata)
(myeff <- setEffect(myeff,from.w.ind, name='mynet2',
				interaction1='mynet1', interaction2='mynet3', parameter=-1))
(ans3 <- siena07(mymodel, data=mydata, effects=myeff))
ans3$targets
ind3 <- diag(1/(colSums(s503) + 1))
tpw3 <- s501 %*% ind3 %*% t(s501)
diag(tpw3) <- 0
sum(tpw3*s503) # OK


myeff <- getEffects(mydata)
(myeff <- setEffect(myeff,from.w.ind, name='mynet2',
				interaction1='mynet1', interaction2='mynet3', parameter=+1))
(ans3 <- siena07(mymodel, data=mydata, effects=myeff))
ans3$targets
ind3 <- diag(colSums(s503))
tpw3 <- s501 %*% ind3 %*% t(s501)
diag(tpw3) <- 0
sum(tpw3*s503) # OK


############################################
### check transtripX
############################################

mynet <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))

# construct actor covariate
(five <- rep(1:5,10))
five <- coCovar(five, centered=FALSE)

mydata <- sienaDataCreate(mynet, five)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,transTrip)
myeff <- includeEffects(myeff,transTripX,interaction1='five')
myeff
mymodel <- sienaModelCreate(projname=NULL, seed=1234)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(s502)
sum(s502*t(s502))
twop <- s502 %*% s502
diag(twop) <- 0
sum(twop*s502) # OK
twopc <- s502 %*% diag(rep(1:5,10)) %*% (s502)
diag(twopc) <- 0
sum(twopc*s502) # OK

############################################
############################################