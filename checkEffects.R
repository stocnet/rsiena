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
### Note: this does not check the effect itself, only the statistic.
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
### check homXTransTrip, homXTransRecTrip
############################################

mynet <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))

# construct actor covariate
(two <- rep(1:2,25))
two <- coCovar(two, centered=FALSE)

mydata <- sienaDataCreate(mynet, two)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,transTrip)
myeff <- includeEffects(myeff,homXTransTrip,interaction1='two')
myeff
mymodel <- sienaModelCreate(projname=NULL, seed=1234)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(s502)
sum(s502*t(s502))
twop <- s502 %*% s502
diag(twop) <- 0
sum(twop*s502) # OK

eq <- 1*outer(two, two, "==") # same two
twop.eq <- (s502*eq) %*% s502
tt.eq <- s502 * eq * twop.eq
sum(tt.eq) # OK


myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,transRecTrip)
myeff <- includeEffects(myeff,homXTransRecTrip,interaction1='two')
myeff
mymodel <- sienaModelCreate(projname=NULL, seed=1234)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(s502)
sum(s502*t(s502))
twop <- s502 %*% s502
diag(twop) <- 0
sum(twop*(s502 * t(s502))) # OK

eq <- 1*outer(two, two, "==") # same two
twop.eq <- (s502*eq) %*% s502
tt.eq.r <- s502 * t(s502) * eq * twop.eq
sum(tt.eq.r) # OK

############################################
### check to, toU
############################################

advice <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
trust <- sienaDependent(array(c(s503, s501), dim=c(50, 50, 2)))
# dyadic covariate
suppressWarnings(mat <- matrix(c(0,1,0,0,4,0,2,0,0), 50,50))
dcov <- coDyadCovar(mat, centered=FALSE)
mydata <- sienaDataCreate(advice,trust,dcov)

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,to,name="trust",interaction1="advice")
myeff
mymodel <- sienaModelCreate(projname = NULL, seed=123)
(myans <- siena07(mymodel, data = mydata, effects = myeff))
myans$targets

# for to effect:
# W=advice , X=trust
WX <- s501 %*% s501
diag(WX) <- 0
WXX <- WX * s501
sum(diag(WXX)) # 0 OK
sum(WXX) # 86 OK

myeff <- getEffects(mydata)
myeff <- setEffect(myeff,to,name="trust",interaction1="advice", parameter=2)
myeff
mymodel <- sienaModelCreate(projname = NULL, seed=123)
(myans <- siena07(mymodel, data = mydata, effects = myeff))
myans$targets

# for to effect, parameter=2:
# W=advice , X=trust
WX <- s501 %*% s501
diag(WX) <- 0
SWXX <- sqrt(WX) * s501
sum(diag(SWXX)) # 0 OK
sum(SWXX) # 78.38478 OK

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,toU,name="trust",interaction1="advice", interaction2="dcov")
(myans <- siena07(mymodel, data = mydata, effects = myeff))
myans$targets
# for toU effect:
# W=advice, X=trust, U=dcov
WUX <- (dcov*s501) %*% s501
diag(WUX) <- 0
WUXX <- WUX * s501
sum(diag(WUXX)) # 0 OK
sum(WUXX) # 90 OK

############################################
### check toBack, MixedInXW
############################################

mynet1 <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaNet(array(c(s503, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet1, mynet2)

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,mixedInXW, toBack, name='mynet2',
                interaction1='mynet1')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets

# mixedInXW
twop <- (s502) %*% (s502)
sum(twop * s501) # 65 OK
# toBack
inst <- (s502) %*% t(s502)
sum(inst * s501) # 60 OK


myeff <- getEffects(mydata)
myeff <- setEffect(myeff,mixedInXW, parameter=2, name='mynet2',
                interaction1='mynet1')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
# mixedInXW, parameter 2
mixos <- t(s502) %*% s501
sum(s502 * sqrt(mixos)) #  58.45997 OK


############################################
### check cl.XWX, cl.XWX1, cl.XWX2
############################################

advice <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
trust <- sienaDependent(array(c(s503, s501), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(advice,trust)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, cl.XWX, name="trust", interaction1="advice")
mymodel <- sienaModelCreate(projname = NULL, seed=123)
myans <- siena07(mymodel, data = mydata, effects = myeff)
myans
myans$targets

# for cl.XWX effect:
XW <- s501 %*% s501
diag(XW) <- 0
XWX <- XW * s501
diag(XWX)
2*sum(XWX) # 172 OK

# for cl.XWX1 effect:
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, cl.XWX1, name="trust", interaction1="advice")
myans <- siena07(mymodel, data = mydata, effects = myeff)
myans
myans$targets
sum(XWX) # 86 OK

# for cl.XWX2 effect:
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, cl.XWX2, name="trust", interaction1="advice")
myans <- siena07(mymodel, data = mydata, effects = myeff)
myans
myans$targets
sum(XWX) # 86 OK

############################################
############################################

