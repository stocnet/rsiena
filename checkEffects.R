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
### Checks of later effects are appended at the end.
################################################################################

# library(RSiena)

################################################################################
### check from.w.ind
################################################################################

mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
mynet3 <- sienaDependent(array(c(s503, s501), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet1, mynet2, mynet3)
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=1234)
mymodel <- getEffects(mydata)

mymodel <- includeEffects(mymodel,from.w.ind, name='mynet2',
				interaction1='mynet1', interaction2='mynet1')
ans2 <- siena07(mycontrols, data=mydata, effects=mymodel)
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

(mymodel <- setEffect(mymodel,from.w.ind, name='mynet2',
				interaction1='mynet1', interaction2='mynet1', parameter=2))
ans2 <- siena07(mycontrols, data=mydata, effects=mymodel)
ans2
ans2$targets
ind3 <- diag(sqrt(colSums(s501)))
tpw3 <- s501 %*% ind3 %*% t(s501)
diag(tpw3) <- 0
sum(tpw3*s503) # OK

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,from.w.ind, name='mynet2',
				interaction1='mynet1', interaction2='mynet1', parameter=-1))
(ans2 <- siena07(mycontrols, data=mydata, effects=mymodel))
ans2$targets
ind3 <- diag(1/(colSums(s501) + 1))
tpw3 <- s501 %*% ind3 %*% t(s501)
diag(tpw3) <- 0
sum(tpw3*s503) # OK

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,from.w.ind, name='mynet2',
				interaction1='mynet1', interaction2='mynet3', parameter=-1))
(ans3 <- siena07(mycontrols, data=mydata, effects=mymodel))
ans3$targets
ind3 <- diag(1/(colSums(s503) + 1))
tpw3 <- s501 %*% ind3 %*% t(s501)
diag(tpw3) <- 0
sum(tpw3*s503) # OK

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,from.w.ind, name='mynet2',
				interaction1='mynet1', interaction2='mynet3', parameter=+1))
(ans3 <- siena07(mycontrols, data=mydata, effects=mymodel))
ans3$targets
ind3 <- diag(colSums(s503))
tpw3 <- s501 %*% ind3 %*% t(s501)
diag(tpw3) <- 0
sum(tpw3*s503) # OK

################################################################################
### check transtripX
################################################################################

mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))

# construct actor covariate
(five <- rep(1:5,10))
five <- coCovar(five, centered=FALSE)

mydata <- sienaDataCreate(mynet, five)
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,transTrip)
mymodel <- includeEffects(mymodel,transTripX,interaction1='five')
mymodel
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=1234)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
sum(s502)
sum(s502*t(s502))
twop <- s502 %*% s502
diag(twop) <- 0
sum(twop*s502) # OK
twopc <- s502 %*% diag(rep(1:5,10)) %*% (s502)
diag(twopc) <- 0
sum(twopc*s502) # OK

################################################################################
### check homXTransTrip, homXTransRecTrip
################################################################################

mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))

# construct actor covariate
(two <- rep(1:2,25))
two <- coCovar(two, centered=FALSE)

mydata <- sienaDataCreate(mynet, two)
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,transTrip)
mymodel <- includeEffects(mymodel,homXTransTrip,interaction1='two')
mymodel
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=1234)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
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

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,transRecTrip)
mymodel <- includeEffects(mymodel,homXTransRecTrip,interaction1='two')
mymodel
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=1234)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
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

################################################################################
### check to, toU
################################################################################

advice <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
trust <- sienaDependent(array(c(s503, s501), dim=c(50, 50, 2)))
# dyadic covariate
suppressWarnings(mat <- matrix(c(0,1,0,0,4,0,2,0,0), 50,50))
dcov <- coDyadCovar(mat, centered=FALSE)
mydata <- sienaDataCreate(advice,trust,dcov)

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,to,name="trust",interaction1="advice")
mymodel
mycontrols <- sienaAlgorithmCreate(projname = NULL, seed=123)
(myans <- siena07(mycontrols, data = mydata, effects = mymodel))
myans$targets

# for to effect:
# W=advice , X=trust
WX <- s501 %*% s501
diag(WX) <- 0
WXX <- WX * s501
sum(diag(WXX)) # 0 OK
sum(WXX) # 86 OK

mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel,to,name="trust",interaction1="advice", parameter=2)
mymodel
mycontrols <- sienaAlgorithmCreate(projname = NULL, seed=123)
(myans <- siena07(mycontrols, data = mydata, effects = mymodel))
myans$targets

# for to effect, parameter=2:
# W=advice , X=trust
WX <- s501 %*% s501
diag(WX) <- 0
SWXX <- sqrt(WX) * s501
sum(diag(SWXX)) # 0 OK
sum(SWXX) # 78.38478 OK

mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel,to,name="trust",interaction1="advice", parameter=3)
mymodel
mycontrols <- sienaAlgorithmCreate(projname = NULL, seed=123)
(myans <- siena07(mycontrols, data = mydata, effects = mymodel))
myans$targets

# for to effect, parameter=3:
# W=advice , X=trust
WX <- s501 %*% s501
diag(WX) <- 0
tWXX <- 1*(WX >= 1) * s501
sum(diag(tWXX)) # 0 OK
sum(tWXX) # 73 OK

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,toU,name="trust",interaction1="advice", interaction2="dcov")
(myans <- siena07(mycontrols, data = mydata, effects = mymodel))
myans$targets
# for toU effect:
# W=advice, X=trust, U=dcov
WUX <- (dcov*s501) %*% s501
diag(WUX) <- 0
WUXX <- WUX * s501
sum(diag(WUXX)) # 0 OK
sum(WUXX) # 90 OK

################################################################################
### check toBack, MixedInXW
################################################################################

mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaDependent(array(c(s503, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet1, mynet2)

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,mixedInXW, toBack, name='mynet2',
                interaction1='mynet1')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets

# mixedInXW
twop <- (s502) %*% (s502)
sum(twop * s501) # 65 OK
# toBack
inst <- (s502) %*% t(s502)
sum(inst * s501) # 60 OK


mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel,mixedInXW, parameter=2, name='mynet2',
                interaction1='mynet1')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# mixedInXW, parameter 2
mixos <- t(s502) %*% s501
sum(s502 * sqrt(mixos)) #  58.45997 OK

################################################################################
### check cl.XWX, cl.XWX1, cl.XWX2
################################################################################

advice <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
trust <- sienaDependent(array(c(s503, s501), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(advice,trust)
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, cl.XWX, name="trust", interaction1="advice")
mycontrols <- sienaAlgorithmCreate(projname = NULL, seed=123)
myans <- siena07(mycontrols, data = mydata, effects = mymodel)
myans
myans$targets

# for cl.XWX effect:
XW <- s501 %*% s501
diag(XW) <- 0
XWX <- XW * s501
diag(XWX)
2*sum(XWX) # 172 OK

# for cl.XWX1 effect:
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, cl.XWX1, name="trust", interaction1="advice")
(myans <- siena07(mycontrols, data = mydata, effects = mymodel))
myans$targets
sum(XWX) # 86 OK

# for cl.XWX2 effect:
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, cl.XWX2, name="trust", interaction1="advice")
(myans <- siena07(mycontrols, data = mydata, effects = mymodel))
myans$targets
sum(XWX) # 86 OK


################################################################################
### check avSim, totSim, avInSim, totInSim, avInSimPopAlt, totInSimPopAlt
################################################################################

mynet <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[,2:3], type="behavior")
mydata <- sienaDataCreate(mynet, mybeh)
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,avSim,name='mybeh',interaction1='mynet')
mymodel
mycontrols <- sienaAlgorithmCreate(projname=NULL)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# linear shape and quadratic shape
(mbh <- mean(mybeh))
sum(mybeh[,,2] - mbh) # OK linear shape
sum((mybeh[,,2] - mbh)^2) # OK quadratic shape

# for avSim effect:
# behavior at wave 2:
zz <- mybeh[,,2]
rangez <- attr(mydata$depvars[[2]], 'range')
# similarity at wave 2:
simi2 <- outer(zz, zz, function(x,y){1 - abs(x-y)/rangez})
simi2 <- simi2 - attr(mydata$depvars[[2]], 'simMean')
diag(simi2) <- 0
divi <- function(a,b){ifelse(b==0, 0, a/b)}
outdeg <- rowSums(s502)
outdegs <- outer(outdeg, outdeg, function(x,y){divi(1,x)})
diag(outdegs) <- 0
sum(s502 * simi2*outdegs) # OK avSim
# This can also be calculated by pre-multiplication by a diagonal matrix:
sum(diag(divi(1,outdeg)) %*% s502 * simi2 ) # OK avSim

# totSim effect:
mymodel <- getEffects(mydata)
# It turns out that for this data set, to get an estimate for totSim
# you first have to include the outdeg effect.
mymodel <- includeEffects(mymodel,outdeg,name='mybeh',interaction1='mynet')
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
mymodel <- includeEffects(mymodel,totSim,name='mybeh',interaction1='mynet')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel, prevAns=ans))
ans$targets
sum(s502 * simi2 ) # OK totSim
# Now try to get a model including totSim but without outdeg
mymodel <- includeEffects(mymodel,outdeg,name='mybeh',interaction1='mynet', include=FALSE)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel, prevAns=ans))
mymodel <- updateTheta(mymodel, ans)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))

# for avSimPopAlt effect:
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,avSimPopAlt,name='mybeh',interaction1='mynet')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets

indeg <- colSums(s502)
sum(diag(divi(1,outdeg)) %*% (s502 * simi2) %*% diag(indeg)) # avSimPopAlt OK

# for avSimPopEgo effect:
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,avSimPopEgo,name='mybeh',interaction1='mynet')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
indeg <- colSums(s502)
sum(diag(divi(indeg,outdeg)) %*% s502 * simi2 ) #  OK avSimPopEgo

# for totSimPopAlt effect:
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,totSimPopAlt,name='mybeh',interaction1='mynet')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
indeg <- colSums(s502)
sum((s502 * simi2) %*% diag(indeg)) # avSimPopAlt OK

# for avInSim effect:
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,avInSim,name='mybeh',interaction1='mynet')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets

indeg <- colSums(s502)
indegs <- outer(indeg, indeg, function(x,y){divi(1,x)})
diag(indegs) <- 0
sum(t(s502) * simi2*indegs) # OK avInSim
# alternative:
sum(diag(divi(1,indeg)) %*% t(s502) * simi2 ) # OK avInSim

# for avInSimPopAlt effect:
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,avInSimPopAlt,name='mybeh',interaction1='mynet')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
sum(diag(divi(1,indeg)) %*% t(s502) * simi2  %*% diag(indeg)) # OK avInSimPopAlt

# for totInSim effect:
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,totInSim,name='mybeh',interaction1='mynet')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
sum(t(s502) * simi2 ) # OK totInSim

# for totInSimPopAlt effect:
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,totInSimPopAlt,name='mybeh',interaction1='mynet')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
sum( t(s502) * simi2  %*% diag(indeg)) # OK totInSimPopAlt


################################################################################
### check avAttHigher, avAttLower, totAttHigher, totAttLower
################################################################################

mynet <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[,2:3], type="behavior")
mydata <- sienaDataCreate(mynet, mybeh)
mymodel <- getEffects(mydata)
mycontrols <- sienaAlgorithmCreate(projname=NULL)
# for avAttHigher effect:
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,avAttHigher,name='mybeh',interaction1='mynet')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets

# behavior at wave 2:
zz <- mybeh[,,2]
rangez <- attr(mydata$depvars[[2]], 'range')
# higher similarity at wave 2:
simim <- outer(zz, zz, function(x,y){1 - pmax(x-y,0)/rangez})
diag(simim) <- 0
divi <- function(a,b){ifelse(b==0, 0, a/b)}
outdeg <- rowSums(s502)
sum(diag(divi(1,outdeg)) %*% s502 * t(simim) ) # OK avAttHigher

# for avAttLower effect:
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,avAttLower,name='mybeh',interaction1='mynet')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
sum(diag(divi(1,outdeg)) %*% s502 * (simim) ) # OK avAttLower

# for totAttHigher effect:
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,totAttHigher,name='mybeh',interaction1='mynet')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
sum(s502 * t(simim) ) # OK totAttHigher

# for totAttLower effect:
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,totAttLower,name='mybeh',interaction1='mynet')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
sum(s502 * (simim) ) # OK totAttLower


################################################################################
### check crprodInActIntn
################################################################################

mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaDependent(array(c(s503, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet1, mynet2)

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, crprodInActIntn, name='mynet2',
                        interaction1='mynet1')
mymodel

# for crprodInActIntn effect (parameter = 2):
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
ctr <- mean(c(colSums(s501), colSums(s502)))
sum(rowSums(s502 * s501) * (sqrt(colSums(s501)) - sqrt(ctr))) # -4.937096 OK

# for crprodInActIntn effect (parameter = 1):
mymodel <- setEffect(mymodel, crprodInActIntn, name='mynet2',
                   interaction1='mynet1', parameter=1)
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
sum(rowSums(s502 * s501) * (colSums(s501) - ctr)) # -3.53 OK


################################################################################
### check rateX
################################################################################


mynet <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[,2:3], type="behavior")
mycov <- coCovar(s50s[,2])
mydata <- sienaDataCreate(mynet, mybeh, mycov)
mymodel <- getEffects(mydata)
(mymodel <- includeEffects(mymodel,RateX,type='rate',
							name='mybeh',interaction1='mycov'))
mycontrols <- sienaAlgorithmCreate(projname=NULL)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
sum(abs(mybeh[,,1]-mybeh[,,2])) # OK rate
mco <- mean(mycov)
sum((mycov-mco)*abs(mybeh[,,1]-mybeh[,,2])) # OK RateX

mymodel <- getEffects(mydata)
(mymodel <- includeEffects(mymodel,RateX,type='rate',
							name='mybeh',interaction1='mybeh'))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
sum(abs(mybeh[,,1]-mybeh[,,2])) # OK rate
mbh <- mean(mybeh)
sum((mybeh[,,1])*abs(mybeh[,,1]-mybeh[,,2])) # OK RateX non-centered

################################################################################
### check totExposure and other diffusion effects
################################################################################


mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
sm1 <- 1*(s50s[,2] >= 2)
sm2 <- 1*(s50s[,3] >= 2)
sm2 <- pmax(sm1,sm2)
sm2[c(33,28,29,44)] <- 1
sminfl1 <- apply(s501,1,function(x){sum(x*sm1)})
table(sminfl1*((sm2-sm1))) # no values larger than 2.
sm1r <- sm1
sm1r[11] <- 0
table(sminfl1*((sm2-sm1r))) # one value larger than 2.
mybeh <- sienaDependent(cbind(sm1,sm2), type="behavior")
mybeh.r <- sienaDependent(cbind(sm1r,sm2), type="behavior")
mycov <- coCovar(s50a[,1])
(mydata <- sienaDataCreate(mynet, mybeh, mycov))
(mydata.r <- sienaDataCreate(mynet, mybeh.r, mycov))
# mydata.r for testing negative parameters

mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=1234, firstg=0.01)
mycontrols2 <- sienaModelCreate(projname=NULL, seed=1234, firstg=0.1)
mycontrols4 <- sienaModelCreate(projname=NULL, seed=1234, firstg=0.001)

mymodel <- getEffects(mydata)
(ans <- siena07(mycontrols2, data=mydata, effects=mymodel))
mymodel0 <- updateTheta(mymodel, ans)

mymodel.r <- getEffects(mydata.r)
(ans.r <- siena07(mycontrols2, data=mydata.r, effects=mymodel.r))
mymodel0.r <- updateTheta(mymodel.r, ans.r)

mymodel <- includeEffects(mymodel0,avExposure,type='rate',
							name='mybeh',interaction1='mynet')
(ans1 <- siena07(mycontrols4, data=mydata, effects=mymodel))
ans1$targets
sum(sm2-sm1) # OK rate

divi <- function(x,y){ifelse(y==0, 0, x/y)}
amean <- function(x,z){ifelse((sum(x)>0), sum(x*z)/sum(x), 0)}
eff <- apply(s501,1,function(x){amean(x,sm2)})
eff00 <- apply(s501,1,function(x){amean(x,sm1)})
sum((sm2-sm1)*eff) #  this is not the target statistic for avExposure
sum((sm2-sm1)*eff00) # OK  avExposure

(mymodel5 <- setEffect(mymodel0,avExposure,type='rate', parameter=2,
							name='mybeh',interaction1='mynet'))
(ans5 <- siena07(mycontrols4, data=mydata, effects=mymodel5, prevAns=ans1))
ans5$targets
eff3 <- apply(s501,1,function(x){sum(x*sm1)})
eff35 <- eff3
eff35[eff35==1] <- 0
outd <- rowSums(s501)
sum((sm2-sm1)*divi(eff35, outd)) #  OK avExposure p=2

(mymodel6.r <- setEffect(mymodel.r, avExposure,type='rate', parameter=-2,
							name='mybeh.r',interaction1='mynet'))
(ans6.r <- siena07(mycontrols4, data=mydata.r, effects=mymodel6.r, prevAns=ans5))
ans6.r$targets
eff3.r <- apply(s501,1,function(x){sum(x*sm1r)})
eff56 <- eff3.r
eff56[eff56==1] <- 0
eff56[eff56>2] <- 2
sum((sm2-sm1r)*divi(eff56, outd)) #   OK avExposure p=-2

(mymodel2 <- setEffect(mymodel0,totExposure,type='rate',
							name='mybeh',interaction1='mynet'))
(ans2 <- siena07(mycontrols4, data=mydata, effects=mymodel2))
ans2$targets
sum(sm2-sm1) # OK rate
eff3 <- apply(s501,1,function(x){sum(x*sm1)})
sum((sm2-sm1)*eff3) # OK totExposure

(mymodel3 <- setEffect(mymodel0,totExposure,type='rate', parameter=2,
							name='mybeh',interaction1='mynet'))
(ans3 <- siena07(mycontrols4, data=mydata, effects=mymodel3, prevAns=ans2))
ans3$targets
eff3 <- apply(s501,1,function(x){sum(x*sm1)})
eff32 <- eff3
eff32[eff32==1] <- 0
sum((sm2-sm1)*eff32) # OK totExposure p=2

(mymodel4.r <- setEffect(mymodel0.r,totExposure,type='rate', parameter=-2,
							name='mybeh.r',interaction1='mynet'))
(ans4.r <- siena07(mycontrols4, data=mydata.r, effects=mymodel4.r, prevAns=ans3))
ans4.r$targets
eff3.r <- apply(s501,1,function(x){sum(x*sm1r)})
eff34.r <- eff3.r
eff34.r[eff34.r==1] <- 0
eff34.r[eff34.r>2] <- 2
sum((sm2-sm1r)*eff34.r) # OK totExposure p=-2

(mymodel6 <- setEffect(mymodel,infectIn,type='rate',
							name='mybeh',interaction1='mynet'))
(ans6 <- siena07(mycontrols4, data=mydata, effects=mymodel6))
ans6$targets
ind <- colSums(s501)
eff3 <- apply(s501,1,function(x){sum(x*sm1)})
eff7 <- vapply(1:50, function(i){sum(s501[i,]*ind*sm1)}, FUN.VALUE=1)
sum((sm2-sm1)*eff7) #  OK infectIn

(mymodel7 <- setEffect(mymodel,infectIn,type='rate', parameter=2,
							name='mybeh',interaction1='mynet'))
(ans7 <- siena07(mycontrols4, data=mydata, effects=mymodel7, prevAns=ans6))
ans7$targets
eff72 <- eff7
eff72[eff3<2] <- 0
sum((sm2-sm1)*eff72) #  OK infectIn p=2

(mymodel8 <- setEffect(mymodel0,infectOut,type='rate', parameter=2,
							name='mybeh', interaction1='mynet'))
(ans8 <- siena07(mycontrols4, data=mydata, effects=mymodel8))
ans8$targets
outd <- rowSums(s501)
eff3 <- apply(s501,1,function(x){sum(x*sm1)})
eff8 <- vapply(1:50, function(i){sum(s501[i,]*outd*sm1)}, FUN.VALUE=1)
eff82 <- eff8
eff82[eff3<2] <- 0
sum((sm2-sm1)*eff82) #  OK infectOut p=2

(mymodel9 <- setEffect(mymodel0,infectCovar,type='rate',
							name='mybeh', interaction1='mynet',
							interaction2='mycov'))
(ans9 <- siena07(mycontrols4, data=mydata, effects=mymodel9))
ans9$targets
(mcm <- mean(mycov))
eff9 <- vapply(1:50, function(i){sum(s501[i,]*(mycov-mcm)*sm1)}, FUN.VALUE=1)
sum((sm2-sm1)*eff9) #   infectCovar OK

(mymodelA <- setEffect(mymodel.r,infectCovar,type='rate', parameter=2,
							name='mybeh.r', interaction1='mynet',
							interaction2='mycov'))
(ansA <- siena07(mycontrols4, data=mydata.r, effects=mymodelA, prevAns=ans9))
ansA$targets #
eff3 <- apply(s501,1,function(x){sum(x*sm1r)})
(mcm <- mean(mycov))
effA <- vapply(1:50, function(i){sum(s501[i,]*(mycov-mcm)*sm1r)}, FUN.VALUE=1)
effA2 <- effA
effA2[eff3==1] <- 0
sum((sm2-sm1r)*effA2) # OK infectCovar p=2

(mymodelB <- setEffect(mymodel0.r,susceptAvIn,type='rate',
							name='mybeh.r',interaction1='mynet'))
(ansB <- siena07(mycontrols4, data=mydata.r, effects=mymodelB))
ansB$targets
eff3 <- apply(s501,1,function(x){sum(x*sm1r)})
outd <- rowSums(s501)
sum((sm2-sm1r)*colSums(s501)*divi(eff3, outd)) #  OK susceptAvIn

(mymodelC <- setEffect(mymodel0.r, susceptAvIn,type='rate', parameter=2,
							name='mybeh.r',interaction1='mynet'))
(ansC <- siena07(mycontrols4, data=mydata.r, effects=mymodelC, prevAns=ansB))
ansC$targets
eff35 <- eff3
eff35[eff3==1] <- 0
sum((sm2-sm1r)*colSums(s501)*divi(eff35, outd)) #  OK susceptAvIn p=2

(mymodelD <- setEffect(mymodel0.r,susceptAvIn,type='rate', parameter=-2,
							name='mybeh.r',interaction1='mynet'))
(ansD <- siena07(mycontrols4, data=mydata.r, effects=mymodelD, prevAns=ansC))
ansD$targets
eff35[eff3>2] <- 2
outd <- rowSums(s501)
sum((sm2-sm1r)*colSums(s501)*divi(eff35, outd)) #  OK susceptAvIn p=-2

(mymodelE <- setEffect(mymodel0.r,susceptAvCovar,type='rate', parameter=0,
							name='mybeh.r',interaction1='mynet',
							interaction2='mycov'))
(ansE <- siena07(mycontrols4, data=mydata.r, effects=mymodelE))
ansE$targets
eff3 <- apply(s501,1,function(x){sum(x*sm1r)})
eff35 <- eff3
(mcm <- mean(mycov))
outd <- rowSums(s501)
sum((sm2-sm1r)*(mycov-mcm)*divi(eff35, outd)) # susceptAvCovar OK

(mymodelE1 <- setEffect(mymodel0.r,susceptAvCovar,type='rate', parameter=2,
							name='mybeh.r',interaction1='mynet',
							interaction2='mycov'))
(ansE1 <- siena07(mycontrols4, data=mydata.r, effects=mymodelE1, prevAns=ansE))
ansE1$targets
eff35[eff3==1] <- 0
sum((sm2-sm1r)*(mycov-mcm)*divi(eff35, outd)) # susceptAvCovar p=2 OK.

(mymodelF <- setEffect(mymodel0.r,susceptAvCovar,type='rate', parameter=-2,
							name='mybeh.r',interaction1='mynet',
							interaction2='mycov'))
(ansF <- siena07(mycontrols4, data=mydata.r, effects=mymodelF, prevAns=ansE1))
ansF$targets
eff35[eff3>2] <- 2
sum((sm2-sm1r)*(mycov-mcm)*divi(eff35, outd)) # susceptAvCovar p=-2 OK


################################################################################
### check outOutActIntn and outOutAvIntn
################################################################################

mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
# construct actor covariate
in1 <- colSums(s501)
out1 <- rowSums(s501)
center <- in1 + out1
center <- 1*(center >= 5)
central <- coCovar(center)
mydata <- sienaDataCreate(mynet1, mynet2, central)
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=1234)

mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, outOutActIntn, name="mynet1", interaction1="mynet2", parameter=1)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
sum(rowSums(s502)) # ok outdegree
sum(rowSums(s502 * t(s502))) # OK recip
(avdeg <- mean(rowSums(s502) + rowSums(s501))/2)
sum(rowSums(s502)* (s502 %*% (rowSums(s502) - avdeg))) # OK outOutActIntn

mymodel2 <- getEffects(mydata)
mymodel2 <- setEffect(mymodel2, outOutActIntn, name="mynet1", interaction1="mynet2", parameter=2)
(ans2 <- siena07(mycontrols, data=mydata, effects=mymodel2))
ans2$targets
sum(rowSums(s502)* (s502 %*% (sqrt(rowSums(s502)) - sqrt(avdeg)))) # OK outOutActIntn p=2

mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, outOutAvIntn, name="mynet1", interaction1="mynet2", parameter=1)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
(avdeg <- mean(rowSums(s502) + rowSums(s501))/2)
sum((s502 %*% (rowSums(s502) - avdeg)) ) # OK outOutAvIntn, note that rowSums(s502) cancel

mymodel2 <- getEffects(mydata)
mymodel2 <- setEffect(mymodel2, outOutAvIntn, name="mynet1", interaction1="mynet2", parameter=2)
(ans2 <- siena07(mycontrols, data=mydata, effects=mymodel2))
ans2$targets
sum( (s502 %*% (sqrt(rowSums(s502)) - sqrt(avdeg)))) # OK outOutAvIntn p=2


################################################################################
### check inRateInv and inRateLog
################################################################################


mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet)
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=534, cond=FALSE, firstg=0.01)
mycontrols2 <- sienaModelCreate(projname=NULL, seed=534, cond=FALSE, firstg=0.02)
mymodel <- getEffects(mydata)
effectsDocumentation(mymodel)
mymodel <- includeEffects(mymodel, inRate, type='rate')
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
(ans1 <- siena07(mycontrols2, data=mydata, effects=mymodel, prevAns=ans))
ans$targets
sum(abs(s501-s502)) ## OK Rate
sum(colSums(s501)* rowSums(abs(s501-s502))) # OK inRate
sum(s502) ## OK density
sum(s502*t(s502)) ## OK recip

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, inRateInv, type='rate')
mycontrols0 <- sienaModelCreate(projname=NULL, seed=534, cond=FALSE)
(ans <- siena07(mycontrols0, data=mydata, effects=mymodel))
ans$targets
sum(abs(s501-s502)) ## OK Rate
sum((1/(colSums(s501)+1))* rowSums(abs(s501-s502))) # OK inRateInv

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, inRateLog, type='rate')
(ans <- siena07(mycontrols0, data=mydata, effects=mymodel))
ans$targets
sum(abs(s501-s502)) ## OK Rate
sum((log(colSums(s501)+1))* rowSums(abs(s501-s502))) # OK inRateLog


################################################################################
### check absOutDiffIntn
################################################################################

mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet1, mynet2)
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=1234)

mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, absOutDiffIntn, name="mynet2", interaction1="mynet1", parameter=1)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
sum(rowSums(s503)) # OK outdegree
sum(rowSums(s503 * t(s503))) # OK recip
ad <- abs(outer(rowSums(s501),rowSums(s501),"-"))
sum(s503*ad) # OK absOutDiffIntn

mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, absOutDiffIntn, name="mynet2", interaction1="mynet1", parameter=2)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
ad2 <- abs(outer(sqrt(rowSums(s501)),sqrt(rowSums(s501)),"-"))
sum(s503*ad2) # OK absOutDiffIntn

################################################################################
### check recipRateInv and recipRateLog
################################################################################

mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet)
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=534, cond=FALSE, firstg=0.01)
mycontrols2 <- sienaModelCreate(projname=NULL, seed=534, cond=FALSE, firstg=0.02)
mymodel <- getEffects(mydata)
effectsDocumentation(mymodel)
mymodel <- includeEffects(mymodel, recipRate, type='rate')
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
(ans1 <- siena07(mycontrols2, data=mydata, effects=mymodel, prevAns=ans))
ans$targets
sum(abs(s501-s502)) ## OK Rate
sum(colSums(s501*t(s501))* rowSums(abs(s501-s502))) # OK recipRate
sum(s502) ## OK density
sum(s502*t(s502)) ## OK recip

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, recipRateInv, type='rate')
mycontrols0 <- sienaModelCreate(projname=NULL, seed=534, cond=FALSE)
(ans <- siena07(mycontrols0, data=mydata, effects=mymodel))
ans$targets
sum(abs(s501-s502)) ## OK Rate
sum((1/(colSums(s501*t(s501))+1))* rowSums(abs(s501-s502))) # OK recipRateInv

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, recipRateLog, type='rate')
(ans <- siena07(mycontrols0, data=mydata, effects=mymodel))
(ans <- siena07(mycontrols0, data=mydata, effects=mymodel, prevAns=ans))
ans$targets
sum(abs(s501-s502)) ## OK Rate
sum((log(colSums(s501*t(s501))+1))* rowSums(abs(s501-s502))) # OK recipRateLog

##########################################################################
### check  SimAllNear,SimAllFar
##########################################################################

mynet <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
z1 <- 1+.29*(1:50)
set.seed(143)
z2 <- z1 + rnorm(50) +0.2
z1 <- trunc(z1)
z2 <- trunc(z2)
table(z2-z1)
trunc(cbind(z1,z2))
z2 <- pmax(1,pmin(z2,15))
table(z2,z1)
mybeh <- sienaDependent(cbind(z1,z2), type="behavior") # range 1:15
(mydata <- sienaDataCreate(mynet, mybeh))
mymodel <- getEffects(mydata)
effectsDocumentation(mymodel)
# for SimAllNear effect:
(mymodel <- setEffect(mymodel,simAllNear,name='mybeh', parameter=2))
mymodel
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=514)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets

(mbh <- mean(mybeh))
sum(abs(z2-z1)) # OK rate
sum((z2-mbh)) # OK linear
sum((z2-mbh)^2) # OK quad

p <- 2
aot <- abs(outer(z2,z2,"-"))
NN <- 1*(aot <= p)
diag(NN) <- 0
sum(NN*(p-aot)) # OK simAllNear

# for SimAllFar effect:
mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,simAllFar,name='mybeh', parameter=4))
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=514, firstg=0.05, diagonalize=0.5)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets

p <- 4
aot <- abs(outer(z2,z2,"-"))
FF <- 1*(aot >= p)
diag(FF) <- 0
sum(FF*(p-aot)) # OK simAllFar

################################################################################
### check avDegIntn
################################################################################

mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- sienaDependent(array(c(s503, s502, s501), dim=c(50, 50, 3)))
(mydata <- sienaDataCreate(mynet1, mynet2))

mymodel <- getEffects(mydata)
(mymodel <- includeEffects(mymodel, avDegIntn, name='mynet2',
                interaction1='mynet1'))
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=514)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel, prevAns=ans))
# convergence is very slow, which is natural for this effect with only 3 waves
ans$targets
################################################################################
### check avDeg
################################################################################

mynet <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
(mydata <- sienaDataCreate(mynet))

mymodel <- getEffects(mydata)
(mymodel <- includeEffects(mymodel, avDeg))
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=514, diagonalize=0.6)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel, prevAns=ans))
# convergence is very slow, which is natural for this effect with only 3 waves
ans$targets

(ad2 <- sum(s502)/50)
(ad3 <- sum(s503)/50)
p <- 2
sum((ad2-p)*s502 + (ad3-p)*s503) # OK avDeg

sum(s502 + s501) # OK density
sum(s502*t(s502) + s501*t(s501)) # OK recip
(ad1 <- sum(s501)/50)
(ad2 <- sum(s502)/50)
p <- 2
sum((ad1-p)*s502 + (ad2-p)*s501) # OK avDegIntn

################################################################################
### check avDeg
################################################################################

mynet <- sienaDependent(array(c(s501, s502, s503, s502, s501), dim=c(50, 50, 5)))
(mydata <- sienaDataCreate(mynet))

mymodel <- getEffects(mydata)
(mymodel <- includeEffects(mymodel, avDeg))
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=514, diagonalize=0.6, nsub=5)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel, prevAns=ans))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel, prevAns=ans))
# convergence is very slow, which is natural for this effect with only 5 waves
ans$targets

(ad2 <- sum(s502)/50)
(ad3 <- sum(s503)/50)
(ad1 <- sum(s501)/50)
p <- 2
sum(2*(ad2-p)*s502 + (ad3-p)*s503 + (ad1-p)*s501) # OK avDeg

mymodel <- setEffect(mymodel, density, initialValue=ans$theta[1], fix=TRUE, test=TRUE)
mymodel <- setEffect(mymodel, recip, initialValue=ans$theta[2], fix=TRUE, test=TRUE)
(ans1 <- siena07(mycontrols, data=mydata, effects=mymodel))
mymodel$fix <- FALSE
mymodel$test <- FALSE
(ans2 <- siena07(mycontrols, data=mydata, effects=mymodel, prevAns=ans1))
(ans3 <- siena07(mycontrols, data=mydata, effects=mymodel, prevAns=ans2))

################################################################################
### check sharedTo
################################################################################

mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet1, mynet2)
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=1234)
mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, sharedTo, name="mynet2", interaction1="mynet1", parameter=1)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 115 116  70 106 122  90 154
intwostars1 <- (s501) %*% t(s501)
intwostars3 <- (s503) %*% t(s503)
diag(intwostars1) <- 0
diag(intwostars3) <- 0
sum(intwostars1 * intwostars3) # OK

###
mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, sharedTo, name="mynet2", interaction1="mynet1", parameter=2)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 115 116  70 106 122  90 139.7484
sum(sqrt(intwostars1) * intwostars3) # OK

###
mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, sharedTo, name="mynet2", interaction1="mynet1", parameter=3)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 115 116  70 106 122  90 139.7484
n <- 50
intwostars2 <- (s502) %*% t(s502)
diag(intwostars2) <- 0
(cc <- (n/(2*(n-1)))*mean(intwostars1 + intwostars2))
sum((intwostars1-cc) * intwostars3) # OK

###
mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, sharedTo, name="mynet2", interaction1="mynet1", parameter=4)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 115 116  70 106 122  90 44.01341
sum((sqrt(intwostars1)-sqrt(cc)) * intwostars3) # OK

################################################################################
### check totInAltW, avInAltW, totWInAlt, avWInAlt
################################################################################

beh2 <- s50a[,2] - mean(s50a[,1:2])
set.seed(123)
dcova <- matrix(runif(2500),50,50)
s3 <- dcova
diag(s3) <- NA
dcov_c <- dcova - mean(s3,na.rm = TRUE)

mynet <- sienaDependent(array(c(s501, s502), dim = c(50, 50, 2)))
beh <- sienaDependent(s50a[,1:2],type = 'behavior')
dcov <- coDyadCovar(dcova)
mydata <- sienaDataCreate(mynet,dcov,beh)
mycontrols <- sienaAlgorithmCreate(projname = NULL, seed = 1234)

# avInAltW p1
mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel,avInAltW, name = 'beh',
                   interaction1 = 'mynet', interaction2 = 'dcov',
                   parameter = 1)
(ans1 <- siena07(mycontrols, data = mydata, effects = mymodel))
ans1$targets

(mbh <- mean(beh))
divi <- function(a,b){ifelse(b==0, 0, a/b)}
mdc <- attr(mydata$dycCovars[[1]],'mean')
avinwalt <-  divi( ((dcova - mdc)*t(s501)) %*% (beh[,,2] - mbh), colSums(s501))
sum( (beh[,,2] - mbh) * avinwalt ) # avInAltW OK

# avInAltW p2
mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel,avInAltW, name = 'beh',
                   interaction1 = 'mynet', interaction2 = 'dcov',
                   parameter = 2)
(ans1 <- siena07(mycontrols, data = mydata, effects = mymodel))
ans1$targets
avinwalt <-  divi( ((dcova - mdc)*t(s501)) %*% (beh[,,2] - mbh),
					rowSums(t(s501)*(dcova-mdc)))
sum( (beh[,,2] - mbh) * avinwalt ) # avInAltW OK

# avWInAlt
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,avWInAlt, name = 'beh',
                        interaction1 = 'mynet', interaction2 = 'dcov')
(ans1 <- siena07(mycontrols, data = mydata, effects = mymodel))
ans1$targets
sum(divi(rowSums(t(s501) * dcov_c), colSums(s501)) * beh2) # avWInAlt OK

# totInAltW
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,totInAltW, name = 'beh',
                        interaction1 = 'mynet', interaction2 = 'dcov')
(ans1 <- siena07(mycontrols, data = mydata, effects = mymodel))
ans1$targets
sum(rowSums(t(s501) * dcov_c %*% diag(beh2)) * beh2) # totInAltW OK

# totWInAlt
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,totWInAlt, name = 'beh',
                        interaction1 = 'mynet', interaction2 = 'dcov')
(ans1 <- siena07(mycontrols, data = mydata, effects = mymodel))
ans1$targets
sum(rowSums(t(s501) * dcov_c) * beh2) # totWInAlt OK


################################################################################
### check toAny
################################################################################

advice <- sienaDependent(array(c(s502, s501), dim=c(50, 50, 2)))
trust <- sienaDependent(array(c(s503, s501), dim=c(50, 50, 2)))
# dyadic covariate
suppressWarnings(mat <- matrix(c(0,1,0,0,4,0,2,0,0), 50,50))
dcov <- coDyadCovar(mat, centered=FALSE)
mydata <- sienaDataCreate(advice,trust,dcov)

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,toAny,name="trust",interaction1="advice")
mymodel
mycontrols <- sienaAlgorithmCreate(projname = NULL, seed=123)
(myans <- siena07(mycontrols, data = mydata, effects = mymodel))
myans$targets

# for toAny effect:
# W=advice , X=trust
XX <- s501 %*% t(s501)
diag(XX) <- 0
tXX <- 1*(XX >= 1)
sum(tXX * s502) # 56 OK


################################################################################
### check outAct_ego, reciAct_ego, inPop_dya
################################################################################

mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet)
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=534)
mymodel <- getEffects(mydata)
effectsDocumentation(mymodel)
mymodel <- includeEffects(mymodel, outAct)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
sum(rowSums(s502)^2) # OK
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, outAct_ego)
(ans1 <- siena07(mycontrols, data=mydata, effects=mymodel))
ans1$targets
(xbar <- (mean(rowSums(s501)) + mean(rowSums(s502)))/2)
sum(rowSums(s502)*(rowSums(s502)-xbar)) # OK
mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, outAct_ego, parameter=2)
(ans2 <- siena07(mycontrols, data=mydata, effects=mymodel))
ans2$targets
sum(rowSums(s502)*(sqrt(rowSums(s502))-sqrt(xbar))) # OK

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, reciAct_ego)
(ans1 <- siena07(mycontrols, data=mydata, effects=mymodel))
ans1$targets
(xbar <- (mean(rowSums(s501*t(s501))) + mean(rowSums(s502*t(s502))))/2)
sum(rowSums(s502)*(rowSums(s502*t(s502))-xbar)) # OK
mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, reciAct_ego, parameter=2)
(ans2 <- siena07(mycontrols, data=mydata, effects=mymodel))
ans2$targets
sum(rowSums(s502)*(sqrt(rowSums(s502*t(s502)))-sqrt(xbar))) # OK

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, inPop_dya)
(ans1 <- siena07(mycontrols, data=mydata, effects=mymodel))
ans1$targets
(xbar <- (mean(colSums(s501)) + mean(colSums(s502)))/2)
sum(colSums(s502)*(colSums(s502)-xbar)) # OK
mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, inPop_dya, parameter=2)
(ans2 <- siena07(mycontrols, data=mydata, effects=mymodel))
ans2$targets
sum(colSums(s502)*(sqrt(colSums(s502))-sqrt(xbar))) # OK



################################################################################
### check sameXInPop, diffXInPop
################################################################################

mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))

# construct actor covariate
(five <- rep(1:5,10))
five <- coCovar(five, centered=FALSE)

mydata <- sienaDataCreate(mynet, five)
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,inPop)
mymodel <- includeEffects(mymodel,sameXInPop,interaction1='five')
mymodel
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=1234)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# check target statistics: indegree popularity
sum((colSums(s502))^2) # 386 OK
# check target statistics: sameXInPop
samef <- 1*outer(five,five,FUN="==")
sum(diag(t(s502) %*% samef %*% s502))
sum(s502 * (samef %*% s502))

mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel,sameXInPop,interaction1='five', parameter=2)
mymodel
(ans1 <- siena07(mycontrols, data=mydata, effects=mymodel))
ans1$targets
sum(s502 * sqrt(samef %*% s502)) # OK


mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,diffXInPop,interaction1='five')
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
difff <- 1*outer(five,five,FUN="!=")
sum(s502 * (difff %*% s502)) # OK

mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel,diffXInPop,interaction1='five', parameter=2)
mymodel
(ans1 <- siena07(mycontrols, data=mydata, effects=mymodel))
ans1$targets
sum(s502 * sqrt(difff %*% s502)) # OK


################################################################################
### check altInDist2, totInDist2
################################################################################

mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaDependent(array(c(s503, s502), dim=c(50, 50, 2)))
# construct actor covariate
in1 <- colSums(s501)
out1 <- rowSums(s501)
center <- in1 + out1
table(center)
center <- round(center/5)
table(center)

central <- coCovar(center)
mydata <- sienaDataCreate(mynet1, mynet2, central)

mymodel <- getEffects(mydata)
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=1234, nsub=2, n3=200)
mymodel <- includeEffects(mymodel,altInDist2,name='mynet2',
                 interaction1='central')

ans <- siena07(mycontrols, data=mydata, effects=mymodel)
ans$targets
# 115.00 116.00  70.00 106.00 116.00  70.00  11.64

sum(mynet1[,,2]) # OK
sum(mynet1[,,2]*t(mynet1[,,2])) # OK
sum(mynet2[,,2]) # OK
sum(mynet2[,,2]*t(mynet2[,,2])) # OK
ind <- colSums(mynet2[,,2])

# for altInDist2
cova <- central - mean(central)
divi <- function(x,y){ifelse(y==0, 0, x/y)}
vv <- matrix(NA, length(central),length(central))
for (i in seq_along(central)) {
  for(j in seq_along(central)){vv[i,j] <-
      divi((sum(mynet2[,j,2]*cova) - mynet2[i,j,2]*cova[i]),(ind[j] - mynet2[i,j,2]))}}

sum(mynet2[,,2]*vv) # OK

# for totInDist2
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,totInDist2,name='mynet2',
                 interaction1='central')
ans <- siena07(mycontrols, data=mydata, effects=mymodel)
ans$targets
# [1] 115.0 116.0  70.0 106.0 116.0  70.0  35.6
for (i in seq_along(central)) {
  for(j in seq_along(central)){vv[i,j] <-
      (sum(mynet2[,j,2]*cova) - mynet2[i,j,2]*cova[i])}}

sum(mynet2[,,2]*vv) # OK

################################################################################

# Now for two-mode

# Create fake network data
set.seed(12321)
wave1 <- matrix(0,12,10)
for(i in 1:120){
  wave1[i] <- sample(c(0,1),1, prob=c(0.7,0.3))
}
wave2 <- wave1
for (i in 1:120){
  wave2[i] <- abs(wave1[i] - 0.5 + sample(c(-0.5,0.5),1, prob=c(0.3,0.7)))
}

# Create fake ego and alter covariate data
covego1 <- sample(0:8, 12, replace=T)
covego2 <- sample(0:8, 12, replace=T)
covalt <- sample(0:8, 10, replace=T)
covego3 <- c(1,1,1,1,2,2,2,2,3,3,3,3)

# Identify nodesets
senders <- sienaNodeSet(12, nodeSetName="senders")
recipients <- sienaNodeSet(10, nodeSetName="recipients")

# Make dependent networks and behaviour
network <- sienaDependent(array(c(wave1, wave2), dim=c(12,10,2)),
                          type="bipartite", nodeSet=c("senders","recipients"),
						  allowOnly=FALSE)

# Make covariates
covaralt <- coCovar(covalt, nodeSet="recipients")
covarego <- coCovar(covego1, nodeSet="senders")
covarego3 <- coCovar(covego3, nodeSet="senders")

# Put it all together
nbdata <- sienaDataCreate(network, covaralt, covarego,
                nodeSets=list(senders,recipients))
nbdata

nbEffects <- getEffects(nbdata)
nbEffects <- setEffect(nbEffects , altInDist2, interaction1="covarego")
(ans <- siena07(mycontrols, data=nbdata, effects=nbEffects))
ans$targets
length(covarego)
ind <- colSums(wave2)
cova <- covarego - mean(covarego)
divi <- function(x,y){ifelse(y==0, 0, x/y)}
vv <- matrix(NA, 12,10)
for (i in 1:12) {
  for(j in 1:10){vv[i,j] <-
      divi((sum(wave2[,j]*cova) - wave2[i,j]*cova[i]),(ind[j] - wave2[i,j]))}}

sum(wave2*vv) # OK

################################################################################
### check altDist2, totDist2
################################################################################

mynet <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
myvar <- coCovar(s50a[,2])
mydata <- sienaDataCreate(mynet, myvar)

myalgorithm <- sienaAlgorithmCreate(projname=NULL, nsub=2, n3=100, seed=12345)
myalgorithm <- sienaAlgorithmCreate(projname=NULL, seed=12345)

# outdist
mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel, altDist2, interaction1="myvar", parameter=1)
(ans <- siena07(myalgorithm, data=mydata, effects=mymodel))
ans$targets

# calculation of target for altDist2:
divi <- function(x,y){ifelse(y==0, 0, x/y)}
mat <- s503
diag(mat) <- 0
table(mat, useNA='always')
cova <- as.double( myvar - mean(myvar))
#mat
totv <- mat %*% cova
avv <- divi(totv,rowSums(mat))
ef <- mat %*% (avv) # effect
sum(ef) # 12.35 OK

################################################################################
### check unequalX, sameWXClosure, diffWXClosure, sameXWClosure, diffCWClosure,
###       sameWWClosure, diffWWClosure
################################################################################

mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaDependent(array(c(s503, s502), dim=c(50, 50, 2)))
# construct actor covariate
in1 <- colSums(s501)
out1 <- rowSums(s501)
center <- in1 + out1
table(center)
center <- round(center/5)
table(center)
eq <- 1*outer(center, center, "==")# matrix I{center[i] = center[j]}
diag(eq) <- 0
table(eq)
sum(center) # 46

central <- coCovar(center)
mydata <- sienaDataCreate(mynet1, mynet2, central)
myalgorithm <- sienaAlgorithmCreate(projname=NULL, nsub=2, n3=100, seed=12345)
myalgorithm <- sienaAlgorithmCreate(projname=NULL, seed=12345)


mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,unequalX,name='mynet1', interaction1='central')
(ans <- siena07(myalgorithm, data=mydata, effects=mymodel))
ans$targets
sum((1-eq)*s502) # 49 OK

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,sameWXClosure,name='mynet2',
                interaction1='mynet1', interaction2='central')
(ans <- siena07(myalgorithm, data=mydata, effects=mymodel))
ans$targets
# sameWXClosure
mtwop <- (s501*eq) %*% s502
diag(mtwop) <- 0
sum(mtwop * s502)
# 36 OK

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,diffWXClosure,name='mynet2',
                interaction1='mynet1', interaction2='central')
(ans <- siena07(myalgorithm, data=mydata, effects=mymodel))
ans$targets
# diffWXClosure
mtwop <- (s501*(1-eq)) %*% s502
diag(mtwop) <- 0
sum(mtwop * s502)
# 24 OK

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,sameXWClosure,name='mynet1',
                interaction1='mynet2', interaction2='central')
(ans <- siena07(myalgorithm, data=mydata, effects=mymodel))
ans$targets
# sameXWClosure
mtwop <- (s502*eq) %*% s503
diag(mtwop) <- 0
sum(mtwop * s502)
# 54 OK

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,diffXWClosure,name='mynet1',
                interaction1='mynet2', interaction2='central')
(ans <- siena07(myalgorithm, data=mydata, effects=mymodel))
ans$targets
# diffXWClosure
mtwop <- (s502*(1-eq)) %*% s503
diag(mtwop) <- 0
sum(mtwop * s502)
# 38 OK

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,sameWWClosure,name='mynet2',
                interaction1='mynet1', interaction2='central')
(ans <- siena07(myalgorithm, data=mydata, effects=mymodel))
ans$targets
# sameWWClosure
mtwop <- (s501*eq) %*% s501
diag(mtwop) <- 0
sum(mtwop * s502)
# 36 OK

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel,diffWWClosure,name='mynet2',
                interaction1='mynet1', interaction2='central')
(ans <- siena07(myalgorithm, data=mydata, effects=mymodel))
ans$targets
# diffWWClosure
mtwop <- (s501*(1-eq)) %*% s501
diag(mtwop) <- 0
sum(mtwop * s502)
# 20 OK


################################################################################
### check outdeg for continuous behavior and two-mode network
################################################################################

# Create fake network data
set.seed(12321)
wave1 <- matrix(0,30,10)
for(i in 1:300){
  wave1[i] <- sample(c(0,1),1, prob=c(0.7,0.3))
}
wave2 <- wave1
for (i in 1:300){
  wave2[i] <- abs(wave1[i] - 0.5 + sample(c(-0.5,0.5),1, prob=c(0.3,0.7)))
}

# Identify nodesets
senders <- sienaNodeSet(30, nodeSetName="senders")
recipients <- sienaNodeSet(10, nodeSetName="recipients")

# Make dependent network
network <- sienaDependent(array(c(wave1, wave2), dim=c(30,10,2)),
                          type="bipartite", nodeSet=c("senders","recipients"),
						  allowOnly=FALSE)
# Make behavior
w1 <- 1:30
w12 <- 1:30 + runif(30, 0, 0.5) - 0.25
beh <- sienaDependent(cbind(w1,w12), type="continuous", nodeSet="senders")

# Make data set
(nbdata <- sienaDataCreate(network, beh, nodeSets=list(senders,recipients)))

myalgorithm <- sienaAlgorithmCreate(projname=NULL, seed=12345)
nbEffects <- getEffects(nbdata, onePeriodSde=TRUE)
nbEffects <- includeEffects(nbEffects, outdeg, name="beh", interaction1="network")
(ans <- siena07(myalgorithm, data=nbdata, effects=nbEffects))
ans$targets
sum((beh[,1,2] - mean(beh)) * rowSums(wave1))  # This would be it for non-continuous behavior
sum((beh[,1,2] ) * rowSums(wave1)) # 1625.097, OK

# Note: for continuous behavior effects, behavior is not centered.


################################################################################
### check sameXInPop, diffXInPop, sameXOutAct, diffXOutAct parameter 1 and 2
################################################################################


mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
sets <- 1*((1:50) < 22)
binary <- coCovar(sets)
(mydata <- sienaDataCreate(mynet, binary))
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=844)

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, sameXInPop, interaction1="binary")
mymodel <- includeEffects(mymodel, inPop)
mymodel
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# check target statistics: indegree popularity
sum((colSums(s502))^2) # 386 OK
# check target statistics: sameXInPop
mat <- outer(1:50,1:50,function(i,j){sets[i]==sets[j]})
mat1 <- t(s502) %*% mat %*% s502
sum(diag((mat1))) # 276 OK
sum((rowSums(s502*mat))^2)


mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel, sameXInPop, interaction1="binary", parameter=2))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# check target statistics: sameXInPop parameter=2
mat <- outer(1:50,1:50,function(i,j){sets[i]==sets[j]})
mat1 <- sqrt(mat %*% s502)
sum(s502*mat1) # 174.3673 OK

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel, diffXInPop, interaction1="binary", parameter=2))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# check target statistics: diffXInPop parameter=2
mat <- outer(1:50,1:50,function(i,j){sets[i]!=sets[j]})
mat1 <- sqrt(mat %*% s502)
sum(s502*mat1) # 82.72646   OK


mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, sameXOutAct, interaction1="binary")
(mymodel <- includeEffects(mymodel, outAct))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# check target statistics: outdegree activity
sum(rowSums(s502)) # 116 OK
sum((rowSums(s502))^2) # 350 OK
# check target statistics: sameXOutAct
mat <- outer(1:50,1:50,function(i,j){sets[i]==sets[j]})
diag(mat) <- 0
sum((rowSums(s502 * mat))^2) # 194 OK


mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel, sameXOutAct, interaction1="binary", parameter=2))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# check target statistics: sameXOutAct parameter=2
mat <- outer(1:50,1:50,function(i,j){sets[i]==sets[j]})
diag(mat) <- 0
sum((rowSums(s502 * mat))*sqrt(rowSums(s502 * mat))) # 121.9956  OK

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, diffXOutAct, interaction1="binary")
(mymodel <- includeEffects(mymodel, outAct))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# check target statistics: diffXOutAct
dmat <- outer(1:50,1:50,function(i,j){sets[i]!=sets[j]})
sum((rowSums(s502 * dmat))^2) # 64 OK

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel, diffXOutAct, interaction1="binary", parameter=2))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# check target statistics: diffXOutAct parameter=2
sum((rowSums(s502 * dmat))*sqrt(rowSums(s502 * dmat))) # 46.73059   OK

################################################################################
### check simEgoInDist2 and avInSimDist2
################################################################################


mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[,1:2], type="behavior")
mydata <- sienaDataCreate(mynet, mybeh)

mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=842)

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,simEgoInDist2, interaction1='mybeh'))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# [1] 115.00000 116.00000  70.00000  11.58412  27.00000   5.50000  71.10500

sum(s502) # OK
sum(s502*t(s502)) # OK
sum(mybeh[,,2]-mean(mybeh)) # 5.5 OK
sum((mybeh[,,2]-mean(mybeh))^2) # 71.105  OK

divi <- function(x,y){ifelse(y==0, 0, x/y)}

cova <- mybeh[,1,1] - mean(mybeh)
ind <- colSums(mynet[,,2])
vv <- matrix(NA, length(cova),length(cova))
for (i in seq_along(cova)) {
  for(j in seq_along(cova)){vv[j,i] <-
      divi((sum(mynet[,j,2]*cova) - mynet[i,j,2]*cova[i]),(ind[j] - mynet[i,j,2]))}}
(range.c <- max(mybeh) - min(mybeh))

simi0 <- matrix(NA, length(cova),length(cova))
for (i in seq_along(cova)) {
  for(j in seq_along(cova)){simi0[i,j] <- 1 - abs(cova[i] - vv[j,i])/range.c}}
simi <- simi0 - attr(mydata$depvars$mybeh, 'simMean')
sum(mynet[,,2]*simi) #  OK simEgoInDist2

# Note: if ind[j] - mynet[i,j,2] = 0, vv[j,i] will be 0,
# and since cova is centered, this is the mean, as stated in the definition of the effect.

(mymodel <- setEffect(mymodel,avInSimDist2, name='mybeh', interaction1='mynet'))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 115.000000 116.000000  70.000000  11.584116  27.000000   5.500000  71.105000   3.675713

depva <- mybeh[,1,2] - mean(mybeh)
ind <- colSums(mynet[,,1])
outd <- rowSums(mynet[,,1])
vv <- matrix(NA, length(depva),length(depva))
for (i in seq_along(depva)) {
  for(j in seq_along(depva)){vv[j,i] <-
      divi((sum(mynet[,j,1]*depva) - mynet[i,j,1]*depva[i]),(ind[j] - mynet[i,j,1]))}}
(range.c <- max(mybeh) - min(mybeh))

(sme <- attr(mydata$depvars$mybeh, 'simMean'))
simi0 <- matrix(NA, length(depva),length(depva))
for (i in seq_along(depva)) {
  for(j in seq_along(depva)){simi0[i,j] <-
		ifelse((ind[j] - mynet[i,j,1])>0, 1 - (abs(depva[i] - vv[j,i])/range.c) - sme, 0)}}
sum(divi(rowSums(mynet[,,1]*simi0), outd)) # OK avInSimDist2

# The treatment of cases with ind[j] - mynet[i,j,2] = 0 is a bit different:
# for simEgoInDist2 vv[j,i] is the mean, i.e., in this case, 0;
# for avInSimDist2 simi0[i,j] is 0, i.e., the mean similarity.

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,totInSimDist2, name='mybeh', interaction1='mynet'))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# [1] 115.000000 116.000000  70.000000  27.000000   5.500000  71.105000   9.769558

depva <- mybeh[,1,2] - mean(mybeh)
ind <- colSums(mynet[,,1])
outd <- rowSums(mynet[,,1])
vv <- matrix(NA, length(depva),length(depva))
for (i in seq_along(depva)) {
  for(j in seq_along(depva)){vv[j,i] <-
      divi((sum(mynet[,j,1]*depva) - mynet[i,j,1]*depva[i]),(ind[j] - mynet[i,j,1]))}}
(range.c <- max(mybeh) - min(mybeh))

(sme <- attr(mydata$depvars$mybeh, 'simMean'))
simi0 <- matrix(NA, length(depva),length(depva))
for (i in seq_along(depva)) {
  for(j in seq_along(depva)){simi0[i,j] <-
		ifelse((ind[j] - mynet[i,j,1])>0, 1 - (abs(depva[i] - vv[j,i])/range.c) - sme, 0)}}
sum(mynet[,,1]*simi0) # OK totInSimDist2


################################################################################
### check divOutEgoIntn and interaction outActIntnX * divOutEgoIntn
### also divInEgoIntn, divOutAltIntn, divInAltIntn.
################################################################################

mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
alc <- coCovar(s50a[,2])
mydata <- sienaDataCreate(mynet1, mynet2, alc)

mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=842)

divi <- function(x,y){ifelse(y==0, 0, x/y)}

# divOutEgoIntn:

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,divOutEgoIntn, name='mynet2', interaction1='mynet1', parameter=1))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 115.00000 116.00000  70.00000 106.00000 122.00000  90.00000  54.53333
sum(divi(rowSums(s503), rowSums(s501))) # divOutEgoIntn OK

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,divOutEgoIntn, name='mynet2', interaction1='mynet1', parameter=2))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 115.00000 116.00000  70.00000 106.00000 122.00000  90.00000  78.98977
sum(divi(rowSums(s503), sqrt(rowSums(s501)))) # divOutEgoIntn OK

# and now the interaction:

mymodel <- getEffects(mydata)
mymodel <- setEffect(mymodel,outActIntnX, name='mynet2', interaction1='mynet1',
					interaction2='alc', parameter=1, include=FALSE)
mymodel <- setEffect(mymodel,divOutEgoIntn, name='mynet2', interaction1='mynet1', parameter=1, include=FALSE)
(mymodel <- includeInteraction(mymodel,outActIntnX,divOutEgoIntn, name='mynet2',
					interaction1=c('mynet1', 'mynet1'), interaction2=c('alc', '')))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
#  115.0 116.0  70.0 106.0 122.0  90.0  30.5

alc.c <- as.vector(alc) - mean(alc)
sumalc.c <- apply(s501, 1, function(x){sum(x*alc.c)})
sum(rowSums(s503) * divi(sumalc.c, rowSums(s501))) # outActIntnX * divOutEgoIntn parameter 1 OK

mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
alc <- coCovar(s50a[,2]-1, centered=FALSE)
mydata <- sienaDataCreate(mynet1, mynet2, alc)
mymodel <- setEffect(mymodel,outActIntnX, name='mynet2', interaction1='mynet1',
					interaction2='alc', parameter=2, include=FALSE)
mymodel <- setEffect(mymodel,divOutEgoIntn, name='mynet2', interaction1='mynet1', parameter=2, include=FALSE)
(mymodel <- includeInteraction(mymodel,outActIntnX,divOutEgoIntn, name='mynet2',
					interaction1=c('mynet1', 'mynet1'), interaction2=c('alc', '')))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targetss
#  115.0000 116.0000  70.0000 106.0000 122.0000  90.0000 179.2528

alc.v <- as.vector(alc)
sumalc.v <- apply(s501, 1, function(x){sum(x*alc.v)})
sum(rowSums(s503) * sqrt(divi(sumalc.v, rowSums(s501)))) # outActIntnX * divOutEgoIntn parameter 2 OK

# divInEgoIntn:

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,divInEgoIntn, name='mynet2', interaction1='mynet1', parameter=1))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 115.00000 116.00000  70.00000 106.00000 122.00000  90.00000  59.625
sum(divi(rowSums(s503), colSums(s501))) # divInEgoIntn OK

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,divInEgoIntn, name='mynet2', interaction1='mynet1', parameter=2))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 115.00000 116.00000  70.00000 106.00000 122.00000  81.40552
sum(divi(rowSums(s503), sqrt(colSums(s501)))) # divInEgoIntn OK

# divOutAltIntn:

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,divOutAltIntn, name='mynet2', interaction1='mynet1', parameter=1))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 115.00000 116.00000  70.00000 106.00000 122.00000  90.00000  52.76667
sum(divi(colSums(s503), rowSums(s501))) # divOutAltIntn OK

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,divOutAltIntn, name='mynet2', interaction1='mynet1', parameter=2))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 115.00000 116.00000  70.00000 106.00000  76.49343
sum(divi(colSums(s503), sqrt(rowSums(s501)))) # divOutAltIntn OK

# divInAltIntn:

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,divInAltIntn, name='mynet2', interaction1='mynet1', parameter=1))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 115.00000 116.00000  70.00000 106.00000 122.00000  90.00000  59.13333
sum(divi(colSums(s503), colSums(s501))) # divInAltIntn OK

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,divInAltIntn, name='mynet2', interaction1='mynet1', parameter=2))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 115.00000 116.00000  70.00000 106.00000  80.37862
sum(divi(colSums(s503), sqrt(colSums(s501)))) # divInAltIntn OK


################################################################################
### check sameEgoDist2 and sameEgoInDist2
################################################################################


mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mycova <- coCovar(s50a[,1])
mydata <- sienaDataCreate(mynet, mycova)

mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=842)

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,sameEgoInDist2, interaction1='mycova'))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 116.0  70.0  36.6

sum(s502) # OK
sum(s502*t(s502)) # OK
eqalc <- outer(mycova, mycova, "==")
sum(diag(s502)) # 0
divi <- function(x,y){ifelse(y==0, 0, x/y)}

a <- divi((eqalc %*% s502 - s502), (matrix(colSums(s502), 50, 50, byrow=TRUE) - s502))
sum(a * s502) # OK sameEgoInDist2


mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,sameEgoInDist2, interaction1='mycova', parameter=0))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
#  116  70  64

eqalc <- outer(mycova, mycova, "==")
a <- (eqalc %*% s502) - s502
sum(s502*pmin(a,1)) # OK


################################################################################
### check threshold, threshold2, threshold3, threshold4
################################################################################

(mydata <- sienaDataCreate(
	mybeh = sienaDependent(s50a[,1:2], type="behavior")))

mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, linear, quad, include=FALSE)
mymodel <- setEffect(mymodel, threshold, parameter=2)
mymodel <- setEffect(mymodel, threshold2, parameter=3)
mymodel <- setEffect(mymodel, threshold3, parameter=4)
(mymodel <- setEffect(mymodel, threshold4, parameter=5))

mycontrols <- sienaModelCreate(projname=NULL, seed=842)
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets
# 47 31 19 8

vapply(2:5, function(t){return(table(s50a[, 2] >= t))}, 1:2)[2, ]
# gives same values; OK

################################################################################
### check WWX, OutWWX
################################################################################

mynet <- sienaNet(array(c(s502, s503), dim=c(50, 50, 2)))
firstnet <- coDyadCovar(s501) # + t(s501)) # to get not just a 0-1 covariate
mydata <- sienaDataCreate(mynet, firstnet)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,WWX,interaction1='firstnet')
myeff <- includeEffects(myeff,OutWWX,interaction1='firstnet')
myeff <- includeEffects(myeff,X,interaction1='firstnet')
myeff
mymodel <- sienaModelCreate(projname=NULL, seed=842)
ans <- siena07(mymodel, data=mydata, effects=myeff)
ans
ans$targets

m1 <- sum(firstnet/(50*49))
sum(s503*(s501-m1)) # X OK

# Note: for WWX etc., no centering is applied internally.

tp <- s501 %*% s501
diag(tp) <- 0
sum(tp*s503) # WWX OK

tp <- s501 %*% t(s501)
diag(tp) <- 0
sum(tp*s503) # OutWWX OK

# Now with a dyadic covariate that has values outside {0,1}

firstnet <- coDyadCovar(s501 + t(s501))
mydata <- sienaDataCreate(mynet, firstnet)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,WWX,interaction1='firstnet')
myeff <- includeEffects(myeff,X,interaction1='firstnet')
myeff
ans <- siena07(mymodel, data=mydata, effects=myeff)
ans
ans$targets

m1 <- sum(firstnet/(50*49))
sum(s503*(firstnet-m1)) # X OK
sum(diag(firstnet))
tp <- firstnet %*% firstnet
diag(tp) <- 0
sum(tp*s503) # WWX OK

################################################################################
### check outXMore
################################################################################

mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mycova <- coCovar(s50a[,1])
mydata <- sienaDataCreate(mynet, mycova)
poscov <- (mydata$cCovars$mycova > 0)
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=138)

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,outXMore, interaction1='mycova', parameter=2))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets

poscovdeg <- apply(s502, 1, function(x){sum(x[poscov]) })
sum(pmax(poscovdeg-2,0)) # 15 OK


################################################################################
### check outMore
################################################################################

mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet)
poscov <- (mydata$cCovars$mycova > 0)
mycontrols <- sienaAlgorithmCreate(projname=NULL, seed=138)

rowSums(s502)

mymodel <- getEffects(mydata)
(mymodel <- setEffect(mymodel,outMore, parameter=3))
(ans <- siena07(mycontrols, data=mydata, effects=mymodel))
ans$targets

sum(pmax(rowSums(s502)-3,0)) # 11 OK
