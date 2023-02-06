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

library(RSiena)

################################################################################
### check from.w.ind
################################################################################

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

################################################################################
### check transtripX
################################################################################

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

################################################################################
### check homXTransTrip, homXTransRecTrip
################################################################################

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

################################################################################
### check to, toU
################################################################################

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
myeff <- setEffect(myeff,to,name="trust",interaction1="advice", parameter=3)
myeff
mymodel <- sienaModelCreate(projname = NULL, seed=123)
(myans <- siena07(mymodel, data = mydata, effects = myeff))
myans$targets

# for to effect, parameter=3:
# W=advice , X=trust
WX <- s501 %*% s501
diag(WX) <- 0
tWXX <- 1*(WX >= 1) * s501
sum(diag(tWXX)) # 0 OK
sum(tWXX) # 73 OK

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

################################################################################
### check toBack, MixedInXW
################################################################################

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

################################################################################
### check cl.XWX, cl.XWX1, cl.XWX2
################################################################################

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
(myans <- siena07(mymodel, data = mydata, effects = myeff))
myans$targets
sum(XWX) # 86 OK

# for cl.XWX2 effect:
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, cl.XWX2, name="trust", interaction1="advice")
(myans <- siena07(mymodel, data = mydata, effects = myeff))
myans$targets
sum(XWX) # 86 OK

################################################################################
### check XWX
################################################################################

a1 <- as.matrix(read.table("aadv_t1.txt",header=F))
a2 <- as.matrix(read.table("aadv_t2.txt",header=F))
t1 <- as.matrix(read.table("ttr_t1.txt",header=F))
t2 <- as.matrix(read.table("ttr_t2.txt",header=F))
hi <- as.matrix(read.table("hier.txt"))

#a1 <- a1[1:20,1:20]
#a2 <- a2[1:20,1:20]
#t1 <- t1[1:20,1:20]
#t2 <- t2[1:20,1:20]
a1 <- a1[30:49,30:49]
a2 <- a2[30:49,30:49]
t1 <- t1[30:49,30:49]
t2 <- t2[30:49,30:49]
a2[8:12,10:15] <- a1[8:12,10:15]
t2[8:12,10:15] <- t1[8:12,10:15]
hi <- hi[30:49,30:49]

behav <- cbind(rowSums(a1) + colSums(t1), rowSums(a2) + colSums(t2))
(behav <- behav %/% 3 + 1)
behav[behav >= 3] <- 3
table(behav)
behav

#### Defining RSiena variables
)

library(RSienaTest)
advice <- sienaNet(array(c(a1,a2), dim=c(20,20,2)))
trust  <- sienaNet(array(c(t1,t2), dim=c(20,20,2)))
trust1 <- coDyadCovar(t1)
hierarchy <- coDyadCovar(hi)
beh <- sienaDependent(behav, type="behavior")

# We just give the two networks to sienaDataCreate.
# They have the same dimensions; this is necessary,
# as they have the same nodeSet.

mydata <- sienaDataCreate(advice,trust1)

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, X, XWX, interaction1='trust1')
myeff

mymodel <- sienaModelCreate(projname="XWX", seed=1234)
ans <- siena07(mymodel, data=mydata, effects=myeff)
ans
ans$targets
diag(a1) <- 0
diag(a2) <- 0
sum(abs(a1-a2))
sum(a2)
sum(a2*t(a2))
diag(t1) <- 0
mt1 <- sum(t1)/(20*19)
sum(a2*(t1-mt1))
XW <-
XW <- a2 %*% t1
diag(XW) <- 0
XWX <- XW * a2
diag(XWX)
sum(XWX) # OK
# In other words: XWX uses uncentered W.

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


################################################################################
### check avSim, totSim, avInSim, totInSim, avInSimPopAlt, totInSimPopAlt
################################################################################

mynet <- sienaNet(array(c(s502, s503), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[,2:3], type="behavior")
mydata <- sienaDataCreate(mynet, mybeh)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,avSim,name='mybeh',interaction1='mynet')
myeff
mymodel <- sienaModelCreate(projname=NULL)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
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
myeff <- getEffects(mydata)
# It turns out that for this data set, to get an estimate for totSim
# you first have to include the outdeg effect.
myeff <- includeEffects(myeff,outdeg,name='mybeh',interaction1='mynet')
(ans <- siena07(mymodel, data=mydata, effects=myeff))
myeff <- includeEffects(myeff,totSim,name='mybeh',interaction1='mynet')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff, prevAns=ans))
ans$targets
sum(s502 * simi2 ) # OK totSim
# Now try to get a model including totSim but without outdeg
myeff <- includeEffects(myeff,outdeg,name='mybeh',interaction1='mynet', include=FALSE)
(ans <- siena07(mymodel, data=mydata, effects=myeff, prevAns=ans))
myeff <- updateTheta(myeff, ans)
(ans <- siena07(mymodel, data=mydata, effects=myeff))

# for avSimPopAlt effect:
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,avSimPopAlt,name='mybeh',interaction1='mynet')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets

indeg <- colSums(s502)
sum(diag(divi(1,outdeg)) %*% (s502 * simi2) %*% diag(indeg)) # avSimPopAlt OK

# for avSimPopEgo effect:
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,avSimPopEgo,name='mybeh',interaction1='mynet')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
indeg <- colSums(s502)
sum(diag(divi(indeg,outdeg)) %*% s502 * simi2 ) #  OK avSimPopEgo

# for totSimPopAlt effect:
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,totSimPopAlt,name='mybeh',interaction1='mynet')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
indeg <- colSums(s502)
sum((s502 * simi2) %*% diag(indeg)) # avSimPopAlt OK

# for avInSim effect:
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,avInSim,name='mybeh',interaction1='mynet')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets

indeg <- colSums(s502)
indegs <- outer(indeg, indeg, function(x,y){divi(1,x)})
diag(indegs) <- 0
sum(t(s502) * simi2*indegs) # OK avInSim
# alternative:
sum(diag(divi(1,indeg)) %*% t(s502) * simi2 ) # OK avInSim

# for avInSimPopAlt effect:
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,avInSimPopAlt,name='mybeh',interaction1='mynet')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(diag(divi(1,indeg)) %*% t(s502) * simi2  %*% diag(indeg)) # OK avInSimPopAlt

# for totInSim effect:
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,totInSim,name='mybeh',interaction1='mynet')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(t(s502) * simi2 ) # OK totInSim

# for totInSimPopAlt effect:
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,totInSimPopAlt,name='mybeh',interaction1='mynet')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum( t(s502) * simi2  %*% diag(indeg)) # OK totInSimPopAlt


################################################################################
### check avAttHigher, avAttLower, totAttHigher, totAttLower
################################################################################

mynet <- sienaNet(array(c(s502, s503), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[,2:3], type="behavior")
mydata <- sienaDataCreate(mynet, mybeh)
myeff <- getEffects(mydata)

# for avAttHigher effect:
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,avAttHigher,name='mybeh',interaction1='mynet')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
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
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,avAttLower,name='mybeh',interaction1='mynet')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(diag(divi(1,outdeg)) %*% s502 * (simim) ) # OK avAttLower

# for totAttHigher effect:
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,totAttHigher,name='mybeh',interaction1='mynet')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(s502 * t(simim) ) # OK totAttHigher

# for totAttLower effect:
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,totAttLower,name='mybeh',interaction1='mynet')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(s502 * (simim) ) # OK totAttLower


################################################################################
### check crprodInActIntn
################################################################################

mynet1 <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaNet(array(c(s503, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet1, mynet2)

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, crprodInActIntn, name='mynet2',
                        interaction1='mynet1')
myeff

# for crprodInActIntn effect (parameter = 2):
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
ctr <- mean(c(colSums(s501), colSums(s502)))
sum(rowSums(s502 * s501) * (sqrt(colSums(s501)) - sqrt(ctr))) # -4.937096 OK

# for crprodInActIntn effect (parameter = 1):
myeff <- setEffect(myeff, crprodInActIntn, name='mynet2',
                   interaction1='mynet1', parameter=1)
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(rowSums(s502 * s501) * (colSums(s501) - ctr)) # -3.53 OK


################################################################################
### check rateX
################################################################################

library(RSiena)

mynet <- sienaNet(array(c(s502, s503), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[,2:3], type="behavior")
mycov <- coCovar(s50s[,2])
mydata <- sienaDataCreate(mynet, mybeh, mycov)
myeff <- getEffects(mydata)
(myeff <- includeEffects(myeff,RateX,type='rate',
							name='mybeh',interaction1='mycov'))
mymodel <- sienaModelCreate(projname=NULL)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(abs(mybeh[,,1]-mybeh[,,2])) # OK rate
mco <- mean(mycov)
sum((mycov-mco)*abs(mybeh[,,1]-mybeh[,,2])) # OK RateX

myeff <- getEffects(mydata)
(myeff <- includeEffects(myeff,RateX,type='rate',
							name='mybeh',interaction1='mybeh'))
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(abs(mybeh[,,1]-mybeh[,,2])) # OK rate
mbh <- mean(mybeh)
sum((mybeh[,,1])*abs(mybeh[,,1]-mybeh[,,2])) # OK RateX non-centered

################################################################################
### check totExposure and other diffusion effects
################################################################################

library(RSiena)

mynet <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
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

mymodel <- sienaModelCreate(projname=NULL, seed=1234, firstg=0.01)
mymodel2 <- sienaModelCreate(projname=NULL, seed=1234, firstg=0.1)
mymodel4 <- sienaModelCreate(projname=NULL, seed=1234, firstg=0.001)

myeff <- getEffects(mydata)
(ans <- siena07(mymodel2, data=mydata, effects=myeff))
myeff0 <- updateTheta(myeff, ans)

myeff.r <- getEffects(mydata.r)
(ans.r <- siena07(mymodel2, data=mydata.r, effects=myeff.r))
myeff0.r <- updateTheta(myeff.r, ans.r)

myeff <- includeEffects(myeff0,avExposure,type='rate',
							name='mybeh',interaction1='mynet')
(ans1 <- siena07(mymodel4, data=mydata, effects=myeff))
ans1$targets
sum(sm2-sm1) # OK rate

divi <- function(x,y){ifelse(y==0, 0, x/y)}
amean <- function(x,z){ifelse((sum(x)>0), sum(x*z)/sum(x), 0)}
eff <- apply(s501,1,function(x){amean(x,sm2)})
eff00 <- apply(s501,1,function(x){amean(x,sm1)})
sum((sm2-sm1)*eff) #  this is not the target statistic for avExposure
sum((sm2-sm1)*eff00) # OK  avExposure

(myeff5 <- setEffect(myeff0,avExposure,type='rate', parameter=2,
							name='mybeh',interaction1='mynet'))
(ans5 <- siena07(mymodel4, data=mydata, effects=myeff5, prevAns=ans1))
ans5$targets
eff3 <- apply(s501,1,function(x){sum(x*sm1)})
eff35 <- eff3
eff35[eff35==1] <- 0
outd <- rowSums(s501)
sum((sm2-sm1)*divi(eff35, outd)) #  OK avExposure p=2

(myeff6.r <- setEffect(myeff.r, avExposure,type='rate', parameter=-2,
							name='mybeh.r',interaction1='mynet'))
(ans6.r <- siena07(mymodel4, data=mydata.r, effects=myeff6.r, prevAns=ans5))
ans6.r$targets
eff3.r <- apply(s501,1,function(x){sum(x*sm1r)})
eff56 <- eff3.r
eff56[eff56==1] <- 0
eff56[eff56>2] <- 2
sum((sm2-sm1r)*divi(eff56, outd)) #   OK avExposure p=-2

(myeff2 <- setEffect(myeff0,totExposure,type='rate',
							name='mybeh',interaction1='mynet'))
(ans2 <- siena07(mymodel4, data=mydata, effects=myeff2))
ans2$targets
sum(sm2-sm1) # OK rate
eff3 <- apply(s501,1,function(x){sum(x*sm1)})
sum((sm2-sm1)*eff3) # OK totExposure

(myeff3 <- setEffect(myeff0,totExposure,type='rate', parameter=2,
							name='mybeh',interaction1='mynet'))
(ans3 <- siena07(mymodel4, data=mydata, effects=myeff3, prevAns=ans2))
ans3$targets
eff3 <- apply(s501,1,function(x){sum(x*sm1)})
eff32 <- eff3
eff32[eff32==1] <- 0
sum((sm2-sm1)*eff32) # OK totExposure p=2

(myeff4.r <- setEffect(myeff0.r,totExposure,type='rate', parameter=-2,
							name='mybeh.r',interaction1='mynet'))
(ans4.r <- siena07(mymodel4, data=mydata.r, effects=myeff4.r, prevAns=ans3))
ans4.r$targets
eff3.r <- apply(s501,1,function(x){sum(x*sm1r)})
eff34.r <- eff3.r
eff34.r[eff34.r==1] <- 0
eff34.r[eff34.r>2] <- 2
sum((sm2-sm1r)*eff34.r) # OK totExposure p=-2

(myeff6 <- setEffect(myeff,infectIn,type='rate',
							name='mybeh',interaction1='mynet'))
(ans6 <- siena07(mymodel4, data=mydata, effects=myeff6))
ans6$targets
ind <- colSums(s501)
eff3 <- apply(s501,1,function(x){sum(x*sm1)})
eff7 <- sapply(1:50, function(i){sum(s501[i,]*ind*sm1)})
sum((sm2-sm1)*eff7) #  OK infectIn

(myeff7 <- setEffect(myeff,infectIn,type='rate', parameter=2,
							name='mybeh',interaction1='mynet'))
(ans7 <- siena07(mymodel4, data=mydata, effects=myeff7, prevAns=ans6))
ans7$targets
eff72 <- eff7
eff72[eff3<2] <- 0
sum((sm2-sm1)*eff72) #  OK infectIn p=2

(myeff8 <- setEffect(myeff0,infectOut,type='rate', parameter=2,
							name='mybeh', interaction1='mynet'))
(ans8 <- siena07(mymodel4, data=mydata, effects=myeff8))
ans8$targets
outd <- rowSums(s501)
eff3 <- apply(s501,1,function(x){sum(x*sm1)})
eff8 <- sapply(1:50, function(i){sum(s501[i,]*outd*sm1)})
eff82 <- eff8
eff82[eff3<2] <- 0
sum((sm2-sm1)*eff82) #  OK infectOut p=2

(myeff9 <- setEffect(myeff0,infectCovar,type='rate',
							name='mybeh', interaction1='mynet',
							interaction2='mycov'))
(ans9 <- siena07(mymodel4, data=mydata, effects=myeff9))
ans9$targets
(mcm <- mean(mycov))
eff9 <- sapply(1:50, function(i){sum(s501[i,]*(mycov-mcm)*sm1)})
sum((sm2-sm1)*eff9) #   infectCovar OK

(myeffA <- setEffect(myeff.r,infectCovar,type='rate', parameter=2,
							name='mybeh.r', interaction1='mynet',
							interaction2='mycov'))
(ansA <- siena07(mymodel4, data=mydata.r, effects=myeffA, prevAns=ans9))
ansA$targets #
eff3 <- apply(s501,1,function(x){sum(x*sm1r)})
(mcm <- mean(mycov))
effA <- sapply(1:50, function(i){sum(s501[i,]*(mycov-mcm)*sm1r)})
effA2 <- effA
effA2[eff3==1] <- 0
sum((sm2-sm1r)*effA2) # OK infectCovar p=2

(myeffB <- setEffect(myeff0.r,susceptAvIn,type='rate',
							name='mybeh.r',interaction1='mynet'))
(ansB <- siena07(mymodel4, data=mydata.r, effects=myeffB))
ansB$targets
eff3 <- apply(s501,1,function(x){sum(x*sm1r)})
outd <- rowSums(s501)
sum((sm2-sm1r)*colSums(s501)*divi(eff3, outd)) #  OK susceptAvIn

(myeffC <- setEffect(myeff0.r, susceptAvIn,type='rate', parameter=2,
							name='mybeh.r',interaction1='mynet'))
(ansC <- siena07(mymodel4, data=mydata.r, effects=myeffC, prevAns=ansB))
ansC$targets
eff35 <- eff3
eff35[eff3==1] <- 0
sum((sm2-sm1r)*colSums(s501)*divi(eff35, outd)) #  OK susceptAvIn p=2

(myeffD <- setEffect(myeff0.r,susceptAvIn,type='rate', parameter=-2,
							name='mybeh.r',interaction1='mynet'))
(ansD <- siena07(mymodel4, data=mydata.r, effects=myeffD, prevAns=ansC))
ansD$targets
eff35[eff3>2] <- 2
outd <- rowSums(s501)
sum((sm2-sm1r)*colSums(s501)*divi(eff35, outd)) #  OK susceptAvIn p=-2

(myeffE <- setEffect(myeff0.r,susceptAvCovar,type='rate', parameter=0,
							name='mybeh.r',interaction1='mynet',
							interaction2='mycov'))
(ansE <- siena07(mymodel4, data=mydata.r, effects=myeffE))
ansE$targets
eff3 <- apply(s501,1,function(x){sum(x*sm1r)})
eff35 <- eff3
(mcm <- mean(mycov))
outd <- rowSums(s501)
sum((sm2-sm1r)*(mycov-mcm)*divi(eff35, outd)) # susceptAvCovar OK

(myeffE1 <- setEffect(myeff0.r,susceptAvCovar,type='rate', parameter=2,
							name='mybeh.r',interaction1='mynet',
							interaction2='mycov'))
(ansE1 <- siena07(mymodel4, data=mydata.r, effects=myeffE1, prevAns=ansE))
ansE1$targets
eff35[eff3==1] <- 0
sum((sm2-sm1r)*(mycov-mcm)*divi(eff35, outd)) # susceptAvCovar p=2 OK.

(myeffF <- setEffect(myeff0.r,susceptAvCovar,type='rate', parameter=-2,
							name='mybeh.r',interaction1='mynet',
							interaction2='mycov'))
(ansF <- siena07(mymodel4, data=mydata.r, effects=myeffF, prevAns=ansE1))
ansF$targets
eff35[eff3>2] <- 2
sum((sm2-sm1r)*(mycov-mcm)*divi(eff35, outd)) # susceptAvCovar p=-2 OK


################################################################################
### check outOutActIntn and outOutAvIntn
################################################################################

mynet1 <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaNet(array(c(s502, s503), dim=c(50, 50, 2)))
# construct actor covariate
in1 <- colSums(s501)
out1 <- rowSums(s501)
center <- in1 + out1
center <- 1*(center >= 5)
central <- coCovar(center)
mydata <- sienaDataCreate(mynet1, mynet2, central)
mymodel <- sienaModelCreate(projname=NULL, seed=1234)

myeff <- getEffects(mydata)
myeff <- setEffect(myeff, outOutActIntn, name="mynet1", interaction1="mynet2", parameter=1)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(rowSums(s502)) # ok outdegree
sum(rowSums(s502 * t(s502))) # OK recip
(avdeg <- mean(rowSums(s502) + rowSums(s501))/2)
sum(rowSums(s502)* (s502 %*% (rowSums(s502) - avdeg))) # OK outOutActIntn

myeff2 <- getEffects(mydata)
myeff2 <- setEffect(myeff2, outOutActIntn, name="mynet1", interaction1="mynet2", parameter=2)
(ans2 <- siena07(mymodel, data=mydata, effects=myeff2))
ans2$targets
sum(rowSums(s502)* (s502 %*% (sqrt(rowSums(s502)) - sqrt(avdeg)))) # OK outOutActIntn p=2

myeff <- getEffects(mydata)
myeff <- setEffect(myeff, outOutAvIntn, name="mynet1", interaction1="mynet2", parameter=1)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
(avdeg <- mean(rowSums(s502) + rowSums(s501))/2)
sum((s502 %*% (rowSums(s502) - avdeg)) ) # OK outOutAvIntn, note that rowSums(s502) cancel

myeff2 <- getEffects(mydata)
myeff2 <- setEffect(myeff2, outOutAvIntn, name="mynet1", interaction1="mynet2", parameter=2)
(ans2 <- siena07(mymodel, data=mydata, effects=myeff2))
ans2$targets
sum( (s502 %*% (sqrt(rowSums(s502)) - sqrt(avdeg)))) # OK outOutAvIntn p=2


################################################################################
### check inRateInv and inRateLog
################################################################################


mynet <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet)
mymodel <- sienaModelCreate(projname=NULL, seed=534, cond=FALSE, firstg=0.01)
mymodel2 <- sienaModelCreate(projname=NULL, seed=534, cond=FALSE, firstg=0.02)
myeff <- getEffects(mydata)
effectsDocumentation(myeff)
myeff <- includeEffects(myeff, inRate, type='rate')
(ans <- siena07(mymodel, data=mydata, effects=myeff))
(ans1 <- siena07(mymodel2, data=mydata, effects=myeff, prevAns=ans))
ans$targets
sum(abs(s501-s502)) ## OK Rate
sum(colSums(s501)* rowSums(abs(s501-s502))) # OK inRate
sum(s502) ## OK density
sum(s502*t(s502)) ## OK recip

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, inRateInv, type='rate')
mymodel0 <- sienaModelCreate(projname=NULL, seed=534, cond=FALSE)
(ans <- siena07(mymodel0, data=mydata, effects=myeff))
ans$targets
sum(abs(s501-s502)) ## OK Rate
sum((1/(colSums(s501)+1))* rowSums(abs(s501-s502))) # OK inRateInv

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, inRateLog, type='rate')
(ans <- siena07(mymodel0, data=mydata, effects=myeff))
ans$targets
sum(abs(s501-s502)) ## OK Rate
sum((log(colSums(s501)+1))* rowSums(abs(s501-s502))) # OK inRateLog


################################################################################
### check absOutDiffIntn
################################################################################

mynet1 <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaNet(array(c(s502, s503), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet1, mynet2)
mymodel <- sienaModelCreate(projname=NULL, seed=1234)

myeff <- getEffects(mydata)
myeff <- setEffect(myeff, absOutDiffIntn, name="mynet2", interaction1="mynet1", parameter=1)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(rowSums(s503)) # OK outdegree
sum(rowSums(s503 * t(s503))) # OK recip
ad <- abs(outer(rowSums(s501),rowSums(s501),"-"))
sum(s503*ad) # OK absOutDiffIntn

myeff <- getEffects(mydata)
myeff <- setEffect(myeff, absOutDiffIntn, name="mynet2", interaction1="mynet1", parameter=2)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
ad2 <- abs(outer(sqrt(rowSums(s501)),sqrt(rowSums(s501)),"-"))
sum(s503*ad2) # OK absOutDiffIntn

################################################################################
### check recipRateInv and recipRateLog
################################################################################

mynet <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet)
mymodel <- sienaModelCreate(projname=NULL, seed=534, cond=FALSE, firstg=0.01)
mymodel2 <- sienaModelCreate(projname=NULL, seed=534, cond=FALSE, firstg=0.02)
myeff <- getEffects(mydata)
effectsDocumentation(myeff)
myeff <- includeEffects(myeff, recipRate, type='rate')
(ans <- siena07(mymodel, data=mydata, effects=myeff))
(ans1 <- siena07(mymodel2, data=mydata, effects=myeff, prevAns=ans))
ans$targets
sum(abs(s501-s502)) ## OK Rate
sum(colSums(s501*t(s501))* rowSums(abs(s501-s502))) # OK recipRate
sum(s502) ## OK density
sum(s502*t(s502)) ## OK recip

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, recipRateInv, type='rate')
mymodel0 <- sienaModelCreate(projname=NULL, seed=534, cond=FALSE)
(ans <- siena07(mymodel0, data=mydata, effects=myeff))
ans$targets
sum(abs(s501-s502)) ## OK Rate
sum((1/(colSums(s501*t(s501))+1))* rowSums(abs(s501-s502))) # OK recipRateInv

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, recipRateLog, type='rate')
(ans <- siena07(mymodel0, data=mydata, effects=myeff))
(ans <- siena07(mymodel0, data=mydata, effects=myeff, prevAns=ans))
ans$targets
sum(abs(s501-s502)) ## OK Rate
sum((log(colSums(s501*t(s501))+1))* rowSums(abs(s501-s502))) # OK recipRateLog

##########################################################################
### check  SimAllNear,SimAllFar
##########################################################################

mynet <- sienaNet(array(c(s502, s503), dim=c(50, 50, 2)))
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
myeff <- getEffects(mydata)
effectsDocumentation(myeff)
# for SimAllNear effect:
(myeff <- setEffect(myeff,simAllNear,name='mybeh', parameter=2))
myeff
mymodel <- sienaModelCreate(projname=NULL, seed=514)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
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
myeff <- getEffects(mydata)
(myeff <- setEffect(myeff,simAllFar,name='mybeh', parameter=4))
mymodel <- sienaModelCreate(projname=NULL, seed=514, firstg=0.05, diagonalize=0.5)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets

p <- 4
aot <- abs(outer(z2,z2,"-"))
FF <- 1*(aot >= p)
diag(FF) <- 0
sum(FF*(p-aot)) # OK simAllFar

################################################################################
### check avDegIntn
################################################################################

mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- sienaNet(array(c(s503, s502, s501), dim=c(50, 50, 3)))
(mydata <- sienaDataCreate(mynet1, mynet2))

myeff <- getEffects(mydata)
(myeff <- includeEffects(myeff, avDegIntn, name='mynet2',
                interaction1='mynet1'))
mymodel <- sienaModelCreate(projname=NULL, seed=514)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
(ans <- siena07(mymodel, data=mydata, effects=myeff, prevAns=ans))
# convergence is very slow, which is natural for this effect with only 3 waves
ans$targets
################################################################################
### check avDeg
################################################################################

mynet <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
(mydata <- sienaDataCreate(mynet))

myeff <- getEffects(mydata)
(myeff <- includeEffects(myeff, avDeg))
mymodel <- sienaModelCreate(projname=NULL, seed=514, diagonalize=0.6)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
(ans <- siena07(mymodel, data=mydata, effects=myeff, prevAns=ans))
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

mynet <- sienaNet(array(c(s501, s502, s503, s502, s501), dim=c(50, 50, 5)))
(mydata <- sienaDataCreate(mynet))

myeff <- getEffects(mydata)
(myeff <- includeEffects(myeff, avDeg))
mymodel <- sienaModelCreate(projname=NULL, seed=514, diagonalize=0.6, nsub=5)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
(ans <- siena07(mymodel, data=mydata, effects=myeff, prevAns=ans))
(ans <- siena07(mymodel, data=mydata, effects=myeff, prevAns=ans))
# convergence is very slow, which is natural for this effect with only 5 waves
ans$targets

(ad2 <- sum(s502)/50)
(ad3 <- sum(s503)/50)
(ad1 <- sum(s501)/50)
p <- 2
sum(2*(ad2-p)*s502 + (ad3-p)*s503 + (ad1-p)*s501) # OK avDeg

myeff <- setEffect(myeff, density, initialValue=ans$theta[1], fix=TRUE, test=TRUE)
myeff <- setEffect(myeff, recip, initialValue=ans$theta[2], fix=TRUE, test=TRUE)
(ans1 <- siena07(mymodel, data=mydata, effects=myeff))
myeff$fix <- FALSE
myeff$test <- FALSE
(ans2 <- siena07(mymodel, data=mydata, effects=myeff, prevAns=ans1))
(ans3 <- siena07(mymodel, data=mydata, effects=myeff, prevAns=ans2))

################################################################################
### check sharedTo
################################################################################

mynet1 <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaNet(array(c(s502, s503), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet1, mynet2)
mymodel <- sienaModelCreate(projname=NULL, seed=1234)
myeff <- getEffects(mydata)
myeff <- setEffect(myeff, sharedTo, name="mynet2", interaction1="mynet1", parameter=1)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
# 115 116  70 106 122  90 154
intwostars1 <- (s501) %*% t(s501)
intwostars3 <- (s503) %*% t(s503)
diag(intwostars1) <- 0
diag(intwostars3) <- 0
sum(intwostars1 * intwostars3) # OK

###
myeff <- getEffects(mydata)
myeff <- setEffect(myeff, sharedTo, name="mynet2", interaction1="mynet1", parameter=2)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
# 115 116  70 106 122  90 139.7484
sum(sqrt(intwostars1) * intwostars3) # OK

###
myeff <- getEffects(mydata)
myeff <- setEffect(myeff, sharedTo, name="mynet2", interaction1="mynet1", parameter=3)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
# 115 116  70 106 122  90 139.7484
n <- 50
intwostars2 <- (s502) %*% t(s502)
diag(intwostars2) <- 0
(cc <- (n/(2*(n-1)))*mean(intwostars1 + intwostars2))
sum((intwostars1-cc) * intwostars3) # OK

###
myeff <- getEffects(mydata)
myeff <- setEffect(myeff, sharedTo, name="mynet2", interaction1="mynet1", parameter=4)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
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

mynet <- sienaNet(array(c(s501, s502), dim = c(50, 50, 2)))
beh <- sienaDependent(s50a[,1:2],type = 'behavior')
dcov <- coDyadCovar(dcova)
mydata <- sienaDataCreate(mynet,dcov,beh)
mymodel <- sienaModelCreate(projname = NULL, seed = 1234)

# avInAltW p1
myeff <- getEffects(mydata)
myeff <- setEffect(myeff,avInAltW, name = 'beh',
                   interaction1 = 'mynet', interaction2 = 'dcov',
                   parameter = 1)
(ans1 <- siena07(mymodel, data = mydata, effects = myeff))
ans1$targets

(mbh <- mean(beh))
divi <- function(a,b){ifelse(b==0, 0, a/b)}
mdc <- attr(mydata$dycCovars[[1]],'mean')
avinwalt <-  divi( ((dcova - mdc)*t(s501)) %*% (beh[,,2] - mbh), colSums(s501))
sum( (beh[,,2] - mbh) * avinwalt ) # avInAltW OK

# avInAltW p2
myeff <- getEffects(mydata)
myeff <- setEffect(myeff,avInAltW, name = 'beh',
                   interaction1 = 'mynet', interaction2 = 'dcov',
                   parameter = 2)
(ans1 <- siena07(mymodel, data = mydata, effects = myeff))
ans1$targets
avinwalt <-  divi( ((dcova - mdc)*t(s501)) %*% (beh[,,2] - mbh),
					rowSums(t(s501)*(dcova-mdc)))
sum( (beh[,,2] - mbh) * avinwalt ) # avInAltW OK

# avWInAlt
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,avWInAlt, name = 'beh',
                        interaction1 = 'mynet', interaction2 = 'dcov')
(ans1 <- siena07(mymodel, data = mydata, effects = myeff))
ans1$targets
sum(divi(rowSums(t(s501) * dcov_c), colSums(s501)) * beh2) # avWInAlt OK

# totInAltW
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,totInAltW, name = 'beh',
                        interaction1 = 'mynet', interaction2 = 'dcov')
(ans1 <- siena07(mymodel, data = mydata, effects = myeff))
ans1$targets
sum(rowSums(t(s501) * dcov_c %*% diag(beh2)) * beh2) # totInAltW OK

# totWInAlt
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,totWInAlt, name = 'beh',
                        interaction1 = 'mynet', interaction2 = 'dcov')
(ans1 <- siena07(mymodel, data = mydata, effects = myeff))
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

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,toAny,name="trust",interaction1="advice")
myeff
mymodel <- sienaModelCreate(projname = NULL, seed=123)
(myans <- siena07(mymodel, data = mydata, effects = myeff))
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

mynet <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet)
mymodel <- sienaModelCreate(projname=NULL, seed=534)
myeff <- getEffects(mydata)
effectsDocumentation(myeff)
myeff <- includeEffects(myeff, outAct)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
sum(rowSums(s502)^2) # OK
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, outAct_ego)
(ans1 <- siena07(mymodel, data=mydata, effects=myeff))
ans1$targets
(xbar <- (mean(rowSums(s501)) + mean(rowSums(s502)))/2)
sum(rowSums(s502)*(rowSums(s502)-xbar)) # OK
myeff <- getEffects(mydata)
myeff <- setEffect(myeff, outAct_ego, parameter=2)
(ans2 <- siena07(mymodel, data=mydata, effects=myeff))
ans2$targets
sum(rowSums(s502)*(sqrt(rowSums(s502))-sqrt(xbar))) # OK

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, reciAct_ego)
(ans1 <- siena07(mymodel, data=mydata, effects=myeff))
ans1$targets
(xbar <- (mean(rowSums(s501*t(s501))) + mean(rowSums(s502*t(s502))))/2)
sum(rowSums(s502)*(rowSums(s502*t(s502))-xbar)) # OK
myeff <- getEffects(mydata)
myeff <- setEffect(myeff, reciAct_ego, parameter=2)
(ans2 <- siena07(mymodel, data=mydata, effects=myeff))
ans2$targets
sum(rowSums(s502)*(sqrt(rowSums(s502*t(s502)))-sqrt(xbar))) # OK

myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, inPop_dya)
(ans1 <- siena07(mymodel, data=mydata, effects=myeff))
ans1$targets
(xbar <- (mean(colSums(s501)) + mean(colSums(s502)))/2)
sum(colSums(s502)*(colSums(s502)-xbar)) # OK
myeff <- getEffects(mydata)
myeff <- setEffect(myeff, inPop_dya, parameter=2)
(ans2 <- siena07(mymodel, data=mydata, effects=myeff))
ans2$targets
sum(colSums(s502)*(sqrt(colSums(s502))-sqrt(xbar))) # OK



################################################################################
### check sameXInPop, diffXInPop
################################################################################

mynet <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))

# construct actor covariate
(five <- rep(1:5,10))
five <- coCovar(five, centered=FALSE)

mydata <- sienaDataCreate(mynet, five)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,inPop)
myeff <- includeEffects(myeff,sameXInPop,interaction1='five')
myeff
mymodel <- sienaModelCreate(projname=NULL, seed=1234)
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
# check target statistics: indegree popularity
sum((colSums(s502))^2) # 386 OK
# check target statistics: sameXInPop
samef <- 1*outer(five,five,FUN="==")
sum(diag(t(s502) %*% samef %*% s502))
sum(s502 * (samef %*% s502))

myeff <- getEffects(mydata)
myeff <- setEffect(myeff,sameXInPop,interaction1='five', parameter=2)
myeff
(ans1 <- siena07(mymodel, data=mydata, effects=myeff))
ans1$targets
sum(s502 * sqrt(samef %*% s502)) # OK


myeff <- getEffects(mydata)
myeff <- includeEffects(myeff,diffXInPop,interaction1='five')
myeff
(ans <- siena07(mymodel, data=mydata, effects=myeff))
ans$targets
difff <- 1*outer(five,five,FUN="!=")
sum(s502 * (difff %*% s502)) # OK

myeff <- getEffects(mydata)
myeff <- setEffect(myeff,diffXInPop,interaction1='five', parameter=2)
myeff
(ans1 <- siena07(mymodel, data=mydata, effects=myeff))
ans1$targets
sum(s502 * sqrt(difff %*% s502)) # OK




