library(RSiena)


##test3
mynet1 <- sienaDependent(array(c(tmp3, tmp4),dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff<- getEffects(mydata)
mymodel<- model.create(findiff=TRUE, fn = simstats0c,
                       cond=FALSE, nsub=2, n3=50, seed=3)
print('test3')
ans<- siena07(mymodel, data=mydata, effects=myeff,
     batch=TRUE, parallelTesting=TRUE, silent=TRUE)
#,dll='../siena/src/RSiena.dll')
ans
(myeff <- includeEffects(myeff, transTrip, cycle4))
(myeff <- includeEffects(myeff, cycle4, include=FALSE))
##test4
mymodel$cconditional <- TRUE
mymodel$condvarno<- 1
print('test4')
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE,
              parallelTesting=TRUE, silent=TRUE)
##, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
ans
##test5
mynet1 <- sienaDependent(array(c(tmp3,tmp4),dim=c(32,32,2)))
mydata <- sienaDataCreate(mynet1)
myeff<- getEffects(mydata)
mymodel<- model.create(fn = simstats0c,  nsub=2, n3=50,
                       cond=FALSE, seed=5)
print('test5')
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE,
              parallelTesting=TRUE, silent=TRUE)
ans
(myeff <- includeEffects(myeff, recip, inPop))
(myeff <- includeEffects(myeff, outAct, fix=TRUE, test=TRUE))
(myeff <- includeInteraction(myeff, recip, inPop, fix=TRUE, test=TRUE))
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE,
              parallelTesting=TRUE, silent=TRUE)
ans
score.Test(ans)
##test6
mynet1 <- sienaDependent(array(c(tmp3,tmp4),dim=c(32,32,2)))
mydata <- sienaDataCreate(mynet1)
myeff<- getEffects(mydata)
mymodel<- model.create(fn = simstats0c,  nsub=2, n3=50,
                       cond=FALSE, doubleAveraging=0,seed=5)
print('test6')
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE,
              parallelTesting=TRUE, silent=TRUE)
ans
myeff <- includeEffects(myeff, recip, include=FALSE)
myeff <- includeEffects(myeff, recip, type='endow')
myeff <- includeEffects(myeff, recip, type='creation')
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE,
              parallelTesting=TRUE, silent=TRUE)
ans
testSame.RSiena(ans, 3, 4)
##test7
mynet1 <- sienaDependent(array(c(tmp3,tmp4),dim=c(32,32,2)))
mydata <- sienaDataCreate(mynet1)
myeff<- getEffects(mydata)
mymodel<- model.create(fn = simstats0c,  nsub=2, n3=50,
                       cond=FALSE,  diagonalize=0.5, seed=5)
print('test7')
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE,
              parallelTesting=TRUE, silent=TRUE)
##, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
ans
##test8
mymodel<- model.create(fn = simstats0c,  nsub=1, n3=50,
                       cond=TRUE, condvarno=1, seed=5)
print('test8')
ans <- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE,
              parallelTesting=TRUE, silent=TRUE)
##, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
ans
##test9
mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- sienaDependent(s50a,type='behavior')
mydata <- sienaDataCreate(mynet1, mynet2)
myeff <- getEffects(mydata)
myeff <- setEffect(myeff, linear, initialValue=0.34699930338, name="mynet2")
myeff <- setEffect(myeff, avAlt, name="mynet2", interaction1="mynet1")
##myeff$initialValue[98] <- 0.34699930338 ## siena3 starting values differ
##test10
print('test10')
mymodel$cconditional <- TRUE
mymodel$condvarno<- 1
ans <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE,
               parallelTesting=TRUE, silent=TRUE)
##, verbose=TRUE)
ans
##test11
print('test11')
mymodel<- model.create(fn = simstats0c,  nsub=1, n3=50,
                       behModelType=c(mynet2=2), seed=6)
(ans <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE,
               parallelTesting=TRUE, silent=TRUE))
##test12
print('test12')
use<- 1:30
mynet1 <- sienaDependent(array(c(s501[use,], s502[use,], s503[use,]),
                         dim=c(length(use), 50,3)), type='bipartite',
                         nodeSet=c('Senders','receivers'))
receivers <- sienaNodeSet(50,'receivers')
senders <- sienaNodeSet(30,'Senders')
myvar1 <- coCovar(s50a[1:30,2], nodeSet='Senders')
mydata <- sienaDataCreate(mynet1, myvar1, nodeSets=list(senders, receivers))
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, inPop)
myeff <- setEffect(myeff, altInDist2, interaction1="myvar1", parameter=1)
ans <- siena07(sienaModelCreate(n3=50, nsub=2, seed=1),
               data=mydata, effects=myeff, batch=TRUE, silent=TRUE)
ans
tt <- sienaTimeTest(ans)
summary(tt)
##test13
print('test13')
use<- 1:30
mynet1 <- sienaDependent(array(c(s502[,use], s503[,use]),
                         dim=c(50, length(use), 2)), type='bipartite',
                         nodeSet=c('Senders','receivers'))
receivers <- sienaNodeSet(30,'receivers')
senders <- sienaNodeSet(50,'Senders')
myvar1 <- coCovar(s50a[1:50,2], nodeSet='Senders')
mydata <- sienaDataCreate(mynet1, myvar1, nodeSets=list(senders, receivers))
myeff <- getEffects(mydata)
myeff <- setEffect(myeff, altInDist2, interaction1="myvar1", parameter=1)
myeff <- setEffect(myeff, egoX, interaction1="myvar1")
(ans <- siena07(sienaModelCreate(n3=50, nsub=2, seed=1),
               data=mydata, effects=myeff, batch=TRUE, silent=TRUE))
##test14
print('test14')
net <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
dataset <- sienaDataCreate(net)
myeff <- getEffects(dataset)
myeff <- includeEffects(myeff, inPop)
algo <- sienaAlgorithmCreate(nsub=1, n3=20, maxlike=TRUE, seed=15, mult=1, prML=1)
(ans <- siena07(algo, data=dataset, effects=myeff, batch=TRUE, silent=TRUE))
##test 15
print('test15')
mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- sienaDependent(s50a,type='behavior')
mydata <- sienaDataCreate(mynet1, mynet2)
myeff <- getEffects(mydata)
(myeff <- includeEffects(myeff, transTrip))
(myeff <- includeEffects(myeff, egoX, simX, interaction1="mynet2"))
(myeff <- includeEffects(myeff, avSim, name="mynet2", interaction1="mynet1"))
(myeff <- includeGMoMStatistics(myeff, simX_gmm, interaction1="mynet2"))
algo <- sienaAlgorithmCreate(nsub=1, n3=100, gmm=TRUE, seed=6)
(ans <- siena07(algo, data=mydata, effects=myeff, batch=TRUE,
                parallelTesting=TRUE, silent=TRUE))
##test16
print('test16')
set.seed(123) # simulate behavior data according to dZ(t) = [-0.1 Z + 1] dt + 1 dW(t)
y1 <- rnorm(50, 0,3)
y2 <- exp(-0.1) * y1 + (1-exp(-0.1)) * 1/ -0.1 + rnorm(50, 0, (exp(-0.2)- 1) / -0.2 * 1^2)
friend <- sienaDependent(array(c(s501, s502), dim = c(50,50,2)))
behavior <- sienaDependent(matrix(c(y1,y2), 50,2), type = "continuous")
(mydata <- sienaDataCreate(friend, behavior))
(myeff <- getEffects(mydata, onePeriodSde = TRUE))
algorithmMoM <- sienaAlgorithmCreate(nsub=1, n3=20, seed=321)
(ans <- siena07(algorithmMoM, data = mydata, effects = myeff, batch=TRUE,
                                                silent=TRUE))
##test17
print('test17')
mynet <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
sm1 <- 1*(s50s[,2] >= 2)
sm2 <- 1*(s50s[,3] >= 2)
sm2 <- pmax(sm1,sm2)
sm2[c(33,28,29,44)] <- 1
mybeh <- sienaDependent(cbind(sm1,sm2), type="behavior")
(mydata <- sienaDataCreate(mynet, mybeh))
mymodel <- sienaModelCreate(projname=NULL, seed=1234, firstg=0.001, nsub=1, n3=10)
myeff <- getEffects(mydata)
(myeff <- setEffect(myeff,avExposure,type='rate',parameter=2,
                            name='mybeh',interaction1='mynet'))
(ans <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE, silent=TRUE))
##test18
print('test18')
myalgorithm <- sienaAlgorithmCreate(nsub=1, n3=10, seed=1293)
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
cova <- coCovar(1:32)
cova2 <- coCovar(rep(0,32), warn=FALSE)
mydata <- sienaDataCreate(mynet1, cova)
mydata2 <- sienaDataCreate(mynet1, cova=cova2)
mygroup <- sienaGroupCreate(list(mydata,mydata2))
myeff <- getEffects(mygroup)
myeff <- setEffect(myeff, simX, interaction1='cova')
(ans <- siena07(myalgorithm, data=mygroup, effects=myeff, batch=TRUE,
                              silent=TRUE))
##test19
print('test19')
myalgorithm <- sienaAlgorithmCreate(nsub=1, n3=20, seed=1943)
mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
alc <- coCovar(s50a[,1])
smoke <- coCovar(s50s[,1])
mydata <- sienaDataCreate(mynet, alc, smoke)
myeff <- getEffects(mydata)
myeff <- setEffect(myeff, gwespFF)
myeff <- setEffect(myeff, gwespFF, parameter=20)
myeff <- setEffect(myeff, outTrunc, parameter=2, include=FALSE)
myeff <- includeInteraction(myeff, outTrunc, egoX, egoX, 
						interaction1=c("","smoke","alc"))
(ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE,
                              silent=TRUE))
## delete output file
if (file.exists('Siena.txt')){unlink('Siena.txt')}