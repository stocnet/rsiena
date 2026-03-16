library(RSiena)

# When new effects are added, the numbering of effects changes.
# This will have consequences for the output of set_interaction,
# and will require adaptation of parallel.Rout.save.

##test3
mynet1 <- as_dependent_rsiena(array(c(tmp3, tmp4),dim=c(32, 32, 2)))
mydata <- make_data_rsiena(mynet1)
myeff<- make_specification(mydata)
print('test3')
alg_alg <- set_algorithm_saom(cond=FALSE, seed=3, n3=50, nsub=2, findiff=TRUE)
ans <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg)
ans
(myeff <- set_effect(myeff, list(transTrip, cycle4)))
(myeff <- set_effect(myeff, cycle4, include=FALSE))
##test4
print('test4')
alg_alg <- set_algorithm_saom(cond=TRUE, condvarno=1, seed=3, n3=50, nsub=2, findiff=TRUE)
ans <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg)
ans
##test5
mynet1 <- as_dependent_rsiena(array(c(tmp3,tmp4),dim=c(32,32,2)))
mydata <- make_data_rsiena(mynet1)
myeff<- make_specification(mydata)
print('test5')
alg_alg <- set_algorithm_saom(cond=FALSE, seed=5, n3=50, nsub=2)
ans <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg)
ans
(myeff <- set_effect(myeff, list(recip, inPop)))
(myeff <- set_effect(myeff, outAct, fix=TRUE, test=TRUE))
(myeff <- set_interaction(myeff, list(recip, inPop), fix=TRUE, test=TRUE))
alg_alg <- set_algorithm_saom(cond=FALSE, seed=5, n3=50, nsub=2)
ans <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg)
ans
test_parameter(ans, method='score')
##test6
mynet1 <- as_dependent_rsiena(array(c(tmp3,tmp4),dim=c(32,32,2)))
mydata <- make_data_rsiena(mynet1)
myeff<- make_specification(mydata)
print('test6')
alg_alg <- set_algorithm_saom(cond=FALSE, seed=5, n3=50, nsub=2)
ans <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg)
ans
myeff <- set_effect(myeff, recip, include=FALSE)
myeff <- set_effect(myeff, recip, type="endow")
myeff <- set_effect(myeff, recip, type="creation")
alg_alg <- set_algorithm_saom(cond=FALSE, seed=5, n3=50, nsub=2)
ans <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg)
ans
test_parameter(ans, method='same', tested=3, tested2=4)
##test7
mynet1 <- as_dependent_rsiena(array(c(tmp3,tmp4),dim=c(32,32,2)))
mydata <- make_data_rsiena(mynet1)
myeff<- make_specification(mydata)
print('test7')
alg_alg <- set_algorithm_saom(cond=FALSE, seed=5, n3=50, nsub=2,
         diagonalize=0.5)
ans <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg)
ans
##test8
print('test8')
alg_alg <- set_algorithm_saom(cond=TRUE, condvarno=1, seed=5, n3=50, nsub=1)
ans  <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg)
ans
##test9
mynet1 <- as_dependent_rsiena(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- as_dependent_rsiena(s50a,type='behavior')
mydata <- make_data_rsiena(mynet1, mynet2)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, linear, depvar="mynet2", initialValue=0.34699930338)
myeff <- set_effect(myeff, avAlt, depvar="mynet2", covar1="mynet1")
##test10
print('test10')
alg_alg <- set_algorithm_saom(cond=TRUE, condvarno=1, seed=5, n3=50, nsub=1)
ans  <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg)
##, verbose=TRUE)
ans
##test11
print('test11')
alg_model <- set_model_saom(behModelType=c(mynet2=2))
alg_alg <- set_algorithm_saom(seed=6, n3=50, nsub=1)
(ans  <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_model=alg_model, control_algo=alg_alg))
##test12
print('test12')
use<- 1:30
mynet1 <- as_dependent_rsiena(array(c(s501[use,], s502[use,], s503[use,]),
                         dim=c(length(use), 50,3)), type='bipartite',
                         nodeSet=c('Senders','receivers'))
receivers <- as_nodeset_rsiena(50,'receivers')
senders <- as_nodeset_rsiena(30,'Senders')
myvar1 <- as_covariate_rsiena(s50a[1:30,2], nodeSet='Senders')
mydata <- make_data_rsiena(mynet1, myvar1, nodeSets=list(senders, receivers))
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, inPop)
myeff <- set_effect(myeff, altInDist2, covar1="myvar1", parameter=1)
alg_alg <- set_algorithm_saom(seed=1, n3=50, nsub=2)
ans  <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg)
ans
tt <- test_time(ans)
summary(tt)
##test13
print('test13')
use<- 1:30
mynet1 <- as_dependent_rsiena(array(c(s502[,use], s503[,use]),
                         dim=c(50, length(use), 2)), type='bipartite',
                         nodeSet=c('Senders','receivers'))
receivers <- as_nodeset_rsiena(30,'receivers')
senders <- as_nodeset_rsiena(50,'Senders')
myvar1 <- as_covariate_rsiena(s50a[1:50,2], nodeSet='Senders')
mydata <- make_data_rsiena(mynet1, myvar1, nodeSets=list(senders, receivers))
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, altInDist2, covar1="myvar1", parameter=1)
myeff <- set_effect(myeff, egoX, covar1="myvar1")
alg_alg <- set_algorithm_saom(seed=1, n3=50, nsub=2)
(ans  <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg))
##test14
print('test14')
net <- as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
dataset <- make_data_rsiena(net)
myeff <- make_specification(dataset)
myeff <- set_effect(myeff, inPop)
alg_alg <- set_algorithm_saom(maxlike=TRUE, cond=FALSE, seed=15, n3=20, nsub=1,
         diagonalize=0, mult=1)
(ans  <- siena(data=dataset, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg))
##test 15
print('test15')
mynet1 <- as_dependent_rsiena(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- as_dependent_rsiena(s50a,type='behavior')
mydata <- make_data_rsiena(mynet1, mynet2)
myeff <- make_specification(mydata)
(myeff <- set_effect(myeff, transTrip))
(myeff <- set_effect(myeff, list(egoX, simX), covar1="mynet2"))
(myeff <- set_effect(myeff, avSim, depvar="mynet2", covar1="mynet1"))
(myeff <- includeGMoMStatistics(myeff, simX_gmm, covar1="mynet2"))
alg_alg <- set_algorithm_saom(gmm=TRUE, seed=6, n3=100, nsub=1)
(ans  <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg))
##test16
print('test16')
set.seed(123) # simulate behavior data according to dZ(t) = [-0.1 Z + 1] dt + 1 dW(t)
y1 <- rnorm(50, 0,3)
y2 <- exp(-0.1) * y1 + (1-exp(-0.1)) * 1/ -0.1 + rnorm(50, 0, (exp(-0.2)- 1) / -0.2 * 1^2)
friend <- as_dependent_rsiena(array(c(s501, s502), dim = c(50,50,2)))
behavior <- as_dependent_rsiena(matrix(c(y1,y2), 50,2), type = "continuous")
(mydata <- make_data_rsiena(friend, behavior))
(myeff <- make_specification(mydata, onePeriodSde = TRUE))
alg_alg <- set_algorithm_saom(seed=321, n3=20, nsub=1)
(ans  <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg))
##test17
print('test17')
mynet <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
sm1 <- 1*(s50s[,2] >= 2)
sm2 <- 1*(s50s[,3] >= 2)
sm2 <- pmax(sm1,sm2)
sm2[c(33,28,29,44)] <- 1
mybeh <- as_dependent_rsiena(cbind(sm1,sm2), type="behavior")
(mydata <- make_data_rsiena(mynet, mybeh))
myeff <- make_specification(mydata)
(myeff <- set_effect(myeff, avExposure, type="rate", depvar="mybeh",
         covar1="mynet", parameter=2))
alg_alg <- set_algorithm_saom(seed=1234, n3=10, nsub=1, firstg=0.001)
(ans  <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg))
##test18
print('test18')
mynet1 <- as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
cova <- as_covariate_rsiena(1:32)
cova2 <- as_covariate_rsiena(rep(0,32), warn=FALSE)
mydata <- make_data_rsiena(mynet1, cova)
mydata2 <- make_data_rsiena(mynet1, cova=cova2)
mygroup <- make_group_rsiena(list(mydata,mydata2))
myeff <- make_specification(mygroup)
myeff <- set_effect(myeff, simX, covar1="cova")
alg_alg <- set_algorithm_saom(seed=1293, n3=10, nsub=1)
(ans  <- siena(data=mygroup, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg))
##test19
print('test19')
mynet <- as_dependent_rsiena(array(c(s501, s502), dim=c(50, 50, 2)))
alc <- as_covariate_rsiena(s50a[,1])
smoke <- as_covariate_rsiena(s50s[,1])
mydata <- make_data_rsiena(mynet, alc, smoke)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, gwespFF)
myeff <- set_effect(myeff, gwespFF, parameter=20)
myeff <- set_effect(myeff, outTrunc, parameter=2, include=FALSE)
myeff <- set_interaction(myeff, list(outTrunc, egoX, egoX), covar1=c("",
         "smoke", "alc"))
alg_alg <- set_algorithm_saom(seed=1943, n3=20, nsub=1)
(ans  <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg))
##test20
print('test20')
mynet1 <- as_dependent_rsiena(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- as_dependent_rsiena(array(c(s503, s502), dim=c(50, 50, 2)))
mydata <- make_data_rsiena(mynet1, mynet2)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, crprod, depvar="mynet2", covar1="mynet1")
myeff <- set_effect(myeff, from, depvar="mynet1", covar1="mynet2")
(myeff <- includeGMoMStatistics(myeff, from_gmm, depvar='mynet1',
                                covar1='mynet2'))
alg_alg <- set_algorithm_saom(seed=1293, n3=50, nsub=2)
(ans  <- siena(data=mydata, effects=myeff[myeff$type != "gmm", ], batch=TRUE,
         silent=TRUE, control_algo=alg_alg))
alg_alg <- set_algorithm_saom(gmm=TRUE, seed=1293, n3=50, nsub=2)
(ans1  <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         control_algo=alg_alg))
alg_alg <- set_algorithm_saom(gmm=TRUE, seed=1293, n3=50, nsub=2)
(ans2  <- siena(data=mydata, effects=myeff, batch=TRUE, silent=TRUE,
         prevAns=ans1, control_algo=alg_alg))

##test21
# Run simple test model ----
mynet <- as_dependent_rsiena(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata <- make_data_rsiena(mynet)
mymodel <- make_specification(mydata)
## TransitiveTriplets model
mymodel <- set_effect(mymodel, transTrip, depvar="mynet")
# Test returnChangeContributions when running siena directly ----

print('test21')
alg_out <- set_output_saom(returnChangeContributions=TRUE)
alg_alg <- set_algorithm_saom(cond=FALSE, seed=42, n3=60, nsub=1)
ans  <- siena(data=mydata, effects=mymodel, batch=TRUE, silent=TRUE,
         control_algo=alg_alg, control_out=alg_out)
length(ans$changeContributions) # 60 as expected
head(ans$changeContributions[[1]][[1]][[1]])
dim(ans$changeContributions[[1]][[1]][[1]][[1]]) # 3 x 50 as expected

# Test returnChangeContributions when setting nsub = 0 and prevAns = ans, 
#using batch mode
alg_out <- set_output_saom(returnChangeContributions=TRUE)
alg_alg <- set_algorithm_saom(cond=FALSE, seed=42, n3=60, nsub=1)
ans2  <- siena(data=mydata, effects=mymodel, batch=TRUE, prevAns=ans,
         control_algo=alg_alg, control_out=alg_out)
length(ans$changeContributions) # 60 chains as expected
head(ans$changeContributions[[1]][[1]][[1]])
##test22
# Test interpret_size_dynamics ----
## Unconditional Estimation ----
### Use changeContributions from ans ----
print('test22')
RIDynamics1 <- interpret_size_dynamics(data=mydata, 
     ans=ans, 
     useChangeContributions=TRUE, 
     intervalsPerPeriod=10)
RIDynamics1
### Unconditional Estimation with model rerun
mycontrols2 <- sienaAlgorithmCreate(nsub=2, n3=60, cond=FALSE, seed = 84)
mynet2 <- as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata2 <- make_data_rsiena(mynet2)
mymodel2 <- make_specification(mydata2)
mymodel2 <- set_effect(mymodel2, list(density, recip, transTies))
alg_alg <- set_algorithm_saom(cond=FALSE, seed=84, n3=60, nsub=2)
ans2  <- siena(data=mydata2, effects=mymodel2, batch=TRUE, silent=TRUE,
         control_algo=alg_alg)
RIDynamics2 <-  interpret_size_dynamics(mydata2, ans=ans2,
     useChangeContributions=FALSE, algorithm = mycontrols2)
RIDynamics2
### Don't use ans but previously estimated coefficients ----
RIDynamics3 <- interpret_size_dynamics(data=mydata2, theta=c(ans2$rate,ans2$theta),
             algorithm=mycontrols2, effects=mymodel2, intervalsPerPeriod=10)
RIDynamics3

## Conditional Estimation ----
mycontrols3 <- sienaAlgorithmCreate(nsub=2, n3=60, cond=TRUE, seed = 4242)
ans3 <- siena07(mycontrols3, data=mydata2, effects=mymodel2, batch=TRUE, 
     silent=TRUE)

RIDynamics4 <- interpret_size_dynamics(mydata2, ans=ans3, effects = mymodel2, 
     useChangeContributions=FALSE, algorithm=mycontrols3)
RIDynamics4
##test23
## Test interpret_size_dynamics with behavior variable ----
mynet3 <- as_dependent_rsiena(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- as_dependent_rsiena(s50a, type="behavior")
mydata3 <- make_data_rsiena(mynet3, mybeh)
mymodel3 <- make_specification(mydata3)
mymodel3 <- set_effect(mymodel3, list(density, recip, transTies, transTrip))
print('test23')
alg_out <- set_output_saom(returnChangeContributions=TRUE)
alg_alg <- set_algorithm_saom(seed=8484, n3=60, nsub=2)
ans4  <- siena(data=mydata3, effects=mymodel3, batch=TRUE, silent=TRUE,
         control_algo=alg_alg, control_out=alg_out)
length(ans4$changeContributions) # 60 chains as expected
length(ans4$changeContributions[[1]][[1]]) # 2 periods as expected
beh_steps <- sapply(ans4$changeContributions[[1]][[1]][[2]], 
     function(ministep) attr(ministep, "networkName")) %in%  "mybeh"
any(beh_steps) # TRUE as expected
ministeps <- ans4$changeContributions[[1]][[1]][[1]]
getDepvarName <- function(ministep) attr(ministep, "networkName")
beh_steps <- Filter(function(ministep) getDepvarName(ministep)  == "mybeh", 
     ans4$changeContributions[[1]][[1]][[1]]) # 1st period, 1st group, 1st chain, behavior steps
RIDynamics5 <- interpret_size_dynamics(mydata3, ans=ans4, depvar="mybeh",
                               useChangeContributions=TRUE)
RIDynamics5
net_steps <- Filter(function(ministep) getDepvarName(ministep)  == "mynet3", 
     ans4$changeContributions[[1]][[1]][[1]]) # 1st period, 1st group, 1st chain, network steps
RIDynamics6 <- interpret_size_dynamics(mydata3, ans=ans4, depvar="mynet3", useChangeContributions=TRUE)
RIDynamics6
##test24
mynet1 <- as_dependent_rsiena(array(c(tmp3,tmp4),dim=c(32,32,2)))
mydata <- make_data_rsiena(mynet1)
myeff<- make_specification(mydata)
myeff <- set_effect(myeff, list(recip, inPop, outAct))
thv <- matrix(NA, 10, 5)
thv[,1] <- 3 +0.1*(1:10)
thv[,2] <- c(-2, -2.4, -2.3, -2.5, -2 + (1:6)*0.15)
thv[,3] <- rep(2, 10)
thv[,4] <- 0.02 * (1:10)
thv[,5] <- 0.02 * (10:1) 
myalg <- set_algorithm_saom(nsub=2, n3=50, cond=FALSE, seed=5, simOnly=TRUE,
				thetaValue=thv)
print('test24')
(ans <- siena(mymodel, data=mydata, effects=myeff, control_algo=myalg, batch=TRUE, silent=TRUE))
## delete output files
if (file.exists('mydata_report.txt')){unlink('mydata_report.txt')}
if (file.exists('mydata_out.txt')){unlink('mydata_out.txt')}
if (file.exists('mydata2_out.txt')){unlink('mydata2_out.txt')}
if (file.exists('mydata3_out.txt')){unlink('mydata3_out.txt')}
if (file.exists('mygroup_out.txt')){unlink('mygroup_out.txt')}
if (file.exists('dataset_out.txt')){unlink('dataset_out.txt')}
