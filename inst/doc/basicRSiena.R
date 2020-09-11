## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
getwd()

## -----------------------------------------------------------------------------
# setwd()

## -----------------------------------------------------------------------------
# setwd("C:\\Users\\tom.snijders\\Documents\\Siena\\s50_script")

## -----------------------------------------------------------------------------
# friend.data.w1 <- as.matrix(read.table("s50-network1.dat"))
# friend.data.w2 <- as.matrix(read.table("s50-network2.dat"))
# friend.data.w3 <- as.matrix(read.table("s50-network3.dat"))
# drink <- as.matrix(read.table("s50-alcohol.dat"))
# smoke <- as.matrix(read.table("s50-smoke.dat"))

## -----------------------------------------------------------------------------
library(RSiena)
# Now we use the internally available s50 data set.
# Look at its description:
?s50
# 3 waves, 50 actors
# Look at the start and end of the first wave matrix
head(s501)
tail(s501)
# and at the alcohol variable
s50a
# Now define the objects with the same names as above
# (this step is superfluous if you read the data already).
        friend.data.w1 <- s501
        friend.data.w2 <- s502
        friend.data.w3 <- s503
        drink <- s50a
        smoke <- s50s

## -----------------------------------------------------------------------------
?sienaDependent
# First create a 50 * 50 * 3 array composed of the 3 adjacency matrices
friendshipData <- array( c( friend.data.w1, friend.data.w2, friend.data.w3 ),
           dim = c( 50, 50, 3 ) )
# and next give this the role of the dependent variable:
	friendship <- sienaDependent(friendshipData)

## -----------------------------------------------------------------------------
smoke1 <- coCovar( smoke[ , 1 ] )
# A variable actor covariate is defined for drinking:
alcohol <- varCovar( drink )
# (This choice is purely for the purpose of illustration here.)

## -----------------------------------------------------------------------------
?sienaDataCreate
mydata <- sienaDataCreate( friendship, smoke1, alcohol )
# Check what we have
mydata

## -----------------------------------------------------------------------------
print01Report( mydata, modelname="s50")
# For the model specification we need to create the effects object
myeff <- getEffects( mydata )
# All the effects that are available given the structure
# of this data set can be seen from
effectsDocumentation(myeff)
# For a precise description of all effects, see Chapter 12 in the RSiena manual.
# A basic specification of the structural effects:
?includeEffects
myeff <- includeEffects( myeff, transTrip, cycle3)
# and some covariate effects:
myeff <- includeEffects( myeff, egoX, altX, simX, interaction1 = "alcohol" )
myeff <- includeEffects( myeff, simX, interaction1 = "smoke1" )
myeff

## -----------------------------------------------------------------------------
?sienaAlgorithmCreate
myalgorithm <- sienaAlgorithmCreate( projname = 's50' )

## -----------------------------------------------------------------------------
?siena07
ans <- siena07( myalgorithm, data = mydata, effects = myeff)
ans

## -----------------------------------------------------------------------------
ans <- siena07( myalgorithm, data = mydata, effects = myeff, prevAns=ans)
ans
# If convergence is good, you can look at the estimates.
# More extensive results
summary(ans)

## -----------------------------------------------------------------------------
(myeff <- includeEffects( myeff, transRecTrip))
(ans1 <- siena07( myalgorithm, data = mydata, effects = myeff, prevAns=ans))
# If necessary, repeat the estimation with the new result:
(ans1 <- siena07( myalgorithm, data = mydata, effects = myeff, prevAns=ans1))

## -----------------------------------------------------------------------------
# Once again, look at the help file
?sienaDependent
# now paying special attention to the <<type>> parameter.
drinking <- sienaDependent( drink, type = "behavior" )

## -----------------------------------------------------------------------------
NBdata <- sienaDataCreate( friendship, smoke1, drinking )
NBdata
NBeff <- getEffects( NBdata )
effectsDocumentation(NBeff)
NBeff <- includeEffects( NBeff, transTrip, transRecTrip )
NBeff <- includeEffects( NBeff, egoX, egoSqX, altX, altSqX, diffSqX,
                         interaction1 = "drinking" )
NBeff <- includeEffects( NBeff, egoX, altX, simX, interaction1 = "smoke1" )
NBeff
# For including effects also for the dependent behaviour variable, see
?includeEffects
NBeff <- includeEffects( NBeff, avAlt, name="drinking",
                         interaction1 = "friendship" )
NBeff
# Define an algorithm with a new project name
myalgorithm1 <- sienaAlgorithmCreate( projname = 's50_NB' )

# Estimate again, using the second algorithm right from the start.
NBans <- siena07( myalgorithm1, data = NBdata, effects = NBeff)
# You may improve convergence (considering the overall maximum
# convergence ratio) by repeated estimation in the same way as above.

# Look at results
NBans
# Make a nicer listing of the results
siena.table(NBans, type="html", sig=TRUE)

