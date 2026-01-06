##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: https://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: checkImpossibleChanges.r
## *
## * Description: This file contains the function checkImpossibleChanges
## * which checks for impossible changes from structural values to
## * different observed values.
## * Used for maximum likelihood and Bayesian estimation.
## *
## ****************************************************************************/
##args:x: Siena data object
## checkImpossibleChanges returns 1*I{any structural zero changed to 1} +
##                          2*I{any structural one changed to 0} +
##       10*I{any structural zero changed to NA and then changed to 1} +
##       20*I{any structural one  changed to NA and then changed to 0} 


notNetwork <- function(depv){
	ifelse(length(dim(depv))<= 2, TRUE, 			
		(!(attr(depv, 'type') %in% c("oneMode", "bipartite"))))
}

##@checkImpossibleChanges checks for likelihood calculations
checkImpossibleChanges <- function(x)
{
	if (!inherits(x,'sienadata'))
	{
		stop('checkImpossibleChanges can only be applied to siena data objects')
	}
	xd <- x$depvars
	impossibleChangeOne <- function(dv)
	{
		ifelse(notNetwork(dv), 0, 
			1*any(sapply((2:dim(dv)[3]),
					FUN=function(i){any((dv[,,i-1] == 10) &
					                    (dv[,,i] %in% c(1,11)), na.rm=TRUE)})))
	}
	impossibleChangeZero <- function(dv)
	{
		ifelse(notNetwork(dv), 0, 
			1*any(sapply((2:dim(dv)[3]),
					FUN=function(i){any((dv[,,i-1] == 11) &
					                  (dv[,,i] %in% c(0,10)), na.rm=TRUE)})))
	}
	impossibleChangeNAOne <- function(dv)
	{
		ifelse(notNetwork(dv), 0, 
			ifelse((dim(dv)[3] >= 3),
			1*any(sapply((3:dim(dv)[3]),
					FUN=function(i){
		any((dv[,,i-2] == 10) & (is.na(dv[,,i-1])) & (dv[,,i] %in% c(1,11)),
								na.rm=TRUE)})), 0))
					
	}
	impossibleChangeNAZero <- function(dv)
	{
		ifelse(notNetwork(dv), 0, 
			ifelse((dim(dv)[3] >= 3),
			1*any(sapply((3:dim(dv)[3]),
					FUN=function(i){
		any(((dv[,,i-2] == 11) & (is.na(dv[,,i-1])) &
					(dv[,,i] %in% c(0,10))), na.rm=TRUE)})), 0))					
	}
# This is not complete; see the manual;
# it does not check for patterns 10, NA, ..., NA, 1
# where ... are all NA; neither for patterns 11, NA, ..., NA, 0
	max(sapply(xd,impossibleChangeOne)) + 2*max(sapply(xd,impossibleChangeZero)) +
			10*max(sapply(xd,impossibleChangeNAOne)) + 
			20*max(sapply(xd,impossibleChangeNAZero))
}


##@checkZeroChanges checks for likelihood calculations
checkZeroChanges <- function(x)
{
	if (!inherits(x,'sienadata'))
	{
		stop('checkZeroChanges can only be applied to siena data objects')
	}
	zeroChange <- function(dv){
		if (attr(dv, 'type') %in% c("oneMode", "bipartite"))
		{
			dv[dv==10] <- 0
			dv[dv==11] <- 1
		}
		ifelse((length(dim(dv))==3),
			1*any(sapply((2:dim(dv)[3]),
					FUN=function(i){all(dv[,,i-1] == dv[,,i], na.rm = TRUE)})), 0)
	}
	sum(sapply(x$depvars,zeroChange))
}
