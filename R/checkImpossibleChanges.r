##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snijders/siena
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
## returns 1*I{any structural zero changed to 1} +
##                          2*I{any structural one changed to 0} +
##                          10*I{any structural zero changed to NA} +
##@checkImpossibleChanges checks for likelihood calculations
checkImpossibleChanges <- function(x)
{
	if (!inherits(x,'siena'))
	{
		stop('checkImpossibleChanges can only be applied to siena data objects')
	}
	xd <- x$depvars
	impossibleChangeOne <- function(dv)
	{
		ifelse((length(dim(dv))==3),
			1*any(sapply((2:dim(dv)[3]),
					FUN=function(i){any((dv[,,i-1] == 10) &
					                    (dv[,,i]==1), na.rm=TRUE)})), 0)
	}
	impossibleChangeZero <- function(dv)
	{
		ifelse((length(dim(dv))==3),
			1*any(sapply((2:dim(dv)[3]),
					FUN=function(i){any((dv[,,i-1] == 11) &
					                    (dv[,,i]==0), na.rm=TRUE)})), 0)
	}
	impossibleChangeNA <- function(dv)
	{
		ifelse((length(dim(dv))==3),
			ifelse((dim(dv)[3] >= 3),
			1*any(sapply((3:dim(dv)[3]),
					FUN=function(i){
					any(((dv[,,i-2] == 10) & (is.na(dv[,,i])) & (dv[,,i]==1)), na.rm=TRUE)})
					), 0), 0)
	}
	max(sapply(xd,impossibleChangeOne)) + 2*max(sapply(xd,impossibleChangeZero)) +
			10*max(sapply(xd,impossibleChangeNA))
}


##@checkZeroChanges checks for likelihood calculations
checkZeroChanges <- function(x)
{
	if (!inherits(x,'siena'))
	{
		stop('checkZeroChanges can only be applied to siena data objects')
	}
	zeroChange <- function(dv){
		dv[dv==10] <- 0
		dv[dv==11] <- 1
		ifelse((length(dim(dv))==3),
			1*any(sapply((2:dim(dv)[3]),
					FUN=function(i){all(dv[,,i-1] == dv[,,i], na.rm = TRUE)})), 0)
	}
	sum(sapply(x$depvars,zeroChange))
}
