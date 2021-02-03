#/******************************************************************************
#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: sienaeffects.r
# *
# * Description: This module contains utilities for updating an effects object
# *****************************************************************************/
##@includeEffect DataCreate
includeEffects <- function(myeff, ..., include=TRUE, name=myeff$name[1],
						   type="eval", interaction1="", interaction2="",
						   fix=FALSE, test=FALSE,
						   character=FALSE, verbose=TRUE)
{
	if (!inherits(myeff, 'sienaEffects'))
	{
		stop("The first argument is not of class <sienaEffects>.")
	}
	if (character)
	{
		dots <- sapply(list(...), function(x)x)
	}
	else
	{
		dots <- substitute(list(...))[-1] ##first entry is the word 'list'
	}
	if (length(dots) == 0)
	{
		stop("This function needs some effect short names.")
	}
	if (!character)
	{
		effectNames <- sapply(dots, function(x)deparse(x))
	}
	else
	{
		effectNames <- dots
	}
	if ("AltsAvAlt" %in% effectNames)
	{
		stop("Effect AltsAvAlt now is called avXAlt. Use the new name please.")
	}
	use <- myeff$shortName %in% effectNames &
	myeff$type==type &
	myeff$name==name &
	myeff$interaction1 == interaction1 &
	myeff$interaction2 == interaction2
	myeff[use, "include"] <- include
	myeff[use, "test"] <- test
	myeff[use, "fix"] <- fix
  	if (sum(myeff[use, "type"]=="gmm") > 0)
  	{
	    stop("\n To include a GMoM statistic use the function includeGMoMStatistics")
	}
	if (sum(use) <= 0)
	{
		cat(paste("There is no effect with short name "))
		cat(paste(effectNames,", \n", sep=""))
		cat(paste("and with interaction1 = <",interaction1,">, ", sep=""))
		cat(paste("interaction2 = <",interaction2,">, ", sep=""))
		cat(paste("and type = <",type,">, \n", sep=""))
		cat(paste("for dependent variable",name,".\n"))
		cat("See effectsDocumentation() for this effects object.\n")
	}
	else
	{
#		print.data.frame(myeff[use, c("name", "shortName", "type",
#			"interaction1", "interaction2", "include")])
		if (verbose)
		{
			print.sienaEffects(myeff[use,], includeRandoms = 
					any(myeff$random & (myeff$shortName != 'density')))
		}
	}
	if (hasArg('initialValue'))
	{
		warning("Warning: argument 'initialValue' has no effect in includeEffects; use setEffect.\n")
	}
	if (hasArg('parameter'))
	{
		warning("Warning: argument 'parameter' has no effect in includeEffects; use setEffect.\n")
	}
	if (hasArg('random'))
	{
		warning("Warning: argument 'random' has no effect in includeEffects; use setEffect.\n")
	}
	myeff
}

##@includeInteraction DataCreate
includeInteraction <- function(myeff, ...,
				include=TRUE, name=myeff$name[1],
				type="eval", interaction1=rep("", 3), interaction2=rep("", 3),
				fix=FALSE, test=FALSE, parameter=NULL, random=FALSE,
				character=FALSE, verbose=TRUE)
{
	if (character)
	{
		dots <- sapply(list(...), function(x)x)
	}
	else
	{
		## check we have 2 or 3 short names
		dots <- substitute(list(...))[-1] ##first entry is the word 'list'
	}
	if (length(dots) == 0)
	{
		stop("need some effect short names")
	}
	if (length(dots) < 2 || length(dots) > 3)
	{
		 stop("need exactly two or three effect short names")
	}
	if (!character)
	{
		shortNames <- sapply(dots, function(x)deparse(x))
	}
	else
	{
		shortNames <- dots
	}
	## find the first underlying effect
	shortName <- shortNames[1]
	interact1 <- interaction1[1]
	interact2 <- interaction2[1]
	use <- myeff$shortName == shortName &
	myeff$type==type &
	myeff$name==name &
	myeff$interaction1 == interact1 &
	myeff$interaction2 == interact2
	if (sum(use) == 0)
	{
		stop("First effect not found")
	}
	if (sum(use) > 1)
	{
		stop("First effect not unique")
	}
	effect1 <- myeff[use, "effectNumber"]
	## find the second underlying effect
	shortName <- shortNames[2]
	interact1 <- ifelse (length(interaction1) > 1, interaction1[2], "")
	interact2 <- ifelse (length(interaction2) > 1, interaction2[2], "")
	use <- myeff$shortName == shortName &
	myeff$type==type &
	myeff$name==name &
	myeff$interaction1 == interact1 &
	myeff$interaction2 == interact2
	if (sum(use) == 0)
	{
		stop("Second effect not found")
	}
	if (sum(use) > 1)
	{
		stop("Second effect not unique")
	}
	effect2 <- myeff[use, "effectNumber"]
	## find the third underlying effect, if any
	if (length(shortNames) > 2)
	{
		shortName <- shortNames[3]
		interact1 <- ifelse (length(interaction1) > 2, interaction1[3], "")
		interact2 <- ifelse (length(interaction2) > 2, interaction2[3], "")
		use <- myeff$shortName == shortName &
		myeff$type==type &
		myeff$name==name &
		myeff$interaction1 == interact1 &
		myeff$interaction2 == interact2
		if (sum(use) == 0)
		{
			stop("Third effect not found")
		}
		if (sum(use) > 1)
		{
			stop("Third effect not unique")
		}
		effect3 <- myeff[use, "effectNumber"]
	}
	else
	{
		effect3 <- 0
	}
	## does the effect already exist?
	intn <- (myeff$effect1 == effect1) & (myeff$effect2 == effect2)
	if (effect3 > 0)
	{
		intn <- intn & (myeff$effect3 == effect3)
	}
	# intn indicates the rows of the effects object with this interaction
	if (sum(intn) >= 2)
	{
		cat('The given effects object is corrupted:\n')
		cat('It already contains more than one copy of this interaction.\n')
		cat('Make a new effects object from scratch.\n')
		stop('Corrupted effects object')
	}

	## if want to include, check that we have a spare row
	if ((include) && (sum(intn) == 0))
	{
		# The interaction must be created
		ints <- myeff[myeff$name == name & myeff$shortName %in%
			c("unspInt", "behUnspInt") &
			(is.na(myeff$effect1) | myeff$effect1 == 0)&
			myeff$type == type, ]
		if (nrow(ints) == 0)
		{
			baseEffect<- myeff[myeff$name == name, ][1, ]
			if (baseEffect$netType == "behavior")
			{
				stop("Use getEffects() with a larger value for behNintn.")
			}
			else
			{
				stop("Use getEffects() with a larger value for nintn.")
			}
		}
		ints <- ints[1, ]
		intn <- myeff$effectNumber == ints$effectNumber
	}
	if (include)
	{
		myeff[intn, "include"] <- include
		myeff[intn, c("effect1", "effect2", "effect3")] <-
			c(effect1, effect2, effect3)
		myeff[intn, "fix"] <- fix
		myeff[intn, "test"] <- test
		if (!is.null(parameter)) {myeff[intn, "parm"] <- parameter}
		myeff[intn, "randomEffects"] <- random
	}
	else
	{
		if (sum(intn) == 0)
		{
		warning('Note: there was no such interaction in this effects object.')
		}
		else
		{
			myeff[intn, "include"] <- FALSE
		}
	}
	if (verbose)
	{
		print.sienaEffects(myeff[intn,], includeRandoms = 
					any(myeff$random & (myeff$shortName != 'density')))
	}
	myeff
}

##@setEffect DataCreate
setEffect <- function(myeff, shortName, parameter=NULL,
					fix=FALSE, test=FALSE, random=FALSE,
					initialValue=0,
					timeDummy=",",
					include=TRUE, name=myeff$name[1],
					type="eval", interaction1="", interaction2="",
					effect1=0, effect2=0, effect3=0,
					period=1, group=1, character=FALSE, verbose=TRUE)
{
	if (!character)
	{
		shortName <- deparse(substitute(shortName))
	}
	if (shortName=="AltsAvAlt")
	{
		stop("Effect AltsAvAlt renamed to avXAlt.")
	}
	use <- myeff$shortName == shortName &
	myeff$name == name &
	myeff$type == type &
	myeff$interaction1 == interaction1 &
	myeff$interaction2 == interaction2 &
	(is.na(myeff$period) | myeff$period == period) &
	myeff$group == group
	if (shortName %in% c("unspInt", "behUnspInt"))
	{
		use <- use & (myeff$include) & (myeff$effect1 == effect1) &
			(myeff$effect2 == effect2) & (myeff$effect3 == effect3)
	}
	if (sum(use) == 0)
	{
		cat(paste("There is no effect with short name "))
		cat(paste(shortName,", \n", sep=""))
		cat(paste("and with interaction1 = <",interaction1,">, ", sep=""))
		cat(paste("interaction2 = <",interaction2,">, ", sep=""))
		cat(paste("type = <",type,">, ", sep=""))
		cat(paste("period = <",period,">, ", sep=""))
		if (shortName %in% c("unspInt", "behUnspInt"))
		{
		cat(paste("effects1-2-3 = <",effect1, effect2, effect3,">,", sep=" "))
		}
		cat(paste("and group = <",group,">, \n ", sep=""))
		cat(paste("for dependent variable",name,".\n"))
		stop("Effect not found")
	}
	if (sum(use) > 1)
	{
		stop("Effect not unique")
	}
	if (!is.null(parameter)) {myeff[use, "parm"] <- parameter}
	myeff[use, "include"] <- include
	myeff[use, "fix"] <- fix
	myeff[use, "test"] <- test
	myeff[use, "initialValue"] <- initialValue
	myeff[use, "timeDummy"] <- timeDummy
	myeff[use, "randomEffects"] <- random
	# print.data.frame(myeff[use, c("name", "shortName", "type", "interaction1",
	# 		"interaction2", "include", "parm", "fix", "test",
	# 		"initialValue", "timeDummy", "period", "group")])
	if (verbose)
	{
	  if (any(myeff$type[use] %in% "gmm"))
	  {
	    myeff2 <- as.data.frame(myeff[use,])
	    rownames(myeff2) <- 1:nrow(myeff2)
	    print.data.frame(myeff2[, c("name", "shortName", "type","include")])
	  }
	  else
	  {
	    print.sienaEffects(myeff[use,], includeRandoms = random)
	  }
	}
	myeff
}

##@includeGMoMStatistics DataCreate
includeGMoMStatistics <- function(myeff, ..., include=TRUE, name=myeff$name[1],
                                  interaction1="", interaction2="",
                                  character=FALSE, verbose=TRUE)
{
  if (character)
  {
    dots <- sapply(list(...), function(x)x)
  }
  else
  {
    dots <- substitute(list(...))[-1] ##first entry is the word 'list'
  }
  if (length(dots) == 0)
  {
    stop("This function needs some effect short names.")
  }
  if (!character)
  {
    effectNames <- sapply(dots, function(x)deparse(x))
  }
  else
  {
    effectNames <- dots
  }
  use <- myeff$shortName %in% effectNames &
    myeff$type=="gmm" &
    myeff$name==name &
    myeff$interaction1 == interaction1 &
    myeff$interaction2 == interaction2
  myeff[use, "include"] <- include
  myeff[use, "test"] <- FALSE
  myeff[use, "fix"] <- TRUE
  myeff[use, "initialValue"] <- 0
  if (sum(use) <= 0)
  {
    cat(paste("There is no GMoM statistic with short name "))
    cat(paste(effectNames,", \n", sep=""))
    cat(paste("and with interaction1 = <",interaction1,">, ", sep=""))
    cat(paste("interaction2 = <",interaction2,">, ", sep=""))
    cat(paste("for dependent variable",name,".\n"))
  }
  else
  {
    #		print.data.frame(myeff[use, c("name", "shortName", "type",
    #			"interaction1", "interaction2", "include")])
    if (verbose)
    {
      myeff2 <- as.data.frame(myeff[use,])
      rownames(myeff2) <- 1:nrow(myeff2)
      print.data.frame(myeff2[, c("name", "shortName", "type","include")])
    }
  }
  myeff
}
