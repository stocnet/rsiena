#/******************************************************************************
#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena
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
	    stop("\n To include a GMoM statistic use the function includeGMoMStatistics.")
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
			myeff2 <- myeff[use,]		
			print.sienaEffects(myeff2, includeOnly=FALSE, 
								includeRandoms=any(myeff$random & (myeff$shortName != 'density')),
								includeShortNames=TRUE)
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
				fix=FALSE, test=FALSE, random=FALSE,
				initialValue=0,
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
	if (hasArg("parameter"))
	{
		stop("includeInteraction should not mention a parameter; see the help page")
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
    ## interaction effects not yet implemented for continuous behavior
 #   if (any(myeff$netType[c(effect1, effect2, effect3)] == "continuous"))
 #   {
 #       stop("Interaction effects not yet implemented for continuous behavior")
 #   }
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
			c("unspInt", "behUnspInt", "contUnspInt") &
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
		myeff[intn, "randomEffects"] <- random
		myeff[intn, "initialValue"] <- initialValue
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
#
	myeff <- fixUpEffectNames(myeff)
	if (verbose)
	{
		intn2 <- intn
		intn[effect1] <- TRUE
		intn[effect2] <- TRUE
		if (effect3 > 0)
		{
			intn[effect3] <- TRUE
		}
		myeff2 <- myeff[intn,]	
		myeff2$effectName <- paste(which(intn), myeff2$effectName)
		print.sienaEffects(myeff2, includeOnly=FALSE, includeRandoms=random,
								includeShortNames=TRUE)
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
	replace.myeff <- function(meff, use, oldPar, newPar){
		meff[use, "effectName"] <- 
              gsub(oldPar, newPar, meff[use, "effectName"], fixed=TRUE)
		meff[use, "functionName"] <- 
              gsub(oldPar, newPar, meff[use, "functionName"], fixed=TRUE)	
		meff[use,]
	}
	if (!character)
	{
		shortName <- deparse(substitute(shortName))
	}
	if (shortName=="AltsAvAlt")
	{
		stop("Effect AltsAvAlt renamed to avXAlt.")
	}
	if (type=="gmm")
	{
	    stop("\n To include a GMoM statistic use the function includeGMoMStatistics.")
	}
	use <- myeff$shortName == shortName &
			myeff$name == name &
			myeff$type == type &
			myeff$interaction1 == interaction1 &
			myeff$interaction2 == interaction2 &
			(is.na(myeff$period) | myeff$period == period) &
			myeff$group == group
	if (shortName %in% c("unspInt", "behUnspInt", "contUnspInt"))
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
		if (shortName %in% c("unspInt", "behUnspInt", "contUnspInt"))
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
	if (!is.null(parameter)) 
	{
		olderParameter <- myeff[use, "parm"]
		myeff[use, "parm"] <- parameter
	}
	myeff[use, "include"] <- include
	myeff[use, "fix"] <- fix
	myeff[use, "test"] <- test
	myeff[use, "initialValue"] <- initialValue
	myeff[use, "timeDummy"] <- timeDummy
	myeff[use, "randomEffects"] <- random
	if (grepl("#", myeff[use, "effectName"], fixed=TRUE))
	{
# This means the original names in row[use,] in allEffects.csv were not yet changed
		myeff[use, "effectName"] <- 
              gsub("#", myeff[use, "parm"], myeff[use, "effectName"], fixed=TRUE)
		myeff[use, "functionName"] <- 
              gsub("#", myeff[use, "parm"], myeff[use, "functionName"], fixed=TRUE)
	}
	else if (!is.null(parameter))
    {
# replace what originally were strings "1/#" or "(#)" or "= #" or  "+ #" or "- #" or "#-" 
# Other original uses of # will not be replaced!
		oldParameter <- paste("1/", olderParameter, sep="")
		newParameter <- paste("1/", myeff[use, "parm"], sep="")
		myeff[use,] <- replace.myeff(myeff, use, oldParameter, newParameter)
		oldParameter <- paste("(", olderParameter, ")", sep="")
		newParameter <- paste("(", myeff[use, "parm"], ")", sep="")
		myeff[use,] <- replace.myeff(myeff, use, oldParameter, newParameter)
		oldParameter <- paste("= ", olderParameter, sep="")
		newParameter <- paste("= ", myeff[use, "parm"], sep="")
		myeff[use,] <- replace.myeff(myeff, use, oldParameter, newParameter)
		oldParameter <- paste("+ ", olderParameter, sep="")
		newParameter <- paste("+ ", myeff[use, "parm"], sep="")
		myeff[use,] <- replace.myeff(myeff, use, oldParameter, newParameter)
		oldParameter <- paste("- ", olderParameter, sep="")
		newParameter <- paste("- ", myeff[use, "parm"], sep="")
		myeff[use,] <- replace.myeff(myeff, use, oldParameter, newParameter)
		oldParameter <- paste(olderParameter, "-", sep="")
		newParameter <- paste(myeff[use, "parm"], "-", sep="")
		myeff[use,] <- replace.myeff(myeff, use, oldParameter, newParameter)
	}
	if (verbose)
	{		
		myeff2 <- myeff[use,]		
		print.sienaEffects(myeff2, includeOnly=FALSE, includeRandoms=random,
								includeShortNames=TRUE)
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


fixUpEffectNames <- function(effects)
##@fixUpEffectNames siena07 Replace # and construct interaction names
{
    ## replace # by the parm value in function and effect names:
    effects$effectName <-
        sapply(1:nrow(effects), function(x, y)
           {
               y <- y[x, ]
               gsub("#", y$parm, y$effectName)
           }, y=effects)
    effects$functionName <-
        sapply(1:nrow(effects), function(x, y)
           {
               y <- y[x, ]
               gsub("#", y$parm, y$functionName)
           }, y=effects)
    ##validate user-specified network interactions
    interactions <- effects[effects$shortName == "unspInt" & effects$include &
                            effects$effect1 > 0 , ]
    if (nrow(interactions) > 0)
    {
        unspIntNames <-
            sapply(1:nrow(interactions), function(x, y, z)
               {
                   y <- y[x, ] ## get the interaction effect
                   twoway <- y$effect3 == 0
                   ## now get the rows which are to interact
                   inter1 <- z[z$effectNumber == y$effect1, ]
                   if (nrow(inter1) != 1 )
                   {
                       stop("invalid network interaction specification: ",
                            "effect number 1")
                   }
                   inter2 <- z[z$effectNumber == y$effect2, ]
                   if (nrow(inter2) != 1 )
                   {
                       stop("invalid network interaction specification: ",
                            "effect number 2")
                   }
                   if (!twoway)
                   {
                       inter3 <- z[z$effectNumber == y$effect3, ]
                       if (nrow(inter3) != 1)
                       {
                           stop("invalid network interaction specification: ",
                                "effect number 3")
                       }
                   }
                   else
                   {
                       inter3 <- z[is.na(z$effectNumber), ]
                       ## should be empty row
                   }
                   if (twoway)
                   {
                       if (inter1$name != inter2$name)
                       {
                           stop("invalid network interaction specification: ",
                                "must all be same network")
                       }
                       if (inter1$type != inter2$type)
                       {
                           stop("invalid network interaction specification: ",
                                "must all be same type: ",
                                "evaluation, endowment or creation")
                       }
                   }
                   else
                   {
                       if (inter1$name != inter2$name ||
                           inter1$name != inter3$name)
                       {
                           stop("invalid network interaction specification: ",
                                "must all be same network")
                       }
                       if (inter1$type != inter2$type ||
                           inter1$type != inter3$type)
                       {
                           stop("invalid network interaction specification: ",
                                "must all be ",
                                "same type: evaluation, endowment or creation ")
                       }
                   }
                   ## check types
                   inters <- rbind(inter1, inter2, inter3)
                   egos <- which(inters$interactionType == "ego")
                   egoCount <- length(egos)
                   dyads <- which(inters$interactionType == "dyadic")
                   dyadCount <- length(dyads)
                   if (twoway)
                   {
                       if (egoCount < 1 && dyadCount != 2)
                       {
                           stop("invalid network interaction specification: ",
                                "must be at least one ego or both dyadic ",
                                "effects")
                       }
                   }
                   else
                   {
                       if (egoCount < 2 && (egoCount + dyadCount < 3))
                       {
                      stop("invalid network 3-way interaction specification: ",
									"must be at least two ego effects ",
									"or all ego or dyadic effects")
                       }
                   }
                   ## construct a name
                   ## make sure the egos are at the front of inters
                   if (egoCount > 0)
                   {
                       inters <- rbind(inters[egos, ], inters[-egos, ])
                   }
 				   tmpnames <- inters$effectName
				   tmpnames[-1] <- sub(paste(inters$name[1], ": ",
											 sep=""), "", tmpnames[-1])
				   tmpname <- paste(tmpnames, collapse = " x ")
# following lines dropped, might be restored if desired
#                   if (twoway && nchar(tmpname) < 38)
#                   {
#                       tmpname <- paste("int. ", tmpname)
#                   }
#                   if (!twoway)
#                  {
#                      tmpname <- paste("i3.", tmpname)
#                  }
                   tmpname
               }, y=interactions, z=effects)
        effects[effects$shortName == "unspInt" & effects$include &
                !is.na(effects$effect1), c("effectName", "functionName")] <-
                    unspIntNames
    }
    ##validate user-specified behavior interactions
    interactions <- effects[effects$shortName == "behUnspInt" &
                            effects$include &
                            effects$effect1 > 0 , ]
    if (nrow(interactions) > 0)
    {
        unspIntNames <-
            sapply(1:nrow(interactions), function(x, y, z)
               {
                   y <- y[x, ] ## get the interaction effect
                   twoway <- y$effect3 == 0
                   ## now get the rows which are to interact
                   inter1 <- z[z$effectNumber == y$effect1, ]
                   if (nrow(inter1) != 1 )
                   {
                       stop("invalid behavior interaction specification: ",
                            "effect number 1")
                   }
                   inter2 <- z[z$effectNumber == y$effect2, ]
                   if (nrow(inter2) != 1 )
                   {
                       stop("invalid behavior interaction specification: ",
                            "effect number 2")
                   }
                   if (!twoway)
                   {
                       inter3 <- z[z$effectNumber == y$effect3, ]
                       if (nrow(inter3) != 1)
                       {
                           stop("invalid behavior interaction specification: ",
                                "effect number 3")
                       }
                   }
                   else
                   {
                       inter3 <- z[is.na(z$effectNumber), ]
                       ## should be empty row
                   }
                   if (twoway)
                   {
                       if (inter1$name != inter2$name)
                       {
                           stop("invalid behavior interaction specification: ",
                                "must all be same behavior variable")
                       }
                       if (inter1$type != inter2$type)
                       {
                           stop("invalid behavior interaction specification: ",
                                "must be same type: evaluation, endowment ",
                                "or creation")
                       }
                   }
                   else
                   {
                       if (inter1$name != inter2$name ||
                           inter1$name != inter3$name)
                       {
                           stop("invalid behavior interaction specification: ",
                                "must all be same behavior variable")
                       }
                       if (inter1$type != inter2$type ||
                           inter1$type != inter3$type)
                       {
                           stop("invalid behavior interaction specification: ",
                                "must all be ",
                                "same type: evaluation, endowment or creation")
                       }
                   }
                   ## check types - at most one should be not OK here
                   inters <- rbind(inter1, inter2, inter3)
                   if (length(which(inters$interactionType != "OK")) > 1)
				   {
				   	   stop("invalid behavior interaction specification: ",
				   			"at most one effect with interactionType ",
				   			"not OK is allowed")
                   }
                   ##if (any(inters$interactionType != "OK"))
                   ##{
                   ##    stop("invalid behavior interaction specification: ",
                   ##         "only effects with interactionType OK are allowed")
                   ##}
                   ## construct a name
				   tmpnames <- inters$effectName
				   tmpnames[-1] <- sub(paste("behavior ", inters$name[1], " ",
											 sep=""), "", tmpnames[-1])
				   tmpname <- paste(tmpnames, collapse = " x ")
                   if (twoway && nchar(tmpname) < 38)
                   {
                       tmpname <- paste("int. ", tmpname)
                   }
                   if (!twoway)
                   {
                       tmpname <- paste("i3.", tmpname)
                   }
                   tmpname
               }, y=interactions, z=effects)
        effects[effects$shortName == "behUnspInt" & effects$include &
                !is.na(effects$effect1), c("effectName", "functionName")] <-
                    unspIntNames
    }
    ##validate user-specified continuous interactions
    interactions <- effects[effects$shortName == "contUnspInt" &
                            effects$include &
                            effects$effect1 > 0 , ]
    if (nrow(interactions) > 0)
    {
        unspIntNames <-
            sapply(1:nrow(interactions), function(x, y, z)
               {
                   #browser()
				   y <- y[x, ] ## get the interaction effect
                   twoway <- y$effect3 == 0
                   ## now get the rows which are to interact
                   inter1 <- z[z$effectNumber == y$effect1, ]
                   if (nrow(inter1) != 1 )
                   {
                       stop("invalid behavior interaction specification: ",
                            "effect number 1")
                   }
                   inter2 <- z[z$effectNumber == y$effect2, ]
                   if (nrow(inter2) != 1 )
                   {
                       stop("invalid behavior interaction specification: ",
                            "effect number 2")
                   }
                   if (!twoway)
                   {
                       inter3 <- z[z$effectNumber == y$effect3, ]
                       if (nrow(inter3) != 1)
                       {
                           stop("invalid behavior interaction specification: ",
                                "effect number 3")
                       }
                   }
                   else
                   {
                       inter3 <- z[is.na(z$effectNumber), ]
                       ## should be empty row
                   }
                   if (twoway)
                   {
                       if (inter1$name != inter2$name)
                       {
                           stop("invalid behavior interaction specification: ",
                                "must all be same behavior variable")
                       }
                       if (inter1$type != inter2$type)
                       {
                           stop("invalid behavior interaction specification: ",
                                "must be same type: evaluation")
                       }
                   }
                   else
                   {
                       if (inter1$name != inter2$name ||
                           inter1$name != inter3$name)
                       {
                           stop("invalid behavior interaction specification: ",
                                "must all be same behavior variable")
                       }
                       if (inter1$type != inter2$type ||
                           inter1$type != inter3$type)
                       {
                           stop("invalid behavior interaction specification: ",
                                "must all be same type: evaluation")
                       }
                   }
                   ## check types - all should be OK here
                   inters <- rbind(inter1, inter2, inter3)
                   if (any(inters$interactionType != "OK"))
				   {
				   	   stop("invalid behavior interaction specification: ",
				   			"only effects with interactionType OK are allowed")
                   }
                   ## construct a name
				   tmpnames <- inters$effectName
				   tmpnames[-1] <- sub(paste("behavior ", inters$name[1], " ",
											 sep=""), "", tmpnames[-1])
				   tmpname <- paste(tmpnames, collapse = " x ")
                   if (twoway && nchar(tmpname) < 38)
                   {
                       tmpname <- paste("int. ", tmpname)
                   }
                   if (!twoway)
                   {
                       tmpname <- paste("i3.", tmpname)
                   }
                   tmpname
               }, y=interactions, z=effects)
        effects[effects$shortName == "contUnspInt" & effects$include &
                !is.na(effects$effect1), c("effectName", "functionName")] <-
                    unspIntNames
    }
    effects
}
