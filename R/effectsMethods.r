##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: https://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: effectsMethods.r
## *
## * Description: This file contains the print and edit methods for the
## * class sienaEffects
## *
## ****************************************************************************/
##@print.sienaEffects Methods
print.sienaEffects <- function(x, fileName=NULL, includeOnly=TRUE, 
    expandDummies=FALSE, includeRandoms=FALSE, dropRates=FALSE, includeShortNames=FALSE, ...)
{
    if (!inherits(x, "sienaEffects"))
        stop("not a legitimate Siena effects object")

    if (typeof(fileName)=="character")
    {
        sink(fileName, split=TRUE)
    }

    interactions <- x[x$shortName %in% c("unspInt", "behUnspInt", "contUnspInt") & x$include &
        x$effect1 > 0, ]
    if (expandDummies)
    {
        if(includeOnly && !all(x[x$include, "timeDummy"] == ",")
            || !all(x[, "timeDummy"] == "," ))
        {
            x <- sienaTimeFix(x)$effects
#            x <- fixUpEffectNames(x)
        }
        else
        {
            if (nrow(interactions) > 0)
            {
#                x <- fixUpEffectNames(x)
            }
        }
    }
    if (nrow(x) > 0)
    {
        nDependents <- length(unique(x$name))
        userSpecifieds <- x$shortName[x$include] %in%
            c("unspInt", "behUnspInt", "contUnspInt")
        endowments <- !x$type[x$include] %in% c("rate", "eval", "gmm")
        # includes creations and gmm
        gmm <- any(x$type[x$include] %in% "gmm")
        timeDummies <- !x$timeDummy[x$include] == ","
		if (includeShortNames)
		{
			specs <- as.data.frame(x[, c("name", "effectName", "shortName", "include", "fix",
                "test", "initialValue", "parm")])		
		}
		else
		{
			specs <- as.data.frame(x[, c("name", "effectName", "include", "fix",
                "test", "initialValue", "parm")])
		}
        if (includeOnly)
        {
            included <- x$include
        }
        else
        {
            included <- rep(TRUE, nrow(x))
        }
        if (dropRates)
        {
            included <- included & !x$basicRate
        }
        specs <- specs[included, ]
        if (nrow(specs) > 0)
        {
            row.names(specs) <- 1:nrow(specs)
        }
        if (nDependents == 1)
        {
            specs <- specs[, -1]  # drop name of dependent variable
        }
        if (any(endowments))
        {
            specs <- cbind(specs, type=x[included, "type"])
        }
        if (gmm)
        {
          specs <- cbind(specs, type=x[included, "type"])
        }
        if (any(timeDummies))
        {
            specs <- cbind(specs, timeDummy=x[included, "timeDummy"])
        }
        if (any(userSpecifieds))
        {
            specs <- cbind(specs, x[included, c("effect1", "effect2")])
            if (any (x$effect3[included] > 0))
            {
                specs <- cbind(specs, effect3=x[included, "effect3"])
            }
        }
        if (includeRandoms)
        {
            specs <- cbind(specs, randomEffects=x[included, "randomEffects"])
        }
        specs[, "initialValue"] <- format(round(specs$initialValue,digits=5),
            width=10)
        if (nrow(specs) > 0)
        {
          if (gmm)
              {
                cat('\n Effects and statistics for estimation by the Generalized Method of Moments\n')
                if (any(specs$type!='gmm'))
                {
                  cat('\n Effects\n')
                  specs1 <- specs[specs$type!='gmm',]
                  row.names(specs1) <- 1:nrow(specs1)
                  print(as.matrix(specs1), quote=FALSE)
                }
                if (any(specs$type=='gmm'))
                {
                  cat('\n Regular and GMoM statistics\n')
                  specs2 <- specs[specs$type=='gmm',]
                  specs3 <- rbind(specs1, specs2)
                  specs3$Statistic <- c(rep("Regular",nrow(specs1)),
                                        rep("GMoM",nrow(specs2)))
                  row.names(specs3) <- 1:nrow(specs3)
                  print(as.matrix(specs3[,c(1:2,9)]), quote=FALSE)
                }
              }
            else
            {
                print(as.matrix(specs), quote=FALSE)
            }
        }
        else
        {
            print.data.frame(specs)
        }
    }
    else
    {
        print.data.frame(x)
    }
    if (includeRandoms)
    {
        nfeff <- sum((!x$randomEffects) & x$include & (!x$fix) & (!x$basicRate))
        nreff <- sum(x$randomEffects & x$include & (!x$fix) & (!x$basicRate))
        nrate <- sum(x$basicRate & x$include & (x$group==2) & (!x$fix))
        if (sum(x$group != 1) > 0) # else there is only one group, and counting should be different.
        {
            cat('Length of priorSigEta for sienaBayes, if used, should be',
                nfeff, '.\n')
            cat('Dimensions of priorMu and priorSigma for sienaBayes should be',
                nreff, '+', nrate, '=', nreff+nrate,'.\n')
        }
        else if (dim(x)[1] > 1)
        # to avoid the following if called from includeInteraction or setEffect
        {
            cat('There are', nreff, 'random effects.\n')
        }
    }
    if (typeof(fileName)=="character")
    {
        sink()
    }
    invisible(x)
}

##@summary.sienaEffects Methods
summary.sienaEffects <- function(object, fileName=NULL, includeOnly=TRUE,
    expandDummies=FALSE, ...)
{
    if (!inherits(object, "sienaEffects"))
        stop("not a legitimate Siena effects object")
    if (expandDummies && (includeOnly && !all(object[object$include,
                "timeDummy"] == ",")
            || !all(object[, "timeDummy"] == ",")))
    {
        object <- sienaTimeFix(object)$effects
    }
    if (includeOnly)
    {
        object <- object[object$include, ]
    }
    class(object) <- c("summary.sienaEffects", class(object))
    object
}

##@print.summary.sienaEffects Methods
print.summary.sienaEffects <- function(x, fileName=NULL, ...)
{
    if (!inherits(x, "summary.sienaEffects"))
        stop("not a legitimate summary of a Siena effects object")
    ## find out if any columns need removing because they will not print

    problem <- sapply(x, function(x) {
        inherits(try(unlist(x[[1]]), silent=TRUE), "try-error")
    })
    if (any(problem))
    {
        x1 <- x[, problem]
        x <- x[, !problem]
    }
    if (typeof(fileName)=="character")
    {
        sink(fileName, split=TRUE)
    }
    print.data.frame(x)
    if (any(problem))
    {
        lapply(x1, print)
    }
    if (typeof(fileName)=="character")
    {
        sink()
    }
    invisible(x)
}

edit.sienaEffects <- function(name, ...)
{
    if (!interactive())
    {
        return(name)
    }
    ## store the original column order
    originalNames <- names(name)
    ## move function name and other things out of the way
    priorityColumns <- c("name", "effectName", "type", "include", "fix",
        "test", "initialValue", "parm", "shortName")
    priorityX <- name[, priorityColumns]
    notPriorityX <- name[, -c(match(priorityColumns, names(name)))]
    name <- cbind(priorityX, notPriorityX)

    ## get edit.data.frame to do the actual edit
    tmp <- NextMethod(, , edit.row.names=FALSE)

    ## re-sort the columns
    tmp <- tmp[, match(originalNames, names(tmp))]
    class(tmp) <- c("sienaEffects", class(tmp))
	 attr(tmp, "version") <- packageDescription(pkgname, fields = "Version")
    tmp
}

##@updateSpecification Methods add specified effects from other effects object
updateSpecification <- function(effects.to, effects.from,
									effects.extra=NULL, name.to=NULL, name.from=NULL)
{
    if (!inherits(effects.to, "sienaEffects"))
    {
        stop("effects.to is not an effects object")
    }
    if (!inherits(effects.from, "sienaEffects"))
    {
        stop("effects.from is not an effects object")
    }
    if (is.null(name.from))
    {
        prevEffects <-
            effects.from[which((effects.from$type != 'gmm')&((effects.from$include))),]
    }
    else
    {
        if (!is.character(name.from))
        {
            stop("name.from should be a string")
        }
        prevEffects <-
            effects.from[which((effects.from$type != 'gmm')&((effects.from$include))&
                (effects.from$name == name.from)),]
    }
    oldlist <- apply(prevEffects, 1, function(x)
                     paste(x[c("name", "shortName",
                               "type", "groupName",
                               "interaction1", "interaction2",
                               "period", "effect1", "effect2", "effect3")],
                           collapse="|"))
    efflist <- apply(effects.to, 1, function(x)
                     paste(x[c("name", "shortName",
                               "type", "groupName",
                               "interaction1", "interaction2",
                               "period", "effect1", "effect2", "effect3")],
                            collapse="|"))
    if (!(is.null(name.to)))
    {
        if (!is.character(name.to))
        {
            stop("name.to should be a string")
        }
        else if (name.to == "all") # omit matching on "name"
        {
            oldlist <- apply(prevEffects, 1, function(x)
                     paste(x[c("shortName",
                               "type", "groupName",
                               "interaction1", "interaction2",
                               "period", "effect1", "effect2", "effect3")],
                           collapse="|"))
            efflist <- apply(effects.to, 1, function(x)
                     paste(x[c("shortName",
                               "type", "groupName",
                               "interaction1", "interaction2",
                               "period", "effect1", "effect2", "effect3")],
                            collapse="|"))
        }
    }
    use <- (efflist %in% oldlist)
    correspondence <- match(efflist, oldlist)
    if (!(is.null(name.to)))
    {
        if (name.to != "all")
        {
            use <- (use & (effects.to$name == name.to))
        }
    }
    effects.to$include[use] <- TRUE
    effects.to$fix[use] <- prevEffects$fix[correspondence][use]
    effects.to$test[use] <- prevEffects$test[correspondence][use]
    effects.to$parm[use] <- prevEffects$parm[correspondence][use]
    effects.to$randomEffects[use] <- prevEffects$randomEffects[correspondence][use]
    effects.to$initialValue[use] <- prevEffects$initialValue[correspondence][use]
# the above does not transfer interaction effects.
# A lot of work is needed to get the information about the interacting effects.
    inter <- which(prevEffects$include &
						(prevEffects$shortName %in% c("unspInt","behUnspInt","contUnspInt")))
    if (length(inter) >= 1)
	{
	# look up the interacting main effects, and try to get information
	# about which are the interacting effects.
	    cat("Handling interactions ... ")
		efn1 <- prevEffects$effect1[inter]
		efn2 <- prevEffects$effect2[inter]
		efn3 <- prevEffects$effect3[inter]
	# prepare a stopmessage; this will possibly be used at various places.
		if (is.null(effects.extra))
		{
			stopMessage <- paste("Effects object ", deparse(substitute(effects.from)),
					" contains some interactions \n",
					" but there is no information ",
					"for the corresponding main effects.", sep="")
		}
		else
		{
			stopMessage <- paste("Effects object ", deparse(substitute(effects.from)),
					" contains some interactions \n",
					"  and neither this nor ", deparse(substitute(effects.extra)),
					"\n contains information for the corresponding main effects.", sep="")
		}
	# Note that some of efn3 may be 0,
	# and effects.from$effectNumber starts counting from 1.
	# First try to find the corresponding main effects in effects.from
		three <- (efn3 > 0)
		mefn1 <- match(efn1, effects.from$effectNumber)
		mefn2 <- match(efn2, effects.from$effectNumber)
		mefn3 <- ifelse(three, match(efn3, effects.from$effectNumber), 0)
		effects.fr <- effects.from
		if (any(is.na(c(mefn1,mefn2,mefn3))))
		{
			if (is.null(effects.extra))
			{
				stop("Effects object ", deparse(substitute(effects.from)),
					" contains some interactions \n",
					"  without information for the corresponding main effects,\n",
					"  and there is no effects.extra.")
			}
			else
			{
				if (!inherits(effects.extra, "sienaEffects"))
				{
					stop("effects.extra is needed, and is not an effects object")
				}
	# Now try to find the corresponding main effects in effects.extra
				version1 <- attr(effects.from, "version")
				version2 <- attr(effects.extra, "version")
				sameversion <- FALSE
				if ((!is.null(version1)) & (!is.null(version2)))
				{
					sameversion <- (version1==version2)
				}
				if (!sameversion)
				{
					warning("RSiena versions of effects.from and effects.extra",
						" are different;\n",
						"  check that the effects object",
						" generated by updateSpecification is correct.")
				}
				mefn1 <- match(efn1, effects.extra$effectNumber)
				mefn2 <- match(efn2, effects.extra$effectNumber)
				mefn3 <- ifelse(three, match(efn3, effects.extra$effectNumber), 0)
				if (any(is.na(c(mefn1,mefn2,mefn3))))
				{
					stop(stopMessage)
				}
				effects.fr <- effects.extra
			}
		}
        shn1  <- effects.fr[mefn1,"shortName"]
        shn2  <- effects.fr[mefn2,"shortName"]
		shn3  <- rep('', length(three))
        shn3[three]  <- effects.fr[mefn3[three],"shortName"]
        int11  <- effects.fr[mefn1,"interaction1"]
        int12  <- effects.fr[mefn2,"interaction1"]
		int13  <- rep('', length(three))
        int13[three]  <- effects.fr[mefn3[three],"interaction1"]
        int21  <- effects.fr[mefn1,"interaction2"]
        int22  <- effects.fr[mefn2,"interaction2"]
		int23  <- rep('', length(three))
        int23[three]  <- effects.fr[mefn3[three],"interaction2"]
        nam  <- prevEffects[inter,"name"]
        typ   <- prevEffects[inter,"type"]
        fixx  <- prevEffects[inter,"fix"]
        tests <- prevEffects[inter,"test"]
        rand  <- prevEffects[inter,"randomEffects"]
		initv <- prevEffects[inter,"initialValue"]
#browser()
        for (k in seq_along(inter)){
			cat(k," ")
			flush.console()
			if (three[k])
			{
				if (inherits(try(
						effects.to <- includeInteraction(effects.to, shn1[k], shn2[k], shn3[k],
						name=nam[k],
						interaction1=c(int11[k],int12[k],int13[k]),
						interaction2=c(int21[k],int22[k],int23[k]),
						initialValue=initv[k],
						type=typ[k], fix=fixx[k], test=tests[k], random=rand[k],
						verbose=FALSE, character=TRUE), 
					silent=TRUE), "try-error"))
				{
					stop(stopMessage)
				}
			}
			else
			{
				if (inherits(try(
					effects.to <- includeInteraction(effects.to, shn1[k], shn2[k],
						name=nam[k],
						interaction1=c(int11[k],int12[k]),
						interaction2=c(int21[k],int22[k]),
						initialValue=initv[k],
						type=typ[k], fix=fixx[k], test=tests[k], random=rand[k],
						verbose=FALSE, character=TRUE), 
					silent=TRUE), "try-error"))
				{
					stop(stopMessage)
				}
			}
        }
		cat("\n")
    }
    effects.to
}
