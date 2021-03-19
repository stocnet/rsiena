##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: effectsMethods.r
## *
## * Description: This file contains the print and edit methods for the
## * class sienaEffects
## *
## ****************************************************************************/
##@print.sienaEffects Methods
print.sienaEffects <- function(x, fileName=NULL, includeOnly=TRUE,
    expandDummies=FALSE, includeRandoms=FALSE, dropRates=FALSE, ...)
{
    if (!inherits(x, "sienaEffects"))
        stop("not a legitimate Siena effects object")

    if (typeof(fileName)=="character")
    {
        sink(fileName, split=TRUE)
    }

    interactions <- x[x$shortName %in% c("unspInt", "behUnspInt") & x$include &
        x$effect1 > 0, ]
    if (expandDummies)
    {
        if(includeOnly && !all(x[x$include, "timeDummy"] == ",")
            || !all(x[, "timeDummy"] == "," ))
        {
            x <- sienaTimeFix(x)$effects
            x <- fixUpEffectNames(x)
        }
        else
        {
            if (nrow(interactions) > 0)
            {
                x <- fixUpEffectNames(x)
            }
        }
    }
    if (nrow(x) > 0)
    {
        nDependents <- length(unique(x$name))
        userSpecifieds <- x$shortName[x$include] %in%
            c("unspInt", "behUnspInt")
        endowments <- !x$type[x$include] %in% c("rate", "eval", "gmm")
        # includes creations and gmm
        gmm <- any(x$type[x$include] %in% "gmm")
        timeDummies <- !x$timeDummy[x$include] == ","
        specs <- as.data.frame(x[, c("name", "effectName", "include", "fix",
                "test", "initialValue", "parm")])
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
    tmp
}

##@updateSpecification Methods add specified effects from other effects object
updateSpecification <- function(effects.to, effects.from, name.to=NULL, name.from=NULL)
{
    if (!inherits(effects.to, "data.frame"))
    {
        stop("effects.to is not a data.frame")
    }
    if (!inherits(effects.from, "data.frame"))
    {
        stop("effects.from is not a data.frame")
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
    effects.to$parameter[use] <- prevEffects$parameter[correspondence][use]
    effects.to$randomEffects[use] <- prevEffects$randomEffects[correspondence][use]
    effects.to
}
