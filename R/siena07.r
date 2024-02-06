##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: https://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: siena07.r
## *
## * Description: This file contains the main controlling module for the
## * estimation.
## * Also contains utility functions used within siena07
## ****************************************************************************/

##@siena07 siena07
siena07 <- function(x, batch = FALSE, verbose = FALSE, silent=FALSE,
	useCluster = FALSE, nbrNodes = 2,
	thetaValues = NULL,
    returnThetas = FALSE,
	targets = NULL,
	initC=TRUE,
	clusterString=rep("localhost", nbrNodes), tt=NULL,
	parallelTesting=FALSE, clusterIter=!x$maxlike,
	clusterType=c("PSOCK", "FORK"), cl=NULL, ...)
{
	exitfn <- function()
	{
		if (!is.batch())
		{
		  tcltk::tkdestroy(tkvars$tt)
		}
		## close the report file
		Report(closefiles=TRUE)
		RNGkind("default")
	}
	on.exit(exitfn())

	# If the user is passing clusters through -cl- then change the
	# useCluster to TRUE, and assign the -nbrNodes- to number of nodes
	if (!useCluster & length(cl))
	{
		useCluster <- TRUE
		nbrNodes   <- length(cl)
	}

	time0 <-  proc.time()['elapsed']
	z <- NULL ## z is the object for all control information which may change.
	## x is designed to be readonly. Only z is returned.
	z$x <- x
	z$returnThetas <- returnThetas
	if (useCluster)
	{
		if (parallelTesting)
		{
			stop("cannot parallel test with multiple processes")
		}
		if (!length(cl))
		{
			clusterType <- match.arg(clusterType)
			if (.Platform$OS.type == "windows" && clusterType != "PSOCK")
			{
				stop("cannot use forking processes on Windows")
			}
		}
		# The possibility to use snow now has been dropped
		# because RSiena requires R >= 2.15.0
		# and snow is superseded.
		# if (getRversion() < "2.14.0")
		## fake this to recreate old results
		## if (TRUE)
		#{
		#	require(snow, warn.conflicts=FALSE)
		#	require(rlecuyer)
		#	clusterType <- "SOCK"
		#}
		#else
		#{
		## require(parallel)
		#}
		if (clusterIter)
		{
			x$firstg <- x$firstg * sqrt(nbrNodes)
			z$int <- nbrNodes
			z$int2 <- 1
		}
		else
		{
			z$int <- 1
			z$int2 <- nbrNodes
		}
	}
	else
	{
		z$int <- 1
		z$int2 <- 1
	}
	if (parallelTesting)
	{
		set.seed(1, kind='Wich')
		## randomseed2 is for second generator needed only for parallel testing
		randomseed2 <- .Random.seed
		## .Random.seed[2:4] <- as.integer(c(1,2,3))
		## randomseed2[2:4] <- as.integer(c(3,2,1))
		randomseed2[2:4] <- as.integer(c(1, 2, 3))
		seed <- 1
		newseed <- 1
		z$parallelTesting <- TRUE
	}
	else
	{
		randomseed2 <- NULL
		## x$randomSeed is the user seed, if any
		if (!is.null(x$randomSeed))
		{
			set.seed(x$randomSeed, kind="default")
			seed <- x$randomSeed
		}
		else
		{
			if (!nzchar(Sys.getenv("RSIENA_TESTING")))
			{
				if (exists(".Random.seed"))
				{
					rm(.Random.seed, pos=1)
				}
				newseed <- trunc(runif(1) * 1000000)
				set.seed(newseed)  ## get R to create a random number seed for me.
				seed <- NULL
			}
			else
			{
				newseed <- NULL
				seed <- NULL
			}
		}
	}
	z$randomseed2 <- randomseed2

	x$targets <- targets

	## set the global is.batch
	batchUse <- batch
	if (!batch)
	{
		if (!requireNamespace("tcltk", quietly = TRUE))
		{
				batchUse <- TRUE
				message("Package tcltk not available, forcing use of batch mode")
		}
		else
		{
			if (.Platform$OS.type != "windows")
			{
				if (!capabilities("X11"))
				{
					batchUse <- TRUE
					message("No X11 device available, forcing use of batch mode")
				}
			}
		}
		if(nzchar(Sys.getenv("RSIENA_TESTING")))
		{
			silent <- TRUE
		}
	}
	is.batch(batchUse)

	## open the output file
	Report(openfiles=TRUE, projname=x$projname, verbose=verbose, silent=silent)
	z <- InitReports(z, seed, newseed)

	## reset the globals for interrupts
	NullChecks()

	## create the screen
	if (!is.batch())
	{
		tkvars <- siena07Gui(tt=tt)
		z$tkvars <- tkvars
		z$pb <- list(pb=tkvars$pb, pbval=0, pbmax=1)
	}
	else
	{
		z$pb <- list(pb=NULL, pbval=0, pbmax=1)
	}

	## create theta values for phase 3, if necessary
	if (is.null(thetaValues))
	{
		z$thetaValues <- NULL
		z$thetaFromFile <- FALSE
	}
	else
	{
		if (!x$simOnly)
		{
			cat('The thetaValues parameter was given\n')
			cat('but simOnly was not specified in the algorithm object.\n')
			cat('This is inconsistent.\n')
			stop('Inconsistent combination simOnly - thetaValues')
		}
		if (!is.matrix(thetaValues))
		{
			stop('thetaValues should be a matrix.')
		}
		z$thetaValues <- thetaValues
		z$thetaFromFile <- TRUE
	}

	z <- robmon(z, x, useCluster, nbrNodes, initC, clusterString,
		clusterIter, clusterType, cl, ...)

	time1 <-  proc.time()['elapsed']
	Report(c("Total computation time", round(time1 - time0, digits=2),
			"seconds.\n"), outf)
z$compTime <- round(time1 - time0, digits=2)

	if (useCluster)
	{
		# Only stop cluster if it wasn't provided by the user
		#	  if (!length(cl))
		#		  stopCluster(z$cl)

		## need to reset the random number type to the normal one
		assign(".Random.seed", z$oldRandomNumbers, pos=1)
	}

	class(z) <- "sienaFit"
	z$tkvars <- NULL
	z$pb <- NULL
	z
}

##@InitReports siena07 Print report
InitReports <- function(z, seed, newseed)
{
	Report("\n\n-----------------------------------\n", outf)
	Report("New Analysis started.\n", outf)
	Report(c("Date and time:", format(Sys.time(),"%d/%m/%Y %H:%M:%S")), outf)
	Report("\nNew results follow.\n", outf)
	Report("-----------------------------------\n", outf)
	version <- packageDescription(pkgname, fields = "Version")
	Report(c(paste("\n", pkgname, " version ", sep = ""), version, " (",
		format(as.Date(packageDescription(pkgname, fields = "Date")), "%d %b %y"),
		")", "\n\n"), sep = "", outf)
	if (z$x$simOnly)
	{
		Heading(1, outf, "Simulations.")
	}
	else
	{
		Heading(1, outf, "Estimation by stochastic approximation algorithm.")
	}
	if (is.null(seed))
	{
		Report("Random initialization of random number stream.\n", outf)
		if (!is.null(newseed))
		{
			Report(sprintf("Current random number seed is %d.\n", newseed), outf)
		}
	}
	else
	{
		Report(sprintf("Current random number seed is %d.\n", seed), outf)
	}
	z$version <- version
	z$startingDate <- date()
	z
}

##@AnnouncePhase siena07 Progress reporting
AnnouncePhase <- function(z, x, subphase=NULL)
{
	## require(tcltk)
	if (!is.batch())
	{
		tcltk::tkdelete(z$tkvars$phase, 0, "end")
		tcltk::tkinsert(z$tkvars$phase, 0, paste(" ", z$Phase))
		tcltk::tkdelete(z$tkvars$subphase, 0, "end")
		tcltk::tkdelete(z$tkvars$iteration, 0, "end")
		tcltk::tkinsert(z$tkvars$iteration, 0, format(0, width=6))
	}
	if (missing(subphase))
	{
		if (is.batch())
		{
			##Report(c("\nStart phase", z$Phase, "\n"), cf)
			Report(c("\nStart phase", z$Phase, "\n"))
		}
	}
	else
	{
		if (!is.batch())
		{
			tcltk::tkinsert(z$tkvars$subphase, 0, paste(" ", subphase))
		}
		else
		{
			##Report(c("\nStart phase ", z$Phase, ".", subphase, "\n"),
			##sep="", cf)
			Report(c("\nStart phase ", z$Phase, ".", subphase, "\n"), sep="")
		}
	}
	if (z$Phase == 0)
	{
		if (!is.batch())
		{
			tcltk::tkconfigure(z$tkvars$current, height=min(z$pp, 30))
			tcltk::tkconfigure(z$tkvars$deviation, height=min(z$pp, 30))
			tcltk::tkconfigure(z$tkvars$quasi, height=min(z$pp, 30))
		}
		n1pos <- z$n1 * (z$pp + 1)
		z$n2min0 <- 7 + z$pp
		z$n2min0 <- max(5, z$n2min0 / z$int)
		z$n2minimum<- rep(0, x$nsub)
		z$n2maximum<- rep(0, x$nsub)
		## 2.5198421 = 2^(4/3); this gives a gain parameter of order n^(-3/4) ##
		if (x$nsub > 0)
		{
			z$n2minimum[1] <-
				ifelse(is.null(x$n2start), trunc(z$n2min0 * 2.52), x$n2start)
			z$n2maximum[1] <- z$n2minimum[1] + 200
			if (x$nsub > 1)
			{
				for (i in 2:x$nsub)
				{
					z$n2minimum[i] <- trunc(z$n2minimum[i-1] * 2.52)
					z$n2maximum[i] <- z$n2minimum[i] + 200
				}
			}
		}
		z$n2partsum <- c(0, cumsum(z$n2maximum))
		n2sum <- sum(z$n2maximum)
		##Progress bar
		pbmax <- n1pos + n2sum + x$n3
		z$n1pos<- n1pos
		if (!x$maxlike && z$FinDiff.method)
			pbmax <- pbmax + x$n3 * z$pp
		z$pb$pbval <- 0
		z$pb <- createProgressBar(z$pb, maxvalue=pbmax)
		z$pb$pbmax <- pbmax
	}
	if (z$Phase==2)
	{
		propo <- z$n1pos + z$n2partsum[max(subphase,1)]
		if (propo> getProgressBar(z$pb))
			z$pb <-setProgressBar(z$pb,propo)
	}
	if (z$Phase ==3)
	{
		propo <- z$n1pos + z$n2partsum[x$nsub + 1]
		if (!z$AllUserFixed)
			z$pb <- setProgressBar(z$pb,propo)
		else
		{
			max <- x$n3
			z$pb <-createProgressBar(z$pb,max)
		}
	}
	z
}
##@roundfreq siena07 Prettify interval between progress reports
roundfreq <- function(w)
{
	vec1 <- c(1, 2, 3, 4, 31, 66, 101, 300, 500)
	vec2 <- c(1, 2, 3, 20, 50, 100, 200, 500)
	if (is.batch())
		w <- max(10,vec2[findInterval(w, vec1, all.inside=TRUE)])
	else
		w <- vec2[findInterval(w, vec1[1:7], all.inside=TRUE)]
	w
}

##@WriteOutTheta siena07 Progress reporting
WriteOutTheta <- function(z)
{
	if (!is.batch())
	{
		DisplayTheta(z)
	}
	else
	{
		Report(c("theta:", format(z$theta, digits=3), "\n"))
	}
	Report("Current parameter values:\n", cf)
	Report(format(z$theta), cf, fill=80)
}

##@DisplayThetaAutocor siena07 Progress reporting
DisplayThetaAutocor <- function(z)
{
	if (!is.batch())
	{
		DisplayTheta(z)
		tcltk::tkdelete(z$tkvars$quasi, "1.0", "end")
		tcltk::tkinsert(z$tkvars$quasi, "1.0", FormatString(z$pp, z$ac))
	}
	else
	{
		Report(c("theta", format(z$theta, digits=3),"\n"))
		Report(c("ac", format(z$ac, digits=3), "\n"))
	}

}
##@DisplayandWriteTheta siena07 Progress reporting
DisplayandWritetheta <- function(z)
{
	if (!is.batch())
	{
		DisplayTheta(z)
	}
	else
	{
		Report(c("theta", format(z$theta, digits=3), "\n"))
	}
}
##@DisplayTheta siena07 Progress reporting
DisplayTheta <- function(z)
{
	if (!is.batch())
	{
		tcltk::tkdelete(z$tkvars$current, "1.0", "end")
		tcltk::tkinsert(z$tkvars$current, "1.0", FormatString(z$pp, z$theta))
	}

}
##@FormatString siena07 Progress Reporting
FormatString <- function(pp, value)
{
	ppuse <- min(30, pp)
	nbrs <- format(1:ppuse)
	nch <- nchar(nbrs[1])
	formatstr <- paste("%", nch, "d.%", (13 - nch), ".4f\n", sep="",
		collapse="")
	paste(sprintf(formatstr, 1:ppuse, value[1:ppuse]), collapse="")
}
##@DisplayDeviations siena07 Progress reporting
DisplayDeviations <- function(z, fra)
{
	if (!is.batch())
	{
		tcltk::tkdelete(z$tkvars$deviations, "1.0", "end")
		tcltk::tkinsert(z$tkvars$deviations, "1.0", FormatString(z$pp, fra))
	}
}
##@DisplayIteration siena07 Progress reporting
DisplayIteration <- function(z)
{
	if (!is.batch())
	{
		tcltk::tkdelete(z$tkvars$iteration, 0, "end")
		tcltk::tkinsert(z$tkvars$iteration, 0, format(z$nit, width=6))
		tcltk::tcl("update")
	}
}
##@Root siena07 Safe square root for compatibility with siena3. Probably not necessary in R.
Root<- function(x)
{
	ifelse(abs(x) > 1e-36, sqrt(abs(x)), 1e-18)
}

##@getProgressBar siena07 Progress reporting
getProgressBar <- function(pb)
{
	if (is.batch())
		val <- pb$pbval
	else
		val <- as.numeric(tcltk::tclvalue(tcltk::tkcget(pb$pb, "-value")))
	val
}

##@setProgressBarProgress siena07 reporting
setProgressBar <- function(pb, val)
{
	if (is.batch())
	{
		pb$pbval <- val
	}
	else
	{
		tcltk::tkconfigure(pb$pb, value=val)
		tcltk::tcl("update")
	}
	pb
}
##@createProgressBar siena07 Progress reporting
createProgressBar <- function(pb, maxvalue)
{
	if (is.batch())
		pb$pbmax <- maxvalue
	else
		tcltk::tkconfigure(pb$pb, maximum=maxvalue)
	pb
}
##@tkErrorMessage Miscellaneous Not used
tkErrorMessage <- function()
{
	tcltk::tkmessageBox(geterrmessage(), icon="error")
}

##@errorHandler Miscellaneous Not used
errorHandler <- function()
{
	## opts <- options()
	if (!is.batch())
	{
		options(show.error.messages=FALSE)
		options(error=tkErrorMessage)
	}
}
