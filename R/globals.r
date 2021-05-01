##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: globals.r
## *
## * Description: This file contains the code to create and use global objects
## *
## ****************************************************************************/
##@outf Objects/File project .txt file
outf <- NULL
##@lf Objects/File suppressed or to console
lf <- NULL
##@bof Objects/File suppressed or to console
bof <- NULL
##@cf Objects/File suppressed or to console
cf <- NULL

##@Reportfun Reporting Part of global mechanism
Reportfun<- function(x, verbose = FALSE, silent=FALSE)
{
    x <- x
    beverbose <- verbose
    besilent <- silent
    noReportFile <- FALSE
    function(txt, dest, fill=FALSE, sep=" ", hdest, openfiles=FALSE,
             closefiles=FALSE, type=c("a", "w", "n"),  projname="Siena" ,
             verbose=FALSE, silent=FALSE)
    {
        if (openfiles)
        {
            type <- match.arg(type)
            beverbose <<- verbose
            besilent <<- silent
			noReportFile <<- FALSE
            if (type =='w')
            {
                x$outf <<- file(paste(projname, ".txt", sep=""), open="w")
            }
            else if (type =="a")
            {
                x$outf <<- file(paste(projname, ".txt", sep=""), open="a")
            }
            else if (type == "n")
            {
                noReportFile <<- TRUE
            }

        }
        else if (closefiles)
        {
			close(x[["outf"]])
			x$outf <<- NULL
		}
        else
        {
            if (missing(dest) && missing(hdest))
            {
                if (!besilent)
                {
                    cat(txt, fill = fill, sep = sep)
                }
            }
            else
            {
                if (missing(dest))
                {
                    if (hdest  %in% c("cf", "lf", "bof"))
                    {
                        if (beverbose)
                        {
                            cat(txt, fill=fill, sep=sep)
                        }
                    }
                    else
                    {
                        if (!noReportFile)
                        {
                            cat(txt, file = x[[hdest]], fill = fill, sep = sep)
                        }
                    }
                }
                else
                {
                    if (deparse(substitute(dest)) %in% c("cf", "lf", "bof"))
                    {
                        if (beverbose)
                        {
                            cat(txt, fill=fill, sep=sep)
                        }
                    }
                    else
                    {
                        if (is.null(x[[deparse(substitute(dest))]]))
                        {
                            if (!besilent)
                            {
                                cat(txt, fill=fill, sep=sep)
                            }
                        }
                        else
                        {
                            if (!noReportFile)
                            {
                                cat(txt, file=x[[deparse(substitute(dest))]],
                                    fill=fill, sep=sep)
                            }
                        }
                    }
                }
            }
        }
	  }
}


##@Report Globals
Report <- local({verbose <-  NULL; silent <- NULL;
                 Reportfun(list(outf=outf, lf=lf, cf=cf, bof=bof), verbose,
                           silent)})
##@UserInterrupt Siena07/GlobalFunctions Global (within siena07)
UserInterrupt <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@EarlyEndPhase2 siena07/GlobalFunctions
EarlyEndPhase2 <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@UserRestart siena07/GlobalFunctions Global (within siena07)
UserRestart <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@UserInterruptFlag siena07/GlobalFunctions Global (within siena07)
UserInterruptFlag <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@EarlyEndPhase2Flag siena07/GlobalFunctions Global (within siena07)
EarlyEndPhase2Flag <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@UserRestartFlag siena07/GlobalFunctions Global (within siena07)
UserRestartFlag <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@is.batch siena07/GlobalFunctions Global (within siena07)
is.batch <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@DONE siena01/GlobalFunctions Used to communicate with siena.exe and sienaScript
DONE <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;invisible(A)}})
##@FRANstore siena07/GlobalFunctions Used to pass data to other processes
FRANstore <- local({A <-  NULL;function(x){if (!missing(x)) A<<-x;A}})

##@Heading Reporting Global function
Heading<- function(level=1, dest, text, fill=FALSE)
{
	ch <- c("=", "-", " ")[level]
	if (missing(dest))
	{
        Report(c("\n", "@", level, "\n", text, "\n"), sep="", fill=fill)
        Report(rep(ch, sum(nchar(text)) + 3), sep="", fill=fill)
		Report("\n\n")
	}
	else
	{
		dest <- deparse(substitute(dest))
		Report(c("\n", "@", level, "\n", text, "\n"), hdest=dest, sep="", fill=fill)
        Report(rep(ch, sum(nchar(text))), hdest=dest, sep="", fill=fill)
		if (level < 3)
			{
			Report("\n\n", hdest = dest)
		}
		else
		{
			Report("\n", hdest = dest)
		}
	}
}

##@PrtOutMat Reporting
PrtOutMat<- function(mat, dest)
	{
	if (is.null(mat))
		{
			return()
		}
		testing <- Sys.getenv("RSIENATESTING")
		testing <- testing != ""
		if (missing(dest))
		{
			Report(format(t(mat), scientific=testing),
				sep=c(rep.int(" ", ncol(mat) - 1), "\n"))
		}
		else
		{
			Report(format(t(mat), scientific=testing),
				sep=c(rep.int(" ", ncol(mat) - 1), "\n"),
				hdest=deparse(substitute(dest)))
			Report("\n", hdest=deparse(substitute(dest)))
		}
}
##@NullChecks siena07/GlobalFunctions Resets global flags
NullChecks <- function()
{
    UserInterrupt(FALSE)
    EarlyEndPhase2(FALSE)
    UserRestart(FALSE)
    UserInterruptFlag(FALSE)
    EarlyEndPhase2Flag(FALSE)
    UserRestartFlag(FALSE)
}

##@CheckBreaks siena07/GlobalFunctions Reads global flags
CheckBreaks <- function()
{
    UserInterruptFlag(UserInterrupt())
    EarlyEndPhase2Flag(EarlyEndPhase2())
    UserRestartFlag(UserRestart())
}


##@gmm siena07/GlobalFunctions checks gmm in a list, if any
# purpose: backward compatibility
gmm <- function(x)
{
    ifelse(is.null(x$gmm), FALSE, x$gmm)
}
