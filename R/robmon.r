#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: robmon.r
# *
# * Description: This module contains the function robmon which controls
# * the phases of the robbins-monro stochastic approximation algorithm.
# *****************************************************************************/
##args:x: model/algorithm object - intended to be read only
##     z: model fitting object
## returns updated z
##@robmon siena07 Controls MOM process
robmon <- function(z, x, useCluster, nbrNodes, initC, clusterString,
                   clusterIter, clusterType, cl, ...)
{
    z$FinDiff.method<- x$FinDiff.method
    z$n <- 0
    z$OK <-  TRUE
    z$error <- FALSE
    z$restarted <- FALSE
    z$DerivativeProblem <- FALSE
    z$ForceFinDifPhase1 <- FALSE
    z$Phase3Interrupt <- FALSE
    z$repeatsForEpsilon <- 0
    z$maxrepeatsubphase <- 4
    z$gain <- x$firstg
    z$haveDfra <- FALSE
    z$maxlike <- x$maxlike
    z$gmm <- gmm(x)	
	if (is.null(x$sf2.byIteration)) # keep compatible
	{
		z$sf2.byIteration <- TRUE
	}
	else
	{
		z$sf2.byIteration <- x$sf2.byIteration
	}
	if (z$maxlike && !is.batch())
	{
		tcltk::tkconfigure(z$tkvars$phaselabel, text="MCMC Burnin")
	}
    #######################################################
    ##do initial setup call of FRAN
    #######################################################
    if (!is.function(x$FRAN))
    {
        x$FRAN <- getFromNamespace(x$FRANname, pkgname)
    }
    z <- initializeFRAN(z, x, initC=FALSE, ...)
   ## If gmm=TRUE and no gmm effects are included, the algorithm stops and 
    ## asks the user to select the regular MoM estimation
    if (gmm(x) & sum(z$gmmEffects)==0)
    {
        stop ("\n No gmm effects are selected. Use the regular Method of Moments estimation. \n")
    }
    ## If gmm=FALSE and gmm effects are included, the algorithm stops and 
    ## asks the user to select use the GMoM estimation
    if (!gmm(x) & sum(z$gmmEffects)>0)
    {
      stop ("\n gmm effects are selected. Use the Generalized Methods of Moments estimation. \n")
    }

	if (z$maxlike && !is.batch())
	{
		tcltk::tkconfigure(z$tkvars$phaselabel, text="Phase")
	}
    ##
    ##if conditional, FRAN changes z$theta etc
    #######################################################
    if (useCluster)
    {
		## check the version has not been updated since we started R: possible
		## on linux and Mac. Not very stringent check, as environment
		## variables or path could lead to different R or different RSiena.
		packageVersion <- packageDescription(pkgname,
											 fields=c("Version", "Date"))
		if (packageVersion$Version != pkgvers$Version ||
			packageVersion$Date != pkgvers$Date)

		{
			stop(c("Incompatible RSiena versions for sub process: exit R ",
				 "and restart"))
		}
       # if (!is.null(x$FRANname) && x$FRANname != "simstats0c")
       # {
       #     stop("Multiple processors only for simstats0c at present")
       # }
		if (z$returnChains)
		{
			stop("returnChains and useCluster are incompatible")
		}
        if (!clusterIter && nbrNodes >= z$observations)
        {
            stop("Not enough observations to use the nodes")
        }
		if (!length(cl)) {
			unlink("cluster.out")
			if (clusterType == "FORK")
			{
				cl <- makeCluster(nbrNodes, type = clusterType,
					outfile = "cluster.out")
			}
			else
			{
				cl <- makeCluster(clusterString, type = clusterType,
					outfile = "cluster.out")
			}
		}
        clusterCall(cl, library, pkgname, character.only = TRUE)
		##parLapply(cl, c('f1','f2'), sink)
		z$oldRandomNumbers <- .Random.seed
		## The possibility to use snow now has been dropped
		## because RSiena requires R >= 2.15.0
		## and snow is superseded.
		## Therefore the call to clusterSetupRNG was dropped.
		#if (getRversion() < "2.14.0")
			## fake this to recreate old results
			##	if (TRUE)
		#{
		#	clusterSetupRNG(cl, seed = as.integer(runif(6,
		#						max=.Machine$integer.max)))
		#}
		#else
		#{
		clusterSetRNGStream(cl, iseed = as.integer(runif(1,
							max=.Machine$integer.max)))
		#}
        clusterCall(cl, storeinFRANstore,  FRANstore())
        if (initC)
        {
			clusterCall(cl, initializeFRAN, z, x, initC=initC)
        }
        z$cl <- cl
    }
    else
    {
        z$cl  <-  NULL
    }

	# check dimensionality of thetaValues
	if (z$thetaFromFile)
	{
		if (dim(z$thetaValues)[2] != z$pp)
		{
			stop(paste('Matrix thetaValues, if given, should have', z$pp, 'columns.'))
		}
	}

    z$newFixed <- rep(FALSE, z$pp)
    z$AllNowFixed <- FALSE
    if (!z$haveDfra)
    {
        z$dinv <- matrix(NA, nrow = z$pp, ncol = z$pp)
    }
    z$scale <- rep(0.1, z$pp)
    Report('\n', outf)
    Report('\nStochastic approximation algorithm.\n', cf)
    if (x$firstg <= 0)
    {
        Report(c('Initial value of the gain parameter is ', x$firstg,
                 '.\n'), outf)
        Report('This is not allowed; changed to 0.0001.\n', outf)
        z$gain <- 0.0001
    }
    if (x$reduceg <= 0)
    {
        Report(c('Reduction factor for the gain parameter is ', x$reduceg,
                 '.\n'), outf)
        Report('This is not allowed; changed to 0.2.\n', outf)
        x$reduceg <- 0.2
    }
    Report(c('Initial value for gain parameter = ', format(z$gain),
			 '.\nReduction factor for gain parameter = ', format(x$reduceg),
             '.\nStart of the algorithm.\n'), cf, sep='')
    Report('Observed function values are \n', cf)
	targets <- if (!z$maxlike) z$targets else z$maxlikeTargets
    ftargets <- format(targets, width = 10, nsmall = 4)
    fnum <- format(1 : z$pp, width = 3)
    Report(c(paste(fnum, '. ', ftargets, sep = '')), cf, fill=80)
    z$epsilon<- pmin(0.1,z$scale)
    z$epsilon[z$posj] <- 0.1 * z$theta[z$posj]
    z$theta0 <- z$theta
	## store starting value without any conditioning variables
    z$anyposj <- any(z$posj)
    if (!gmm(x))
    {
      z$n1 <- 7 + 3 * z$pp
    }
    else
    {
      z$n1 <- 100 + 7 * z$pp
    }
	if (x$dolby){z$n1 <- max(z$n1, 50)}
    if (any(!z$fixed))
    {
        z$AllUserFixed<- FALSE
    }
    else
    {
        z$AllUserFixed <- TRUE
    }

    repeat  ##this is startagain:
    {
      z$epsilonProblem <- FALSE
      repeat ## this loop is simply to break out of! only intend to do it once
        {
            if (any(!z$fixed))
                z$AllNowFixed<- FALSE
            else
                z$AllNowFixed <- TRUE
            phase3Only <- FALSE
            if (!is.batch())
            {
                tcltk::tkdelete(z$tkvars$subphase, 0, "end")
                tcltk::tkconfigure(z$tkvars$earlyEndPhase2, state="disabled")
                tcltk::tkconfigure(z$tkvars$subphase, state="disabled")
                tcltk::tkconfigure(z$tkvars$subphaselabel, state="disabled")
            }
            if (z$AllNowFixed || x$nsub == 0)
            {
                if (z$AllNowFixed)
                {
                    Report('All parameters are fixed.\n', outf)
                    Report(c('All parameters fixed; no estimation;',
                             'only Phase 3.\n'), lf)
                }
                else
                {
                    Report('Number of subphases is specified as 0.\n', outf)
                    Report('0 subphases; no estimation: only phase 3.\n', lf)
                }
                Report(c('Therefore the estimation phase is skipped\n',
                         'and the program passes on immediately to phase 3\n',
                         'for checking the current parameter values and',
                         'calculating standard errors.\n'),
                       outf)
                phase3Only <- TRUE
            }
            NullChecks() ## reset user interrupt variables and flags
            z$Phase <- 0
            z <- AnnouncePhase(z, x)
            WriteOutTheta(z)
            z$fixed[z$newFixed] <- FALSE
            z$newFixed <- rep(FALSE,z$pp) ##not clear what we should do here
            z$restart <- FALSE
            z$OK <- TRUE
            if (!phase3Only)
            {
                if (!z$haveDfra)
                {
                    ##start phase1 and do 10 (approx) iterations,
                    z <- phase1.1(z, x, ...)
                    ##check epsilon
                    if (!z$OK || UserInterruptFlag() || UserRestartFlag())
                    {
                        break
                    }
                    if (z$epsilonProblem && z$repeatsForEpsilon<4)
                    {
                        z$repeatsForEpsilon<- z$repeatsForEpsilon+1
                        z$restart<- !z$restarted
                        break
                    }
                    z<- phase1.2(z, x, ...)
                    if (!z$OK || z$DerivativeProblem ||
                        UserInterruptFlag() || UserRestartFlag())
                    {
                        break
                    }
                }
                if (x$nsub > 0)
                {
                    z <- phase2.1(z, x, ...)
                }
                if (!z$OK || UserInterruptFlag() || UserRestartFlag())
                {
                    if (!is.batch())
                    {
                        tcltk::tkdelete(z$tkvars$subphase,0,'end')
                        tcltk::tkconfigure(z$tkvars$earlyEndPhase2,state='disabled')
                        tcltk::tkconfigure(z$tkvars$subphase,state='disabled')
                        tcltk::tkconfigure(z$tkvars$subphaselabel,state='disabled')
                    }
                    break
                }
                if (x$nsub>1 && !EarlyEndPhase2Flag())
                {
                    z <- proc2subphase(z, x, 2, ...)
                }
                if (!z$OK || z$restart||
                    UserInterruptFlag() || UserRestartFlag() )
                {
                    if (!is.batch())
                    {
                        tcltk::tkdelete(z$tkvars$subphase,0,'end')
                        tcltk::tkconfigure(z$tkvars$earlyEndPhase2, state='disabled')
                        tcltk::tkconfigure(z$tkvars$subphase, state='disabled')
                        tcltk::tkconfigure(z$tkvars$subphaselabel, state='disabled')
                    }
                    break
                }
                if (x$nsub > 2 && !EarlyEndPhase2Flag())
                {
                    for (i in 3 : x$nsub)
                    {
                        z <- proc2subphase(z, x, i, ...)
                        if (!z$OK || UserInterruptFlag() ||
                            UserRestartFlag() || EarlyEndPhase2Flag())
                        {
                            if (!is.batch())
                            {
                                tcltk::tkdelete(z$tkvars$subphase, 0, 'end')
                                tcltk::tkconfigure(z$tkvars$earlyEndPhase2,
                                            state='disabled')
                                tcltk::tkconfigure(z$tkvars$subphase, state='disabled')
                                tcltk::tkconfigure(z$tkvars$subphaselabel,
                                            state='disabled')
                            }
                            break
                        }
                    }
                    if (!z$OK || UserInterruptFlag() ||
                        UserRestartFlag() )
                    {
                        if (!is.batch())
                        {
                            tcltk::tkdelete(z$tkvars$subphase, 0, 'end')
                            tcltk::tkconfigure(z$tkvars$subphase, state='disabled')
                            tcltk::tkconfigure(z$tkvars$subphaselabel,
                                        state='disabled')
                            tcltk::tkconfigure(z$tkvars$earlyEndPhase2,
                                        state='disabled')
                        }
                        break
                    }
                }
            }
            if (!is.batch())
            {
                tcltk::tkdelete(z$tkvars$subphase, 0, 'end')
                tcltk::tkconfigure(z$tkvars$subphase, state='disabled')
                tcltk::tkconfigure(z$tkvars$subphaselabel, state='disabled')
                tcltk::tkconfigure(z$tkvars$earlyEndPhase2, state='disabled')
            }
            z<- phase3(z, x, ...)
            break
        }
        ##stop if not OK or user has asked to
        if (!z$OK || UserInterruptFlag())
        {
            break
        }
        ##break unless we want to start at phase 1 again.
        if (!z$DerivativeProblem && !z$restart && !UserRestartFlag())
        {
            break
        }
    }
    if (!z$OK || UserInterruptFlag())
    {
        if (z$OK)
        {
            if (!z$Phase3Interrupt)
                z$termination <- 'UserInterrupt'
            else
                z$termination <- 'OK'
        }
        else
        {
            z$termination <- 'Error'
        }
        if (!z$OK || !z$Phase3Interrupt)
		{
			if (useCluster)
			{
				stopCluster(cl)
			}
			useCluster <- FALSE
            return(z)
		}
    }
    ## #####################################################
    ## do final call of FRAN
    ## #####################################################
    z <- terminateFRAN(z, x)
    ## #####################################################
    ## call to FRAN changes covariance matrix for conditional estimation
	if (x$simOnly)
	{
		if (z$thetaFromFile)
		{
			z$theta <- colMeans(z$thetaUsed)
			z$covtheta <- cov(z$thetaUsed)
			z$se <- sqrt(diag(z$covtheta))
			z$thetaValues <- NULL # not needed any longer, superseded by z$thetaUsed
		}
		else if (!gmm(x))
		{
			z$covtheta <- matrix(0, z$pp, z$pp )
			z$se <- rep(0, z$pp)
		}
		else if (gmm(x))
	  	{
	    	z$covtheta <- matrix(0, z$pp - sum(z$gmmEffects), z$pp- sum(z$gmmEffects))
	    	z$se <- rep(0, z$pp- sum(z$gmmEffects))
	  	}
	}
	else if (!gmm(x))
	{
		z$diver<- (z$fixed | z$diver | diag(z$covtheta) < 1e-9) & (!z$AllUserFixed)
		z$covtheta[z$diver, ] <- NA # was Root(diag(z$covtheta)) * 33
		##not sure this does not use very small vals
		z$covtheta[, z$diver] <- NA # was Root(diag(z$covtheta)) * 33
		diag(z$covtheta)[z$diver] <- NA # was 999
		z$se <- sqrt(diag(z$covtheta))
	}
#	z$gmm <- FALSE
    z$termination <- 'OK'
    if (useCluster)
    {
		stopCluster(cl)
	}
    z
}

