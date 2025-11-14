##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: https://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: phase2.r
## *
## * Description: This module contains the functions phase2.1, proc2subphase
## * and doIterations which together perform a robbins-monro stochastic
## * approximation algorithm.
## * phase2.1 and proc2subphase are called from robmon in robmon.r.
## ****************************************************************************/
## args: z: internal control object
##       x: model object (readonly as not returned)
## Hidden options exist that are activated 
## if the algorithm object x has a component x$moreUpdates or x$phase2imp

##@storeinFRANstore siena07 Used to avoid Namespace problems with multiple processes
storeinFRANstore <- function(...)
{
    FRANstore(...)
}
##@phase2.1 siena07 Start phase 2
phase2.1<- function(z, x, ...)
{
		if (ifelse(is.null(x$phase2imp),FALSE,x$phase2imp))
		{
			z$addChainToStore <- TRUE
			z$nbrNodes <- 1
			z$storedChains <- TRUE
			z$nGroup <- 1
			z$thetaList <- list()
			z$nCount <- 1
		}
	
    #initialise phase2
    if (x$maxlike)
    {
        z$phase2fras <- array(0, dim=c(4, z$pp, 1000))
        z$rejectprops <- array(0, dim=c(4, 4, 1000))
    }
    z$Phase <- 2
    z$writefreq <- 1
    if (!is.batch())
    {
        tcltk::tkconfigure(z$tkvars$earlyEndPhase2,state='normal')
        tcltk::tkconfigure(z$tkvars$subphase,state='normal')
        tcltk::tkconfigure(z$tkvars$subphaselabel,state='normal')
    }
    z$Deriv <- FALSE
	msf <- as.matrix(cov(z$sf)) # as.matrix() just in case z$pp = 1
#	z$sd <- sqrt(diag(msf))
# Instead of the preceding line,
# the following is used for equality with earlier versions.
	z$sd <- sqrt(pmax(apply(z$sf, 2, function(x) sum(x^2) / nrow(z$sf) - mean(x)^2),0))
    Report("standardization : \n", cf)
	if (!z$gmm) 
	{
   	z$sd[z$fixed] <- 0
	  z$standardization <-
						1/sqrt(pmax(diag(as.matrix(z$dinvv %*% msf %*% t(z$dinvv))),0))
	PrtOutMat( cbind(z$requestedEffects[,"effectName"], round(z$standardization,4)), cf)
	}
	else
	{
	  z$sd[which(z$fixed & !z$gmmEffects)] <- 0
	  z$standardization <- 1/sqrt(pmax(diag(as.matrix(z$dinvv %*% msf %*% t(z$dinvv))),0))
	  Report(round(z$standardization,4), cf)
	}
	if (sum(z$fixed) < z$pp)
	{
		z$sf.invcov <-
		 solve(msf[!z$fixed, !z$fixed] + 0.0001 * diag(z$pp - sum(z$fixed)))
	}
    Report(paste('\nPhase 2 has', x$nsub, 'subphases.\n'), cf)
    z$gain <- x$firstg
    z$reduceg <- x$reduceg
    if (x$nsub <= 0)
    {
        Report('With 0 subphases, there is no phase 2.\n', cf)
    }
    else
    {
        if (z$maxrepeatsubphase >= 2)
        {
            Report(c('Each subphase can be repeated up to', z$maxrepeatsubphase,
                   'times\n'), cf)
        }
        z <- proc2subphase(z, x, 1, ...)
    }
    z
}
##@proc2subphase siena07 Do one subphase of phase 2
proc2subphase <- function(z, x, subphase, ...)
{
	## init subphase of phase 2
	z <- AnnouncePhase(z, x, subphase)
	Report(paste("\nStart phase ", z$Phase, ".", subphase, "\n", sep=""), cf)
	if (subphase <= 0)
	{
		z$n2min <- 5
		z$n2max <- 100
	}
	else
	{
		z$n2min <- z$n2minimum[subphase]
		z$n2max <- z$n2maximum[subphase]
	}
	Report(paste("Number of iterations minimum ", z$n2min, ", maximum", z$n2max, ".\n", sep=""), cf)
	z$repeatsubphase <- 0
	repeat
	{
		z$repeatsubphase <- z$repeatsubphase + 1
		z$truncated <- rep(FALSE, z$n2max)
		z$positivized <- rep(0, z$pp)
		z$ctime <- proc.time()[3]
		z$time1 <- proc.time()[3]
		z$thav <- z$theta
		z$thavn <- 1
		## cat(z$thav, z$theta, '\n')
		z$prod0 <- rep(0, z$pp)
		z$prod1 <- rep(0, z$pp)
		functionAv <- 0
		## ###############################################
		## do the iterations for this repeat of this subphase
		## ##############################################
		z <- doIterations(z, x, subphase, ...)
		if (!z$OK || UserInterruptFlag() || UserRestartFlag() ||
			EarlyEndPhase2Flag())
		{
			break
		}
		##
		## end processing for this repeat of this subphase
		##
		## report truncations and positivizations
		if (any(z$truncated))
		{
			msg<- paste('Intervention 2.', subphase,
				'.1: changes truncated, iterations: ',sep='')
			Report(msg, cf)
			Report(which(z$truncated), cf, fill=80)
		}
		if (any(z$positivized >= 1))
		{
			msg <- paste('Intervention 2.',subphase,
				'.2: positivity restriction:\n ',sep='')
			Report(msg,cf)
			subs <- which((z$positivized) > 0)
			msg <- paste('Positivized: ',
				paste((1:z$pp)[subs], ': ', z$positivized[subs], '; ',
					collapse=' '))
			#            msg<- sapply(subs, function(i, y)
			#                         paste('Observation:', i, 'Coordinate(s):',
			#                               paste((1:z$pp)[y[i,]], collapse = ' ')),
			#                         y = z$positivized)
			Report(msg, cf, fill=80)
		}
		if ((subphase >= 1) && (z$maxacor >= sqrt(2.0 / (z$nit + 1))))
		{
			Report('Note: an autocorrelation is positive at the end',cf)
			Report(' of this subphase.\n',cf)
			Report('Autocorrelations:\n',cf)
			prtmat <- z$prod1 / z$prod0
			PrtOutMat(as.matrix(prtmat), cf)
		}
		if ((z$nit >= z$n2max)
			|| (subphase <= 0)
			|| (z$maxacor < 1e-10)
			|| (z$repeatsubphase >= z$maxrepeatsubphase))
		{
			break
		}
	}
	##finalize the subphase
	if (!z$OK || UserInterruptFlag() || UserRestartFlag())
	{
		return(z)
	}
	if (EarlyEndPhase2Flag())
	{
		Report('The user asked for early end of phase 2.\n', outf)
	}
	##    cat('it',z$nit,'\n')
	##recalculate autocor using -1 instead of -2 as error
	ac <- ifelse (z$prod0 > 1e-12, z$prod1 / z$prod0, -1)
	maxacor <- max(-99, ac[!z$fixed]) ##note -1 > -99
	z$theta <- z$thav / z$thavn   # z$thavn = (z$nit + 1)

	DisplayThetaAutocor(z)
	Report(paste('Phase ', z$Phase,'.', subphase, ' ended after ', z$nit,
			' iterations.\n', sep = ''), cf)
	if (x$checktime)
	{
		time1 <- proc.time()[[3]] - z$ctime
		subphaseTime <- time1 / z$nit
		Report(paste('Time per iteration in phase ', z$Phase, '.', subphase,
				' = ', format(subphaseTime, nsmall=4, digits=4),
				'\n', sep=''), lf)
	}
	if ((maxacor >= sqrt(2 / (z$nit + 1))) && (subphase >= 1))
	{
		Report('Note. Autocorrelation criterion not satisfied.\n', cf)
	}
	WriteOutTheta(z)
	if (EarlyEndPhase2Flag())
	{
		return(z)
	}
	if (subphase == 2 && z$restart) ## this means we restarted in phase 1
		##because of epsilon and need to restart the whole thing again now
	{
		Report('Restart after subphase 2.2 from current parameter values\n', cf)
		Report('because the initial values used ', cf)
		Report('led to questionable epsilon values\n', cf)
		z$fixed[z$newfixed] <- FALSE
		z$restarted <- TRUE
	}
	## For the next subphase:
	if (x$maxlike)
	{
		z$gain <- z$gain * 0.25
	}
	else
	{
		z$gain <- z$gain * z$reduceg
	}
	z
} ##end of this subphase

##@doIterations siena07 Do all iterations for 1 repeat of 1 subphase of phase 2
doIterations<- function(z, x, subphase,...)
{
	z$nit <- 0
	ac <- 0
	z$maxacor <- -1
	z$minacor <- 1
	xsmall <- NULL
	zsmall <- makeZsmall(z)
	z$returnDeps <- FALSE
	sumfra <- 0.0
	if (z$returnThetas)
	{
		thetas <- NULL
		sfs <- NULL
	}
	repeat
	{
		z$n <- z$n+1
		z$nit <- z$nit + 1
		if (z$returnThetas)
		{
			thetas <- rbind(thetas, c(subphase, z$theta))
		}
		if (subphase == 1 && z$nit == 2)
			z$time1 <- proc.time()[[3]]
		if (subphase == 1 && z$nit == 11)
		{
			time1 <- proc.time()[[3]] - z$time1
			if (time1 > 1e-5)
			{
				z$writefreq <- max(1, round(20.0 / time1))
			}
			else
			{
				z$writefreq <- 20
			}
			z$writefreq <- roundfreq(z$writefreq)
			if (is.batch())
			{
				z$writefreq <-  z$writefreq * 10 ##compensation for it
				## running faster with no tcl/tk
			}
		}
		zsmall$nit <- z$nit
		if (x$dolby) {zsmall$Deriv <- TRUE} ## include scores in FRAN
		

#Report(paste("z$theta: ", "\n"), cf)
#PrtOutMat(as.matrix(z$targets), cf)
#PrtOutMat(as.matrix(z$theta), cf)
#PrtOutMat(as.matrix(zsmall$theta), cf)


		if (z$int == 1) ## then no parallel runs at this level
		{
			if (ifelse(is.null(x$phase2imp),FALSE,x$phase2imp))
			{
				zz <- x$FRAN(zsmall, xsmall,  returnLoglik=TRUE, returnChains=TRUE)
				z$myloglik2 <- c(z$myloglik2, sum(zz$loglik))
			} else {
				zz <- x$FRAN(zsmall, xsmall)
			}
			fra <- colSums(zz$fra) - z$targets
			
			
#Report(paste("fra (1): ", "\n"), cf)
#PrtOutMat(as.matrix(fra), cf)

			if (!zz$OK)
			{
				z$OK <- zz$OK
				break
			}
			if (z$returnThetas)
			{
				sfs <- rbind(sfs, c(subphase, fra))
			}
			if (x$dolby)
				## subtract regression on scores;
				## permitted because fra is used in phase2 only for updates
			{
				ssc <- zz$sc # scores;  periods by variables
				ssc <- colSums(ssc) # add over periods
				fra <- fra - (z$regrCoef * ssc)
			}
		}
		else
		{
			zz <- clusterCall(z$cl, simstats0c, zsmall, xsmall)
			if (x$dolby)
				## subtract regression on scores;
			{
				fra <- sapply(zz, function(x){
					ssc <- x$sc # scores;  periods by variables
					ssc <- colSums(ssc) # add over periods
					colSums(x$fra)- z$targets - (z$regrCoef * ssc)
			})
			}
			else
			{
				fra <- sapply(zz, function(x) colSums(x$fra)- z$targets)
			}
			dim(fra) <- c(z$pp, z$int)
			fra <- rowMeans(fra)
			zz$OK <- sapply(zz, function(x) x$OK)
			if (!all(zz$OK))
			{
				z$OK <- FALSE
				break
			}

		}
		
		if ((z$nit <= 10) || (z$nit %% z$writefreq ==0))
		{
			DisplayIteration(z)
			if (is.batch())
			{
				val <- getProgressBar(z$pb)
				increment <- ifelse(z$nit <= 10, 1, z$writefreq)
				Report(paste('Phase ', z$Phase, ' Subphase ', subphase,
						' Iteration ', z$nit,' Progress: ',
						round((increment + val) / z$pb$pbmax * 100),
						'%\n', sep = ''))
				z$pb <- setProgressBar(z$pb, val + increment)
			}
			else
			{
				if (z$nit>1)
				{
					DisplayDeviations(z, fra)
				}
				if  (z$nit %% z$writefreq == 0)
				{
					val <- getProgressBar(z$pb)
					z$pb <- setProgressBar(z$pb, val + z$writefreq)
				}
			}
		}
		## setup check for end of phase. either store or calculate
		if (z$nit %% 2 == 1)
		{
			prev.fra <- fra
		}
		else
		{
			z$prod0 <- z$prod0 + fra * fra
			z$prod1 <- z$prod1 + fra * prev.fra
			ac <- ifelse (z$prod0 > 1e-12, z$prod1 / z$prod0, -2)
			z$maxacor <- max(-99, ac[!z$fixed]) ##note -2 > -99
			z$minacor <- min(1, ac[(!z$fixed) & ac > -1.0])
			z$ac <- ac
			if  (z$nit %% z$writefreq == 0)
			{
				DisplayThetaAutocor(z)
			}
		}
		## limit change.  Reporting is delayed to end of phase.
		# The truncation has been different from version 1.1-227 to 1.1-243,
		# due to a misunderstanding.
		# In version 1.1-244 it was changed back to the old Siena 3 way.
		# Except now the threshold is 5 instead of 10.
		# In version 1.1-285 the threshold is made a parameter from sienaModelCreate,
		# still with default 5;
		# and for the case !x$diagg, a multivariate truncation is used.

		if (x$diagg)
		{
			maxRatio <- max(ifelse(z$fixed, 1.0, abs(fra)/ z$sd), na.rm=TRUE)
			# note: this is a number, not a vector
		}
		else
		{
			if (z$pp > sum(z$fixed))
			{
			  if(!z$gmm)
			  {
			    maxRatio <-
			      max(sqrt((t(fra[!z$fixed]) %*% z$sf.invcov %*% fra[!z$fixed]) /
			                 (z$pp-sum(z$fixed))))			    
			  }
			  else
			  {
			    maxRatio <-
			      max(sqrt((t(fra[!z$fixed & !z$gmmEffects]) %*% z$sf.invcov %*% fra[!z$fixed & !z$gmmEffects]) /
			                 (z$pp-sum(!z$fixed & !z$gmmEffects)))) 
			   }
				# the max() is just to turn the (1,1) matrix into a number.
			}
			else
			{
				maxRatio <- 1.0
			}
		}
		if ((is.na(maxRatio)) || (is.nan(maxRatio)))
		{
			maxRatio <- 1.0
		}
		if (maxRatio > x$truncation)
		{
			fra <- x$truncation*fra/maxRatio
			z$truncated[z$nit] <- TRUE
		}
		if (subphase > x$doubleAveraging)
		{
			sumfra <- sumfra + fra
			fra <- sumfra
		}

#Report(paste("fra: ", "\n"), cf)
#PrtOutMat(as.matrix(fra), cf)

		if (x$standardizeVar)
		{
			if (x$diagg)
			{
				changestep <- fra / z$sd
			}
			else
			{
				if (x$standardizeWithTruncation)
				{
				  if (!z$gmm)
				  {
				    changestep <- (as.vector(z$dinvv %*% fra))*pmin(z$standardization, 1)
				  }
				  else
				  {
				    change <- (as.vector(z$dinvv %*% fra))*pmin(z$standardization, 1)
				    changestep <- rep(0,z$pp)
				    changestep[!z$gmmEffect] <- change
				  }
				}
				else
				{
				  if (!z$gmm)
				  {
				    changestep <- (as.vector(z$dinvv %*% fra))*z$standardization
				  }
				  else
				  {
				    changestep <- (as.vector(z$dinvv %*% fra))*z$standardization
				    changestep <- rep(0,z$pp)
				    changestep[!z$gmmEffect] <- change
		
				  }
				}
			}
		}
		else
		{
			if (x$diagg & !z$gmm)
			{
				changestep <- fra / diag(z$dfra)
			}
		  	else if (z$gmm)
		  	{
			    change <- (as.vector(z$dinvv %*% fra))
			    changestep <- rep(0,z$pp)
			    changestep[!z$gmmEffect] <- change 
			}
			else
			{
				changestep <- as.vector(z$dinvv %*% fra)
			}
		}
		changestep[z$fixed] <- 0.0
		fchange <- as.vector(z$gain * changestep)

		## check positivity restriction
		fchange[is.na(fchange)] <- 0
		z$positivized[fchange > z$theta] <- z$positivized[fchange > z$theta] +1
		z$positivized[!z$posj] <- 0
		fchange <- ifelse(z$posj & (fchange > z$theta), z$theta * 0.5, fchange)
		# make update step
		if (subphase > x$doubleAveraging)
		{
			zsmall$theta <- (z$thav/z$thavn) - fchange
		}
		else
		{
				zsmall$theta <- zsmall$theta - fchange
		}
		z$theta <- zsmall$theta
		z$thav <- z$thav + zsmall$theta
		z$thavn <- z$thavn + 1
 
#Report(paste("thavs: ", round(z$thavn), "\n"), cf)
#PrtOutMat(as.matrix(z$thav), cf)
		
		if (any(!z$fixed))
		{
			if (max(abs(z$theta[!z$fixed])) > z$thetaBound)
			{
				cat("The update steps have led to a parameter with maximum absolute value", 
						max(abs(z$theta[!z$fixed])), 
						",\nwhich is larger than thetaBound =", z$thetaBound, ".\n")
				larger <- rep("", length(z$theta))
				larger[((!z$fixed)&(abs(z$theta > z$thetaBound)))] <- " *****"
				print(cbind(z$requestedEffects$effectName, round(z$theta, 4), larger), quote=FALSE)
				if (interactive())
				{
					cat("If you wish to continue estimation in this session,")
					cat("\ngive a higher value for thetaBound.\n")
					thetaBound0 <- z$thetaBound
					z$thetaBound <- as.numeric(readline(prompt="Give a number: "))
					if (is.na(z$thetaBound)) 
					{
						stop("You gave no numeric answer, therefore the estimation stops.")
					}
					if (z$thetaBound > thetaBound0)
					{
						cat("OK, estimation continues.\n")
						flush.console()
					}
					else
					{
						stop("You gave a lower number, therefore the estimation stops.")
					}
				}
				else
				{
					stop("thetaBound should be set higher.")
				}
			}
		}

# This is a hidden option: 
# it is activated if the algorithm object x has a component x$moreUpdates
		if (x$maxlike && !is.null(x$moreUpdates) && x$moreUpdates > 0)
		{
			z <- doMoreUpdates(z, x, x$moreUpdates * subphase)
			zsmall$theta <- z$theta
		}
	
		if (ifelse(is.null(x$phase2imp),FALSE,x$phase2imp))
		{
			z$thetaList[[z$nCount]] <- z$theta
			z$nCount <- z$nCount + 1
		}
		
		##check for user interrupt
		CheckBreaks()
		if (UserInterruptFlag() || UserRestartFlag() || EarlyEndPhase2Flag())
		{
			break
		}
		## do we stop?
		if (!(is.na(z$minacor) || is.na(z$maxacor)))
		{
			if ((z$nit >= z$n2min && z$maxacor < 1e-10)||
				(z$nit >= z$n2max)
			|| ((z$nit >= 50) && (z$minacor < -0.8) &&
				(z$repeatsubphase < z$maxrepeatsubphase)))
			{
				break
			}
		}
	}
	if (ifelse(is.null(x$phase2imp),FALSE,x$phase2imp))
	{
		cat(z$nit)
		z$nAll <- c(z$nAll, z$nit)
		if (subphase==x$nsub)
		{
			z <- stdError(z,x,subphase=4, ...)
		}
	}
	
	if (z$returnThetas)
	{
		z$thetas <- rbind(z$thetas, thetas)
		z$sfs <- rbind(z$sfs, sfs)
	}
	z
}

