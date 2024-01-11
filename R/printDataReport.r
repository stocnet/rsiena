# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: printDatareport.r
# *
# * Description: This module contains the function to produce the data
# * report from siena07
# *
# *****************************************************************************/
##@DataReport siena07 Print report
DataReport <- function(z, x, f)
{
	if (!is.null(z$effectsName))
	{
		Report(c("Effects object used:", z$effectsName, "\n"), outf)
	}

	if (z$maxlike)
	{
		Report(c(z$nrunMH,
				"MCMC steps per RM step (multiplication factor =",
				x$mult), outf)
		Report(")\n", outf)
	}
	## f could be a group, but has attributes like a group even if not!
	oneMode <- attr(f, "types") == "oneMode"
	bipartite <- attr(f, "types") == "bipartite"
	behavior <- attr(f, "types") == "behavior"
	nOneMode <- sum(oneMode)
	nBehavior <- sum(behavior)
	nBipartites <- sum(bipartite)
	oneModeNames <- attr(f, "netnames")[oneMode]
	behaviorNames <- attr(f, "netnames")[behavior]
	##  behaviorNames <- attr(f, "netnames")[behavior]
	##  bipartiteNames <- attr(f, "netnames")[bipartite]
	nDepVars <- nOneMode + nBehavior + nBipartites
	observations <- attr(f, "observations") ##note this is total number of
	##  periods to process
	exogenous <- attr(f, 'compositionChange')
	exoOptions <- attr(f, 'exooptions')
	if (nOneMode > 0)
	{
		Report("Model Type:\n", outf)
		for (i in 1:nOneMode)
		{
			if (nOneMode > 1)
			{
				Report(sprintf("Network %d %s :", i, oneModeNames[i]), outf)
			}
			Report(sprintf(" %s \n",ModelTypeStrings(z$modelType[i])), outf)
		}
		if (any(z$modelType != x$modelType))
		{
			Report(c("Note that the model type requested has been",
					"over-ridden\n"), outf)
		}
	}
	if (nBehavior > 0)
	{
		Report("Behavioral Model Type:\n",outf)
		for (i in 1:nBehavior)
		{
			if (nBehavior > 1)
			{
				Report(sprintf("Behavior %d %s :", i, behaviorNames[i]), outf)
			}
			Report(sprintf(" %s\n",BehaviorModelTypeStrings(z$behModelType[i])), outf)
		}
	}
	if (x$cconditional  != z$cconditional)
	{
		Report("\nNB. Request for conditional estimation has been over-ridden.\n\n",
			outf)
	}
	Report("Estimation method: ", outf)
	if (z$cconditional)
	{
		Report("conditional moment estimation", outf)
		if (gmm(x))
		{
			Report(' by the Generalized Method of Moments.\n', outf)
		}
		else
		{
			Report('.\n', outf)
		}
		Report(c('Conditioning variable is the total number of observed',
				'changes ("distance") \n'), outf)
		if (z$condtype == 'oneMode')
		{
			if (nOneMode == 1)
			{
				Report('in the network variable.\n', outf)
			}
			else
			{
				Report(c("in network variable", x$condname, '\n'), outf)
			}
		}
		else
		{
			Report(c("in behavioral variable", x$condname, '\n'), outf)
		}
		if (gmm(x))
		{
			Report('in the network variable.\n', outf)
		}
		if (observations == 1)
		{
			## note can only be one group
			Report(c("Distance for simulations is",
					format(f[[1]]$distances, width=4), '.\n'), outf)
		}
		else
		{
			Report("Distances for simulations are\n", outf)
			#            Report(c("period   :", format(1:observations, width=4),
			Report(c("period   :", sapply(f[1], function(x)
						sapply(x$depvars, function(y)
							format(1:observations, width=4))),
					'\n'), sep='', outf)
			Report(c("distance : ", format(sapply(f[1],
							function(x)x$distances),
						width = 4), '.\n'),sep='',  outf)
		}
	}
	else if (!z$maxlike)
	{
		Report("unconditional moment estimation\n", outf)
		if (gmm(x))
		{
			Report(' by the Generalized Method of Moments.\n', outf)
		}
		else
		{
			Report('.\n', outf)
		}

		#if (exogenous)
		#{
		#    Report("Changing composition: no conditional moment estimation.\n",
		#           outf)
		#    ## TODO report separately for different node sets if appropriate
		#    Report(c("Joiners/leavers option:", exoOptions, "\n"), sep="",
		#           outf)
		#}
		Report('\nTime duration for simulations ', outf)
		if (observations >= 2)
		{
			Report('in each period ', outf)
		}
		Report('is 1.0.\n', outf)
	}
	else
	{
		Report("Maximum likelihood estimation\n", outf)
	}
	if (exogenous)
	{
		Report("Changing composition.\n",  outf)
		## TODO report separately for different node sets if appropriate
		Report(c("Joiners/leavers option: ", exoOptions, "\n"), sep="",
			outf)
	}
	if (z$FinDiff.method)
	{
		if ((x$nsub == 0)&(x$simOnly))
		{
			Report(c('Derivatives are estimated with the finite difference',
					'method.\n'), outf)
		}
		else
		{
			Report(c('Standard errors are estimated with the finite difference',
					'method.\n'), outf)
		}
		if (!x$FinDiff.method)
		{
			Report("Note that the option requested has been over-ridden\n",
				outf)
		}
	}
	else
	{
		if ((x$nsub == 0)&(x$simOnly))
		{
			Report(c('Derivatives are estimated with the likelihood ratio',
					'method.\n'), outf)
		}
		else
		{
			Report(c('Standard errors are estimated with the likelihood ratio',
					'method.\n'), outf)
		}
	}
	if (x$dolby)
	{
		Report('Dolby method (regression on scores) is used.\n', outf)
	}
	## any max degree restrictions?
	if (any(x$MaxDegree > 0))
	{
		for (i in 1:length(x$MaxDegree))
		{
			if (x$MaxDegree[i])
			{
				mdnet <- names(x$MaxDegree)[i]
				Report(c("Dependent network variable", mdnet, ':\n'), outf)
				maxod <- max(
					attr(f$Data1$depvars[[match(mdnet, attr(f, "netnames"))]],
						"maxObsOutDegree"))
				if (maxod > x$MaxDegree[i])
				{
					Report(c("The algorithm object requires outdegrees not",
							"larger than", x$MaxDegree[i], '\n',
							"but the maximum observed outdegree is", maxod,
							".\n"), outf)
					Report("This is incompatible.\n", outf)
					cat('Some observed outdegrees higher than MaxDegree.\n')
					stop("Incompatibility between data and MaxDegree in algorithm object.")
				}
				if (attr(f, 'symmetric')[match(mdnet, attr(f, "netnames"))])
				{
					Report(c("All graphs are constrained to having degrees not",
							"larger than", x$MaxDegree[i], '.\n'), outf)
				}
				else
				{
					Report(c("All digraphs are constrained to having",
							"outdegrees not larger than", x$MaxDegree[i],
							'.\n'), outf)
				}

			}
		}
	}
	## any universal offsets?
	if (any(x$UniversalOffset != 0))
	{
		Report("Offset for universal setting or symmetric type-3 model \n", outf)
		for (i in 1:length(x$UniversalOffset))
		{
			if (x$UniversalOffset[i] != 0)
			{
				mdnet <- names(x$UniversalOffset)[i]
				Report(c("Dependent network variable", mdnet, ': '), outf)
				Report(c(sprintf("%9.4f",x$UniversalOffset[i]),'\n'), outf)
			}
		}
		Report("\n", outf)
	}
	if ((x$nsub == 0)&(x$simOnly))
	{
		Report('Parameter values are \n', outf)
	}
	else
	{
		Report(sprintf("Initial value of gain parameter is %10.7f.\n",
				x$firstg), outf)
		Report(sprintf("Reduction factor for gain parameter is %10.7f.\n",
				x$reduceg), outf)
		Report(sprintf("Number of subphases in Phase 2 is %d.\n\n", x$nsub),
			outf)
		Report('Initial parameter values are \n', outf)
	}
	if (z$cconditional)
	{
		if (observations == 1)
		{
			if (nDepVars == 1)
			{
				Report(format('  0. Rate parameter', width = 43), outf)
			}
			else
			{
				Report(format('  0. Rate parameter of conditioning variable',
						width = 43), outf)
			}
			Report(sprintf("%9.4f\n", attr(f, "condEffects")$initialValue[1]),
				outf)
		}
		else ## observations > 1
		{
			for (i in 1:observations)
			{
				if (nDepVars == 1)
				{
					Report(format(paste('  0.', format(format(i), width=2),
								' Rate parameter', sep=''),
							width = 43), outf)
				}
				else
				{
					Report(format(paste('  0.', format(format(i), width=2),
								' Rate parameter cond.',
								'variable period', i, sep=''),
							width = 43), outf)
				}
				Report(sprintf("%9.4f\n",
						attr(f, "condEffects")$initialValue[i]), outf)
			}
		}
	}

	if (!gmm(x))
	{
	  fixed <- ifelse(z$fixed, '  (fixed) ', '')
	}
	else
	{
	  if (!z$cconditional)
	  {
	    fixed <- ifelse((z$fixed & !z$gmmEffects), '  (fixed) ', '')
	  }
	  else
	  {
	    fixed <- ifelse((z$fixed & !(z$requestedEffects$type=="gmm")), '  (fixed) ', '')
	  }
	}
	if (!gmm(x))
	{	tmp <- paste(sprintf("%3d",1:length(z$requestedEffects$effectName)), '. ',
	               format(paste(z$requestedEffects$type, ':  ',
	                            z$requestedEffects$effectName,
	                            sep = ''), width = 52),
	               sprintf("%9.4f", z$requestedEffects$initialValue), fixed, '\n',
	               sep = '', collapse = '')
	}
	else
	{
	  tmp <- paste(sprintf("%3d",1:(length(z$requestedEffects$effectName)-sum(z$gmmEffects))), '. ',
	               format(paste(z$requestedEffects$type[!z$gmmEffects], ':  ',
	                            z$requestedEffects$effectName[!z$gmmEffects],
	                            sep = ''), width = 52),
	               sprintf("%9.4f", z$requestedEffects$initialValue[!z$gmmEffects]), fixed[!z$gmmEffects], '\n',
	               sep = '', collapse = '')
	}
	Report(tmp, outf)
	## targets:
	Report("\n\nValues of target statistics are\n", outf)
	if (z$maxlike)
	{
		targets <- z$maxlikeTargets
	}
	else
	{
		targets <- z$targets
	}
	tmp <- paste(sprintf("%3d",1:length(z$requestedEffects$effectName)), '. ',
		format(z$requestedEffects$functionName, width = 66),
		sprintf("%9.4f",
			ifelse(z$requestedEffects$type=='endow', -targets,
				targets)),
		'\n', sep = '', collapse = '')
	Report(tmp, outf)
	if (attr(z$targets, "fromData"))
	{
		Report("These were calculated from the data.\n", outf)
	}
	else
	{
		Report("These were given by the user in the call of siena07.\n", outf)
	}
	Report(c('\n', nrow(z$requestedEffects)-sum(z$requestedEffects$type=="gmm"), 'parameters,',
			nrow(z$requestedEffects),
			'statistics\n'),outf)
}
