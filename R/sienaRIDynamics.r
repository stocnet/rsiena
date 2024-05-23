#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena/
# *
# * File: sienaRIDynamics.r
# *
# * Description: Used to determine, print, and plots relative importances of effects
# * in sequences of simulated micro-steps
# *****************************************************************************/

##@sienaRIDynamics
sienaRIDynamics <- function(data, ans=NULL, theta=NULL, algorithm=NULL, effects=NULL, depvar=NULL, intervalsPerPeriod=NULL)
{
	if(length(data$depvars)>1){
		if(is.null(depvar)){
			currentNetName <- attr(data$depvars,"names")[1]
			warning("If the models contains more than one dependent variables, \n it should be specified by variable 'depvar' for which dependent variable the relative importances should be calculated. \n\n")
			warning(paste("As 'depvar = NULL', relative importances are calculated for variable ", currentNetName, sep=""))
		}else if(!(depvar %in% attr(data$depvars,"names"))){
			stop("'depvar' is not a name of a dependent variable")
		}else{
			currentNetName <- depvar
		}
	}else{
		currentNetName <- attr(data$depvars,"names")[1]
	}
	if (!inherits(data, "siena"))
	{
		stop("data is not a legitimate Siena data specification")
	}
	if(!is.null(ans))
	{
		if (!inherits(ans, "sienaFit"))
		{
			stop("ans is not a legitimate Siena fit object")
		}
		if(ans$cconditional){
			if(is.null(effects))
			{
				stop("effects = NULL! In case of conditional estimation, the 'sienaEffects' object has to be given to function sienRIDynamics directly.")
			}
			if(!inherits(effects, "sienaEffects"))
			{
				stop("effects is not a legitimate Siena effects object")
			}
			effs <- effects
		}else{
			if(!is.null(algorithm)||!is.null(theta)||!is.null(effects))
			{
				warning("some information are multiply defined \n results will be based on 'theta', 'algorithm', and 'effects' stored in 'ans' (as 'ans$theta', 'ans$x', 'ans$effects')")
			}
			effs <- ans$effects
		}
		if(!is.null(intervalsPerPeriod))
		{
			if(is.numeric(intervalsPerPeriod))
			{
				intervalsPerPeriod <- as.integer(intervalsPerPeriod)
			}else{
				intervalsPerPeriod <- NULL
				warning("'intervalsPerPeriod' has to be of type 'numeric' \n used default settings")
			}
		}
		if(is.null(intervalsPerPeriod))
		{
			intervalsPerPeriod <- 10
		}
		RIValues <- calculateRIDynamics(data = data, theta= c(ans$rate,ans$theta), algorithm = ans$x,  effects = effs, depvar = currentNetName, intervalsPerPeriod=intervalsPerPeriod)
	}else{
		if (!inherits(algorithm, "sienaAlgorithm"))
		{
			stop("algorithm is not a legitimate Siena algorithm specification")
		}
		algo <- algorithm
		if (!inherits(effects, "sienaEffects"))
		{
			stop("effects is not a legitimate Siena effects object")
		}
		effs <- effects
		if(!is.numeric(theta))
		{
			stop("theta is not a legitimate parameter vector")
		}
		if(length(theta) != sum(effs$include==TRUE))
		{
			if(length(theta) == sum(effs$include==TRUE & effs$type!="rate"))
			{
				stop("vector of model parameters has wrong dimension, maybe rate parameters are missing")
			}
			stop("theta is not a legitimate parameter vector \n number of parameters has to match number of effects")
		}
		paras <- theta
		if(!is.null(intervalsPerPeriod))
		{
			if(is.numeric(intervalsPerPeriod))
			{
				intervalsPerPeriod <- as.integer(intervalsPerPeriod)
			}else{
				intervalsPerPeriod <- NULL
				warning("'intervalsPerPeriod' has to be of type 'numeric' \n used default settings")
			}
		}
		if(is.null(intervalsPerPeriod))
		{
			intervalsPerPeriod <- 10
		}
		## all necessary information available
		RIValues <- calculateRIDynamics(data = data, theta= paras, algorithm = algo,  effects = effs, depvar = currentNetName, intervalsPerPeriod=intervalsPerPeriod)
	}
	RIValues
}

##@calculateRIDynamics calculateRIDynamics simulates sequences of micro-steps, and aggregates the relative importances of effects over micro-steps of same time intervals
calculateRIDynamics <- function(data=NULL, theta=NULL, algorithm=NULL, effects=NULL, depvar=NULL, intervalsPerPeriod=10, returnActorStatistics=NULL)
{
	x <- algorithm
	currentNetName <- depvar
	z  <-  NULL
	z$FRAN <- getFromNamespace(x$FRANname, pkgname)
	x$cconditional <-  FALSE
	z$print <- FALSE
	z$Phase <- 3
	z <- initializeFRAN(z, x, data, effects, prevAns=NULL, initC=FALSE, returnDeps=FALSE)
	z$returnChangeContributions <- TRUE
	z$theta <- theta
	if (!is.null(x$randomSeed))
	{
		set.seed(x$randomSeed, kind="default")
	}
	else
	{
		if (exists(".Random.seed"))
		{
			rm(.Random.seed, pos=1)
			RNGkind(kind="default")
		}
	}
	chains <- x$n3
	periods <- data$observation-1
	effects <- effects[effects$include==TRUE,]
	noRate <- effects$type != "rate"
	thetaNoRate <- theta[noRate]
	effectNames <- effects$shortName[noRate]
	effectTypes <- effects$type[noRate]
#	networkName <- effects$name[noRate]
	networkInteraction <- effects$interaction1[noRate]
	effectIds <- paste(effectNames,effectTypes,networkInteraction, sep = ".")
	currentNetObjEffs <- effects$name[noRate] == currentNetName
# The old code leading to a crash...
	RIintervalValues <- list()
# ...threfore, the following line is added
	blanks <- matrix(0, sum(currentNetObjEffs),intervalsPerPeriod)
	for(period in 1:periods)
	{
		RIintervalValues[[period]] <- data.frame(blanks, row.names = effectIds[currentNetObjEffs])
	}
	for (chain in (1:chains))
	{
#cat("The following line leads to an error\n")
#browser()
		ans <- z$FRAN(z, x)
		for(period in 1:periods)
		{
			microSteps <- length(ans$changeContributions[[1]][[period]])
			stepsPerInterval <- microSteps/intervalsPerPeriod
# Tom's new code TS (7 lines).... Natalie slightly modified the next 4 lines
# and moved them one loop higher (out of the "chains"-loop)
#			blanks <- matrix(NA, sum(currentNetObjEffs),
#									ceiling(microSteps/stepsPerInterval))
#			RIintervalValues <-
#						as.list(rep(data.frame(blanks,
#							row.names = effectIds[currentNetObjEffs]), periods))
			blanks <- matrix(0, sum(currentNetObjEffs), intervalsPerPeriod)
			RItmp <- data.frame(blanks, row.names = effectIds[currentNetObjEffs])
			interval <- 1
			stepsInIntervalCounter <- 0
			for(microStep in 1:microSteps)
			{
				#currentNetName <- attr(ans$changeContributions[[1]][[period]][[microStep]],"networkName")
				if(attr(ans$changeContributions[[1]][[period]][[microStep]],
												"networkName")==currentNetName)
				{
					stepsInIntervalCounter <- stepsInIntervalCounter + 1
					## distributions[1,] contains
					## the probabilities of the available choices in this micro-step
					## distributions[2,],distributions[3,],distributions[4,], ...contains
					## the probabilities of the available choices
					## if the parameter of the first, second, third ... effects is set to zero.
					distributions <- calculateDistributions(
						ans$changeContributions[[1]][[period]][[microStep]],
						thetaNoRate[currentNetObjEffs])
					## If one wishes another measure
					## than the L^1-difference between distributions,
					## here is the right place
					## to call some new function instead of "L1D".
					RItmp[,stepsInIntervalCounter] <-
						standardize(L1D(distributions[1,],
									distributions[2:dim(distributions)[1],]))
				}
				if (microStep >
						interval * stepsPerInterval || microStep == microSteps)
				{
					if(chain == 1)
					{
						RIintervalValues[[period]][,interval] <-
										rowSums(RItmp)/length(RItmp)
					}
					else
					{
						RIintervalValues[[period]][,interval] <-
							RIintervalValues[[period]][,interval] +
											rowSums(RItmp)/length(RItmp)
					}
					interval <- interval + 1
					stepsInIntervalCounter <- 0
				}
			}
		}
	}
	for(period in 1:periods)
	{
		RIintervalValues[[period]] <- RIintervalValues[[period]]/chains
	}
	RIDynamics <- NULL
	RIDynamics$intervalValues <-RIintervalValues
	RIDynamics$dependentVariable <- currentNetName
	RIDynamics$effectNames <- paste(effectTypes[currentNetObjEffs], " ",
						(effects$effectName[noRate])[currentNetObjEffs], sep="")
	class(RIDynamics) <- "sienaRIDynamics"
	attr(RIDynamics, "version") <- packageDescription(pkgname, fields = "Version")
	RIDynamics
}

standardize <- function(values = NULL)
{
	newValues <- values/sum(values)
	newValues
}


##@print.sienaRIDynamics Methods
print.sienaRIDynamics <- function(x, ...){
	if (!inherits(x, "sienaRIDynamics"))
	{
		stop("x is not of class 'sienaRIDynamics' ")
	}
	periods <- length(x$intervalValues)
	intervals <- length(x$intervalValues[[1]])
	effs <- length(x$effectNames)
	cat(paste("\n  Relative importance of effects in micro-steps of dependent variable '", x$dependentVariable,"'. \n \n",sep=""))
	cat(paste("  Periods between observations are devided into ", intervals, " intervals. \n \n",sep=""))
	cat(paste("  Displayed results are aggregations over intervals:\n", sep=""))
	for(p in 1:periods){
		cat(paste("\n\nPeriod ",p,":\n", sep=""))
		colNames = paste("int ", 1:intervals, sep="")
		line1 <- paste(format("", width =52), sep="")
		line2 <- paste(format(1:effs,width=3), '. ', format(x$effectNames, width = 45),sep="")
		col <- 0
		for(i in 1:length(colNames))
		{
			col <- col + 1
			line1 <- paste(line1, format(colNames[i], width=8),"  ", sep = "")
			line2 <- paste(line2, format(round(x$intervalValues[[p]][[i]], 4), width=8, nsmall=4),"  ",sep="")
			if(col == 5)
			{
				line2 <- paste(line2, rep('\n',effs), sep="")
				cat(as.matrix(line1),'\n \n', sep='')
				cat(as.matrix(line2),'\n', sep='')
				line1 <- paste(format("", width =52), sep="")
				line2 <- paste(format(1:effs,width=3), '. ', format(x$effectNames, width = 45),sep="")
				col<-0
			}
		}
		if(col>0)
		{
			line2 <- paste(line2, rep('\n',effs), sep="")
			cat(as.matrix(line1),'\n \n', sep='')
			cat(as.matrix(line2),'\n', sep='')
		}
	}
	invisible(x)
}
#
##@summary.sienaRIDynamics Methods
summary.sienaRIDynamics <- function(object, ...)
{
	if (!inherits(x, "sienaRIDynamics"))
	{
		stop("x is not of class 'sienaRIDynamics' ")
	}
	class(x) <- c("summary.sienaRIDynamics", class(x))
	x
}
##@print.summary.sienaRIDynamics Methods
print.summary.sienaRIDynamics <- function(x, ...)
{
	if (!inherits(x, "summary.sienaRIDynamics"))
	{
		stop("not a legitimate summary of a 'sienaRIDynamics' object")
	}
	print.sienaRIDynamics(x)
	invisible(x)
}


##@plot.sienaRIDynamics Methods
plot.sienaRIDynamics <- function(x, staticRI = NULL, col = NULL,
			ylim=NULL, width = NULL, height = NULL, legend = TRUE,
			legendColumns = NULL, legendHeight = NULL, cex.legend = NULL, ...)
{
	if(!inherits(x, "sienaRIDynamics"))
	{
		stop("x is not of class 'sienaRIDynamics' ")
	}
	if(!is.null(staticRI))
	{
		if(!inherits(staticRI, "sienaRI"))
		{
			warning("staticRI is not of class 'sienaRI' and is therefore ignored")
			staticValues <- NULL
		}
		else
		{
			if(staticRI$dependentVariable != x$dependentVariable)
			{
				warning("staticRI does not correspond to x and is therefore ignored,\n staticRI and x do not refer to the same dependent variable")
				staticValues <- NULL
			}
			else if(length(staticRI$effectNames) != length(x$effectNames) || sum(x$effectNames==staticRI$effectNames)!=length(x$effectNames))
			{
				warning("staticRI does not correspond to x and is therefore ignored,\n staticRI and x do not refer to the same effects")
				staticValues <- NULL
			}
			else if(length(staticRI$expectedRI)-1 != length(x$intervalValues))
			{
				warning("staticRI does not correspond to x and is therefore ignored,\n staticRI and x do not refer to the same number of observation moments")
				staticValues <- NULL
			}
			else
			{
				staticValues <- staticRI$expectedRI
			}
		}
	}
	else
	{
		staticValues <- NULL
	}
	if(legend)
	{
		if(!is.null(legendColumns))
		{
			if(is.numeric(legendColumns))
			{
				legendColumns <- as.integer(legendColumns)
			}else{
				legendColumns <- NULL
				warning("legendColumns has to be of type 'numeric' \n used default settings")
			}
		}
		if(is.null(legendColumns))
		{
			legendColumns <- 2
		}
		if(!is.null(legendHeight))
		{
			if(is.numeric(legendHeight))
			{
				legendHeight <- legendHeight
			}else{
				legendHeight <- NULL
				warning("legendHeight has to be of type 'numeric' \n used default settings")
			}
		}
		if(is.null(legendHeight))
		{
			legendHeight <- 1
		}
	}
	if(!is.null(height))
	{
		if(is.numeric(height))
		{
			height <- height
		}else{
			height <- NULL
			warning("height has to be of type 'numeric' \n used default settings")
		}
	}
	if(is.null(height))
	{
		if(legend)
		{
			height <- 4
		}else{
			height <- 3
		}
	}

	if(!is.null(width))
	{
		if(is.numeric(width))
		{
			width <- width
		}else{
			width <- NULL
			warning("width has to be of type 'numeric' \n used default settings")
		}
	}
	if(is.null(width))
	{
		width <- 8
	}

	if(!is.null(cex.legend))
	{
		if(is.numeric(cex.legend))
		{
			cex.legend <- cex.legend
		}else{
			cex.legend <- NULL
			warning("cex.legend has to be of type 'numeric' \n used default settings")
		}
	}
	if(is.null(cex.legend))
	{
		cex.legend <- 1
	}

	if(!is.null(col))
	{
		cl <- col
	}
	if(is.null(col))
	{
		alph <- 175
		green <- rgb(127, 201, 127,alph, maxColorValue = 255)
		lila <-rgb(190, 174, 212,alph, maxColorValue = 255)
		orange <- rgb(253, 192, 134,alph, maxColorValue = 255)
		yellow <- rgb(255, 255, 153,alph, maxColorValue = 255)
		blue <- rgb(56, 108, 176,alph, maxColorValue = 255)
		lightgray <- rgb(184,184,184,alph, maxColorValue = 255)
		darkgray <- rgb(56,56,56,alph, maxColorValue = 255)
		gray <- rgb(120,120,120,alph, maxColorValue = 255)
		pink <- rgb(240,2,127,alph, maxColorValue = 255)
		brown <- rgb(191,91,23,alph, maxColorValue = 255)
		cl <- c(green,lila,orange,yellow,blue,lightgray,darkgray,gray,pink,brown)
		while(length(cl)<length(x$effectNames)){
			alph <- (alph+75)%%255
			green <- rgb(127, 201, 127,alph, maxColorValue = 255)
			lila <-rgb(190, 174, 212,alph, maxColorValue = 255)
			orange <- rgb(253, 192, 134,alph, maxColorValue = 255)
			yellow <- rgb(255, 255, 153,alph, maxColorValue = 255)
			blue <- rgb(56, 108, 176,alph, maxColorValue = 255)
			lightgray <- rgb(184,184,184,alph, maxColorValue = 255)
			darkgray <- rgb(56,56,56,alph, maxColorValue = 255)
			gray <- rgb(120,120,120,alph, maxColorValue = 255)
			pink <- rgb(240,2,127,alph, maxColorValue = 255)
			brown <- rgb(191,91,23,alph, maxColorValue = 255)
			cl <- c(cl,green,lila,orange,yellow,blue,lightgray,darkgray,gray,pink,brown)
		}
	}
#	bordergrey <-"gray25"
	values <- x$intervalValues
	periods <- length(values)
	effectNames <- x$effectNames
	effectNumber <- length(effectNames)
	lineTypes <- rep("solid",effectNumber)
	legendNames <- effectNames
	if(is.null(ylim))
	{
		ylim <-c(0,ceiling(max(unlist(lapply(values, max)))*10)*0.1)
	}else{
		ylim <-ylim
	}
	if(legend)
	{
		layout(rbind(1:periods, rep(periods+1,periods)),
				widths=rep(4, periods),heights=c(3,legendHeight))
	}else{
		layout(rbind(1:periods),widths=rep(4, periods),heights=c(3))
	}
#	par( oma = c( 1, 3, 1, 3 ),mar = par()$mar+c(-5,-4.1,-4,-2.1), xpd=T )
	par( oma = c( 1, 3, 1, 3 ),mar = c(0.1, 0.1, 0.1, 0.1), xpd=T )
	for(period in 1:periods){
		timeseries<-ts(t(values[[period]]))
		plot.ts(timeseries, plot.type = "single",  col = cl,
				lty = lineTypes, lwd = rep(1.5,effectNumber), bty = "n",
				xaxt = "n",yaxt = "n", ylab ="", xlab = "", ylim = ylim)
		for(eff in 1:effectNumber)
		{
			points(ts(t(values[[period]]))[,eff], col = cl[eff], type = "p", pch = 20)
			if(!is.null(staticValues))
			{
				points(xy.coords(1,staticValues[[period]][eff]),col = cl[eff],
						type = "p", pch = 1, cex = 1.75)
			}
		}
		ax <- ((ylim[1]*10):(ylim[2]*10))/10
		if(period==1){
			if(period == periods){
				axis(4, ax)
				axis(2, ax)
			}else{
				axis(4, ax, labels = FALSE)
				axis(2, ax)
			}
		}else if(period == periods){
			axis(4, ax)
			axis(2, ax, labels = FALSE)
		}else{
			axis(2, ax, labels = FALSE)
			axis(4, ax, labels = FALSE)
		}
	}
	if(legend)
	{
		plot(c(0,1), c(0,1), col=rgb(0,0,0,0),axes=FALSE,  ylab = "", xlab = "")
		legend(0, 1, legendNames, col = cl[1:effectNumber], lwd = 2,
				lty = lineTypes, ncol = legendColumns, bty = "n",cex=cex.legend)
	}
	invisible(cl)
}


