#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena/
# *
# * File: sienaRI.r
# *
# * Description: Used to determine, print, and plots various effect sizes
# * for potential decisions of actors at observation moments.
# *****************************************************************************/

##@interpret_size
sienaRI <- function(data, ans=NULL, theta=NULL, effects=NULL,
	getChangeStats=FALSE)
{
	if (!inherits(data, "sienadata") && !inherits(data, "siena"))
	{
		stop("not a legitimate Siena data specification")
	}
	datatypes <- sapply(data$depvars, function(x){attr(x,"type")})
	if (any(datatypes == "bipartite"))
	{
#		stop("interpret_size works only for dependent variables of type 'oneMode' or 'behavior'")
	}
		if (!inherits(ans, "sienaFit"))
		{
			stop(paste("ans is not a legitimate Siena fit object", sep=""))
		}
		if (!is.null(theta)||!is.null(effects))
		{
			warning(paste("some informations are multiply defined \n",
					"results will be based on 'theta' and 'effects'\n",
					"stored in 'ans' (as 'ans$theta','ans$effects')\n", sep=""))
		}
		if (sum(ans$effects$include==TRUE &
				(ans$effects$type =="endow"|ans$effects$type =="creation")) > 0)
		{
			stop("interpret_size does not yet work for models containing endowment or creation effects")
		}
		if (any(ans$effects$shortName %in% c("unspInt", "behUnspInt", "contUnspInt"))){
			stop("interpret_size does not work for models containing interaction effects")
		}
		contributions <- getChangeContributions(data = data,
			effects = ans$effects)
		# contributions[[1]] is periods by effects by actors by actors
		RI <- expectedRelativeImportance(conts = contributions,
			effects = ans$effects, theta =ans$theta, thedata=data,
			getChangeStatistics=getChangeStats)
	RI$theta <- ans$theta[!ans$requestedEffects$basicRate]
	RI$se <- ans$se[!ans$requestedEffects$basicRate]
	RI$effectNames <- ans$requestedEffects$effectName[!ans$requestedEffects$basicRate]
	RI$depName <- ans$requestedEffects$name[!ans$requestedEffects$basicRate]
	class(RI) <- "siena_sizes"
	attr(RI, "version") <- packageDescription(pkgname, fields = "Version")
	RI
}



interpret_size <- function(x, ...) UseMethod("interpret_size", x)

##@interpret_size  method for sienaEffects objects
interpret_size.sienaEffects <- function(x, data,
		theta, getChangeStats=FALSE, ...)
{
    sienaRI(data, theta=theta, effects=x,
					getChangeStats=getChangeStats)
}

##@interpret_size  method for sienaFit objects
interpret_size.sienaFit <- function(x, data, getChangeStats=FALSE, ...)
{
    sienaRI(data=data, ans=x, getChangeStats=getChangeStats)
}

##@getChangeContributions. Use as RSiena:::getChangeContributions
getChangeContributions <- function(data, effects)
{
	# How to check for settings without algorithm? Is it relevant?
    setup <- sienaSetupForCpp(data, effects,
                              includeBehavior = TRUE,
                              includeBipartite = TRUE,
                              returnActorStatistics = FALSE,
                              returnStaticChangeContributions = TRUE,
                              parallelrun = FALSE)

    ans <- .Call(C_getTargets, PACKAGE=pkgname,
                 setup$pData, setup$pModel, setup$myeffects,
                 parallelrun = FALSE,
                 returnActorStatistics = FALSE,
                 returnStaticChangeContributions = TRUE)
    ans
}

expectedRelativeImportance <- function(conts, effects, theta, thedata=NULL,
	getChangeStatistics=FALSE, effectNames = NULL)
{
	rms <- function(xx){sqrt((1/dim(xx)[2])*rowSums(xx^2,na.rm=TRUE))}
	periods <- length(conts[[1]])
	effects <- effects[effects$include == TRUE,]
	noRate <- effects$type != "rate"
	effects <- effects[noRate,]
	if(sum(noRate)!=length(theta))
	{
		theta <- theta[noRate]
	}
	effectNa <- attr(conts,"effectNames")
	effectTypes <- attr(conts,"effectTypes")
	networkNames <- attr(conts,"networkNames")
	networkTypes <- attr(conts,"networkTypes")
	networkInteraction <- effects$interaction1
	effectIds <- paste(effectNa,effectTypes,networkInteraction, sep = ".")
	currentDepName <- ""
	depNumber <- 0
	RI <- list()
	for(eff in 1:length(effectIds))
	{
		if(networkNames[eff] != currentDepName)
		{
			currentDepName <- networkNames[eff]
			actors <- length(conts[[1]][[1]][[1]])
			depNumber <- depNumber + 1
			currentDepEffs <- effects$name == currentDepName
			effNumber <- sum(currentDepEffs)
			depNetwork <- thedata$depvars[[depNumber]]
			if (networkTypes[eff] == "oneMode")
			{
				choices <- actors
			}
			else if (networkTypes[eff] == "behavior")
			{
				choices <- 3
			}
			else if (networkTypes[eff] == "bipartite")
			{
				choices <- dim(depNetwork)[2] + 1
			}
			else
			{
				stop("interpret_size does not work for dependent variables of type
					'continuous'")
			}

			# impute for period 1
			if (networkTypes[eff] %in% c("oneMode", "bipartite"))
			{
				depNetwork[,,1][is.na(depNetwork[,,1])] <- 0
			}
			else
			{
				depNetwork[,,1][is.na(depNetwork[,,1])] <- attr(depNetwork, 'modes')[1]
			}
			# impute for next periods;
			# this may be undesirable for structurals immediately followed by NA...
			for (m in 2:dim(depNetwork)[3]){depNetwork[,,m][is.na(depNetwork[,,m])] <-
				depNetwork[,,m-1][is.na(depNetwork[,,m])]}
				# Make sure the diagonals are not treated as structurals
				if (networkTypes[eff] == "oneMode")
				{
					for (m in 1:(dim(depNetwork)[3])){diag(depNetwork[,,m]) <- 0}
				}
			structurals <- (depNetwork >= 10)
			if (networkTypes[eff] == "oneMode"){
				if (attr(depNetwork, 'symmetric')){
					message('\nNote that for symmetric networks,
						effect sizes are for modelType 2
						(forcing).')
					}
				}

			# currentDepObjEffsNames <- paste(effects$shortName[currentDepEffs],
			# 	effects$type[currentDepEffs],effects$interaction1[currentDepEffs],sep=".")
			# otherObjEffsNames <- paste(effects$shortName[!currentDepEffs],
			# 	effects$type[!currentDepEffs],effects$interaction1[!currentDepEffs],sep=".")

			entropies <- vector(mode="numeric", length = actors)
			expectedRI <- list()
			expectedI <- list()
			RIActors <- list()
			IActors <- list()
			absoluteSumActors <- list()
			RHActors <-list()
			changeStats <-list()
			sigma <- list()
			sigmas <- matrix(NA, sum(currentDepEffs), periods)
			if (networkTypes[eff] == "behavior")
			{
				toggleProbabilities <- array(0, dim=c(actors, 3, periods))
			}
			else
			{
				toggleProbabilities <- array(0, dim=c(actors, choices, periods))
			}
			for(w in 1:periods)
			{
				currentDepEffectContributions <- conts[[1]][[w]][currentDepEffs]
				if (networkTypes[eff] == "bipartite")
				{
					currentDepEffectContributions <- lapply(
							currentDepEffectContributions,
							function(x){lapply(x,function(xx){xx[1:choices]})})
				}
				# conts[[1]] is periods by effects by actors by actors
				currentDepEffectContributions <-
					sapply(lapply(currentDepEffectContributions, unlist),
						matrix, nrow=actors, ncol=choices, byrow=TRUE,
						simplify="array")
				cdec <- apply(currentDepEffectContributions, c(2,1), as.matrix)
				# cdec is effects by actors (alters) by actors (egos)
				if (dim(currentDepEffectContributions)[3] <= 1) # only one effect
				{
					cdec <- array(cdec, dim=c(1,dim(cdec)))
				}
				rownames(cdec) <- effectNa[currentDepEffs]
				if (getChangeStatistics)
				{
					changeStats[[w]] <- cdec
				}
				# replace structural 0s and 1s by NA,
				# so they are omitted from calculation of RI, R_H, sigma
				if (networkTypes[eff] == "oneMode")
				{
					#	structuralsw <- structurals[,,w]
					for (ff in 1:(dim(cdec)[1])){cdec[ff,,][t(structurals[,,w])] <- NA}
				}
				distributions <- apply(cdec, 3,
					calculateDistributions, theta[which(currentDepEffs)])

				distributions <-
					lapply(apply(distributions, 2, list),
						function(x){matrix(x[[1]], nrow=effNumber+1,
							ncol=choices, byrow=F)})
				# distributions is a list, length = number of actors
				# distributions[[i]] is for actor i, a matrix of dim (effects + 1)
				# * (actors as alters) giving the probability of toggling the
				# tie variable to the alters;
				# the first row is for the unchanged parameter vector theta,
				# each of the following has put one element of theta to 0.
				# for behavior it is a matrix of dim (effects + 1) * 3
				if (networkTypes[eff] == "behavior")
				{
					toggleProbabilities[,,w] <-
						t(vapply(distributions, function(x){x[1,]}, rep(0,3)))
				}
				else
				{
					toggleProbabilities[,,w] <-
						t(vapply(distributions, function(x){x[1,]}, rep(0,choices)))
				}
				# toggleProbabilities is an array referring to ego * choices *
				# period, giving the probability of ego in a ministep at the
				# start of the period to make this choice.
				# For oneMode networks choices are alters;
				# for  bipartite choices are second mode nodes,
				# and the last is "no change";
				# for behavior choices are to add -1, 0, +1.
				# The function "entropy" uses the property that non-choices
				# have probability NA; 
				# this is changed below to 0 for toggleProbabilities,
				# so this function should be applied to "distributions"
				# but not to "toggleProbabilities":
				entropy_vector <- unlist(lapply(distributions,
						function(x){entropy(x[1,])}))
				## If one wishes in the following another measure than the
				## L^1-difference between distributions, here is
				## the right place to call some new function instead of "L1D".
				RIs_list <- lapply(distributions,function(x){L1D(x[1,],
						x[2:dim(x)[1],])})
				RIs_matrix <-(matrix(unlist(RIs_list),nrow=effNumber,
						ncol=actors, byrow=F))
				entropies <- entropy_vector
				# divide by column sums:
				RIActors[[w]] <- t(t(RIs_matrix)/rowSums(t(RIs_matrix), na.rm=TRUE))
				rownames(RIActors[[w]]) <- effectNa[currentDepEffs]
				absoluteSumActors[[w]] <- colSums(RIs_matrix, na.rm=TRUE)
				RHActors[[w]] <- entropies
				expectedRI[[w]] <- rowMeans(RIActors[[w]], na.rm=TRUE)
				IActors[[w]] <- RIs_matrix
				rownames(IActors[[w]]) <- effectNa[currentDepEffs]
				expectedI[[w]] <- rowMeans(RIs_matrix, na.rm=TRUE)
				sigma[[w]] <- apply(cdec, c(1,3), sd, na.rm=TRUE)
				rownames(sigma[[w]]) <- effectNa[currentDepEffs]
				sigmas[,w] <- rms(sigma[[w]])
			}
			rownames(sigmas) <- effectNa[currentDepEffs]
			meansigmas <- rms(sigmas)
			names(meansigmas) <- effectNa[currentDepEffs]
			RItmp <- NULL
			RItmp$dependentVariable <- currentDepName
			RItmp$expectedRI <- expectedRI
			RItmp$RIActors <- RIActors
			RItmp$expectedI <- expectedI
			RItmp$IActors <- IActors
			RItmp$absoluteSumActors <- absoluteSumActors
			RItmp$RHActors <- RHActors
			RItmp$sigma <- sigma
			RItmp$sigmas <- sigmas
			RItmp$meansigmas <- meansigmas
			if(!is.null(effectNames))
			{
				RItmp$effectNames <- effectNames[currentDepEffs]
			}else{
				RItmp$effectNames <-
					paste(effectTypes[currentDepEffs], " ",
						effects$effectName[currentDepEffs], sep="")
			}
			if (getChangeStatistics){
				RItmp$changeStatistics <- changeStats
				toggleProbabilities[is.na(toggleProbabilities)] <- 0
				RItmp$toggleProbabilities <- toggleProbabilities
			}
				RI[[depNumber]]<-RItmp
			}
		}
	names(RI) <- sapply(RI, function(rr){rr$dependentVariable})
	RI
}

calculateDistributions <- function(effectContributions = NULL, theta = NULL)
{
	neffects <- dim(effectContributions)[1]
	nchoices <- dim(effectContributions)[2]
	distributions <- array(NA, dim = c(neffects+1,nchoices))
	the.choices <- !is.na(colSums(effectContributions))
	if (sum(the.choices) >= 2)
	{
		distributions[1,the.choices] <-
			exp(colSums(theta*effectContributions[,the.choices,drop=FALSE], na.rm=TRUE))/
			sum(exp(colSums(theta*effectContributions[,the.choices,drop=FALSE], na.rm=TRUE)))
		for(eff in 1:neffects)
		{
			th <- theta
			th[eff] <- 0
			distributions[eff+1,the.choices] <-
				exp(colSums(th*effectContributions[,the.choices,drop=FALSE], na.rm=TRUE))/
				sum(exp(colSums(th*effectContributions[,the.choices,drop=FALSE], na.rm=TRUE)))
		}
	}
	distributions
}

entropy <- function(distribution = NULL)
{
	if (sum(!is.na((distribution))) <= 1) # only constant choice
	{
		certainty <- NA
	}
	else
	{
		entropy <- -1*(sum(distribution * log(distribution), na.rm=TRUE)/
			log(sum(!is.na((distribution)))))
		certainty <- 1 - entropy
	}
	certainty
}

KLD <- function(referenz = NULL, distributions = NULL)
{
	if (sum(!is.na((referenz))) <= 1) # only constant choice
	{
		kld <- rep(NA, dim(distributions)[1])
	}
	else
	{
		if(is.vector(distributions))
		{
			kld <- sum(referenz * (log(referenz)-log(distributions)),
				na.rm=TRUE)/log(sum(!is.na((referenz))))
		}
		else
		{
			kld <- colSums(referenz * (log(referenz)-t(log(distributions))),
				na.rm=TRUE)/log(sum(!is.na((referenz))))
		}
	}
	kld
}

## calculates the L^1-differenz between distribution "reference"
## (which is a vector of length n)
## and each row of distributions (which is a matrix with n columns)
L1D <- function(referenz = NULL, distributions = NULL)
{
	if (sum(!is.na((referenz))) <= 1) # only constant choice
	{
		l1d <- rep(NA, dim(distributions)[1])
	}
	else
	{
		if(is.vector(distributions))
		{
			l1d <- sum(abs(referenz-distributions), na.rm=TRUE)
		}
		else
		{
			l1d <- colSums(abs(referenz-t(distributions)), na.rm=TRUE)
		}
	}
	l1d
}


##@print.siena_sizes Methods
print.siena_sizes <- function(x, ...){
	makenice <- function(r)
	{
		paste(format(round(r, 4), width=8, nsmall=4), "   ", sep="")
	}
	makenice <- Vectorize(makenice)
	divi <- function(a,b){ifelse(b==0, NA, a/b)}
	present.p <- function(p){
		if (is.na(p))
		{
			pval <- "  NA"
		}
		else
	{
			if (p < 0.001)
			{
				pval <- "<0.001"
			}
			else if (p < 0.10)
			{
				pval <- format(round(p, 3), width=5, nsmall=3)
		}
			else
			{
				pval <- paste(format(round(p, 2), width=4, nsmall=2), " ", sep="")
			}
		}
		pval
	}
	present.p <- Vectorize(present.p)
	pris0 <- function(xx)
	{
# This internal function handles results for each dependent variable
		cat("\n****    Dependent variable ", format(xx$dependentVariable, width=wid+2), "**** \n")
		effs <- dim(xx$sigmas)[1]
		periods <- dim(xx$sigmas)[2]
	colNames = paste("period ", 1:periods, sep="")
		line1 <- format("", width = wid+7)
	line2 <- paste(format(1:effs,width=3), '. ',
			format(xx$effectNames, width = wid),sep="")
	line3 <- line2
		line4 <- format("R_H ('degree of certainty')", width = wid+5)
	line5 <- line2
		line6 <- line2
		xth <- x$theta[x$depName == xx$dependentVariable]
		seth <- x$theta[x$depName == xx$dependentVariable]
		fixed <- is.na(seth)
		
		xth[fixed] <- NA
		RH_overall <- mean(vapply(xx$RHActors, mean, FUN.VALUE=1))
		sigma_overall  <- xx$meansigmas
		semisth <- sigma_overall * xth
		line3 <- paste(line3, makenice(cbind(xth, seth)), sep="")
		for (w in 1:periods)
	{
		line1 <- paste(line1, format(colNames[w], width=8),"  ", sep = "")
			line5 <- paste(line5, makenice(xx$sigmas[,w]), sep="")
			line6 <- paste(line6, makenice(xth * xx$sigmas[,w]), sep="")
			line4 <- paste(line4, makenice(mean(xx$RHActors[[w]], na.rm=TRUE)), sep="")
		}
		line1 <- paste(line1, "overall")
		line5 <- paste(line5, makenice(sigma_overall), rep('\n',effs), sep="")
		line6 <- paste(line6, makenice(semisth), rep('\n',effs), sep="")
		line4 <- paste(line4, makenice(RH_overall), sep="")		
		cat(paste("\n Observation moments:\n", sep=""))
		cat(as.matrix(line1),'\n', sep='')
		cat("\nSemi-standardized parameter estimates:\n")
		cat('\n',as.matrix(line6),'\n', sep='')
		cat("sigma (average within-ego standard deviation of change statistics):\n")
		cat('\n',as.matrix(line5),'\n', sep='')
		cat(as.matrix(line4),'\n', sep='')
		}

	if (!inherits(x, "siena_sizes"))
	{
		stop("not a legitimate siena_sizes object")
	}
	neff <- length(x$theta)
	wid  <- max(nchar(x$effectNames) + 4)	
	pvalue <- 1-pchisq((divi(x$theta, x$se))^2,1)
	xdf <- data.frame(format(x$effectNames, width = wid), 
				makenice(x$theta),  
				signiftext(divi(x$theta, x$se)), 
				makenice(x$se), 
				rep("   ", neff),
				present.p(pvalue))
	colnames(xdf) <- c("effect", "   par.", " ", "   s.e.", " ", "  p ")	
	print(xdf, right=FALSE)				
	cat("\nSignif. codes:  '.'  p < 0.1; '*' p < 0.05; '**' p < 0.01; '***' p < 0.001. \n")
	lapply(1:length(unique(x$depName)), function(k){pris0(x[[k]])})
	invisible(x)
}


##@summary.siena_sizes Methods
summary.siena_sizes <- function(object, ...)
{
	if (!inherits(object, "siena_sizes"))
	{
		stop("not a legitimate Siena relative importance of effects object")
	}
	class(object) <- c("summary.siena_sizes", class(object))
	object
}

##@print.summary.siena_sizes Methods
print.summary.siena_sizes <- function(x, ...)
{
	if (!inherits(x, "summary.siena_sizes"))
	{
		stop("not a legitimate summary of a Siena relative importance of effects object")
	}
prsri0 <- function(xx){
# This internal function handles results for each dependent variable
		cat(paste("\n  Expected relative importance of effects for dependent variable '",
				xx$dependentVariable,"' at observation moments:\n\n\n", sep=""))
		periods <- length(xx$expectedRI)
		effs <- length(xx$effectNames)
		colNames = paste("period ", 1:periods, sep="")
		line1 <- format("", width =63)
		line2 <- paste(format(1:effs,width=3), '. ',
				format(xx$effectNames, width = 56),sep="")
		line3 <- line2
		line5 <- line2
		for (w in 1:length(colNames))
		{
			line1 <- paste(line1, format(colNames[w], width=8),"  ", sep = "")
			line2 <- paste(line2, format(round(xx$expectedRI[[w]], 4),
				width=8, nsmall=4),"  ",sep="")
			line3 <- paste(line3, format(round(xx$expectedI[[w]], 4),
				width=8, nsmall=4),"  ",sep="")
		}
		line2 <- paste(line2, rep('\n',effs), sep="")
		line3 <- paste(line3, rep('\n',effs), sep="")
		cat(as.matrix(line1),'\n \n', sep='')
		cat(as.matrix(line2),'\n', sep='')
		cat("\n  Expected importance of effects for this dependent variable:\n\n")
		cat(as.matrix(line3),'\n\n', sep='')
	}
	print.siena_sizes(x)
	lapply(1:length(unique(x$depName)), function(k){prsri0(x[[k]])})
	invisible(x)
}


##@plot.siena_sizes Methods
plot.siena_sizes <- function(x, depvar = NULL, actors = NULL, col = NULL, addPieChart = FALSE,
	radius = 1, width = NULL, height = NULL, legend = TRUE,
	legendColumns = NULL, legendHeight = NULL, cex.legend = NULL,
	cex.names = NULL, ...)
{
	if (!inherits(x, "siena_sizes"))
	{
		stop("not a legitimate Siena relative importance of effects object")
	}
	if (is.null(depvar))
	{
		depvar <- x$depName[1]
	}
	xx <- x[[depvar]]
	periods <- length(xx$expectedRI)
	if (is.null(actors))
	{
		nactors <- dim(xx$RIActors[[1]])[2]
		actors <- (1:nactors)
	}
	else
	{
		if ((!inherits(actors,"integer")) ||
			(min(actors) < 1) || (max(actors) > dim(xx$RIActors[[1]])[2]))
		{
			stop(paste("parameter <actors> must be a set of integers from 1 to",
					dim(xx$RIActors[[1]])[2]))
		}
		nactors <- length(actors)
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
			legendColumns <-floor((nactors+2)/11)
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
			legendHeight <-
				max(0.8,ceiling(length(xx$effectNames)/legendColumns)*0.2)
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
		height <- 1
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
		if(addPieChart)
		{
			width = (nactors/3+4)
		}else{
			width = (nactors/3+3)
		}
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
		cex.legend <- 1.3
	}

	if(!is.null(cex.names))
	{
		if(is.numeric(cex.names))
		{
			cex.names <- cex.names
		}else{
			cex.names <- NULL
warning("cex.names has to be of type 'numeric' \n used default settings")
		}
	}
	if(is.null(cex.names))
	{
		cex.names <- 1
	}

	if(!is.null(radius))
	{
		if(is.numeric(radius))
		{
			rad <- radius
		}else{
			rad <- NULL
		warning("radius has to be of type 'numeric' \n used default settings")
		}
	}
	if(is.null(radius))
	{
		rad <- 1
	}

	if(!is.null(col))
	{
		cl <- col
	}else{
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
		while(length(cl)<length(xx$effectNames)){
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
	bordergrey <-"gray25"

	if(addPieChart)
	{
		if(legend)
		{
			layoutMatrix <- matrix(c(1:(2*periods+1),(2*periods+1)), byrow= TRUE,
				ncol=2, nrow=(periods+1))
			layout(layoutMatrix,widths= c((nactors/6)+10,3.5+2.5*(rad^2)),
				heights=c(rep(height,periods),legendHeight))
		}else{
			layoutMatrix <- matrix((1:(2*periods)), byrow= TRUE,
				ncol=2, nrow=periods)
			layout(layoutMatrix,widths = c((nactors/6)+10,7+2.5*(rad^2)),
				heights=rep(height,periods))
		}
	}else{
		if(legend)
		{
			layoutMatrix <- matrix((1:(periods+1)), byrow= TRUE,
				ncol=1, nrow=(periods+1))
			layout(layoutMatrix)
		}else{
			layoutMatrix <- matrix((1:periods), byrow= TRUE, ncol=1, nrow=periods)
			layout(layoutMatrix, heights=2*rep(height,periods))
			# no widths, because these are only relative numbers,
			# so requiring constant widths is redundant
		}
		par( oma = c( 1, 1, 2, 1 ), xpd=T , cex = 0.75, no.readonly = TRUE )
	}
	par(mar = c(3,3,1,1))
	for(w in 1:periods)
	{
		barplot(cbind(xx$RIActors[[w]][,actors], xx$expectedRI[[w]]),
			space=c(rep(0.1,nactors),1.5),width=c(rep(1,nactors),1),
			beside =FALSE, yaxt = "n", xlab="Actor", cex.names = cex.names,
			ylab=paste("period ", w, sep=""),border=bordergrey,
			col = cl, names.arg=c(actors,"exp. rel. imp."))
		axis(2, at=c(0,0.25,0.5,0.75,1),labels=c("0","","0.5","","1"))
		axis(4, at=c(0,0.25,0.5,0.75,1),labels=c("0","","0.5","","1"))
		if(addPieChart)
		{
			pie(xx$expectedRI[[w]], col = cl, labels=NA, border = bordergrey,
				radius = rad)
			mtext("exp. rel. imp.",side = 1, line = 1, cex=cex.names*0.75)
		}
	}
	if(legend)
	{
		plot(c(0,1), c(0,1), col=rgb(0,0,0,0), axes=FALSE, ylab = "", xlab = "")
		legend(0, 1, xx$effectNames, fill=cl, ncol = legendColumns,
			bty = "n", cex=cex.legend)
	}
	invisible(cl)
}

##@sienaSetupForCpp Use as RSiena:::sienaSetupForCpp
sienaSetupForCpp <- function(data, effects,
                             includeBehavior = TRUE,
                             includeBipartite = TRUE,
                             returnActorStatistics = FALSE,
                             returnStaticChangeContributions = FALSE,
                             parallelrun = FALSE,
							 initC = FALSE) {
    f <- unpackData(data)
    effects <- effects[effects$include,]
    effects$setting <- rep("", nrow(effects))
    pData <- .Call(C_setupData, PACKAGE=pkgname,
                   list(as.integer(f$observations)),
                   list(f$nodeSets))
    reg.finalizer(pData, clearData, onexit = FALSE)
    .Call(C_OneMode, PACKAGE=pkgname, pData, list(f$nets))
    if (includeBipartite && !is.null(f$bipartites)) {
        .Call(C_Bipartite, PACKAGE=pkgname, pData, list(f$bipartites))
    }
    if (includeBehavior && !is.null(f$behavs)) {
        .Call(C_Behavior, PACKAGE=pkgname, pData, list(f$behavs))
    }
    .Call(C_ConstantCovariates, PACKAGE=pkgname, pData, list(f$cCovars))
    .Call(C_ChangingCovariates, PACKAGE=pkgname, pData, list(f$vCovars))
    .Call(C_DyadicCovariates, PACKAGE=pkgname, pData, list(f$dycCovars))
    .Call(C_ChangingDyadicCovariates, PACKAGE=pkgname, pData, list(f$dyvCovars))
    # if (is.null(f$exog) || !is.list(f$exog)) {
	# 	f$exog <- vector("list", length(f$depvars))
	# }
	# ## split the names of the constraints
	# higher <- attr(f, "allHigher")
	# disjoint <- attr(f, "allDisjoint")
	# atLeastOne <- attr(f, "allAtLeastOne")
	# froms <- sapply(strsplit(names(higher), ","), function(x)x[1])
	# tos <- sapply(strsplit(names(higher), ","), function(x)x[2])
	# .Call(C_Constraints, PACKAGE=pkgname,
	# 	pData, froms[higher], tos[higher],
	# 	froms[disjoint], tos[disjoint],
	# 	froms[atLeastOne], tos[atLeastOne])
	# siena07 only does this with !initC
	storage.mode(effects$parm) <- 'integer'
    storage.mode(effects$group) <- 'integer'
    storage.mode(effects$period) <- 'integer'
    effects$effectPtr <- rep(NA, nrow(effects))
    depvarnames <- names(data$depvars)
    tmpeffects <- split(effects, effects$name)
    myeffectsOrder <- match(depvarnames, names(tmpeffects))
    ans <- .Call(C_effects, PACKAGE=pkgname, pData, tmpeffects)
    pModel <- ans[[1]][[1]]
	reg.finalizer(pModel, clearModel, onexit = FALSE)

    for (i in 1:length(ans[[2]])) {
        effectPtr <- ans[[2]][[i]]
        tmpeffects[[i]]$effectPtr <- effectPtr
    }
    myeffects <- tmpeffects
    for (i in 1:length(myeffectsOrder)) {
        myeffects[[i]] <- tmpeffects[[myeffectsOrder[i]]]
    }
    list(pData = pData, pModel = pModel, myeffects = myeffects)
}