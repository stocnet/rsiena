#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: sienaModelCreate.r
# *
# * Description: This module contains the function for creating model objects.
# *
# *****************************************************************************/

##@sienaModelCreate DataCreate
sienaModelCreate <- function(fn,
	projname="Siena", MaxDegree=NULL, Offset=NULL,
	useStdInits=FALSE,
	n3=1000, nsub=4, n2start = NULL, dolby=TRUE,
	maxlike=FALSE, gmm=FALSE, diagonalize=0.2*!maxlike,
	condvarno=0, condname='',
	firstg=0.2, reduceg=0.5, cond=NA, findiff=FALSE,  seed=NULL,
	pridg=0.05, prcdg=0.05, prper=0.2, pripr=0.3, prdpr=0.3,
	prirms=0.05, prdrms=0.05, maximumPermutationLength=40,
	minimumPermutationLength=2, initialPermutationLength=20,
	modelType=NULL, behModelType=NULL, mult=5, simOnly=FALSE, localML=FALSE,
	truncation=5, doubleAveraging=0, standardizeVar=(diagonalize<1),
	lessMem=FALSE)
{
	model <- NULL
	checking <- any(grepl("_R_CHECK", names(Sys.getenv()))) 
	if (is.null(projname) | checking)
	{
		model$projname <- tempfile("Siena")
		if (checking)
		{		
	cat('If you use this algorithm object, siena07 will create/use an output file', 
				paste('Siena','.txt',sep=''),'.\n')
		
		}
		else
		{
	cat('If you use this algorithm object, siena07 will create/use an output file',
				paste(model$projname,'.txt',sep=''),'.\n')
			cat('This is a temporary file for this R session.\n')
		}
	}
	else
	{
		if (is.character(projname))
		{
			model$projname <- projname
			cat('If you use this algorithm object, siena07 will create/use an output file', 
				paste(model$projname,'.txt',sep=''),'.\n')
		}
		else
		{
			stop('projname should be a character string.')
		}
	}
	model$useStdInits <- useStdInits
	model$checktime <- TRUE
	model$n3 <- n3
	model$firstg <- firstg
	model$reduceg <- reduceg
	model$maxrat <- 1.0
	model$normSetRates <- FALSE
	model$maxlike <- maxlike
	model$gmm <- gmm
	model$simOnly <- simOnly
	model$localML <- localML
	model$FRANname <- deparse(substitute(fn))
	if (maxlike)
	{
		if (missing(fn))
		{
			model$FRANname <- "maxlikec"
		}
		if (is.na(cond))
		{
			cond <- FALSE
		}
		if (cond)
		{
			stop("Conditional estimation is not possible with",
				"maximum likelihood estimation")
		}
		if (findiff)
		{
			stop("Finite differences estimation of derivatives",
				"is not possible with maximum likelihood estimation")
		}
	}
	else
	{
		if (missing(fn))
		{
			model$FRANname <- "simstats0c"
		}
	}
	model$cconditional <- cond
	if (!is.na(cond) && cond && condvarno == 0 && condname == "")
	{
		model$condvarno <-  1
		model$condname <- ""
	}
	else
	{
		model$condvarno <-  condvarno
		model$condname <- condname
	}
	model$FinDiff.method <-  findiff
	model$nsub <- nsub
	model$n2start <- n2start
	model$dolby <- (dolby && (!maxlike)&& (!gmm))
	if (diagonalize < 0) {diagonalize <- 0}
	if (diagonalize > 1) {diagonalize <- 1}
	model$diagg <- (diagonalize >= 0.9999)
	model$diagonalize <- diagonalize

	if (!is.null(modelType))
	{
		if (any(!(modelType %in% 1:6)))
		{
			warning('modelType can only have values from 1 to 6; other values changed to 1\n')
			model$modelType[!(modelType %in% 1:6)] <- 1
		}
		if (any(is.null(names(modelType))))
		{
			stop('modelType must have names of dependent networks')
		}
		model$modelType <- modelType
	}
	if (!is.null(behModelType))
	{
		if (any(!(behModelType %in% c(1,2))))
		{
			warning('behModelType can only have values 1 or 2; other values changed to 1\n')
			model$behModelType[!(behModelType %in% c(1,2))] <- 1
		}
		if (any(is.null(names(behModelType))))
		{
			stop('behModelType must have names of dependent behavioral variables')
		}
		model$behModelType <- behModelType
	}
	model$MaxDegree <- MaxDegree
	if (!is.null(Offset))
	{
		if (any(is.null(names(Offset))))
		{
			stop('Offset must have names of dependent behavioral variables')
		}
		model$UniversalOffset <- Offset
	}
	model$randomSeed <- seed
	model$pridg <- pridg
	model$prcdg <- prcdg
	model$prper <- prper
	model$pripr <- pripr
	model$prdpr <- prdpr
	model$prirms <- prirms
	model$prdrms <- prdrms
	model$maximumPermutationLength <- maximumPermutationLength
	model$minimumPermutationLength <- minimumPermutationLength
	model$initialPermutationLength <- initialPermutationLength
	model$mult <- mult
	model$truncation <- truncation
	model$doubleAveraging <- doubleAveraging
	model$sf2.byIteration <- !lessMem
	model$standardizeWithTruncation <- standardizeVar
	model$standardizeVar <- standardizeVar
	# The difference between these two is a hidden, non-documented option,
	# perhaps for being tried out
	# by later modification of the sienaAlgorithm object.
	model$noAggregation <- FALSE
	# This also is a hidden, non-documented option, perhaps for being tried out.
	#  \item{noAggregation}{Logical:
	#   do not replace current parameter value after subphase 1
	#   by the mean over subphase 1, if some quasi-autocorrelation
	#   then is larger than .5. May be helpful if initial value was very far away.
	# The two options model$noAggregation and model$standardizeWithTruncation
	# are used only in phase2.r.
	class(model) <- "sienaAlgorithm"
	model
}

model.create <- sienaModelCreate

##@sienaAlgorithmCreate DataCreate
sienaAlgorithmCreate <- sienaModelCreate

##@ModelTypeStrings DataCreate
ModelTypeStrings <- function(i){
	ifelse(((i %in% 1:6) && (!is.null(i))),
		switch(i,
			"Standard actor-oriented model",
			"Forcing model",
			"Initiative model",
			"Pairwise forcing model",
			"Pairwise mutual model",
			"Pairwise joint model"), "")
}

##@BehaviorModelTypeStrings DataCreate
BehaviorModelTypeStrings <- function(i){
	ifelse(((i %in% 1:2) && (!is.null(i))),
		switch(i,
			"Standard behavior actor-oriented model ('restrict')",
			"Boundary-absorbing behavior model"),
		"")
}
