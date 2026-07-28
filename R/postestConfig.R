#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: postestConfig.R
# *
# * Description: Configuration-object constructors for postestimation
# * (marginalEffects.sienaFit / predict.sienaFit) uncertainty and control
# * settings. These bundle the flat arguments currently passed individually
# * into two validated, printable objects: `sienaPostestUncertainty` and
# * `sienaPostestControl`.
# *
# * This file is purely additive: it does not change marginalEffects.sienaFit
# * or any other existing function. Wiring these objects into the
# * postestimation entry points is done elsewhere.
# *****************************************************************************/


# ---------------------------------------------------------------------------
# Internal validation helpers (not exported).
# ---------------------------------------------------------------------------

# A single non-NA logical (TRUE/FALSE flag).
.checkFlag <- function(x, name)
{
	if (!is.logical(x) || length(x) != 1L || is.na(x))
	{
		stop("'", name, "' must be a single non-NA logical (TRUE or FALSE).")
	}
	invisible(x)
}

# A single finite numeric, optionally required to be > 0.
# Returns the (possibly integer-coerced) value.
.checkSingleNumeric <- function(x, name, positive = TRUE, asInteger = FALSE)
{
	if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
	{
		stop("'", name, "' must be a single finite numeric value.")
	}
	if (positive && x <= 0)
	{
		stop("'", name, "' must be > 0.")
	}
	if (asInteger) as.integer(x) else x
}

# As .checkSingleNumeric, but NULL passes through unchanged.
.checkSingleNumericOrNull <- function(x, name, positive = TRUE, asInteger = FALSE)
{
	if (is.null(x)) return(NULL)
	.checkSingleNumeric(x, name, positive = positive, asInteger = asInteger)
}

# A single non-NA character string.
.checkSingleChar <- function(x, name)
{
	if (!is.character(x) || length(x) != 1L || is.na(x))
	{
		stop("'", name, "' must be a single non-NA character string.")
	}
	invisible(x)
}

# As .checkSingleChar, but NULL passes through unchanged.
.checkSingleCharOrNull <- function(x, name)
{
	if (is.null(x)) return(NULL)
	.checkSingleChar(x, name)
	x
}


# ---------------------------------------------------------------------------
# set_postest_uncertainty_saom
# ---------------------------------------------------------------------------

##@set_postest_uncertainty_saom PostestConfig
set_postest_uncertainty_saom <- function(
	enabled    = TRUE,
	mode       = c("bootstrap", "delta", "deltaFull"),
	nsim       = 1000,
	sd         = TRUE,
	ci         = TRUE,
	ciInterval = c(0.025, 0.975),
	simMean    = FALSE,
	simMedian  = FALSE
)
{
	mode <- match.arg(mode)

	.checkFlag(enabled, "enabled")
	.checkFlag(sd, "sd")
	.checkFlag(ci, "ci")
	.checkFlag(simMean, "simMean")
	.checkFlag(simMedian, "simMedian")

	nsim <- .checkSingleNumeric(nsim, "nsim", positive = TRUE, asInteger = TRUE)

	if (!is.numeric(ciInterval) || length(ciInterval) != 2L ||
		any(!is.finite(ciInterval)))
	{
		stop("'ciInterval' must be a numeric vector of length 2 with ",
			 "finite values.")
	}
	if (any(ciInterval <= 0) || any(ciInterval >= 1))
	{
		stop("'ciInterval' must have both values strictly between 0 and 1.")
	}
	if (!(ciInterval[1] < ciInterval[2]))
	{
		stop("'ciInterval' must be strictly increasing: ",
			 "ciInterval[1] < ciInterval[2].")
	}

	# simMean / simMedian are summaries of the simulation draw distribution;
	# delta and deltaFull are analytic and produce no draws to summarize.
	if ((simMean || simMedian) && mode != "bootstrap")
	{
		offending <- c("simMean", "simMedian")[c(simMean, simMedian)]
		verb <- if (length(offending) > 1L) "require" else "requires"
		stop(paste(paste0(offending, " = TRUE"), collapse = " and "), " ",
			 verb, " mode = 'bootstrap' (delta modes draw no simulations). ",
			 "Current mode = '", mode, "'.")
	}

	obj <- list(
		enabled    = enabled,
		mode       = mode,
		nsim       = nsim,
		sd         = sd,
		ci         = ci,
		ciInterval = ciInterval,
		simMean    = simMean,
		simMedian  = simMedian
	)
	class(obj) <- "sienaPostestUncertainty"
	attr(obj, "version") <- utils::packageDescription("RSiena", fields = "Version")
	obj
}


##@print.sienaPostestUncertainty Methods
print.sienaPostestUncertainty <- function(x, ...)
{
	cat("RSiena postestimation uncertainty\n")

	if (!isTRUE(x$enabled))
	{
		cat("  status  : disabled (point estimates only, no uncertainty computed)\n")
		return(invisible(x))
	}

	modeDesc <- switch(x$mode,
		bootstrap = paste0("bootstrap  (nsim = ", x$nsim, " draws)"),
		delta     = "delta  (analytic; no simulation draws)",
		deltaFull = "deltaFull  (analytic, full covariance; no simulation draws)"
	)
	cat("  mode    : ", modeDesc, "\n", sep = "")

	summaryParts <- character(0)
	if (isTRUE(x$sd)) summaryParts <- c(summaryParts, "SE")
	if (isTRUE(x$ci))
	{
		pct <- diff(x$ciInterval) * 100
		pctStr <- if (isTRUE(all.equal(pct, round(pct))))
			sprintf("%d%%", as.integer(round(pct))) else sprintf("%.3g%%", pct)
		summaryParts <- c(summaryParts,
			sprintf("%s CI [%s, %s]", pctStr, x$ciInterval[1], x$ciInterval[2]))
	}
	if (identical(x$mode, "bootstrap"))
	{
		if (isTRUE(x$simMean))   summaryParts <- c(summaryParts, "mean of draws")
		if (isTRUE(x$simMedian)) summaryParts <- c(summaryParts, "median of draws")
	}
	if (length(summaryParts) == 0L) summaryParts <- "none requested"
	cat("  summary : ", paste(summaryParts, collapse = ", "), "\n", sep = "")

	invisible(x)
}


# ---------------------------------------------------------------------------
# set_postest_algo_saom
# ---------------------------------------------------------------------------

##@set_postest_algo_saom PostestConfig
set_postest_algo_saom <- function(
	dynamic                = FALSE,
	algorithm              = NULL,
	n3                     = 200,
	n3PointEst             = NULL,
	n3BatchSize            = 100L,
	chainStoreMode         = c("auto", "disk", "memory"),
	useChangeContributions = FALSE,
	chainStorePath         = NULL,
	useCluster             = FALSE,
	nbrNodes               = 1,
	clusterType            = c("PSOCK", "FORK"),
	cl                     = NULL,
	batchDir               = "temp",
	prefix                 = "simBatch_b",
	combineBatch           = TRUE,
	batchSize              = NULL,
	keepBatch              = FALSE,
	verbose                = TRUE,
	details                = FALSE,
	memoryScale            = NULL,
	batchUnitBudget        = 2.5e+08,
	dynamicMinistepFactor  = 10,
	saveDir                = NULL,
	gcEachBatch            = FALSE,
	gcEachSim              = FALSE
)
{
	chainStoreMode <- match.arg(chainStoreMode)
	clusterType    <- match.arg(clusterType)

	# ---- logical flags: single non-NA logical ----
	# Checked up front (before these values are used in conditionals below)
	# so that e.g. dynamic = NA or useCluster = NA fail with a clear message
	# naming the offending argument, rather than a cryptic base-R error from
	# `if (NA && ...)`.
	# `verbose` is deliberately excluded here: it is used as a level
	# elsewhere (verbose >= 1, verbose >= 2), so any single non-NA
	# numeric-or-logical is accepted and stored unchanged.
	.checkFlag(dynamic, "dynamic")
	.checkFlag(useChangeContributions, "useChangeContributions")
	.checkFlag(useCluster, "useCluster")
	.checkFlag(combineBatch, "combineBatch")
	.checkFlag(keepBatch, "keepBatch")
	.checkFlag(details, "details")
	.checkFlag(gcEachBatch, "gcEachBatch")
	.checkFlag(gcEachSim, "gcEachSim")

	if (!((is.logical(verbose) || is.numeric(verbose)) &&
		  length(verbose) == 1L && !is.na(verbose)))
	{
		stop("'verbose' must be a single non-NA logical or numeric value.")
	}

	# ---- cluster normalisation (lifted from siena07.r) ----
	# Passing `cl` without `useCluster = TRUE` currently runs single-core
	# silently; make the intent explicit instead.
	if (!useCluster && length(cl) > 0)
	{
		useCluster <- TRUE
		nbrNodes   <- length(cl)
	}

	if (dynamic && is.null(algorithm))
	{
		stop("'algorithm' must be provided when dynamic = TRUE.")
	}

	n3          <- .checkSingleNumeric(n3, "n3", positive = TRUE, asInteger = TRUE)
	nbrNodes    <- .checkSingleNumeric(nbrNodes, "nbrNodes", positive = TRUE,
										asInteger = TRUE)
	n3BatchSize <- .checkSingleNumeric(n3BatchSize, "n3BatchSize", positive = TRUE,
										asInteger = TRUE)
	n3PointEst  <- .checkSingleNumericOrNull(n3PointEst, "n3PointEst",
											  positive = TRUE, asInteger = TRUE)
	batchSize   <- .checkSingleNumericOrNull(batchSize, "batchSize",
											  positive = TRUE, asInteger = TRUE)

	batchUnitBudget <- .checkSingleNumeric(batchUnitBudget, "batchUnitBudget",
											positive = TRUE, asInteger = FALSE)
	dynamicMinistepFactor <- .checkSingleNumeric(dynamicMinistepFactor,
												  "dynamicMinistepFactor",
												  positive = TRUE, asInteger = FALSE)

	memoryScale <- .checkSingleNumericOrNull(memoryScale, "memoryScale",
											  positive = TRUE, asInteger = FALSE)

	.checkSingleChar(batchDir, "batchDir")
	.checkSingleChar(prefix, "prefix")
	chainStorePath <- .checkSingleCharOrNull(chainStorePath, "chainStorePath")
	saveDir        <- .checkSingleCharOrNull(saveDir, "saveDir")

	obj <- list(
		dynamic                = dynamic,
		algorithm              = algorithm,
		n3                     = n3,
		n3PointEst             = n3PointEst,
		n3BatchSize            = n3BatchSize,
		chainStoreMode         = chainStoreMode,
		useChangeContributions = useChangeContributions,
		chainStorePath         = chainStorePath,
		useCluster             = useCluster,
		nbrNodes               = nbrNodes,
		clusterType            = clusterType,
		cl                     = cl,
		batchDir               = batchDir,
		prefix                 = prefix,
		combineBatch           = combineBatch,
		batchSize              = batchSize,
		keepBatch              = keepBatch,
		verbose                = verbose,
		details                = details,
		memoryScale            = memoryScale,
		batchUnitBudget        = batchUnitBudget,
		dynamicMinistepFactor  = dynamicMinistepFactor,
		saveDir                = saveDir,
		gcEachBatch            = gcEachBatch,
		gcEachSim              = gcEachSim
	)
	class(obj) <- "sienaPostestControl"
	attr(obj, "version") <- utils::packageDescription("RSiena", fields = "Version")
	obj
}


##@print.sienaPostestControl Methods
print.sienaPostestControl <- function(x, ...)
{
	kind <- if (isTRUE(x$dynamic)) "dynamic" else "static"
	cat("RSiena postestimation control (", kind, ")\n", sep = "")

	simLine <- paste0("n3 = ", x$n3)
	if (!is.null(x$n3PointEst))
	{
		simLine <- paste0(simLine, ", n3PointEst = ", x$n3PointEst)
	}
	simLine <- paste0(simLine, ", batches of ", x$n3BatchSize)
	cat("  simulation : ", simLine, "\n", sep = "")

	if (isTRUE(x$useCluster))
	{
		cat("  cluster    : ", x$nbrNodes,
			if (x$nbrNodes == 1L) " node (" else " nodes (",
			x$clusterType, ")\n", sep = "")
	}
	else
	{
		cat("  cluster    : single-core\n")
	}

	chainsLine <- x$chainStoreMode
	if (isTRUE(x$useChangeContributions))
	{
		chainsLine <- paste0(chainsLine, "  (reusing pre-computed change contributions)")
	}
	cat("  chains     : ", chainsLine, "\n", sep = "")
	if (!is.null(x$chainStorePath))
	{
		cat("  chain path : ", x$chainStorePath, "\n", sep = "")
	}

	batchSizeStr <- if (is.null(x$batchSize)) "auto" else as.character(x$batchSize)
	dynBit <- if (isTRUE(x$dynamic))
		paste0(", ministepFactor = ", x$dynamicMinistepFactor) else ""
	cat("  batching   : size = ", batchSizeStr,
		", budget = ", format(x$batchUnitBudget, scientific = TRUE),
		" units", dynBit, "\n", sep = "")
	if (!is.null(x$memoryScale))
	{
		cat("  memScale   : ", x$memoryScale, "\n", sep = "")
	}
	cat("  batch files: dir = \"", x$batchDir, "\", prefix = \"", x$prefix,
		"\", combine = ", x$combineBatch, ", keep = ", x$keepBatch, "\n", sep = "")

	vb <- x$verbose
	loggingBits <- if (isTRUE(vb)) {
		"verbose"
	} else if (is.numeric(vb) && vb >= 1) {
		paste0("verbose (level ", vb, ")")
	} else {
		"quiet"
	}
	if (isTRUE(x$details)) loggingBits <- paste0(loggingBits, ", details")
	cat("  logging    : ", loggingBits, "\n", sep = "")

	gcBits <- character(0)
	if (isTRUE(x$gcEachBatch)) gcBits <- c(gcBits, "after each batch")
	if (isTRUE(x$gcEachSim))   gcBits <- c(gcBits, "after each simulation")
	if (length(gcBits) > 0L)
	{
		cat("  gc         : ", paste(gcBits, collapse = ", "), "\n", sep = "")
	}

	if (!is.null(x$saveDir))
	{
		cat("  saveDir    : ", x$saveDir, "\n", sep = "")
	}

	invisible(x)
}
