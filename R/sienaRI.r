#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena/
# *
# * File: sienaRI.r
# *
# * Description: Used to determine, print, and plots relative importances of effects
# * for potential decisions of actors at observation moments.
# *****************************************************************************/




# ---- relativeImportance generic + sienaFit method --------------------------

##@relativeImportance Generic for computing relative importance of effects
relativeImportance <- function(object, ...) {
    UseMethod("relativeImportance", object)
}

##@relativeImportance.sienaFit  Relative importance for sienaFit objects
# Computes the relative importance of effects at observation moments
# (Indlekofer & Brandes 2013) for a fitted SAOM.
#
# Returns a sienaRI object containing per-ego-period RI values
# and within-ego standard deviations of change statistics.
# Optionally includes uncertainty quantification
# via theta draws from the asymptotic distribution.
#
# When dynamic = TRUE, Markov chains are simulated between
# observation moments using siena07, and RI is computed at
# each simulated ministep rather than at observation moments only.
relativeImportance.sienaFit <- function(object, data,
    depvar = NULL,
    effects = NULL,
    dynamic = FALSE,
    distFun = c("L1D", "KLD"),
    algorithm = NULL,
    n3 = 500,
    useChangeContributions = FALSE,
    uncertainty = FALSE,
    level = "ego",
    sum_fun = mean,
    na.rm = TRUE,
    nsim = 1000,
    uncertaintyMode = c("batch", "stream"),
    uncertaintySd = TRUE,
    uncertaintyCi = TRUE,
    uncertaintyProbs = c(0.025, 0.5, 0.975),
    useCluster = FALSE,
    nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    batchDir = "temp",
    batchSize = NULL,
    keepBatch = FALSE,
    verbose = TRUE,
    memoryScale = NULL,
    batchUnitBudget = 2.5e8,
    dynamicMinistepFactor = 10,
    getChangeStats = FALSE,
    ...) {

    uncertaintyMode <- match.arg(uncertaintyMode)
    distFun <- match.arg(distFun)
    if (is.null(effects)) effects <- object$requestedEffects
    if (!inherits(effects, "sienaEffects"))
        stop("effects is not a legitimate Siena effects object")
		# endowment and creation should be possible now
    if (sum(effects$include == TRUE &
            (effects$type == "endow" | effects$type == "creation")) > 0)
        stop("relativeImportance does not yet work for models containing ",
             "endowment or creation effects")
    if (dynamic && is.null(algorithm))
        stop("'algorithm' must be provided when dynamic = TRUE")
    # Check interaction effects have their base effects present
    interactionShortNames <- c("unspInt", "behUnspInt", "contUnspInt")
    interactionEffs <- effects[effects$include &
        effects$shortName %in% interactionShortNames, ]
    if (nrow(interactionEffs) > 0) {
        referenced <- unique(c(interactionEffs$effect1,
            interactionEffs$effect2, interactionEffs$effect3))
        referenced <- referenced[referenced > 0]
        missingBase <- setdiff(referenced, effects$effectNumber)
        if (length(missingBase) > 0)
            stop("Requires full effects object when interaction ",
                 "effects are present without base effects.")
    }
    if (is.null(depvar)) {
        depvars <- names(data[["depvars"]])
    } else {
        depvars <- depvar
    }

    thetaHat <- nameThetaFromEffects(object$theta, effects)

    result_list <- list()
    for (dv in depvars) {

        if (dynamic) {
            # Dynamic path: simulate chains, compute RI at each ministep
            riData <- predictRelativeImportanceDynamic(
                ans = object, data = data, theta = thetaHat,
                algorithm = algorithm, effects = effects,
                depvar = dv, n3 = n3,
                useChangeContributions = useChangeContributions,
                distFun = distFun,
                ...
            )
            effNames <- grep("^RI_", names(riData), value = TRUE)
            effNames <- sub("^RI_", "", effNames)
        } else {
            # Static path: extract contributions at observation moments
            staticContribs <- getStaticChangeContributions(
                ans = object, data = data, effects = effects,
                depvar = dv, returnWide = TRUE
            )
            riData <- predictRelativeImportanceStatic(staticContribs, thetaHat,
                distFun = distFun)
            effNames <- staticContribs$effectNames
        }

        # Build the new-style sienaRI object
        reqEffects <- effects[effects$include == TRUE, ]
        noRate     <- reqEffects$type != "rate"
        reqEffects <- reqEffects[noRate, ]
        dvEffects  <- reqEffects$name == dv
        oldEffNames <- paste(reqEffects$type[dvEffects], " ",
                              reqEffects$effectName[dvEffects], sep = "")

        riObj <- list(
            dependentVariable = dv,
            data              = riData,
            effectNames       = oldEffNames,
            shortEffectNames  = effNames,
            dynamic           = dynamic
        )

        # Attach change statistics if requested (static only)
        if (getChangeStats && !dynamic) {
            riObj$changeStatistics <- staticContribs
        }

        # Uncertainty via sienaPostestimate
        if (uncertainty) {
            if (is.null(batchSize)) {
                batchSize <- planBatch(
                    data = data, depvar = dv, nsim = nsim,
                    nbrNodes = nbrNodes, useCluster = useCluster,
                    dynamic = dynamic, n3 = n3,
                    unitBudget = batchUnitBudget,
                    dynamicMinistepFactor = dynamicMinistepFactor,
                    memoryScale = memoryScale
                )
            }
            riCols <- grep("^RI_", names(riData), value = TRUE)

            if (dynamic) {
                predictFun  <- predictRelativeImportanceDynamic
                predictArgs <- list(
                    ans       = object,
                    data      = data,
                    algorithm = algorithm,
                    effects   = effects,
                    depvar    = dv,
                    n3        = n3,
                    useChangeContributions = useChangeContributions
                )
            } else {
                predictFun  <- predictRelativeImportanceStatic
                predictArgs <- list(staticContributions = staticContribs,
                                    distFun = distFun)
            }

            riObj$uncertainty <- sienaPostestimate(
                predictFun       = predictFun,
                predictArgs      = predictArgs,
                outcomeName      = riCols[1],
                level            = level,
                condition        = NULL,
                sum_fun          = sum_fun,
                na.rm            = na.rm,
                thetaHat         = thetaHat,
                covTheta         = object$covtheta,
                uncertainty      = TRUE,
                uncertaintyMode  = uncertaintyMode,
                nsim             = nsim,
                uncertaintySd    = uncertaintySd,
                uncertaintyCi    = uncertaintyCi,
                uncertaintyProbs = uncertaintyProbs,
                useCluster       = useCluster,
                nbrNodes         = nbrNodes,
                clusterType      = clusterType,
                batchDir         = batchDir,
                batchSize        = batchSize,
                keepBatch        = keepBatch,
                verbose          = verbose
            )
        }

        class(riObj) <- "sienaRI"
        attr(riObj, "version") <- packageDescription(pkgname, fields = "Version")
        result_list[[dv]] <- riObj
    }

    if (length(result_list) == 1L) {
        return(result_list[[1L]])
    } else {
        message("more than one dependent variable\n",
                "return value is therefore not of class 'sienaRI'\n",
                "but a list of objects of class 'sienaRI'.")
        return(result_list)
    }
}


##@relativeImportance.sienaEffects  Point-estimate RI from effects + theta
# Computes RI without a sienaFit object, using a user-supplied
# parameter vector theta and a sienaEffects object.
# Does not support uncertainty quantification or dynamic RI.
relativeImportance.sienaEffects <- function(object, data,
        theta, depvar = NULL, getChangeStats = FALSE,
        distFun = c("L1D", "KLD"), ...) {
    distFun  <- match.arg(distFun)
    effects  <- object
    if (is.null(depvar)) {
        depvars <- names(data[["depvars"]])
    } else {
        depvars <- depvar
    }
    thetaHat <- nameThetaFromEffects(theta, effects)

    result_list <- lapply(depvars, function(dv) {
        contribs <- getStaticChangeContributions(
            ans = NULL, data = data, effects = effects,
            depvar = dv, returnWide = TRUE)
        riData <- computeRelativeImportance(contribs, thetaHat,
                                            distFun = distFun)

        reqEff <- effects[effects$include & effects$type != "rate" &
                          effects$name == dv, ]
        effNames <- sub("^RI_", "", grep("^RI_", names(riData), value = TRUE))
        oldEffNames <- paste(reqEff$type, reqEff$effectName, sep = " ")

        riObj <- list(
            dependentVariable = dv,
            data              = riData,
            effectNames       = oldEffNames,
            shortEffectNames  = effNames,
            dynamic           = FALSE
        )
        if (getChangeStats) riObj$changeStatistics <- contribs
        class(riObj) <- "sienaRI"
        attr(riObj, "version") <- packageDescription(pkgname, fields = "Version")
        riObj
    })

    if (length(result_list) == 1L) {
        result_list[[1L]]
    } else {
        message("more than one dependent variable\n",
                "return value is therefore not of class 'sienaRI'\n",
                "but a list of objects of class 'sienaRI'.")
        result_list
    }
}


# ---- Relative importance (vectorised on wide struct) -----------------------


# Static relative importance: reuses pre-built static contributions.
predictRelativeImportanceStatic <- function(staticContributions, theta,
                                             structurals = NULL,
                                             distFun = c("L1D", "KLD")) {
    computeRelativeImportance(staticContributions, theta,
                               structurals = structurals,
                               distFun = distFun)
}

# Dynamic relative importance: simulates chains then computes RI.
predictRelativeImportanceDynamic <- function(ans, data, theta, algorithm,
                                              effects, depvar, n3 = NULL,
                                              useChangeContributions = FALSE,
                                              batch = TRUE, silent = TRUE,
                                              distFun = c("L1D", "KLD")) {
    dynamicContributions <- getDynamicChangeContributions(ans = ans,
        theta = theta, data = data, algorithm = algorithm,
        effects = effects, depvar = depvar, n3 = n3,
        useChangeContributions = useChangeContributions,
        returnWide = TRUE, batch = batch, silent = silent
    )
    computeRelativeImportance(dynamicContributions, theta, distFun = distFun)
}

## Swappable grouped divergence between full-theta and LOO choice distributions.
## ref     : nRows vector (full-theta choice probabilities)
## loo     : nRows x nEff matrix (LOO probs, col k = effect k zeroed)
## group_id: integer vector of contiguous group labels (one block per ego x period)
## Returns : nGroups x nEff divergence matrix; degenerate groups (<=1 valid
##           choice) get NA rows.  Change the body here to plug in a new measure.
computeGroupedDistance <- function(ref, loo, group_id,
                                    distFun = c("L1D", "KLD")) {
    distFun <- match.arg(distFun)
    gid <- as.integer(group_id)
    if (distFun == "L1D") l1d_grouped(ref, loo, gid)
    else                   kld_grouped(ref, loo, gid)
}

# Core vectorised relative importance computation on any wide struct.
computeRelativeImportance <- function(contributions, theta,
                                       structurals = NULL,
                                       distFun = c("L1D", "KLD")) {
    distFun    <- match.arg(distFun)
    theta_use  <- theta[contributions$effectNames]
    contribMat <- contributions$contribMat
    group_id   <- contributions$group_id
    effNames   <- contributions$effectNames
    nEff       <- length(effNames)

    # --- full-theta probabilities ---
    utility  <- calculateUtility(contribMat, theta_use)
    fullProb <- calculateChangeProb(utility, group_id)
    if (!is.null(structurals)) fullProb[structurals] <- NA

    # --- leave-one-out probabilities (one C++ call replacing the R loop) ---
    # Exploits util_k = util_full - contribMat[,k]*theta[k], so O(n) per effect
    # instead of a full O(n*nEff) matrix multiply for each k.
    looProb <- loo_change_probs(contribMat, theta_use, as.integer(group_id))
    if (!is.null(structurals)) looProb[structurals, ] <- NA_real_

    # --- per-group divergence -> normalised RI ---
    uGroups   <- unique(group_id)
    nGroups   <- length(uGroups)
    firstIdx  <- match(uGroups, group_id)
    coordCols <- groupColsList(contributions, firstIdx)
    coordCols[["choice"]] <- NULL

    divMat    <- computeGroupedDistance(fullProb, looProb, group_id, distFun)
    absSumVec <- rowSums(divMat, na.rm = TRUE)
    degRows   <- is.na(divMat[, 1L])   # groups with <= 1 valid choice
    riMat     <- divMat / pmax(absSumVec, .Machine$double.eps)
    riMat[degRows, ]   <- NA_real_
    absSumVec[degRows] <- NA_real_

    # --- within-ego sigma of change statistics ---
    sigmaMat <- matrix(NA_real_, nrow = nGroups, ncol = nEff)
    for (i in seq_len(nGroups)) {
        rows <- which(group_id == uGroups[i])
        for (k in seq_len(nEff)) {
            sigmaMat[i, k] <- sd(contribMat[rows, k], na.rm = TRUE)
        }
    }

    out <- data.frame(coordCols, stringsAsFactors = FALSE)
    for (k in seq_len(nEff)) {
        out[[paste0("RI_", effNames[k])]] <- riMat[, k]
    }
    out[["absSum"]]    <- absSumVec
    for (k in seq_len(nEff)) {
        out[[paste0("sigma_", effNames[k])]] <- sigmaMat[, k]
    }
    # Sigma-standardised RI: importance expressed in units of within-ego
    # change-statistic variation (NA when sigma is zero or NA).
    for (k in seq_len(nEff)) {
        sigK <- sigmaMat[, k]
        out[[paste0("stdRI_", effNames[k])]] <-
            ifelse(!is.na(sigK) & sigK > 0, riMat[, k] / sigK, NA_real_)
    }
    if (requireNamespace("data.table", quietly = TRUE)) data.table::setDT(out)
    out
}

# expectedRelativeImportance and calculateDistributions have been removed.
# The new path uses computeRelativeImportance on the wide struct produced by
# getStaticChangeContributions(returnWide = TRUE).  calculateDistributions
# duplicated calculateUtility + calculateChangeProb (the C++ path).

## KLD is kept as a utility for alternative divergence measures (not used in the main RI path).
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

##@print.sienaRI Methods
print.sienaRI <- function(x, printSigma = FALSE, ...){
	object <- x
	if (!inherits(object, "sienaRI"))
	{
		if (inherits(object[[1]], "sienaRI"))
		{
			cat("This object is a list, the components of which\n")
			cat("are Siena relative importance of effects objects.\n")
			cat("Apply the print function to the separate components.\n")
		}
		stop("not a legitimate Siena relative importance of effects object")
	}

	# ---- New format (data.frame-based) ----
	if (!is.null(object$data)) {
		rms <- function(xx) sqrt((1/ncol(xx)) * rowSums(xx^2, na.rm = TRUE))
		riDF <- if (requireNamespace("data.table", quietly = TRUE) &&
				data.table::is.data.table(object$data))
			as.data.frame(object$data) else object$data
		effNames <- object$shortEffectNames
		riCols   <- paste0("RI_", effNames)
		sigCols  <- paste0("sigma_", effNames)
		periods  <- sort(unique(riDF$period))
		nPeriods <- length(periods)
		nEff     <- length(object$effectNames)

		if (isTRUE(object$dynamic)) {
			cat(paste("\n  Expected relative importance of effects for dependent variable '",
				object$dependentVariable, "'\n  (dynamic: averaged across simulated ministeps):\n\n\n", sep = ""))
		} else {
			cat(paste("\n  Expected relative importance of effects for dependent variable '",
				object$dependentVariable, "' at observation moments:\n\n\n", sep = ""))
		}

		colNames <- paste("period ", periods, sep = "")
		line1 <- format("", width = 63)
		line2 <- paste(format(1:nEff, width = 3), '. ',
			format(object$effectNames, width = 56), sep = "")
		line3 <- line2
		line5 <- line2

		for (w in seq_len(nPeriods)) {
			pRows <- riDF[riDF$period == periods[w], ]
			eriW  <- colMeans(pRows[, riCols, drop = FALSE], na.rm = TRUE)
			# un-normalised importance: RI * absSum
			absSumW <- pRows[["absSum"]]
			eiW <- numeric(nEff)
			for (k in seq_len(nEff)) {
				eiW[k] <- mean(pRows[[riCols[k]]] * absSumW, na.rm = TRUE)
			}

			line1 <- paste(line1, format(colNames[w], width = 8), "  ", sep = "")
			line2 <- paste(line2, format(round(eriW, 4),
				width = 8, nsmall = 4), "  ", sep = "")
			line3 <- paste(line3, format(round(eiW, 4),
				width = 8, nsmall = 4), "  ", sep = "")

			if (printSigma) {
				# sigma per period: RMS across actors
				sMat <- as.matrix(pRows[, sigCols, drop = FALSE])
				sigW <- sqrt(colMeans(sMat^2, na.rm = TRUE))
				line5 <- paste(line5,
					format(round(sigW, 4), width = 8, nsmall = 4), "  ", sep = "")
			}
		}
		line2 <- paste(line2, rep('\n', nEff), sep = "")
		line3 <- paste(line3, rep('\n', nEff), sep = "")
		cat(as.matrix(line1), '\n \n', sep = '')
		cat(as.matrix(line2), '\n', sep = '')
		cat("\n  Expected importance of effects for this dependent variable:\n\n")
		cat(as.matrix(line3), '\n\n', sep = '')
		if (printSigma) {
			cat("\n sigma (average within-ego standard deviation of change statistics):\n\n")
			line5 <- paste(line5, rep('\n', nEff), sep = "")
			cat('\n', as.matrix(line5), '\n', sep = '')
		}
		if (!is.null(object$uncertainty)) {
			cat("\n  Uncertainty estimates available in object$uncertainty\n")
		}
		invisible(object)
		return(invisible(object))
	}

	# ---- Old format (list-based, from sienaRI() direct call) ----
	cat(paste("\n  Expected relative importance of effects for dependent variable '",
			object$dependentVariable,"' at observation moments:\n\n\n", sep=""))
	periods <- length(object$expectedRI)
	effects <- length(object$effectNames)
	colNames = paste("period ", 1:periods, sep="")
	line1 <- format("", width =63)
	line2 <- paste(format(1:effects,width=3), '. ',
		format(object$effectNames, width = 56),sep="")
	line3 <- line2
	line5 <- line2
	for (w in 1:length(colNames))
	{
		line1 <- paste(line1, format(colNames[w], width=8),"  ", sep = "")
		line2 <- paste(line2, format(round(object$expectedRI[[w]], 4),
				width=8, nsmall=4),"  ",sep="")
		line3 <- paste(line3, format(round(object$expectedI[[w]], 4),
				width=8, nsmall=4),"  ",sep="")
		if (printSigma)
		{
			line5 <- paste(line5,
				format(round(object$sigmas[,w], 4), width=8, nsmall=4),"  ",sep="")
		}
	}
	line2 <- paste(line2, rep('\n',effects), sep="")
	line3 <- paste(line3, rep('\n',effects), sep="")
	cat(as.matrix(line1),'\n \n', sep='')
	cat(as.matrix(line2),'\n', sep='')
	cat("\n  Expected importance of effects for this dependent variable:\n\n")
	cat(as.matrix(line3),'\n\n', sep='')
	if (printSigma)
	{
		cat("\n sigma (average within-ego standard deviation of change statistics):\n\n")
		line5 <- paste(line5, rep('\n',effects), sep="")
		cat('\n',as.matrix(line5),'\n', sep='')
	}
	invisible(object)
}

##@summary.sienaRI Methods
summary.sienaRI <- function(object, ...)
{
	if (!inherits(object, "sienaRI"))
	{
		stop("not a legitimate Siena relative importance of effects object")
	}
	class(object) <- c("summary.sienaRI", class(object))
	object
}
##@print.summary.sienaRI Methods
print.summary.sienaRI <- function(x, ...)
{
	object <- x
	if (!inherits(object, "summary.sienaRI"))
	{
		stop("not a legitimate summary of a Siena relative importance of effects object")
	}
	print.sienaRI(object)
	invisible(object)
}


##@plot.sienaRI Methods
plot.sienaRI <- function(x, actors = NULL, col = NULL, addPieChart = FALSE,
	radius = 1, width = NULL, height = NULL, legend = TRUE,
	legendColumns = NULL, legendHeight = NULL, cex.legend = NULL,
	cex.names = NULL, intervalsPerPeriod = 20, ...)
{
	object <- x
	if (!inherits(object, "sienaRI"))
	{
		stop("not a legitimate Siena relative importance of effects object")
	}
	# Convert new format to old-style fields for plotting
	if (!is.null(object$data)) {
		riDF <- if (requireNamespace("data.table", quietly = TRUE) &&
				data.table::is.data.table(object$data))
			as.data.frame(object$data) else object$data
		effNames <- object$shortEffectNames
		riCols   <- paste0("RI_", effNames)
		periods_vec <- sort(unique(riDF$period))
		nPeriods <- length(periods_vec)
		nEff     <- length(effNames)

		# ---- Dynamic: time-series line plot ----
		if (isTRUE(object$dynamic)) {
			# Generate colour palette
			if (!is.null(col)) {
				cl <- col
			} else {
				alph <- 175
				cl <- c(
					rgb(127,201,127,alph, maxColorValue=255),
					rgb(190,174,212,alph, maxColorValue=255),
					rgb(253,192,134,alph, maxColorValue=255),
					rgb(255,255,153,alph, maxColorValue=255),
					rgb( 56,108,176,alph, maxColorValue=255),
					rgb(184,184,184,alph, maxColorValue=255),
					rgb( 56, 56, 56,alph, maxColorValue=255),
					rgb(120,120,120,alph, maxColorValue=255),
					rgb(240,  2,127,alph, maxColorValue=255),
					rgb(191, 91, 23,alph, maxColorValue=255))
				while (length(cl) < nEff) {
					alph <- (alph + 75) %% 255
					cl <- c(cl,
						rgb(127,201,127,alph, maxColorValue=255),
						rgb(190,174,212,alph, maxColorValue=255),
						rgb(253,192,134,alph, maxColorValue=255),
						rgb(255,255,153,alph, maxColorValue=255),
						rgb( 56,108,176,alph, maxColorValue=255),
						rgb(184,184,184,alph, maxColorValue=255),
						rgb( 56, 56, 56,alph, maxColorValue=255),
						rgb(120,120,120,alph, maxColorValue=255),
						rgb(240,  2,127,alph, maxColorValue=255),
						rgb(191, 91, 23,alph, maxColorValue=255))
				}
			}
			if (is.null(cex.legend)) cex.legend <- 1.3

			# Layout: one panel per period + optional legend
			if (legend) {
				layout(matrix(seq_len(nPeriods + 1L), ncol = 1))
			} else {
				layout(matrix(seq_len(nPeriods), ncol = 1))
			}
			par(mar = c(4, 4, 2, 1), oma = c(1, 1, 2, 1))

			for (w in seq_len(nPeriods)) {
				pRows <- riDF[riDF$period == periods_vec[w], ]
				ms <- pRows$ministep
				if (length(unique(ms)) < 2L) {
					# Only one ministep — degenerate case, skip panel
					plot.new()
					title(main = paste("Period", periods_vec[w]))
					next
				}
				nInt <- min(intervalsPerPeriod, length(unique(ms)))
				breaks <- seq(min(ms), max(ms), length.out = nInt + 1L)
				bins   <- findInterval(ms, breaks, all.inside = TRUE)
				xvals  <- (breaks[-length(breaks)] + breaks[-1]) / 2

				riMeans <- matrix(NA_real_, nrow = nInt, ncol = nEff)
				for (b in seq_len(nInt)) {
					bRows <- pRows[bins == b, , drop = FALSE]
					if (nrow(bRows) > 0L) {
						for (k in seq_len(nEff)) {
							riMeans[b, k] <- mean(bRows[[riCols[k]]],
								na.rm = TRUE)
						}
					}
				}

				ymax <- max(riMeans, na.rm = TRUE)
				if (!is.finite(ymax) || ymax == 0) ymax <- 1
				matplot(xvals, riMeans, type = "l", lty = 1,
					col = cl[seq_len(nEff)], lwd = 2,
					xlab = "Ministep", ylab = "Relative importance",
					main = paste("Period", periods_vec[w]),
					ylim = c(0, ymax * 1.1))
			}

			if (legend) {
				plot.new()
				legend("center", legend = object$effectNames,
					col = cl[seq_len(nEff)], lty = 1, lwd = 2,
					ncol = if (!is.null(legendColumns)) legendColumns else 2L,
					bty = "n", cex = cex.legend)
			}
			return(invisible(cl))
		}

		# ---- Static: convert to old-style for barplots ----
		object$expectedRI <- list()
		object$RIActors   <- list()
		for (w in seq_len(nPeriods)) {
			pRows <- riDF[riDF$period == periods_vec[w], ]
			riMat <- t(as.matrix(pRows[, riCols, drop = FALSE]))
			rownames(riMat) <- effNames
			object$RIActors[[w]] <- riMat
			object$expectedRI[[w]] <- rowMeans(riMat, na.rm = TRUE)
		}
	}
	periods <- length(object$expectedRI)
	if (is.null(actors))
	{
		nactors <- dim(object$RIActors[[1]])[2]
		actors <- (1:nactors)
	}
	else
	{
		if ((!inherits(actors,"integer")) ||
			(min(actors) < 1) || (max(actors) > dim(object$RIActors[[1]])[2]))
		{
			stop(paste("parameter <actors> must be a set of integers from 1 to",
					dim(object$RIActors[[1]])[2]))
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
				max(0.8,ceiling(length(object$effectNames)/legendColumns)*0.2)
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
		while(length(cl)<length(object$effectNames)){
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
		barplot(cbind(object$RIActors[[w]][,actors], object$expectedRI[[w]]),
			space=c(rep(0.1,nactors),1.5),width=c(rep(1,nactors),1),
			beside =FALSE, yaxt = "n", xlab="Actor", cex.names = cex.names,
			ylab=paste("period ", w, sep=""),border=bordergrey,
			col = cl, names.arg=c(actors,"exp. rel. imp."))
		axis(2, at=c(0,0.25,0.5,0.75,1),labels=c("0","","0.5","","1"))
		axis(4, at=c(0,0.25,0.5,0.75,1),labels=c("0","","0.5","","1"))
		if(addPieChart)
		{
			pie(object$expectedRI[[w]], col = cl, labels=NA, border = bordergrey,
				radius = rad)
			mtext("exp. rel. imp.",side = 1, line = 1, cex=cex.names*0.75)
		}
	}
	if(legend)
	{
		plot(c(0,1), c(0,1), col=rgb(0,0,0,0), axes=FALSE, ylab = "", xlab = "")
		legend(0, 1, object$effectNames, fill=cl, ncol = legendColumns,
			bty = "n", cex=cex.legend)
	}
	invisible(cl)
}

# Deprecated -- use interpret_size() instead.
#
# When ans is provided, delegates directly to
# interpret_size(ans, data, ...).
# When only theta + effects are provided (ans = NULL),
# delegates to interpret_size(effects, data, theta = theta, ...).
sienaRI <- function(data, ans = NULL, theta = NULL, effects = NULL,
    getChangeStats = FALSE)
{
    .Deprecated("interpret_size",
        msg = "sienaRI() is deprecated. Use interpret_size() or relativeImportance() instead.")
    if (!inherits(data, "sienadata") && !inherits(data, "siena"))
        stop("not a legitimate Siena data specification")
    if (!is.null(ans)) {
        if (!inherits(ans, "sienaFit"))
            stop("ans is not a legitimate Siena fit object")
        suppressWarnings(
            interpret_size(ans, data = data, effects = effects,
                getChangeStats = getChangeStats))
    } else {
        if (is.null(effects))
            stop("Either 'ans' (sienaFit) or 'effects' (sienaEffects) must be provided.")
        if (!inherits(effects, "sienaEffects"))
            stop("effects is not a legitimate Siena effects object")
        if (is.null(theta))
            stop("'theta' must be provided when ans = NULL.")
        suppressWarnings(
            interpret_size(effects, data = data, theta = theta,
                getChangeStats = getChangeStats))
    }
}
