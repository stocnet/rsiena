##@marginalEffects Generic
marginalEffects <- function(object, ...) UseMethod("marginalEffects", object)

##@marginalEffects.sienaFit Method
marginalEffects.sienaFit <- function(
    object,
    data,
    effects = NULL,
    depvar = NULL,
    effectName1,
    diff1 = NULL,
    contrast1 = NULL,
    interaction1 = FALSE,
    intEffectNames1 = NULL,
    modEffectNames1 = NULL,
    second = FALSE,
    effectName2 = NULL,
    diff2 = NULL,
    contrast2 = NULL,
    interaction2 = FALSE,
    intEffectNames2 = NULL,
    modEffectNames2 = NULL,
    effectList = NULL,
    level = "period",
    condition = NULL,
    type = c("changeProb", "tieProb"),
    mainEffect = "riskDifference", # allow utility later as well
    perturbType1 = NULL,
    perturbType2 = NULL,
    massContrasts = NULL,
    egoNormalize = TRUE,
    returnDecisionDetails = FALSE,
    sum_fun = mean,
    na.rm = TRUE,
    dynamic = FALSE,
    algorithm = NULL,
    n3 = 200,
    n3PointEst = NULL,
    useChangeContributions = FALSE,
    uncertainty = TRUE,
    nsim = 1000,
    uncertaintySd = TRUE,
    uncertaintyCi = TRUE,
    uncertaintyMean = FALSE,
    uncertaintyMedian = FALSE,
    uncertaintyProbs = c(0.025, 0.5, 0.975),
    uncertaintyMcse = FALSE,
    uncertaintymcseBatches = NULL,
    # Multicore
    useCluster = FALSE,
    nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    cl = NULL,
    batchDir = "temp",
    prefix = "simBatch_b",
    combineBatch = TRUE,
    batchSize = NULL,
    keepBatch = FALSE,
    verbose = TRUE,
    details = FALSE,
    memoryScale = NULL,
    batchUnitBudget = 2.5e8,
    dynamicMinistepFactor = 10,
    saveDir = NULL,
    gcEachBatch = FALSE,
    gcEachSim = FALSE,
    ...
) {
    if (inherits(data, "sienaGroup"))
      stop("marginalEffects does not support multi-group data (sienaGroup).")
    type <- match.arg(type)
    if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]

    # ---- behaviour DV guard ----
    dvType <- attr(data[["depvars"]][[depvar]], "type")
    if (!is.null(dvType) && dvType == "behavior")
      stop("marginalEffects is currently only implemented for network ",
           "dependent variables. Behaviour variable '", depvar,
           "' was selected.")

    if (dynamic && is.null(algorithm))
      stop("'algorithm' must be provided when dynamic = TRUE")

    # ================================================================
    # Unified code path: if scalar effectName1 args are given (single
    # effect), wrap them into a one-element effectList so everything
    # flows through the same code.
    # ================================================================
    .single_effect <- FALSE
    if (is.null(effectList) && !missing(effectName1)) {
        effectList <- list(single = list(
            effectName1     = effectName1,
            diff1           = diff1,
            contrast1       = contrast1,
            interaction1    = interaction1,
            intEffectNames1 = intEffectNames1,
            modEffectNames1 = modEffectNames1,
            second          = second,
            effectName2     = effectName2,
            diff2           = diff2,
            contrast2       = contrast2,
            interaction2    = interaction2,
            intEffectNames2 = intEffectNames2,
            modEffectNames2 = modEffectNames2,
            perturbType1    = perturbType1,
            perturbType2    = perturbType2,
            massContrasts   = massContrasts,
            returnDecisionDetails = returnDecisionDetails
        ))
        .single_effect <- TRUE
    }

    if (is.null(effectList))
        stop("Either 'effectName1' or 'effectList' must be provided.")

    # ================================================================
    # Main computation path (handles both single and multi-effect).
    #
    # Static:  getStaticChangeContributions() once (theta-independent),
    #          then N effects x nsim draws evaluated cheaply.
    # Dynamic: getDynamicChangeContributions() once per theta draw
    #          (n3 forward sims shared across all N effects).
    #
    # level, condition, type, mainEffect are fixed for all effects.
    # Per-effect entries only specify effect-identifying params
    # (effectName1, diff1, contrast1, interaction1, etc.).
    # ================================================================
    if (!is.list(effectList) || length(effectList) == 0L)
        stop("'effectList' must be a non-empty named list.")
    if (is.null(names(effectList)) || any(nchar(names(effectList)) == 0L))
        stop("All elements of 'effectList' must be named.")

        if (is.null(effects)) effects <- object$requestedEffects
        thetaHat    <- coef(object)
        covTheta    <- vcov(object)

        # ---- saveDir: check for completed effects --
        if (!is.null(saveDir)) {
            if (!dir.exists(saveDir))
                dir.create(saveDir, recursive = TRUE)
            done <- vapply(names(effectList), function(nm)
                file.exists(file.path(saveDir, paste0(nm, ".rds"))),
                logical(1L))
            if (all(done)) {
                if (verbose) message("All effects found in saveDir, loading.")
                results <- lapply(names(effectList), function(nm)
                    readRDS(file.path(saveDir, paste0(nm, ".rds"))))
                names(results) <- names(effectList)
                return(if (.single_effect) results[["single"]] else results)
            }
        }

        # ---- Build shared contribution function ----
        if (!dynamic) {
            # Strip stored chains before passing object to getStaticChangeContributions
            # (chains are theta-independent and not needed for static contributions).
            # object itself is rm()'d before forking below.
            object$changeContributions <- NULL
            staticContributions <- getStaticChangeContributions(
                ans     = object,
                data    = data,
                effects = effects,
                depvar  = depvar,
                returnWide = TRUE
            )

            staticContributions$changeStats <- contribToChangeStats(
                staticContributions$contribMat,
                staticContributions$effectNames,
                theta = NULL  # theta populated per-draw in predictFirstDiff
            )
            csNames <- staticContributions$changeStats$csNames
            knownEffectNames <- csNames
            # Re-key effectInteractionTypes to changeStats names.
            eitOld <- staticContributions$effectInteractionTypes
            # Map composite names to changeStats: strip the type suffix.
            eitBases <- sub("_[^_]+$", "", names(eitOld))
            eitCanon <- setNames(
              eitOld[match(csNames, eitBases)],
              csNames
            )
            getContribFun <- function(theta) staticContributions
        } else {
            if ("basicRate" %in% names(effects) && !any(effects$basicRate)) {
                stop(
                    "The 'effects' object contains no rate effects. ",
                    "This is required when dynamic = TRUE because siena07 ",
                    "is rerun internally.\n",
                    "Please pass the full effects object, e.g. effects = mymodel."
                )
            }
            knownEffectNames <- getEffectNamesNoRate(effects, depvar)
            # dynArgs: used by uncertainty workers (small n3 for speed).
            dynArgs <- list(
                ans      = object,
                data     = data,
                algorithm = algorithm,
                effects  = effects,
                depvar   = depvar,
                n3       = n3,
                useChangeContributions = useChangeContributions,
                returnWide = TRUE
            )
            # dynArgsHat: used for the point estimate only.
            # n3PointEst=NULL means use algorithm$n3 (estimation-time value),
            # giving a high-quality point estimate independent of the (smaller)
            # n3 used for uncertainty workers.
            dynArgsHat <- dynArgs
            dynArgsHat$n3 <- n3PointEst
            getContribFun <- function(theta) {
                do.call(getDynamicChangeContributions,
                        c(list(theta = theta), dynArgs))
            }
            getContribFunHat <- function(theta) {
                do.call(getDynamicChangeContributions,
                        c(list(theta = theta), dynArgsHat))
            }
        }

        # ---- effectInteractionTypes for resolvePerturbType ----
        if (!dynamic) {
            eiTypes <- eitCanon
        } else if ("interactionType" %in% names(effects)) {
            noRate  <- effects$type != "rate"
            inc     <- effects[
                noRate & effects$include & effects$name == depvar, ]
            eiTypes <- setNames(inc$interactionType, knownEffectNames)
        } else {
            eiTypes <- NULL
        }

        # ---- Uncertainty summary function ----
        uncertainty_summary_fun <- makeUncertaintySummarizer(
            return_sd     = uncertaintySd,
            return_ci     = uncertaintyCi,
            return_mean   = uncertaintyMean,
            return_median = uncertaintyMedian,
            probs         = uncertaintyProbs,
            return_mcse   = uncertaintyMcse,
            mcseBatches   = uncertaintymcseBatches
        )

        # ---- Build validated per-effect spec list ----
        builtSpecList <- vector("list", length(effectList))
        names(builtSpecList) <- names(effectList)

        for (nm in names(effectList)) {
            spec <- effectList[[nm]]
            eff_second <- isTRUE(spec$second)

            if (is.null(spec$effectName1))
                stop("Effect '", nm, "' must specify 'effectName1'.")

            # feells unnecessary to have a helper function here?
            
            resolved <- resolveAMEEffectNames(
                knownEffectNames,
                spec$effectName1,
                spec$intEffectNames1,
                spec$modEffectNames1,
                spec$effectName2,
                spec$intEffectNames2,
                spec$modEffectNames2,
                eff_second
            )

            pt1 <- resolvePerturbType(resolved$effectName1, eiTypes,
                                      spec$perturbType1)
            pt2 <- if (eff_second)
                       resolvePerturbType(resolved$effectName2, eiTypes,
                                          spec$perturbType2)
                   else
                       "alter"

            eff_massC <- if (!is.null(spec$massContrasts)) {
                spec$massContrasts
            } else {
                (pt1 == "ego") || (eff_second && pt2 == "ego")
            }

            eff_diffName <- if (!eff_second) {
                ifelse(mainEffect == "riskDifference", "firstDiff",
                       "firstRiskRatio")
            } else {
                ifelse(mainEffect == "riskDifference", "secondDiff",
                       "secondRiskRatio")
            }

            # Auto-center contrasts for centered covariates
            eff_contrast1 <- spec$contrast1
            eff_contrast2 <- spec$contrast2
            if (!is.null(eff_contrast1)) {
                cm1 <- getCovCenteringMean(resolved$effectName1,
                                           effects, data, depvar)
                if (cm1 != 0) eff_contrast1 <- eff_contrast1 - cm1
            }
            if (eff_second && !is.null(eff_contrast2)) {
                cm2 <- getCovCenteringMean(resolved$effectName2,
                                           effects, data, depvar)
                if (cm2 != 0) eff_contrast2 <- eff_contrast2 - cm2
            }

            eff_retDet    <- isTRUE(spec$returnDecisionDetails)
            needs_attrib  <- !is.null(condition) || eff_retDet

            if (!eff_second) {
                diffFun  <- predictFirstDiff
                diffArgs <- list(
                    type           = type,
                    effectName     = resolved$effectName1,
                    diff           = spec$diff1,
                    contrast       = eff_contrast1,
                    interaction    = isTRUE(spec$interaction1),
                    intEffectNames = resolved$intEffectNames1,
                    modEffectNames = resolved$modEffectNames1,
                    details        = FALSE,
                    calcRiskRatio  = FALSE,
                    mainEffect     = mainEffect,
                    perturbType    = pt1,
                    massContrasts  = eff_massC,
                    attachContribs = needs_attrib
                )
            } else {
                diffFun  <- predictSecondDiff
                diffArgs <- list(
                    type            = type,
                    effectName1     = resolved$effectName1,
                    diff1           = spec$diff1,
                    contrast1       = eff_contrast1,
                    interaction1    = isTRUE(spec$interaction1),
                    intEffectNames1 = resolved$intEffectNames1,
                    modEffectNames1 = resolved$modEffectNames1,
                    effectName2     = resolved$effectName2,
                    diff2           = spec$diff2,
                    contrast2       = eff_contrast2,
                    interaction2    = isTRUE(spec$interaction2),
                    intEffectNames2 = resolved$intEffectNames2,
                    modEffectNames2 = resolved$modEffectNames2,
                    details         = FALSE,
                    calcRiskRatio   = FALSE,
                    mainEffect      = mainEffect,
                    perturbType1    = pt1,
                    perturbType2    = pt2,
                    massContrasts   = eff_massC,
                    attachContribs  = needs_attrib
                )
            }

            builtSpecList[[nm]] <- list(
                diffFun               = diffFun,
                diffArgs              = diffArgs,
                outcomeName           = eff_diffName,
                level                 = level,
                condition             = condition,
                massContrasts         = eff_massC,
                returnDecisionDetails = eff_retDet,
                metadata = list(
                    method      = "marginalEffects",
                    type        = type,
                    effectName1 = resolved$effectName1,
                    contrast1   = eff_contrast1,
                    effectName2 = if (eff_second) resolved$effectName2
                                  else NULL,
                    contrast2   = if (eff_second) eff_contrast2
                                  else NULL,
                    second      = eff_second,
                    level       = level,
                    condition   = condition,
                    depvar      = depvar,
                    dynamic     = dynamic,
                    nsim        = nsim,
                    n3          = if (dynamic) n3 else NULL
                )
            )
        }

        # ---- Default batchSize (batch path uses nbrNodes=1) ----
        if (is.null(batchSize)) {
            batchSize <- planBatch(
                data     = data,
                depvar   = depvar,
                nsim     = nsim,
                dynamic  = dynamic,
                n3       = n3,
                nbrNodes = nbrNodes,
                useCluster = useCluster,
                unitBudget = batchUnitBudget,
                dynamicMinistepFactor = dynamicMinistepFactor,
                memoryScale = memoryScale,
                useChangeContributions = useChangeContributions
            )
        }

        # ---- Point estimates (shared hat contributions) ----
        # Dynamic: use getContribFunHat which uses n3PointEst (default NULL =
        # algorithm$n3 from estimation).  Static: getContribFun is fine (theta-
        # independent, n3 irrelevant).
        ccHat <- if (dynamic) getContribFunHat(thetaHat) else getContribFun(thetaHat)
        baselineHat <- predictProbability(ccHat, thetaHat, type,
                                          returnComponents = TRUE)
        theta_hat_use <- baselineHat$theta_use

        decision_details <- vector("list", length(builtSpecList))
        names(decision_details) <- names(builtSpecList)

        expects <- lapply(names(builtSpecList), function(nm) {
            spec    <- builtSpecList[[nm]]
            pe_args <- spec$diffArgs
            if (isTRUE(spec$returnDecisionDetails))
                pe_args$attachContribs <- TRUE
            unit_pred <- do.call(
                spec$diffFun,
                c(list(changeContributions = ccHat,
                       theta_use = theta_hat_use,
                       baseline = baselineHat),
                  pe_args)
            )
            if (isTRUE(spec$returnDecisionDetails))
                decision_details[[nm]] <<- unit_pred
            main_result <- agg(spec$outcomeName, unit_pred,
                level = spec$level, condition = spec$condition,
                sum_fun = sum_fun, na.rm = na.rm,
                egoNormalize = egoNormalize)
            # Also aggregate mass contrast point estimates alongside main outcome
            massCols_pe <- intersect(c("massCreation", "massDissolution"),
                                     names(unit_pred))
            for (mc in massCols_pe) {
                mc_agg <- agg(mc, unit_pred,
                              level = spec$level, condition = NULL,
                              sum_fun = sum_fun, na.rm = na.rm,
                              egoNormalize = egoNormalize)
                mc_by <- intersect(
                    getGroupVars(level = spec$level, condition = NULL),
                    intersect(names(main_result), names(mc_agg))
                )
                if (length(mc_by) > 0L) {
                    main_result <- merge(main_result, mc_agg, by = mc_by,
                                         all.x = TRUE, sort = FALSE)
                } else {
                    main_result[[mc]] <- mc_agg[[mc]]
                }
            }
            main_result
        })
        names(expects) <- names(builtSpecList)

        if (!uncertainty) {
            results <- mapply(function(nm, spec) {
                res <- stampPostestimate(expects[[nm]], spec$metadata)
                if (isTRUE(spec$returnDecisionDetails) &&
                    !is.null(decision_details[[nm]]))
                    attr(res, "decisionDetails") <- decision_details[[nm]]
                res
            }, names(builtSpecList), builtSpecList, SIMPLIFY = FALSE)
            if (!is.null(saveDir)) {
                for (nm in names(results))
                    saveRDS(results[[nm]],
                            file.path(saveDir, paste0(nm, ".rds")))
            }
            return(if (.single_effect) results[["single"]] else results)
        }

        # ---- Uncertainty via shared simulation loop ----
        rm(ccHat, baselineHat)
        if (dynamic && exists("dynArgs", inherits = FALSE)) {
            # Workers draw their own theta_sim — chains must be generated at
            # that theta, not reused from hat-theta estimation.  Set
            # useChangeContributions=FALSE explicitly (avoids warning fallback)
            # and drop ans entirely: theta/effects/data/algorithm are all
            # explicit in dynArgs, so ans is dead weight on the FALSE path.
            dynArgs$useChangeContributions <- FALSE
            dynArgs$ans <- NULL
        }
        rm(object)
        gc(verbose = FALSE)

        # When saveDir is set, keep batch files until all effects are
        # saved so the simulation phase can resume on crash.
        keepBatch_internal <- if (!is.null(saveDir)) TRUE else keepBatch

        raw_sims_list <- drawSimBatch(
            getContribFun  = getContribFun,
            effectSpecList = builtSpecList,
            thetaHat    = thetaHat,
            covTheta    = covTheta,
            type        = type,
            nsim        = nsim,
            useCluster  = useCluster,
            nbrNodes    = nbrNodes,
            clusterType = clusterType,
            cl     = cl,
            batchSize   = batchSize,
            batchDir    = batchDir,
            prefix      = prefix,
            keepBatch   = keepBatch_internal,
            verbose     = verbose,
            gcEachBatch = gcEachBatch,
            gcEachSim   = gcEachSim,
            egoNormalize = egoNormalize
        )

        # ---- Aggregate uncertainty per effect ----
        results <- vector("list", length(builtSpecList))
        names(results) <- names(builtSpecList)

        for (nm in names(builtSpecList)) {
            # Skip effects already saved from a previous (crashed) run
            if (!is.null(saveDir) &&
                file.exists(file.path(saveDir, paste0(nm, ".rds")))) {
                results[[nm]] <- readRDS(
                    file.path(saveDir, paste0(nm, ".rds")))
                if (verbose) message("  Loading saved result: ", nm)
                next
            }

            spec     <- builtSpecList[[nm]]
            raw_j    <- raw_sims_list[[nm]]
            expect_j <- expects[[nm]]

            uncert_j <- agg(
                spec$outcomeName, raw_j,
                level     = spec$level,
                condition = spec$condition,
                sum_fun   = uncertainty_summary_fun,
                na.rm     = na.rm,
                egoNormalize = egoNormalize
            )
            result_j <- mergeEstimates(expect_j, uncert_j,
                                       level     = spec$level,
                                       condition = spec$condition)

            massCols <- intersect(
                c("massCreation", "massDissolution"), names(raw_j))
            for (mc in massCols) {
                mc_uncert <- agg(mc, raw_j,
                                 level     = spec$level,
                                 condition = NULL,
                                 sum_fun   = uncertainty_summary_fun,
                                 na.rm     = na.rm,
                                 egoNormalize = egoNormalize)
                level_by <- intersect(
                    getGroupVars(level = spec$level, condition = NULL),
                    names(mc_uncert)
                )
                uc_cols <- setdiff(names(mc_uncert), level_by)
                for (uc in uc_cols)
                    names(mc_uncert)[names(mc_uncert) == uc] <-
                        paste0(mc, "_", uc)
                if (length(level_by) > 0L) {
                    result_j <- merge(result_j, mc_uncert,
                                      by = level_by,
                                      all.x = TRUE, sort = FALSE)
                } else {
                    for (col in setdiff(names(mc_uncert), level_by))
                        result_j[[col]] <- mc_uncert[[col]]
                }
            }

            results[[nm]] <- stampPostestimate(result_j, spec$metadata)
            if (isTRUE(spec$returnDecisionDetails) &&
                !is.null(decision_details[[nm]]))
                attr(results[[nm]], "decisionDetails") <-
                    decision_details[[nm]]

            # Save immediately for crash resilience
            if (!is.null(saveDir))
                saveRDS(results[[nm]],
                        file.path(saveDir, paste0(nm, ".rds")))
        }

        # Clean up batch files now that all effects are saved
        if (!is.null(saveDir) && !keepBatch) {
            batch_pattern <- sprintf("^%s\\d{3}\\.rds$", prefix)
            bf <- list.files(batchDir, pattern = batch_pattern,
                             full.names = TRUE)
            for (f in bf) file.remove(f)
        }

    # Return: single data frame for scalar call, named list for effectList
    if (.single_effect) results[["single"]] else results
}



# ---- Thin wrappers (static / dynamic) ----------------------------------------

predictFirstDiffStatic <- function(theta, staticContributions,
    type = "changeProb",
    effectName, diff = NULL, contrast = NULL,
    interaction = FALSE, intEffectNames = NULL, modEffectNames = NULL,
    details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference",
    perturbType = "alter", massContrasts = FALSE, attachContribs = TRUE)
{
  # can theta even still be incorrectly ordered here? maybe safer to align by name just in case
  effectNames <- staticContributions$effectNames
  theta_use   <- theta[effectNames]
  predictFirstDiff(changeContributions = staticContributions, theta_use, type,
      effectName, diff, contrast, interaction, intEffectNames,
      modEffectNames, details, calcRiskRatio, mainEffect,
      perturbType = perturbType, massContrasts = massContrasts,
      attachContribs = attachContribs)
}

predictSecondDiffStatic <- function(theta, staticContributions,
    type = "changeProb",
    effectName1, diff1 = NULL, contrast1 = NULL,
    interaction1 = FALSE, intEffectNames1 = NULL, modEffectNames1 = NULL,
    effectName2 = NULL, diff2 = NULL, contrast2 = NULL,
    interaction2 = FALSE, intEffectNames2 = NULL, modEffectNames2 = NULL,
    mainEffect = "riskDifference", details = FALSE,
    perturbType1 = "alter", perturbType2 = "alter",
    massContrasts = FALSE, attachContribs = TRUE)
{
  effectNames <- staticContributions$effectNames
  theta_use   <- theta[effectNames]
  # can theta even still be incorrectly ordered here? maybe safer to align by name just in case
  predictSecondDiff(changeContributions = staticContributions, theta_use, type,
      effectName1, diff1, contrast1, interaction1, intEffectNames1, modEffectNames1,
      effectName2, diff2, contrast2, interaction2, intEffectNames2, modEffectNames2,
      details, FALSE, mainEffect,
      perturbType1 = perturbType1, perturbType2 = perturbType2,
      massContrasts = massContrasts, attachContribs = attachContribs)
}

predictFirstDiffDynamic <- function(ans, data, theta, effects, algorithm,
    type = "changeProb", depvar = NULL,
    effectName, diff = NULL, contrast = NULL,
    interaction = FALSE, intEffectNames = NULL, modEffectNames = NULL,
    n3 = NULL, useChangeContributions = FALSE,
    details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference",
    perturbType = "alter", massContrasts = FALSE, attachContribs = TRUE)
{
  if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
  dynContrib <- getDynamicChangeContributions(
    ans = ans, theta = theta, data = data, algorithm = algorithm,
    effects = effects, depvar = depvar, n3 = n3,
    useChangeContributions = useChangeContributions, returnWide = TRUE
  )
  dynContrib$changeStats <- contribToChangeStats(dynContrib$contribMat,
                                                dynContrib$effectNames)
  theta_use <- theta[dynContrib$effectNames]
  predictFirstDiff(changeContributions = dynContrib, theta_use, type,
      effectName, diff, contrast, interaction, intEffectNames,
      modEffectNames, details, calcRiskRatio, mainEffect,
      perturbType = perturbType, massContrasts = massContrasts,
      attachContribs = attachContribs)
}

predictSecondDiffDynamic <- function(ans, data, theta, effects, algorithm,
    type = "changeProb", depvar = NULL,
    effectName1, diff1 = NULL, contrast1 = NULL,
    interaction1 = FALSE, intEffectNames1 = NULL, modEffectNames1 = NULL,
    effectName2 = NULL, diff2 = NULL, contrast2 = NULL,
    interaction2 = FALSE, intEffectNames2 = NULL, modEffectNames2 = NULL,
    n3 = NULL, useChangeContributions = FALSE,
    calcRiskRatio = FALSE, mainEffect = "riskDifference", details = FALSE,
    perturbType1 = "alter", perturbType2 = "alter",
    massContrasts = FALSE, attachContribs = TRUE)
{
  if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
  dynContrib <- getDynamicChangeContributions(
    theta = theta, data = data, algorithm = algorithm,
    effects = effects, depvar = depvar, n3 = n3,
    useChangeContributions = useChangeContributions, returnWide = TRUE
  )
  dynContrib$changeStats <- contribToChangeStats(dynContrib$contribMat,
                                                dynContrib$effectNames)
  theta_use <- theta[dynContrib$effectNames]
  predictSecondDiff(changeContributions = dynContrib, theta_use, type,
      effectName1, diff1, contrast1, interaction1, intEffectNames1, modEffectNames1,
      effectName2, diff2, contrast2, interaction2, intEffectNames2, modEffectNames2,
      details, calcRiskRatio, mainEffect,
      perturbType1 = perturbType1, perturbType2 = perturbType2,
      massContrasts = massContrasts, attachContribs = attachContribs)
}

# predictFirstDiff and predictSecondDiff reuse a lot of code, unify?

# Core computation shared by both static and dynamic paths.
# `changeContributions` — unified wide struct (contribMat + coord vectors + effectNames)
# `theta_use`           — named numeric vector aligned to changeContributions$effectNames
predictFirstDiff <- function(changeContributions, theta_use, type,
    effectName, diff, contrast, interaction, intEffectNames,
    modEffectNames, details, calcRiskRatio, mainEffect,
    perturbType = "alter", massContrasts = FALSE, attachContribs = TRUE,
    baseline = NULL)
{
  # Ensure changeStats is cached on the contribution struct.
  if (is.null(changeContributions$changeStats))
    changeContributions$changeStats <- contribToChangeStats(
      changeContributions$contribMat, changeContributions$effectNames)
  cs <- changeContributions$changeStats
  csMat      <- cs$csMat
  csNames    <- cs$csNames
  densityIdx <- cs$densityIdx
  density    <- cs$density

  # Use pre-computed baseline components if supplied, else compute.
  if (!is.null(baseline)) {
    thetaEff   <- baseline$thetaEff
    utility    <- baseline$utility
    changeProb <- baseline$changeProb
    tieProb    <- baseline$tieProb
  } else {
    thetaEff   <- buildThetaEff(theta_use, cs$changeStatsMap)
    if (!is.null(changeContributions$changeUtility) &&
        !all(is.na(changeContributions$changeUtility))) {
      utility    <- changeContributions$changeUtility
      changeProb <- changeContributions$changeProbability
    } else {
      utility    <- calculateUtility(csMat, thetaEff,
                                     changeContributions$permitted, densityIdx)
      changeProb <- calculateChangeProb(utility, changeContributions$group_id)
    }
    tieProb <- if (type == "tieProb")
      calculateTieProb(changeProb, density) else NULL
  }

  fd <- calculateFirstDiff(
    densityValue       = density,
    changeProb         = changeProb,
    changeUtil         = utility,
    effectName         = effectName,
    effectContribution = csMat[, effectName],
    diff               = diff,
    contrast           = contrast,
    interaction        = interaction,
    intEffectNames     = intEffectNames,
    modEffectNames     = modEffectNames,
    modContribution    = if (!is.null(modEffectNames)) csMat[, modEffectNames] else NULL,
    effectNames        = csNames,
    theta              = thetaEff,
    type               = type,
    tieProb            = tieProb,
    details            = details,
    calcRiskRatio      = calcRiskRatio,
    mainEffect         = mainEffect,
    perturbType        = perturbType,
    group_id           = changeContributions$group_id
    # document use of group_id better
  )

  # maybe should be optional to subset to keep?
  keep <- density != 0L
  out  <- groupColsList(changeContributions, keep)
  out[names(fd)] <- lapply(fd, `[`, keep)

  if (details) {
    out[["changeUtil"]] <- utility[keep]
    out[["changeProb"]] <- changeProb[keep]
    if (!is.null(tieProb)) out[["tieProb"]] <- tieProb[keep]
  }

  # Mass contrasts on ego-level (Delta P_i^+ and Delta P_i^-)
  if (massContrasts) {
    diffColName <- intersect(c("firstDiff", "firstRiskRatio"), names(fd))[1L]
    if (!is.na(diffColName) && diffColName == "firstDiff") {
      mc <- computeMassContrasts(
        firstDiff = fd[["firstDiff"]][keep],
        density   = density[keep],
        ego       = changeContributions$ego[keep],
        period    = changeContributions$period[keep],
        group     = if (!is.null(changeContributions$group))
                      changeContributions$group[keep] else rep(1L, sum(keep)),
        type      = type
      )
      out[["massCreation"]]    <- mc[["massCreation"]]
      out[["massDissolution"]] <- mc[["massDissolution"]]
    }
  }

  if (attachContribs) {
    out <- attachContributions(out, csNames,
                               csMat[keep, , drop = FALSE], flip = FALSE)
  }
  attr(out, "row.names") <- .set_row_names(length(out[[1L]]))
  class(out) <- "data.frame"
  out
}


predictSecondDiff <- function(changeContributions, theta_use, type,
    effectName1, diff1, contrast1, interaction1,
    intEffectNames1, modEffectNames1,
    effectName2, diff2, contrast2, interaction2,
    intEffectNames2, modEffectNames2,
    details, calcRiskRatio, mainEffect,
    perturbType1 = "alter", perturbType2 = "alter",
    massContrasts = FALSE, attachContribs = TRUE,
    baseline = NULL)
{
  # Ensure changeStats is cached on the contribution struct.
  if (is.null(changeContributions$changeStats))
    changeContributions$changeStats <- contribToChangeStats(
      changeContributions$contribMat, changeContributions$effectNames)
  cs <- changeContributions$changeStats
  csMat      <- cs$csMat
  csNames    <- cs$csNames
  densityIdx <- cs$densityIdx
  density    <- cs$density

  # Use pre-computed baseline components if supplied, else compute.
  if (!is.null(baseline)) {
    thetaEff   <- baseline$thetaEff
    utility    <- baseline$utility
    changeProb <- baseline$changeProb
    tieProb    <- baseline$tieProb
  } else {
    thetaEff   <- buildThetaEff(theta_use, cs$changeStatsMap)
    if (!is.null(changeContributions$changeUtility) &&
        !all(is.na(changeContributions$changeUtility))) {
      utility    <- changeContributions$changeUtility
      changeProb <- changeContributions$changeProbability
    } else {
      utility    <- calculateUtility(csMat, thetaEff,
                                     changeContributions$permitted, densityIdx)
      changeProb <- calculateChangeProb(utility, changeContributions$group_id)
    }
    tieProb <- if (type == "tieProb")
      calculateTieProb(changeProb, density) else NULL
  }

  sd <- calculateSecondDiff(
    densityValue        = density,
    changeProb          = changeProb,
    changeUtil          = utility,
    effectName1         = effectName1,
    effectContribution1 = csMat[, effectName1],
    diff1               = diff1,
    contrast1           = contrast1,
    interaction1        = interaction1,
    intEffectNames1     = intEffectNames1,
    modEffectNames1     = modEffectNames1,
    modContribution1    = if (!is.null(modEffectNames1)) csMat[, modEffectNames1] else NULL,
    effectName2         = effectName2,
    effectContribution2 = if (!is.null(effectName2)) csMat[, effectName2] else NULL,
    diff2               = diff2,
    contrast2           = contrast2,
    interaction2        = interaction2,
    intEffectNames2     = intEffectNames2,
    modEffectNames2     = modEffectNames2,
    modContribution2    = if (!is.null(modEffectNames2)) csMat[, modEffectNames2] else NULL,
    effectNames         = csNames,
    theta               = thetaEff,
    type                = type,
    tieProb             = tieProb,
    details             = details,
    calcRiskRatio       = calcRiskRatio,
    mainEffect          = mainEffect,
    perturbType1        = perturbType1,
    perturbType2        = perturbType2,
    group_id            = changeContributions$group_id
  )

  keep <- density != 0L
  out  <- groupColsList(changeContributions, keep)
  out[names(sd)] <- lapply(sd, `[`, keep)

  if (details) {
    out[["changeUtil"]] <- utility[keep]
    out[["changeProb"]] <- changeProb[keep]
    if (!is.null(tieProb)) out[["tieProb"]] <- tieProb[keep]
  }

  # Mass contrasts for ego-wide perturbations
  if (massContrasts) {
    diffColName <- intersect(c("secondDiff", "secondRiskRatio"), names(sd))[1L]
    if (!is.na(diffColName) && diffColName == "secondDiff") {
      mc <- computeMassContrasts(
        firstDiff = sd[["secondDiff"]][keep],
        density   = density[keep],
        ego       = changeContributions$ego[keep],
        period    = changeContributions$period[keep],
        group     = if (!is.null(changeContributions$group))
                      changeContributions$group[keep] else rep(1L, sum(keep)),
        type      = type
      )
      out[["massCreation"]]    <- mc[["massCreation"]]
      out[["massDissolution"]] <- mc[["massDissolution"]]
    }
  }

  if (attachContribs) {
    out <- attachContributions(out, csNames,
                               csMat[keep, , drop = FALSE], flip = FALSE)
  }
  # faster than as.data.frame for large nrow, and we don't need factors or automatic names here
  attr(out, "row.names") <- .set_row_names(length(out[[1L]]))
  class(out) <- "data.frame"
  out
}

calculateFirstDiff <- function(densityValue,
                               changeProb,
                               changeUtil,
                               effectName,
                               effectContribution,
                               diff = NULL,
                               contrast = NULL, 
                               interaction = FALSE,
                               intEffectNames = NULL,
                               modEffectNames = NULL,
                               modContribution = NULL,
                               effectNames,
                               theta, 
                               type = "changeProb",
                               tieProb = NULL,
                               details = FALSE,
                               calcRiskRatio = FALSE,
                               mainEffect = "firstDiff",
                               perturbType = "alter",
                               group_id = NULL){

  if (effectName == "density") {
    if((!is.null(diff))) stop("firstDiff for density must be contrast c(-1,1)")
    if(!is.null(contrast)){
      if(any(setdiff(contrast, c(-1,1)))) {stop("firstDiff for density can only be be calculated for c(-1,1)")}
      if(interaction == TRUE) {stop("Interaction with density is not possible")}
      oldChangeStatistic <- densityValue
      newChangeStatistic <- rep(NA, length(oldChangeStatistic))
      newChangeStatistic[oldChangeStatistic == contrast[1]] <- contrast[2]
      newChangeStatistic[oldChangeStatistic == contrast[2]] <- contrast[1]
      utilDiff <- ifelse(is.na(newChangeStatistic), NA, -2 * changeUtil)
      densityValue <- newChangeStatistic
    }
  } else {
    if(!is.null(contrast)){
      # effectContribution is in CS space (delta) — density-independent.
      # Contrast values match delta directly for both ego and alter effects.
      oldChangeStatistic <- effectContribution
      newChangeStatistic <- rep(NA, length(oldChangeStatistic))
      newChangeStatistic[oldChangeStatistic == contrast[1]] <- contrast[2]
      newChangeStatistic[oldChangeStatistic == contrast[2]] <- contrast[1]
      if (sum(!is.na(newChangeStatistic)) == 0L)
        warning("contrast: no rows matched the supplied values ",
                "(after auto-centering, if applicable). Check that the ",
                "contrast values correspond to observed change statistics.",
                call. = FALSE)
      diff <- newChangeStatistic - oldChangeStatistic # this is a vector
    }
    utilDiff <- calculateUtilityDiff(effectName = effectName, diff = diff, 
                                      theta = theta, densityValue = densityValue,
                                      interaction = interaction,
                                      intEffectNames = intEffectNames,
                                      modEffectNames = modEffectNames,
                                      modContribution = modContribution,
                                      effectNames = effectNames)
  }

  changeProb_cf <- mlogit_update_r(changeProb, utilDiff, group_id, perturbType)
  # the counterfactual change probability IS well defined, but also not used for density=0 cases
  # but should probably not use NA anyway
  changeProb_cf[densityValue == 0] <- NA
  if (type == "tieProb") {
    tieProb_cf <- changeProb_cf
    idx <- which(!is.na(densityValue) & densityValue == -1)
    if (length(idx) > 0) tieProb_cf[idx] <- 1 - changeProb_cf[idx]
    firstDiff <- tieProb_cf - tieProb
  } else {
     firstDiff <- changeProb_cf - changeProb
  }

  if(!is.null(contrast)){
    idx_flip <- which(newChangeStatistic == min(contrast))
    firstDiff[idx_flip] <- -firstDiff[idx_flip]
  }

  if (calcRiskRatio || mainEffect == "riskRatio") {
    if (type == "tieProb") {
      firstRiskRatio <- tieProb_cf / tieProb
    } else {
      firstRiskRatio <- changeProb_cf / changeProb
    }
    if (!is.null(contrast)) {
      idx_flip <- which(newChangeStatistic == min(contrast))
      firstRiskRatio[idx_flip] <- 1 / firstRiskRatio[idx_flip]
    }
  }

  if(details){ # mostly for debugging
    out <- data.frame(
      "firstDiff" = firstDiff,
      "utilDiff" = utilDiff,
      "newChangeProb" = changeProb_cf,
      "oldChangeProb" = changeProb
    )
    if(type == "tieProb"){
      out[["newTieProb"]] <- tieProb_cf
      out[["oldTieProb"]] <- tieProb
    }
    if(calcRiskRatio|| mainEffect == "riskRatio"){
      out[["firstRiskRatio"]] <- firstRiskRatio
    }
    return(out)
  } else if (mainEffect == "riskRatio") {
    return(list(firstRiskRatio = firstRiskRatio))
  } else {
    return(list(firstDiff = firstDiff))
  }
}

calculateSecondDiff <- function(densityValue,
                                changeProb,
                                changeUtil,
                                effectName1,
                                effectContribution1,
                                diff1 = NULL, 
                                contrast1 = NULL,
                                interaction1 = FALSE,
                                intEffectNames1 = NULL, # later it would be nice to detect these automatically
                                modEffectNames1 = NULL,
                                modContribution1 = NULL,
                                effectName2, 
                                effectContribution2,
                                diff2 = NULL,
                                contrast2 = NULL,
                                interaction2 = FALSE,
                                intEffectNames2 = NULL, # later it would be nice to detect these automatically
                                modEffectNames2 = NULL,
                                modContribution2 = NULL,
                                effectNames,
                                theta,
                                type = "changeProb",
                                tieProb = NULL,
                                details = FALSE,
                                mainEffect = "riskDifference",
                                calcRiskRatio = FALSE,
                                perturbType1 = "alter",
                                perturbType2 = "alter",
                                group_id = NULL){
  firstDiff <- calculateFirstDiff(
    densityValue = densityValue,
    changeProb = changeProb,
    changeUtil = changeUtil,
    effectName = effectName1,
    effectContribution = effectContribution1,
    diff = diff1,
    contrast = contrast1, 
    interaction = interaction1,
    intEffectNames = intEffectNames1,
    modEffectNames = modEffectNames1,
    modContribution = modContribution1,
    effectNames = effectNames,
    theta = theta, 
    type = type,
    tieProb = tieProb,
    mainEffect = mainEffect,
    details = details,
    perturbType = perturbType1,
    group_id = group_id
  )

  if(!is.null(contrast2)){
    # effectContribution2 is in CS space (delta) — density-independent.
    oldChangeStatistic2 <- effectContribution2
    newChangeStatistic2 <- rep(NA, length(oldChangeStatistic2))
    newChangeStatistic2[oldChangeStatistic2 == contrast2[1]] <- contrast2[2]
    newChangeStatistic2[oldChangeStatistic2 == contrast2[2]] <- contrast2[1]
    if (sum(!is.na(newChangeStatistic2)) == 0L)
      warning("contrast2: no rows matched the supplied values ",
              "(after auto-centering, if applicable). Check that the ",
              "contrast values correspond to observed change statistics.",
              call. = FALSE)
    diff2 <- newChangeStatistic2 - oldChangeStatistic2
  }
  utilDiff21 <- calculateUtilityDiff(effectName = effectName2,
                                      diff = diff2, theta = theta, 
                                      densityValue = densityValue,
                                      interaction = interaction2,
                                      intEffectNames = intEffectNames2,
                                      modEffectNames = modEffectNames2,
                                      modContribution = modContribution2,
                                      effectNames = effectNames)
  changeProb_cf21 <- mlogit_update_r(changeProb, utilDiff21, group_id, perturbType2)
  if (type == "tieProb") {
    tieProb_cf21 <- changeProb_cf21
    idx <- which(!is.na(densityValue) & densityValue == -1)
    if (length(idx) > 0) tieProb_cf21[idx] <- 1 - changeProb_cf21[idx]
  }
  changeUtil21 <- changeUtil + utilDiff21

  # Update moderator for effectName1's interaction after the step-2 shift:
  # if the shifted effect (effectName2) IS one of the moderators, that moderator
  # must reflect the post-shift state; otherwise the interaction contribution
  # to firstDiff2 uses stale values and the utility-level SD is missed.
  # Both diff2 and modContribution1 are in CS space — just add directly.
  if (interaction1 && !is.null(modEffectNames1)) {
    mod_shift <- if (is.null(diff2)) 1 else diff2
    mod_shift[is.na(mod_shift)] <- 0
    match_idx <- which(modEffectNames1 == effectName2)
    if (length(match_idx) > 0L) {
      if (is.matrix(modContribution1)) {
        for (mi in match_idx) modContribution1[, mi] <- modContribution1[, mi] + mod_shift
      } else {
        modContribution1 <- modContribution1 + mod_shift
      }
    }
  }

  # Calculate first difference of changing effect1 if effect2 was already changed
  firstDiff2 <- calculateFirstDiff(
    densityValue = densityValue, # careful if one was density!
    changeProb = changeProb_cf21,
    changeUtil = changeUtil21,
    effectName = effectName1, 
    diff = diff1,
    contrast = contrast1,
    effectContribution = effectContribution1,
    theta = theta, 
    type = type,
    tieProb = tieProb_cf21,
    interaction = interaction1,
    intEffectNames = intEffectNames1,
    modEffectNames = modEffectNames1,
    modContribution = modContribution1,
    effectNames = effectNames,
    mainEffect = mainEffect,
    details = details,
    perturbType = perturbType1,
    group_id = group_id
  )

  secondDiff <- firstDiff2[["firstDiff"]] - firstDiff[["firstDiff"]]
  if(!is.null(contrast2)){
    secondDiff[which(newChangeStatistic2 == min(contrast2))] <- -secondDiff[which(newChangeStatistic2 == min(contrast2))]
  }

  if (mainEffect == "riskRatio") {
      secondRiskRatio <- firstDiff2[["firstRiskRatio"]] / firstDiff[["firstRiskRatio"]]
    if (!is.null(contrast2)) {
      idx_flip <- which(newChangeStatistic2 == min(contrast2))
      secondRiskRatio[idx_flip] <- 1 / secondRiskRatio[idx_flip]
    }
  }
  ## really use data.frame here?
  if(details){
    out <- data.frame(
      "changeProb_base" = changeProb,        
      "changeProb_main" = firstDiff[["newChangeProb"]],        
      "changeProb_mod" = changeProb_cf21,       
      "changeProb_both" = firstDiff2[["newChangeProb"]],         
      "firstDiff1" = firstDiff[["firstDiff"]],
      "firstDiff2" = firstDiff2[["firstDiff"]],
      "secondDiff" = secondDiff
    )
    if (type == "tieProb") {
      out$tieProb_base  <- tieProb
      out$tieProb_main  <- firstDiff[["newTieProb"]]
      out$tieProb_mod   <- firstDiff2[["newTieProb"]]
      out$tieProb_both  <- tieProb_cf21
    }
    if (mainEffect == "riskRatio" || calcRiskRatio) {
      out$firstRiskRatio1 <- firstDiff[["firstRiskRatio"]]
      out$firstRiskRatio2 <- firstDiff2[["firstRiskRatio"]]
      out$secondRiskRatio <- secondRiskRatio
    }
    return(out)
  } else if (mainEffect == "riskRatio") {
    return(list(secondRiskRatio = secondRiskRatio))
  } else {
    return(list(secondDiff = secondDiff))
  }
}

# Compute the utility shift from perturbing effectName by diff.
#
# diff is in "network feature" space (absolute), e.g. +1 = one more unit of
# the structural feature.  The signed change in csMat is d*diff (since
# csMat stores signed Δs).  So: Δu = d * diff * θ_combined.
#
# theta: either a named numeric vector (legacy) or a nEffects x 2 matrix
#   ("thetaEff") with columns "creation" and "dissolution" (changeStats).
#   When matrix, direction-dependent theta is selected per-row via densityValue.
calculateUtilityDiff <- function(effectName, diff, 
                                 theta, densityValue,
                                 interaction = FALSE,
                                 intEffectNames = NULL,
                                 modEffectNames = NULL,
                                 modContribution = NULL,
                                 effectNames = NULL){
  if (is.null(diff)) diff <- 1  # NULL means "+1 unit" perturbation
  effectNum <- which(effectNames == effectName)

  # Helper: extract direction-dependent scalar/vector from theta.
  thetaForEffect <- function(eNum) {
    if (is.matrix(theta)) {
      # Direction-dependent: select per-row by density.
      th <- rep(0, length(densityValue))
      th[densityValue ==  1] <- theta[eNum, "creation"]
      th[densityValue == -1] <- theta[eNum, "dissolution"]
      th
    } else {
      theta[eNum]
    }
  }

  if(interaction == TRUE){
    if (is.null(intEffectNames))
      stop("'intEffectNames' must not be NULL when interaction = TRUE.")
    if (is.null(modEffectNames))
      stop("'modEffectNames' must not be NULL when interaction = TRUE.")
    inner <- thetaForEffect(effectNum)
    nInt <- length(intEffectNames)
    for (k in seq_len(nInt)) {
      mod_k <- if (is.matrix(modContribution)) modContribution[, k] else modContribution
      int_num <- which(effectNames == intEffectNames[k])
      inner <- inner + mod_k * thetaForEffect(int_num)
    }
    util_diff <- densityValue * diff * inner
  } else {
    util_diff <- densityValue * diff * thetaForEffect(effectNum)
  }
  util_diff
 }

# Resolve a user-supplied effect name to the changeStats names used in
# contributions/prediction pipelines.
# changeStats names have no type suffix (e.g. "mynet_recip", not "mynet_recip_eval").
# - Accepts exact changeStats names (e.g. "mynet_recip").
# - Accepts shorthand names (e.g. "recip").
# - Accepts legacy names with type suffix (e.g. "recip_eval", "mynet_recip_eval")
#   — the type suffix is stripped before matching.
resolveEffectName <- function(effectName, effectNames) {
  if (is.null(effectName)) return(NULL)

  resolve_one <- function(nm) {
    if (nm %in% effectNames) return(nm)

    # Strip type suffix if present (legacy names).
    nmStripped <- sub("_(eval|endow|creation)$", "", nm)
    if (nmStripped %in% effectNames) return(nmStripped)

    # Plain shortname: match as suffix of changeStats name.
    # e.g. "recip" matches "mynet_recip", "egoX" matches "mynet_egoX_mybeh"
    if (!grepl("_", nmStripped, fixed = TRUE)) {
      pattern <- paste0("(^|_)", nmStripped, "(_[^_]+)*$")
      m <- grep(pattern, effectNames, perl = TRUE, value = TRUE)
      if (length(m) == 1L) return(m)
      if (length(m) > 1L) {
        stop("Effect '", nm, "' is ambiguous. Matches: ",
             paste(m, collapse = ", "), call. = FALSE)
      }
    } else {
      # Partial qualified name (e.g. "mynet_recip" or "recip_covar").
      m <- grep(paste0("(^|_)", nmStripped, "$"), effectNames, perl = TRUE, value = TRUE)
      if (length(m) == 1L) return(m)
      if (length(m) > 1L) {
        stop("Effect '", nm, "' is ambiguous. Matches: ",
             paste(m, collapse = ", "), call. = FALSE)
      }
    }

    stop("Effect '", nm, "' not found in contribMat columns: ",
         paste(effectNames, collapse = ", "), call. = FALSE)
  }

  vapply(effectName, resolve_one, character(1L), USE.NAMES = FALSE)
}

# Helper: resolve the effect-name arguments of marginalEffects against a known column
# list.  Returns a named list with the same structure used in predictArgs.
resolveAMEEffectNames <- function(effectNames,
                                  effectName1, intEffectNames1, modEffectNames1,
                                  effectName2, intEffectNames2, modEffectNames2,
                                  second) {
  list(
    effectName1     = resolveEffectName(effectName1,      effectNames),
    intEffectNames1 = resolveEffectName(intEffectNames1,  effectNames),
    modEffectNames1 = resolveEffectName(modEffectNames1,  effectNames),
    effectName2     = if (second) resolveEffectName(effectName2,     effectNames) else NULL,
    intEffectNames2 = if (second) resolveEffectName(intEffectNames2, effectNames) else NULL,
    modEffectNames2 = if (second) resolveEffectName(modEffectNames2, effectNames) else NULL
  )
}

# Multinomial-logit probability update with string-based perturbation type.
#
# Thin wrapper around the Rcpp mlogit_update that accepts
# perturbType as "alter" or "ego" (instead of 0L/1L)
# and returns a plain numeric vector.
#
# p:           Numeric vector of baseline probabilities.
# delta_u:     Numeric vector of utility shifts (same length as p).
# group_id:    Integer vector of group identifiers.
# perturbType: Character: "alter" (one-alternative update)
#              or "ego" (ego-wide renormalization).
# Returns numeric vector of updated probabilities.
mlogit_update_r <- function(p, delta_u, group_id, perturbType) {
    if (is.null(group_id)) group_id <- integer(length(p))
    as.vector(mlogit_update(p, delta_u, group_id,
                            perturbTypeToInt(perturbType)))
}