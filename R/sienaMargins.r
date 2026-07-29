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
    mainEffect = "riskDifference", # allow objective later as well
    perturbType1 = NULL,
    perturbType2 = NULL,
    massContrasts = NULL,
    egoNormalize = TRUE,
    accumulated = FALSE,
    rateWeight = FALSE,
    returnDecisionDetails = FALSE,
    dynamic = FALSE,
    algorithm = NULL,
    n3 = 200,
    n3PointEst = NULL,
    n3BatchSize = 100L,
    chainStoreMode = c("auto", "disk", "memory"),
    useChangeContributions = FALSE,
    chainStorePath = NULL,
    uncertainty = TRUE,
    nsim = 1000,
    uncertaintySd = TRUE,
    uncertaintyCi = TRUE,
    uncertaintyMean = FALSE,
    uncertaintyMedian = FALSE,
    ciInterval = c(0.025, 0.975),
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
    uncertaintyMode = c("bootstrap", "delta", "deltaFull"),
    combineSameLevel = TRUE,
    format = c("long", "wide"),
    targets = NULL,
    control_uncertainty = NULL,
    control_algo = NULL,
    control_out = NULL,
    ...
) {
    # ---- configuration objects ---------------------------------------------
    # 'uncertainty' accepts a sienaPostestUncertainty object as well as a
    # logical; 'control' takes a sienaPostestControl.  Both are unpacked into
    # the individual variables the rest of the function already uses, so no
    # code below this block needs to know which form the caller used.
    #
    # Supplying a configuration object AND one of the individual arguments it
    # covers is an error rather than a silent precedence rule: that
    # combination only arises from a partially converted call, where quietly
    # ignoring the stray argument would change results without saying so.
    .mc    <- match.call()
    .given <- function(nms) intersect(nms, names(.mc)[-1L])

    if (!is.null(targets)) {
        .lowered <- .targetsToEffectList(targets)
        .clash <- .given(c("effectList", "effectName1", "effects", "depvar",
                           names(.lowered$defaults)))
        if (length(.clash))
            stop("'targets' was supplied, so these arguments must be set on ",
                 "the targets object instead of passed separately: ",
                 paste(.clash, collapse = ", "), ".", call. = FALSE)
        effectList <- .lowered$effectList
        effects    <- attr(targets, "effects")
        depvar     <- attr(targets, "depvar")
        for (.nm in names(.lowered$defaults))
            assign(.nm, .lowered$defaults[[.nm]])
    }

    if (!is.null(control_uncertainty)) {
        if (!inherits(control_uncertainty, "sienaPostestUncertainty"))
            stop("'control_uncertainty' must be a sienaPostestUncertainty ",
                 "object, as returned by set_postest_uncertainty_saom().",
                 call. = FALSE)
        if ("uncertainty" %in% names(.mc)[-1L])
            stop("Supply either 'control_uncertainty' or 'uncertainty', ",
                 "not both.", call. = FALSE)
        uncertainty <- control_uncertainty
    }

    if (inherits(uncertainty, "sienaPostestUncertainty")) {
        .clash <- .given(c("nsim", "uncertaintySd", "uncertaintyCi",
                           "uncertaintyMean", "uncertaintyMedian",
                           "ciInterval", "uncertaintyMode"))
        if (length(.clash))
            stop("'uncertainty' was given as a sienaPostestUncertainty object, ",
                 "so these arguments must be set inside it instead of passed ",
                 "separately: ", paste(.clash, collapse = ", "), ".",
                 call. = FALSE)
        .u <- uncertainty
        uncertainty       <- .u$enabled
        uncertaintyMode   <- .u$mode
        nsim              <- .u$nsim
        uncertaintySd     <- .u$sd
        uncertaintyCi     <- .u$ci
        ciInterval        <- .u$ciInterval
        uncertaintyMean   <- .u$simMean
        uncertaintyMedian <- .u$simMedian
    }

    if (!is.null(control_algo)) {
        if (!inherits(control_algo, "sienaPostestControl"))
            stop("'control_algo' must be a sienaPostestControl object, as ",
                 "returned by set_postest_algo_saom().", call. = FALSE)
        .clash <- .given(names(control_algo))
        if (length(.clash))
            stop("'control_algo' was supplied, so these arguments must be set ",
                 "inside it instead of passed separately: ",
                 paste(.clash, collapse = ", "), ".", call. = FALSE)
        # Every field of the object is named for the variable it feeds.
        for (.nm in names(control_algo)) assign(.nm, control_algo[[.nm]])
    }

    if (!is.null(control_out)) {
        if (!inherits(control_out, "sienaPostestOutput"))
            stop("'control_out' must be a sienaPostestOutput object, as ",
                 "returned by set_postest_output_saom().", call. = FALSE)
        .clash <- .given(names(control_out))
        if (length(.clash))
            stop("'control_out' was supplied, so these arguments must be set ",
                 "inside it instead of passed separately: ",
                 paste(.clash, collapse = ", "), ".", call. = FALSE)
        for (.nm in names(control_out)) assign(.nm, control_out[[.nm]])
    }

    if (inherits(data, "sienaGroup"))
      stop("marginalEffects does not support multi-group data (sienaGroup).")
    type           <- match.arg(type)
    clusterType    <- match.arg(clusterType)
    uncertaintyMode <- match.arg(uncertaintyMode)
    chainStoreMode <- match.arg(chainStoreMode)
    format         <- match.arg(format)
    if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
    if (!is.null(condition)) condition <- resolveCondition(condition)

    # ---- behaviour DV guard ----
    dvType <- attr(data[["depvars"]][[depvar]], "type")
    if (!is.null(dvType) && dvType == "behavior")
      stop("marginalEffects is currently only implemented for network ",
           "dependent variables. Behaviour variable '", depvar,
           "' was selected.")

    if (dynamic && is.null(algorithm))
      stop("'algorithm' must be provided when dynamic = TRUE")

    # nsim <= 0 means point estimate only for bootstrap uncertainty.
    # Delta and deltaFull compute uncertainty analytically and should remain enabled.
    if (nsim <= 0L && uncertaintyMode == "bootstrap") uncertainty <- FALSE

    # ---- accumulated / rateWeight guards ----
    if (accumulated && !dynamic)
      stop("'accumulated = TRUE' requires 'dynamic = TRUE' ",
           "(accumulated ME sums over simulated ministep chains).")
    if (accumulated && type != "tieProb")
      message("Note: canonical accumulated ME uses type='tieProb' ",
              "(joined tie-presence ME). You have type='", type, "'.")
    if (rateWeight && dynamic)
      message("Note: for dynamic estimation, rate-weighting is absorbed ",
              "by the simulation (actors selected proportional to lambda). ",
              "'rateWeight' has no additional effect and is ignored.")

    if(dynamic && returnDecisionDetails)
      stop("returnDecisionDetails = TRUE is not supported when dynamic = TRUE.")
    # Two resons: dynamic decision details can explode in memory, and are unsafe
    # and batching mode does not support decision details anyway

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

        # ---- Theta / covariance with proper rate handling ----
        # Dynamic path: getDynamicChangeContributions re-runs siena07, which
        # needs rate parameters in effects$initialValue.  Unconditional models
        # store rates inside theta; conditional models store them separately
        # in object$rate with effects rows in attr(object$f, "condEffects").

        # ---- Detect rateWeight before theta block (needed for static branch) ----
        anyRateWeight <- rateWeight ||
            any(vapply(effectList, function(s) isTRUE(s$rateWeight), logical(1L)))

        if (dynamic && isTRUE(object$cconditional)) {
            # Conditional estimation: rates are fixed (not estimated jointly).
            # Splice condEffects back into effects so siena07 has rate rows.
            condEff <- attr(object$f, "condEffects")
            if (!is.null(condEff) && nrow(condEff) > 0L) {
                condEff$initialValue <- object$rate
                effects <- rbind(condEff, effects)
            }
            # theta: non-rate only (rates are in effects$initialValue above).
            thetaHat <- coef(object)
            covTheta <- vcov(object)
            # covTheta stays non-rate-only: uncertainty draws do not perturb
            # the fixed rate parameters.
        } else if (dynamic) {
            # Unconditional: include rate parameters in theta/covTheta so that
            # getDynamicChangeContributions and uncertainty draws get the full
            # parameter vector.  Avoids the fragile backfill fallback.
            thetaHat <- coef(object, dropRates = FALSE)
            covTheta <- vcov(object, dropRates = FALSE)
        } else {
            # Static path: include basic rate parameters in the draw distribution
            # when rateWeight is active on an unconditional model, so that lambda
            # uncertainty is captured jointly with behavioral param uncertainty.
            # For conditional models rates are fixed by design: no change needed.
            # predictProbability / predictFirstDiff use name-based theta lookup so
            # extra rate entries in theta are silently ignored for all other specs.
            if (anyRateWeight && !isTRUE(object$cconditional)) {
                thetaHat <- coef(object, dropRates = FALSE)
                covTheta <- vcov(object, dropRates = FALSE)
            } else {
                thetaHat <- coef(object)
                covTheta <- vcov(object)
            }
        }

        # ---- Rate parameters for rateWeight (static path) ----
        # rateParams: point-estimate rates (for hat draw and conditional fallback).
        # rateIdx:    integer positions of basic rate params within the
        #             dropRates=FALSE theta vector; NULL for conditional models.
        #             Used in one_batch to extract per-draw rates from theta_sim.
        rateParams <- NULL
        rateIdx    <- NULL
        if (anyRateWeight && !dynamic) {
            if (isTRUE(object$cconditional)) {
                # Conditional: rates fixed at point estimate; rateIdx stays NULL.
                rateParams <- object$rate
            } else {
                eff_df  <- as.data.frame(object$requestedEffects)
                eff_inc <- eff_df[eff_df$include, ]
                # Guard: non-constant rates (structural / covariate) make the
                # scalar lambda_p multiplication incorrect.  Those require a
                # per-actor rate evaluation that is not yet implemented.
                hasNonConstantRates <- any(
                    !eff_inc$basicRate & eff_inc$type == "rate")
                if (hasNonConstantRates)
                    stop(
                        "rateWeight = TRUE is not supported when the model ",
                        "includes non-constant rate effects (structural or ",
                        "covariate-dependent). rateWeight multiplies by the ",
                        "basic period rate only and would ignore actor-level ",
                        "rate heterogeneity.")
                rate_idx   <- which(eff_inc$basicRate)
                theta_full <- coef(object, dropRates = FALSE)
                rateParams <- theta_full[rate_idx]
                rateIdx    <- rate_idx
            }
            if (length(rateParams) == 0L)
                stop("rateWeight = TRUE but no basic rate parameters found.")
        }

        # ---- saveDir: check for completed effects ----
        # check if everything below is necessary!
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

            # Keep saveDir reload behavior aligned with fresh computation.
            results <- lapply(results, function(r) {
              class(r) <- c("sienaMarginalEffect", class(r))
              r
            })

            if (isTRUE(combineSameLevel) && !.single_effect) {
              loaded_specs <- lapply(results, function(r)
                list(level = attr(r, "level"),
                   condition = attr(r, "condition")))
              names(loaded_specs) <- names(results)
              results <- combinePostestResults(results, loaded_specs, format = format)
              results <- lapply(results, function(r) {
                if (!inherits(r, "sienaMarginalEffect"))
                  class(r) <- c("sienaMarginalEffect", class(r))
                r
              })
            }

                return(if (.single_effect) results[["single"]] else results)
            }
        }
        # currently returns list with 1 element even if all effects are combineSameLevel

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
                includePermitted = TRUE,
                returnWide = TRUE
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
            effectMetaNoRate <- getEffectMetaNoRate(effects, depvar)
            knownEffectNames <- effectMetaNoRate$base_names

            # ---- Resolve chain source into standalone .chains local ----
            # Chains are kept in a standalone local variable (.chains),
            # NEVER attached to `object`.  This ensures that after
            # rm(.chains), the chains are trivially freed
            
            # Sources (in priority order):
            #   1. chainStorePath (caller pre-serialized to disk)
            #   2. object$changeContributions (in-memory from estimation)
            #   3. NULL (simulate fresh chains later)
            .chains <- NULL
            if (!is.null(chainStorePath) && file.exists(chainStorePath)) {
                .chains <- readRDS(chainStorePath)
                if (verbose) message(sprintf(
                    "Loaded %d chains from chainStorePath.",
                    length(.chains)))
            } else if (!is.null(object$changeContributions)) {
                .chains <- object$changeContributions
            }

            # Strip chains from object immediately — object stays lightweight.
            # The caller's fit may still hold chains (their responsibility),
            # but our local `object` does not.
            object$changeContributions <- NULL

            hasStoredChains <- !is.null(.chains) && useChangeContributions

            # dynArgs: lightweight object for uncertainty path (simulate mode).
            # Never carries chains — each uncertainty draw generates fresh ones.
            dynArgs <- list(
                ans      = object,
                data     = data,
                algorithm = algorithm,
                effects  = effects,
                depvar   = depvar,
                n3       = n3,
                returnWide = TRUE
            )
            getContribFun <- function(theta) {
                do.call(getDynamicChangeContributions,
                        c(list(theta = theta), dynArgs))
            }
        }

        # ---- effectInteractionTypes for resolvePerturbType ----
        if (!dynamic) {
            eiTypes <- eitCanon
        } else if ("interactionType" %in% names(effects)) {
          eiTypes <- effectMetaNoRate$interaction_types[knownEffectNames]
          names(eiTypes) <- knownEffectNames
        } else {
            eiTypes <- NULL
        }

        # ---- Build validated per-effect spec list ----
        builtSpecList <- vector("list", length(effectList))
        names(builtSpecList) <- names(effectList)

        for (nm in names(effectList)) {
            spec <- effectList[[nm]]
            eff_second <- isTRUE(spec$second)

            # Per-spec level/condition overrides (fall back to call-level defaults).
            eff_level     <- if ("level"     %in% names(spec)) spec$level
                             else                               level
            eff_condition <- if ("condition" %in% names(spec)) {
                if (!is.null(spec$condition)) resolveCondition(spec$condition)
                else NULL
            } else {
                condition
            }

            # Per-spec accumulated/rateWeight (OR'd with call-level defaults)
            eff_accumulated <- isTRUE(spec$accumulated) || accumulated
            eff_rateWeight  <- isTRUE(spec$rateWeight)  || rateWeight
            if (eff_accumulated && !dynamic)
                stop("Effect '", nm, "': accumulated = TRUE requires ",
                     "dynamic = TRUE.")
            if (eff_rateWeight && dynamic)
                eff_rateWeight <- FALSE   # silently ignored; message already issued

            if (is.null(spec$effectName1))
                stop("Effect '", nm, "' must specify 'effectName1'.")

            effectName1<- resolveEffectName(spec$effectName1,
                               knownEffectNames)
            intEffectNames1<- resolveEffectName(spec$intEffectNames1,
                                 knownEffectNames)
            modEffectNames1<- resolveEffectName(spec$modEffectNames1,
                                 knownEffectNames)
            effectName2<- if (eff_second)
                resolveEffectName(spec$effectName2,
                        knownEffectNames) else NULL
            intEffectNames2<- if (eff_second)
                resolveEffectName(spec$intEffectNames2,
                        knownEffectNames) else NULL
            modEffectNames2<- if (eff_second)
                resolveEffectName(spec$modEffectNames2,
                        knownEffectNames) else NULL

            pt1 <- resolvePerturbType(effectName1, eiTypes,
                                      spec$perturbType1)
            pt2 <- if (eff_second)
                   resolvePerturbType(effectName2, eiTypes,
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
              cm1 <- getCovCenteringMean(effectName1,
                                           effects, data, depvar)
                if (cm1 != 0) eff_contrast1 <- eff_contrast1 - cm1
            }
            if (eff_second && !is.null(eff_contrast2)) {
              cm2 <- getCovCenteringMean(effectName2,
                                           effects, data, depvar)
                if (cm2 != 0) eff_contrast2 <- eff_contrast2 - cm2
            }

            eff_retDet    <- isTRUE(spec$returnDecisionDetails)
            # attachContribs: TRUE = all columns (returnDecisionDetails),
            # character = only those columns (condition), FALSE = none.
            if (eff_retDet) {
                attachVal <- TRUE
            } else if (!is.null(eff_condition)) {
                attachVal <- eff_condition
            } else {
                attachVal <- FALSE
            }

            if (!eff_second) {
                diffFun  <- predictFirstDiff
                diffArgs <- list(
                    type           = type,
                  effectName     = effectName1,
                    diff           = spec$diff1,
                    contrast       = eff_contrast1,
                    interaction    = isTRUE(spec$interaction1),
                  intEffectNames = intEffectNames1,
                  modEffectNames = modEffectNames1,
                    details        = isTRUE(details),
                    calcRiskRatio  = FALSE,
                    mainEffect     = mainEffect,
                    perturbType    = pt1
                )
            } else {
                diffFun  <- predictSecondDiff
                diffArgs <- list(
                    type            = type,
                  effectName1     = effectName1,
                    diff1           = spec$diff1,
                    contrast1       = eff_contrast1,
                    interaction1    = isTRUE(spec$interaction1),
                  intEffectNames1 = intEffectNames1,
                  modEffectNames1 = modEffectNames1,
                  effectName2     = effectName2,
                    diff2           = spec$diff2,
                    contrast2       = eff_contrast2,
                    interaction2    = isTRUE(spec$interaction2),
                  intEffectNames2 = intEffectNames2,
                  modEffectNames2 = modEffectNames2,
                    details         = isTRUE(details),
                    calcRiskRatio   = FALSE,
                    mainEffect      = mainEffect,
                    perturbType1    = pt1,
                    perturbType2    = pt2
                )
            }

            builtSpecList[[nm]] <- makeSpec(
                predictFun            = diffFun,
                outcomeName           = eff_diffName,
                predictArgs           = diffArgs,
                level                 = eff_level,
                condition             = eff_condition,
                accumulated           = eff_accumulated,
                egoNormalize          = egoNormalize,
                rateWeight            = eff_rateWeight,
                dynamic               = dynamic,
                massContrasts         = eff_massC,
                returnDecisionDetails = eff_retDet,
                jacobianFun           = if (!eff_second) predictFirstDiffJac
                                        else NULL,
                metadata = list(
                    method      = "marginalEffects",
                    type        = type,
                  effectName1 = effectName1,
                    contrast1   = eff_contrast1,
                  effectName2 = if (eff_second) effectName2
                                  else NULL,
                    contrast2   = if (eff_second) eff_contrast2
                                  else NULL,
                    second      = eff_second,
                    level       = eff_level,
                    condition   = eff_condition,
                    depvar      = depvar,
                    dynamic     = dynamic,
                    accumulated = eff_accumulated,
                    rateWeight  = eff_rateWeight,
                    nsim        = nsim,
                    n3          = if (dynamic) n3 else NULL,
                    n3PointEst  = if (dynamic) n3PointEst else NULL
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
                clusterType = clusterType,
                unitBudget = batchUnitBudget,
                dynamicMinistepFactor = dynamicMinistepFactor,
                memoryScale = memoryScale,
                useChangeContributions = useChangeContributions
            )
        }

        # ---- Build contribFun (point estimate) via chainStore ----
        hatStore <- NULL
        if (dynamic) {
            n3Total <- if (!is.null(n3PointEst)) n3PointEst else algorithm$n3
            n3Batch <- if (!is.null(n3BatchSize)) min(n3BatchSize, n3Total) else n3Total

            # OOM-GUARD-ENABLED: pre-flight memory check (revert: remove this block)
            memCheck <- .checkDynMemory(
                data         = data,
                depvar       = depvar,
                effects      = effects,
                n3_per_batch = n3Batch,
                n3_uncert    = n3,
                useCluster   = useCluster,
                nbrNodes     = nbrNodes,
                clusterType  = clusterType,
                uncertainty  = uncertainty,
                verbose      = verbose
            )
            nbrNodes <- memCheck$nbrNodes
            # END OOM-GUARD-ENABLED

            if (hasStoredChains) {
                # Trim or warn about chain count mismatch before building store.
                storedN <- length(.chains)
                if (storedN > n3Total) {
                    .chains <- .chains[seq_len(n3Total)]
                } else if (storedN < n3Total) {
                    warning(sprintf(paste0(
                        "Stored chains (%d) < n3PointEst (%d); ",
                        "using all %d stored chains."),
                        storedN, n3Total, storedN),
                        call. = FALSE, immediate. = TRUE)
                    n3Total <- storedN
                }
                if (verbose) message(sprintf(
                    "Point estimate re-uses %d stored chains.", storedN))
                built <- .buildDynChainStore(.chains, dynArgs, n3Total, n3Batch,
                                             depvar, effects, data, verbose,
                                             chainStoreMode = chainStoreMode,
                                             vmaxGB = memCheck$vmaxGB,
                                             nbrNodes = nbrNodes)
                # .chains is the sole reference to the chain data; free after
                # chainStore_disk has written batch files to disk.
                rm(.chains)
                gc(verbose = FALSE)
                # Remove the caller's original tempfile (now superseded by
                # per-batch RDS files) to avoid 2× disk usage.
                if (!is.null(chainStorePath) && file.exists(chainStorePath))
                    unlink(chainStorePath)
            } else {
                built <- .buildDynChainStore(NULL, dynArgs, n3Total, n3Batch,
                                             depvar, effects, data, verbose,
                                             chainStoreMode = chainStoreMode,
                                             vmaxGB = memCheck$vmaxGB,
                                             nbrNodes = nbrNodes)
            }
            hatStore      <- built$store
            contribFun    <- built$contribFun
            nChainBatches <- built$nChainBatches
            .chainMode    <- built$mode
        } else {
            nChainBatches <- 1L
            contribFun <- makeContribFun(getContribFun = getContribFun)
            .chainMode <- "static"
        }

        if (verbose && dynamic)
            message(sprintf("Chain mode: %s, nChainBatches=%d",
                            .chainMode, nChainBatches))

        # Clean up disk-backed chain store on exit (even on error)
        if (!is.null(hatStore) && hatStore$mode == "disk")
            on.exit(hatStore$cleanup(), add = TRUE)

        # ---- Decision details pre-capture (marginalEffects-specific) ----
        # Only possible for single-batch hat (nChainBatches == 1); with chain
        # batching the raw contribs are never all in memory simultaneously.
        anyReturnDetails <- any(vapply(builtSpecList, function(s)
            isTRUE(s$returnDecisionDetails), logical(1L)))
        if (anyReturnDetails && nChainBatches > 1L)
            warning("returnDecisionDetails is ignored with chain batching.",
                    call. = FALSE)
        decision_details <- setNames(
            vector("list", length(builtSpecList)),
            names(builtSpecList))

        if (anyReturnDetails && nChainBatches == 1L) {
            cc <- contribFun(thetaHat, 1L, 1L)
            baseline <- predictProbability(cc, thetaHat, type,
                                           returnComponents = TRUE)
            for (nm in names(builtSpecList)) {
                spec <- builtSpecList[[nm]]
                if (isTRUE(spec$returnDecisionDetails)) {
                    pe_args <- spec$predictArgs
                    pe_args$attachContribs <- TRUE
                    decision_details[[nm]] <- do.call(spec$predictFun,
                        c(list(changeContributions = cc,
                               theta = thetaHat,
                               baseline = baseline,
                               returnComponents = TRUE), pe_args))
                }
            }
            rm(cc, baseline); gc(verbose = FALSE)
        }

        # ---- Delegate to sienaPostestimate ----
        # sienaPostestimate builds both estimatorFuns internally:
        #   hat:         contribFun (preloaded/per_batch/all_at_once, n3PointEst total)
        #   uncertainty: fresh per_batch from dynArgs, n3 chains per theta draw
        if (verbose && dynamic)
            memReport("pre-postestimate (chains freed)", verbose = TRUE)
        results <- sienaPostestimate(
            contribFun    = contribFun,
            nChainBatches = nChainBatches,
            type          = type,
            rateParams    = rateParams,
            rateIdx       = rateIdx,
            specs         = builtSpecList,
            thetaHat      = thetaHat,
            covTheta      = covTheta,
            dynamic       = dynamic,
            dynArgs       = if (dynamic) dynArgs else NULL,
            n3            = if (dynamic) n3 else NULL,
            n3BatchSize   = if (dynamic) n3BatchSize else NULL,
            useChangeContributions = if (dynamic)
                useChangeContributions else FALSE,
            uncertainty      = uncertainty,
            uncertaintyMode  = uncertaintyMode,
            nsim             = nsim,
            batchSize     = batchSize,
            useCluster    = useCluster,
            nbrNodes      = nbrNodes,
            clusterType   = clusterType,
            cl            = cl,
            batchDir      = batchDir,
            prefix        = prefix,
            keepBatch     = keepBatch,
            verbose       = verbose,
            egoNormalize  = egoNormalize,
            uncertaintySd     = uncertaintySd,
            uncertaintyCi     = uncertaintyCi,
            uncertaintyMean   = uncertaintyMean,
            uncertaintyMedian = uncertaintyMedian,
            ciInterval        = ciInterval,
            decisionDetails = decision_details,
            saveDir         = saveDir,
            gcEachBatch     = gcEachBatch,
            gcEachSim       = gcEachSim
        )

    # Assign S3 class on each result; sienaPostestimate returns plain data.frames.
    # this might actually not be fully correct anymore with the recent changes
    results <- lapply(results, function(r) {
        class(r) <- c("sienaMarginalEffect", class(r))
        r
    })

    if (isTRUE(combineSameLevel) && !.single_effect) {
        results <- combinePostestResults(results, builtSpecList, format = format)
        results <- lapply(results, function(r) {
            if (!inherits(r, "sienaMarginalEffect"))
                class(r) <- c("sienaMarginalEffect", class(r))
            r
        })
    }

    if (.single_effect) return(results[["single"]])
    if (length(results) == 1L) return(results[[1L]])
    results
}

##@print.sienaMarginalEffect S3 print
print.sienaMarginalEffect <- function(x, ...) {
  cat("SAOM Marginal Effect\n")
  cat("  Effect:    ", if (!is.null(attr(x, "effectName1"))) attr(x, "effectName1") else "unknown", "\n")
  if (isTRUE(attr(x, "second"))) {
    cat("  Effect 2:  ", if (!is.null(attr(x, "effectName2"))) attr(x, "effectName2") else "unknown", "\n")
  }
  me <- if (!is.null(attr(x, "mainEffect"))) attr(x, "mainEffect") else "riskDifference"
  cat("  Scale:     ", me, "\n")
  cat("  Type:      ", if (!is.null(attr(x, "type"))) attr(x, "type") else "unknown", "\n")
  cat("  Dep. var.: ", if (!is.null(attr(x, "depvar"))) attr(x, "depvar") else "unknown", "\n")
  cat("  Level:     ", if (!is.null(attr(x, "level"))) attr(x, "level") else "unknown", "\n")
  cat("  Dynamic:   ", if (!is.null(attr(x, "dynamic"))) attr(x, "dynamic") else FALSE, "\n")
  cat("  nsim:      ", if (!is.null(attr(x, "nsim"))) attr(x, "nsim") else NA, "\n")
  cat("\n")

  massCols <- intersect(c("massCreation", "massDissolution"), names(x))
  mass_related <- character(0)
  if (length(massCols) > 0L) {
    mass_related <- grep(
      paste0("^(", paste(massCols, collapse = "|"), ")($|_)"),
      names(x), value = TRUE
    )
  }

  main_cols <- setdiff(names(x), mass_related)
  print.data.frame(x[, main_cols, drop = FALSE], row.names = FALSE, ...)

  for (mc in massCols) {
    mc_uc_pattern <- paste0("^", mc, "_")
    mc_uc_cols <- grep(mc_uc_pattern, names(x), value = TRUE)
    level_cols <- intersect(c("period", "group"), names(x))
    mc_all <- c(level_cols, mc, mc_uc_cols)
    mc_df <- unique(x[, mc_all, drop = FALSE])
    names(mc_df) <- sub(mc_uc_pattern, "", names(mc_df))
    cat("\n  ", mc, ":\n")
    print.data.frame(mc_df, row.names = FALSE, ...)
  }

  invisible(x)
}

##@summary.sienaMarginalEffect S3 summary
summary.sienaMarginalEffect <- function(object, ...) {
  cat("SAOM Marginal Effect Summary\n")
  cat("  Effect:    ", if (!is.null(attr(object, "effectName1"))) attr(object, "effectName1") else "unknown", "\n")
  if (isTRUE(attr(object, "second"))) {
    cat("  Effect 2:  ", if (!is.null(attr(object, "effectName2"))) attr(object, "effectName2") else "unknown", "\n")
  }
  me <- if (!is.null(attr(object, "mainEffect"))) attr(object, "mainEffect") else "riskDifference"
  cat("  Scale:     ", me, "\n")
  cat("  Type:      ", if (!is.null(attr(object, "type"))) attr(object, "type") else "unknown", "\n")
  cat("  Dep. var.: ", if (!is.null(attr(object, "depvar"))) attr(object, "depvar") else "unknown", "\n")
  cat("  Level:     ", if (!is.null(attr(object, "level"))) attr(object, "level") else "unknown", "\n")
  cat("  Dynamic:   ", if (!is.null(attr(object, "dynamic"))) attr(object, "dynamic") else FALSE, "\n")
  cat("  nsim:      ", if (!is.null(attr(object, "nsim"))) attr(object, "nsim") else NA, "\n")
  if (!is.null(attr(object, "condition")))
    cat("  Condition: ", paste(attr(object, "condition"), collapse = ", "), "\n")
  cat("\n")
  cat("  Rows:      ", nrow(object), "\n")
  cat("  Columns:   ", paste(names(object), collapse = ", "), "\n\n")
  summary.data.frame(object, ...)
}


# Core computation shared by both static and dynamic paths.
# `changeContributions` — unified wide struct (contribMat + coord vectors + effectNames)
# `theta`               — named numeric vector (full parameter vector; subsetted internally)
predictFirstDiff <- function(changeContributions, theta, type,
    effectName, diff, contrast, interaction, intEffectNames,
    modEffectNames, details, calcRiskRatio, mainEffect,
    perturbType = "alter", massContrasts = FALSE, attachContribs = TRUE,
    baseline = NULL, outcomesOnly = FALSE)
{
  theta_use <- theta[changeContributions$effectNames]
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
  )

  # ---- outcomesOnly: return raw numeric vectors (no structural cols) ------
  # Used by makeEstimatorFun's one_batch to avoid copying structural columns
  # and let the aggregation layer own sorting/grouping.
  if (outcomesOnly) {
    out <- fd  # named list, e.g. list(firstDiff = numeric(n))
    if (massContrasts) {
      diffColName <- intersect(c("firstDiff", "firstRiskRatio"), names(fd))[1L]
      if (!is.na(diffColName) && diffColName == "firstDiff") {
        mc <- computeMassContrasts(
          firstDiff = fd[["firstDiff"]],
          density   = density,
          ego       = changeContributions$ego,
          period    = changeContributions$period,
          group     = if (!is.null(changeContributions$group))
                        changeContributions$group else
                        rep(1L, length(density)),
          type      = type
        )
        out[["massCreation"]]    <- mc[["massCreation"]]
        out[["massDissolution"]] <- mc[["massDissolution"]]
      }
    }
    return(out)
  }

  # ---- full data.frame path (standalone callers) --------------------------
  # Exclude density=0 (no-change) rows and non-permitted creation alters
  # (at-cap egos where maxDegree blocks tie creation).
  # Dynamic chains never emit such rows (C++ sets density=0, CS=0 for them);
  # static must exclude them explicitly so both paths average over the same
  # effective choice set.
  perm <- changeContributions$permitted
  keep <- density != 0L & (if (is.null(perm)) TRUE else perm)
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

  if (is.character(attachContribs)) {
    # Selective: only attach specified columns (e.g. condition column).
    # Avoids copying the full csMat[keep, ] matrix (can be many GB).
    cols <- resolveEffectName(attachContribs, csNames)
    for (cn in cols) out[[cn]] <- csMat[keep, cn]
  } else if (isTRUE(attachContribs)) {
    out <- attachContributions(out, csNames,
                               csMat[keep, , drop = FALSE], flip = FALSE)
  }
  attr(out, "row.names") <- .set_row_names(length(out[[1L]]))
  class(out) <- "data.frame"
  out
}

# --------------------------------------------------------------------------
# predictFirstDiffJac — analytical Jacobian paired with predictFirstDiff.
#
# Computes the per-row n × K_eff Jacobian for first-difference outcomes.
# Called by evalBatchJacobian via spec$jacobianFun or .resolveBuiltinJac.
#
# Standard Jacobian-fn arguments used:
#   cc         — for cc$effectNames, cc$contribMat, cc$group_id
#   theta      — for theta subsetting to effect names
#   changeProb — baseline choice probabilities (used for softmax Jac + mlogit)
#   density    — ±1/0 density values (tieProb scaling)
#   pa         — spec$predictArgs (effectName, contrast, diff, type, interaction)
#   cs         — changeStats (csMat column lookup, changeStatsMap for ud_jac)
# oc is absorbed by ... (not needed; pa$type carries changeProb/tieProb).
#
# Returns: n × K_eff matrix, or NULL when the analytical path is
# unsupported (interaction effects, density/creation-endow terms, or missing
# change-stats column — caller falls back to finite-difference Jacobian).
# --------------------------------------------------------------------------
predictFirstDiffJac <- function(cc, theta, changeProb, density,
                                 pa, cs, ...) {
  eff_cs_name <- pa$effectName
  contrast    <- pa$contrast
  theta_use   <- theta[cc$effectNames]
  Jp          <- softmax_jac_rcpp(changeProb, cc$contribMat, cc$group_id)

  # ── resolve diff_j and sign_j (mirrors calculateFirstDiff) ───────────────
  if (!is.null(contrast)) {
    col_idx <- match(eff_cs_name, colnames(cs$csMat))
    if (is.na(col_idx)) return(NULL)
    effectContrib <- cs$csMat[, col_idx]

    new_cs <- rep(NA_real_, length(effectContrib))
    new_cs[effectContrib == contrast[1L]] <- contrast[2L]
    new_cs[effectContrib == contrast[2L]] <- contrast[1L]
    diff_j <- new_cs - effectContrib       # NA outside contrast range
    sign_j <- ifelse(!is.na(new_cs) & new_cs == contrast[1L], -1L, 1L)
    sign_j[is.na(diff_j)] <- NA_integer_
  } else {
    diff_pa <- pa$diff
    diff_j  <- if (!is.null(diff_pa)) diff_pa else rep(1, length(density))
    sign_j  <- rep(1L, length(density))
  }

  # ── utility-diff Jacobian (returns NULL for unsupported cases) ────────────
  ud_jac <- calculateUtilityDiffJacobian(
    effectName   = eff_cs_name,
    diff_j       = diff_j,
    densityValue = density,
    theta_use    = theta_use,
    interaction  = isTRUE(pa$interaction),
    csMap        = cs$changeStatsMap
  )
  if (is.null(ud_jac)) return(NULL)

  # ── mlogit update Jacobian ────────────────────────────────────────────────
  d_scale <- if (pa$type == "tieProb") density else 1
  J_diff  <- mlogitUpdateJacobian(Jp, changeProb,
                                   ud_jac$delta_u, ud_jac$A_col,
                                   ud_jac$eff_col)
  (sign_j * d_scale) * J_diff   # NA sign_j propagates
}


predictSecondDiff <- function(changeContributions, theta, type,
    effectName1, diff1, contrast1, interaction1,
    intEffectNames1, modEffectNames1,
    effectName2, diff2, contrast2, interaction2,
    intEffectNames2, modEffectNames2,
    details, calcRiskRatio, mainEffect,
    perturbType1 = "alter", perturbType2 = "alter",
    massContrasts = FALSE, attachContribs = TRUE,
    baseline = NULL, outcomesOnly = FALSE)
{
  theta_use <- theta[changeContributions$effectNames]
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

  # ---- outcomesOnly: return raw numeric vectors (no structural cols) ------
  if (outcomesOnly) {
    out <- sd
    if (massContrasts) {
      diffColName <- intersect(c("secondDiff", "secondRiskRatio"), names(sd))[1L]
      if (!is.na(diffColName) && diffColName == "secondDiff") {
        mc <- computeMassContrasts(
          firstDiff = sd[["secondDiff"]],
          density   = density,
          ego       = changeContributions$ego,
          period    = changeContributions$period,
          group     = if (!is.null(changeContributions$group))
                        changeContributions$group else
                        rep(1L, length(density)),
          type      = type
        )
        out[["massCreation"]]    <- mc[["massCreation"]]
        out[["massDissolution"]] <- mc[["massDissolution"]]
      }
    }
    return(out)
  }

  # ---- full data.frame path (standalone callers) --------------------------
  perm <- changeContributions$permitted
  keep <- density != 0L & (if (is.null(perm)) TRUE else perm)
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

  if (is.character(attachContribs)) {
    cols <- resolveEffectName(attachContribs, csNames)
    for (cn in cols) out[[cn]] <- csMat[keep, cn]
  } else if (isTRUE(attachContribs)) {
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
    if (is.null(contrast))
        stop("firstDiff for density requires contrast = c(-1, 1).")
    if(!is.null(contrast)){
      if (!setequal(contrast, c(-1, 1)))
          stop("firstDiff for density requires contrast = c(-1, 1).")
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
    idx_flip <- which(newChangeStatistic == contrast[1L])
    firstDiff[idx_flip] <- -firstDiff[idx_flip]
  }

  if (calcRiskRatio || mainEffect == "riskRatio") {
    if (type == "tieProb") {
      firstRiskRatio <- tieProb_cf / tieProb
    } else {
      firstRiskRatio <- changeProb_cf / changeProb
    }
    if (!is.null(contrast)) {
      idx_flip <- which(newChangeStatistic == contrast[1L])
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

  # ---- Build intermediates for the "moderated" counterfactual ----
  # Memory note: at n3=1000, every 186M-element double vector is ~1.5 GB.
  # calculateSecondDiff calls calculateFirstDiff twice.  Between the calls
  # we free vectors that are no longer needed and trigger gc() so that the
  # second call does not push peak memory past the R_MAX_VSIZE limit.

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
    # Precompute sign-flip indices; free the full vectors early.
    contrast2_flip_idx <- which(newChangeStatistic2 == contrast2[1L])
    rm(oldChangeStatistic2, newChangeStatistic2)
  }
  # effectContribution2 is consumed; free it.
  rm(effectContribution2)

  utilDiff21 <- calculateUtilityDiff(effectName = effectName2,
                                      diff = diff2, theta = theta, 
                                      densityValue = densityValue,
                                      interaction = interaction2,
                                      intEffectNames = intEffectNames2,
                                      modEffectNames = modEffectNames2,
                                      modContribution = modContribution2,
                                      effectNames = effectNames)
  rm(modContribution2)

  changeProb_cf21 <- mlogit_update_r(changeProb, utilDiff21, group_id, perturbType2)
  if (type == "tieProb") {
    tieProb_cf21 <- changeProb_cf21
    idx <- which(!is.na(densityValue) & densityValue == -1)
    if (length(idx) > 0) tieProb_cf21[idx] <- 1 - changeProb_cf21[idx]
  }
  changeUtil21 <- changeUtil + utilDiff21
  rm(utilDiff21)

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
  rm(diff2)

  # Free intermediates before the second calculateFirstDiff so that
  # the second call's temporaries do not push peak memory over the limit.
  # Cost of gc() is proportional to live objects (~2 sec per 15 GB heap);
  # only worthwhile when live vectors are large enough to matter.
  if (length(densityValue) > 1e7) invisible(gc())

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
    secondDiff[contrast2_flip_idx] <- -secondDiff[contrast2_flip_idx]
  }

  if (mainEffect == "riskRatio") {
      secondRiskRatio <- firstDiff2[["firstRiskRatio"]] / firstDiff[["firstRiskRatio"]]
    if (!is.null(contrast2)) {
      secondRiskRatio[contrast2_flip_idx] <- 1 / secondRiskRatio[contrast2_flip_idx]
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

# Analytical Jacobian components for utility-difference perturbations.
#
# Returns NULL for unsupported analytical cases so caller can fall back to
# finite-difference Jacobians.
calculateUtilityDiffJacobian <- function(effectName, diff_j, densityValue,
                                         theta_use, interaction, csMap) {
  if (interaction) return(NULL)
  if (grepl("density", effectName, fixed = TRUE)) return(NULL)

  members <- which(csMap$bases == effectName)
  if (length(members) == 0L) return(NULL)

  m_types <- csMap$types[members]
  if (any(m_types %in% c("creation", "endow"))) return(NULL)

  eff_col <- members[m_types == "eval"]
  if (length(eff_col) != 1L) return(NULL)

  diff_j_0 <- ifelse(is.na(diff_j), 0, diff_j)
  delta_u  <- densityValue * diff_j_0 * theta_use[eff_col]
  A_col    <- densityValue * diff_j_0

  list(delta_u = delta_u, eff_col = eff_col, A_col = A_col)
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

# Jacobian for first-difference mlogit update with respect to theta columns.
mlogitUpdateJacobian <- function(Jp, cp, delta_u, A_col, eff_col) {
  e_j <- exp(delta_u)
  D_j <- 1 - cp + cp * e_j
  w   <- e_j / D_j^2 - 1
  v   <- cp * e_j * (1 - cp) / D_j^2

  J <- w * Jp
  J[, eff_col] <- J[, eff_col] + v * A_col
  J
}

# Resolve the perturbation type for a named effect.
#
# Used by marginalEffects (and potentially predict) to
# decide whether a counterfactual shift uses the one-alternative
# ("alter") or the ego-wide ("ego") mlogit update.
#
# Auto-detection relies on the interactionType metadata stored
# in the wide struct's effectInteractionTypes vector:
#   "ego"    -> perturbType "ego"
#   "dyadic" -> perturbType "alter"
#   anything else (empty, "OK", structural effects) ->
#             perturbType "alter" (safe default for network effects).
#
# effectName: Character: (resolved) composite effect name.
# effectInteractionTypes: Named character vector; may be NULL.
# override: Optional user override, one of "alter", "ego", or NULL (auto-detect).
# Returns character scalar: "alter" or "ego".
resolvePerturbType <- function(effectName,
                               effectInteractionTypes = NULL,
                               override = NULL) {
    if (!is.null(override)) {
        override <- match.arg(override, c("alter", "ego"))
        return(override)
    }
    if (is.null(effectInteractionTypes) ||
        !(effectName %in% names(effectInteractionTypes))) {
        return("alter")
    }
    iType <- effectInteractionTypes[[effectName]]
    if (identical(iType, "ego")) "ego" else "alter"
}





# Convert perturbType string to the integer code expected by
# the mlogit_update Rcpp function.
# perturbType: "alter" or "ego".
# Returns integer: 0L for "alter", 1L for "ego".
perturbTypeToInt <- function(perturbType) {
    switch(perturbType, alter = 0L, ego = 1L,
           stop("Invalid perturbType: ", perturbType))
}

# Compute ego-level probability-mass contrasts.
#
# For ego-wide perturbations, the natural actor-level QOIs are the
# creation and dissolution probability-mass contrasts
# Delta P_i+ = sum_{j in C_i} (p'_ij - p_ij) and
# Delta P_i- = sum_{j in D_i} (p'_ij - p_ij),
# where the sums run over the creation and dissolution risk sets
# respectively.
#
# firstDiff: Numeric vector of dyad-level first differences
#            (on whichever scale: changeProb or tieProb).
# density:   Integer vector: 1 = creation, -1 = dissolution.
#            Rows with density = 0 should already be removed.
# ego:       Integer vector of ego identifiers.
# period:    Integer vector of period identifiers.
# group:     Integer/character vector of group identifiers.
# type:      Character: "changeProb" or "tieProb".
#            Needed to recover change-probability-scale diffs for tieProb mode.
# Returns a data.frame with columns massCreation and
#   massDissolution, one value per row (broadcast from ego-level aggregate).
computeMassContrasts <- function(firstDiff, density, ego, period, group,
                                 type = "changeProb") {
    n <- length(firstDiff)

    # Convert firstDiff to change-probability scale if needed.
    changeProbDiff <- firstDiff
    if (type == "tieProb") {
        diss <- density == -1L
        changeProbDiff[diss] <- -changeProbDiff[diss]
    }

    # Integer group key: encode (group, period, ego) as a single integer.
    # This avoids paste() on millions of rows.
    nPeriod <- max(period)
    nEgo    <- max(ego)
    grpKey  <- (as.integer(group) - 1L) * (nPeriod * nEgo) +
               (as.integer(period) - 1L) * nEgo +
               as.integer(ego)

    # Replace NA diffs with 0 for rowsum (rowsum has no na.rm).
    cpd <- changeProbDiff
    cpd[is.na(cpd)] <- 0

    # Sum by group × direction using rowsum (compiled C, very fast).
    cre <- density ==  1L
    dis <- density == -1L
    mc_sums <- rowsum(cpd[cre], grpKey[cre], reorder = FALSE)
    md_sums <- rowsum(cpd[dis], grpKey[dis], reorder = FALSE)

    # Broadcast sums back to every row via match().
    mc_keys <- as.integer(rownames(mc_sums))
    md_keys <- as.integer(rownames(md_sums))
    massCreation    <- mc_sums[match(grpKey, mc_keys)]
    massDissolution <- md_sums[match(grpKey, md_keys)]
    massCreation[is.na(massCreation)]       <- 0
    massDissolution[is.na(massDissolution)] <- 0

    data.frame(massCreation = massCreation, massDissolution = massDissolution)
}

# --------------------------------------------------------------------------
# combinePostestResults — opt-in wide merge for specs sharing (level, condition)
#
# When multiple effects in `results` were computed at the same (level, condition)
# and the user opts in via combineSameLevel = TRUE, merge them into a single
# wider data.frame rather than returning N separate ones.
#
# `results`: named list of data.frames (output of aggregatePostEstimation or
#            sienaPostestimate). Each element may carry class-level attributes.
# `specs`:   named list of spec entries (same keys as results).
#
# Returns: named list of data.frames. Each group-of-same-level specs is replaced
# by a single wider data.frame keyed on group columns; elements that are alone at
# their (level, condition) stay as single-element lists (unchanged).
# --------------------------------------------------------------------------
combinePostestResults <- function(results, specs, format = c("long", "wide")) {
  format <- match.arg(format)
  if (length(results) == 0L) return(results)

  group_key <- vapply(names(specs), function(nm) {
    sp  <- specs[[nm]]
    lv  <- if (!is.null(sp$level)) sp$level else "none"
    cnd <- if (!is.null(sp$condition)) paste(sp$condition, collapse = ",") else ""
    paste0(lv, "|", cnd)
  }, character(1L))

  out <- list()
  for (gk in unique(group_key)) {
    nms <- names(group_key)[group_key == gk]
    if (length(nms) == 1L) {
      out[[nms]] <- results[[nms]]
      next
    }

    sp       <- specs[[nms[1L]]]
    lv       <- if (!is.null(sp$level)) sp$level else "none"
    cnd      <- if (!is.null(sp$condition))
                  resolveEffectName(sp$condition, names(results[[nms[1L]]])) else NULL
    key_cols <- getGroupVars(level = lv, condition = cnd)
    combined_nm <- paste(nms, collapse = "+")

    if (format == "long") {
      # Long format stacks targets as rows, so every frame must carry the same
      # columns.  Targets differ in what their quantity column is called --
      # firstDiff, secondDiff, firstRiskRatio, ... -- which would make rbind
      # fail outright (base rbind requires identical names; it does not union
      # and NA-fill).
      #
      # Rather than pad with NAs, the quantity column is renamed to `est` in
      # every frame.  Which quantity a row holds is already carried by the
      # `effect` column, whose default names encode it (`recip_fd`,
      # `recip_transTrip1_sd`), so the information is not lost and the table
      # stays one value per row -- which is what makes it readable as a report.
      frames <- lapply(nms, function(nm) {
        df <- results[[nm]]
        oc <- specs[[nm]]$outcomeName
        if (!is.null(oc) && oc %in% names(df))
          names(df)[match(oc, names(df))] <- "est"
        df[["effect"]] <- nm
        df
      })
      common <- Reduce(intersect, lapply(frames, names))
      dropped <- setdiff(unique(unlist(lapply(frames, names))), common)
      if (length(dropped))
        warning("Combining targets at the same level dropped column(s) not ",
                "shared by all of them: ", paste(dropped, collapse = ", "),
                ". Use combineSameLevel = FALSE to keep each target separate.",
                call. = FALSE)
      frames <- lapply(frames, function(df) df[, common, drop = FALSE])
      merged <- do.call(rbind, frames)
      # Reorder for reading: effect first, then the keys it is broken down by
      # (period, condition, ...), then the quantities.  One row per
      # effect x condition x period, which is the shape a results table wants.
      other_cols <- setdiff(names(merged), c(key_cols, "effect"))
      merged <- merged[, c("effect", key_cols, other_cols), drop = FALSE]

    } else {
      # Wide: cbind non-key columns with effect-name suffix.
      # All frames share identical key columns (same chain/ministep/ego/...
      # rows in the same order) so we can cbind directly without merge.
      # Verify row count consistency first.
      nrows <- vapply(nms, function(nm) nrow(results[[nm]]), integer(1L))
      if (length(unique(nrows)) > 1L)
        stop("combinePostestResults (wide): frames for group '", gk,
             "' have differing row counts: ",
             paste(nms, nrows, sep = "=", collapse = ", "),
             ". Use format='long' instead.")

      # Start with key columns from first frame
      first      <- results[[nms[1L]]]
      merged     <- first[, key_cols, drop = FALSE]

      for (nm in nms) {
        df          <- results[[nm]]
        value_cols  <- setdiff(names(df), key_cols)
        suffix_cols <- paste0(value_cols, "_", nm)
        extra       <- df[, value_cols, drop = FALSE]
        # Rename via setNames to avoid match() fragility with duplicate names
        names(extra) <- suffix_cols
        merged <- cbind(merged, extra)
      }
    }

    attr(merged, "row.names") <- .set_row_names(nrow(merged))
    class(merged) <- "data.frame"
    out[[combined_nm]] <- merged
  }
  out
}