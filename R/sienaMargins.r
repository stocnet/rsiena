##@marginalEffects Generic
marginalEffects <- function(object, ...) UseMethod("marginalEffects", object)

##@marginalEffects.sienaFit Method
marginalEffects.sienaFit <- function(
    object,
    data,
    targets,
    control_uncertainty = set_postest_uncertainty_saom(),
    control_algo        = set_postest_algo_saom(),
    control_out         = set_postest_output_saom()
) {
    # ---- configuration objects ---------------------------------------------
    # Everything the computation needs arrives in one of four objects, each
    # answering a different question -- mirroring siena()'s own inputs:
    #
    #   targets              WHAT do we want        (built from the effects object)
    #   control_uncertainty  how precisely
    #   control_algo         HOW do we get there    (built from the algorithm object)
    #   control_out          what to do with the results
    #
    # Each is unpacked into the individual variables the rest of the function
    # uses, so nothing below needs to know where a setting came from.
    #
    # Unpacking ORDER is load-bearing: control_out goes last, so where two
    # objects name the same variable the output setting wins.  A field living
    # in two objects at once is a bug either way -- it was one for `details`
    # and for `combineSameLevel` -- but the order at least makes the winner
    # predictable rather than accidental.
    #
    # The constructors validate their own arguments (match.arg on type,
    # mainEffect, mode, chainStoreMode, clusterType, format), so there is no
    # re-validation here.
    if (missing(targets))
        stop("'targets' is required: build it with make_marginal_targets().",
             call. = FALSE)
    if (!inherits(targets, "sienaPostestTargets"))
        stop("'targets' must be a sienaPostestTargets object, as returned by ",
             "make_marginal_targets().", call. = FALSE)

    .lowered   <- .targetsToEffectList(targets)
    effectList <- .lowered$effectList
    effects    <- attr(targets, "effects")
    depvar     <- attr(targets, "depvar")
    for (.nm in names(.lowered$defaults))
        assign(.nm, .lowered$defaults[[.nm]])

    if (!inherits(control_uncertainty, "sienaPostestUncertainty"))
        stop("'control_uncertainty' must be a sienaPostestUncertainty object, ",
             "as returned by set_postest_uncertainty_saom().", call. = FALSE)
    .u                <- control_uncertainty
    uncertainty       <- .u$enabled
    uncertaintyMode   <- .u$mode
    nsim              <- .u$nsim
    uncertaintySd     <- .u$sd
    uncertaintyCi     <- .u$ci
    ciInterval        <- .u$ciInterval
    uncertaintyMean   <- .u$simMean
    uncertaintyMedian <- .u$simMedian

    if (!inherits(control_algo, "sienaPostestControl"))
        stop("'control_algo' must be a sienaPostestControl object, as ",
             "returned by set_postest_algo_saom().", call. = FALSE)
    # Every field of the object is named for the variable it feeds.
    for (.nm in names(control_algo)) assign(.nm, control_algo[[.nm]])

    # One seed governs the whole call: the simulations (via the algorithm) and
    # the bootstrap draws.  Falls back to the caller's RNG stream, so a plain
    # set.seed() before this call reproduces everything.
    .seeds    <- .postestSeeds(algorithm, seed)
    algorithm <- .seeds$algorithm
    drawSeed  <- .seeds$drawSeed
    chainSeed <- .seeds$chainSeed

    if (!inherits(control_out, "sienaPostestOutput"))
        stop("'control_out' must be a sienaPostestOutput object, as ",
             "returned by set_postest_output_saom().", call. = FALSE)
    for (.nm in names(control_out)) assign(.nm, control_out[[.nm]])

    # `details` is supplied by NO configuration object: it was withdrawn from
    # both while it is broken (it errors in encodeGroupKeys on every path, a
    # defect that predates this refactor).  Pinned here so the compute
    # functions still receive it, rather than the argument going on accepting
    # a value it cannot honour.  See dev-notes/step5-deliberate-changes.md.
    details <- FALSE

    if (inherits(data, "sienaGroup"))
      stop("marginalEffects does not support multi-group data (sienaGroup).")
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

    # Per-target since the targets object took over the spec, so the guard
    # asks the specs rather than a call-level argument.
    if (dynamic && any(vapply(effectList,
                              function(s) isTRUE(s$returnDecisionDetails),
                              logical(1L))))
      stop("returnDecisionDetails = TRUE is not supported when dynamic = TRUE.")
    # Two resons: dynamic decision details can explode in memory, and are unsafe
    # and batching mode does not support decision details anyway

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
    # (effectName1, diff1, contrast1, effectName2, etc.).
    # ================================================================
    if (!is.list(effectList) || length(effectList) == 0L)
        stop("'effectList' must be a non-empty named list.")
    if (is.null(names(effectList)) || any(nchar(names(effectList)) == 0L))
        stop("All elements of 'effectList' must be named.")

        if (is.null(effects)) effects <- object$requestedEffects

        # ---- Theta / covariance with proper rate handling ----
        .anyRW     <- rateWeight ||
            any(vapply(effectList, function(s) isTRUE(s$rateWeight),
                       logical(1L)))
        .th        <- .postestTheta(object, effects, dynamic, .anyRW)
        thetaHat   <- .th$thetaHat
        covTheta   <- .th$covTheta
        effects    <- .th$effects
        rateParams <- .th$rateParams
        rateIdx    <- .th$rateIdx

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

            if (isTRUE(combineSameLevel)) {
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

                return(if (length(results) == 1L) results[[1L]] else results)
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

            # Per-spec accumulated/rateWeight, override-if-present, exactly
            # as level/condition above.  These used to be OR'd with the
            # model-level value, which let a target switch them ON but never
            # OFF: a target asking for accumulated = FALSE under a model that
            # set TRUE was silently ignored.  Nothing suggests that asymmetry
            # was intended -- the two pairs are the same kind of setting --
            # and it made the model-level value un-overridable in one
            # direction only.  (Step 5d.)
            eff_accumulated <- if ("accumulated" %in% names(spec))
                                   isTRUE(spec$accumulated) else accumulated
            eff_rateWeight  <- if ("rateWeight" %in% names(spec))
                                   isTRUE(spec$rateWeight)  else rateWeight
            if (eff_accumulated && !dynamic)
                stop("Effect '", nm, "': accumulated = TRUE requires ",
                     "dynamic = TRUE.")
            if (eff_rateWeight && dynamic)
                eff_rateWeight <- FALSE   # silently ignored; message already issued

            if (is.null(spec$effectName1))
                stop("Effect '", nm, "' must specify 'effectName1'.")

            effectName1<- resolveEffectName(spec$effectName1,
                               knownEffectNames)
            effectName2<- if (eff_second)
                resolveEffectName(spec$effectName2,
                        knownEffectNames) else NULL

            ## Declared relations name effects by SHORT name; everything
            ## downstream works with resolved names.  Resolved HERE, with the
            ## same resolver the effect names above use, so there is one
            ## notion of what an effect name means rather than two.
            eff_dependencies <- lapply(spec$dependencies, function(d) list(
                target = resolveEffectName(d$target, knownEffectNames),
                terms  = resolveEffectName(d$terms,  knownEffectNames)))

            pt1 <- resolvePerturbType(effectName1, eiTypes,
                                      spec$perturbType1)
            pt2 <- if (eff_second)
                   resolvePerturbType(effectName2, eiTypes,
                                          spec$perturbType2)
                   else
                       "alter"

            # massContrasts is deliberately THREE-STATE (step 5d decided
            # this rather than changing it): NULL means "decide for me" and
            # auto-detects from the perturbation mode, TRUE and FALSE are
            # explicit and both are honoured.  So an explicit FALSE does
            # suppress auto-detection for an ego perturbation -- that is what
            # explicit means, and without it there would be no way to turn
            # the columns off.  This is why the test below is `!is.null` and
            # not `%in% names(spec)` like the other overrides: for them NULL
            # is an override, for this one NULL is a value.
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
                    dependencies   = eff_dependencies,
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
                  effectName2     = effectName2,
                    diff2           = spec$diff2,
                    contrast2       = eff_contrast2,
                    dependencies    = eff_dependencies,
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
                                        else predictSecondDiffJac,
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
            n3PointEst    = n3PointEst,
            n3BatchSize   = if (dynamic) n3BatchSize else NULL,
            useChangeContributions = if (dynamic)
                useChangeContributions else FALSE,
            uncertainty      = uncertainty,
            uncertaintyMode  = uncertaintyMode,
            modeExplicit     = isTRUE(.u$mode_explicit),
            reportConditional = isTRUE(.u$reportConditional),
            nsim             = nsim,
            drawSeed         = drawSeed,
            chainSeed        = chainSeed,
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

    # Assign S3 class on each result; sienaPostestimate returns plain
    # data.frames.  Verified to hold on every exit: single or multiple
    # targets, combineSameLevel on or off, whether the result comes back as a
    # data.frame or as a list of them, every returned frame carries
    # "sienaMarginalEffect" -- combinePostestResults below re-applies it to
    # anything it rebuilds.
    results <- lapply(results, function(r) {
        class(r) <- c("sienaMarginalEffect", class(r))
        r
    })

    if (isTRUE(combineSameLevel)) {
        results <- combinePostestResults(results, builtSpecList, format = format)
        results <- lapply(results, function(r) {
            if (!inherits(r, "sienaMarginalEffect"))
                class(r) <- c("sienaMarginalEffect", class(r))
            r
        })
    }

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
  # How the SE was derived.  The column is called "SE" whatever produced it,
  # so this is where that information lives -- without it the number is not
  # interpretable, since the delta modes differ in what they hold fixed.
  um <- attr(x, "uncertaintyMethod")
  if (!is.null(um)) {
    cat("  SE method: ", switch(um,
        delta      = "delta (conditional)",
        deltaFull  = "deltaFull (conditional + path distribution)",
        bootstrap  = "bootstrap", um), "\n", sep = "")
    if (identical(um, "bootstrap"))
      cat("  nsim:      ", if (!is.null(attr(x, "nsim"))) attr(x, "nsim") else NA, "\n")
  }
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
    effectName, diff, contrast,
    dependencies = NULL,
    details, calcRiskRatio, mainEffect,
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
    ## The raw theta and the declared relations, so the shift-set chain can
    ## build Ddelta -- the same object the Jacobian uses.
    thetaRaw           = theta_use,
    dependencies       = dependencies,
    csMat              = csMat,
    changeStatsMap     = cs$changeStatsMap,
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

  ps <- .perturbationShift(eff_cs_name, pa$diff, contrast, cs$csMat,
                           length(density))
  if (is.null(ps)) return(NULL)
  diff_j <- ps$shift; sign_j <- ps$sign

  # ── the perturbation, as A = d(Du)/d(theta) ───────────────────────────────
  #
  # perturbationJacobian() returns one entry per parameter the perturbation
  # reaches, so a single perturbation may move several: eval always, creation
  # only on rows that create a tie, endowment only on rows that drop one.
  #
  # Still deferred to finite differences:
  #
  #   density        perturbing it changes the ministep's DIRECTION rather than
  #                  a change statistic, so it is not a shift in this
  #                  representation at all and needs its own route.
  #
  # perturbType = "ego" is NOT deferred: firstDiffJacobianEgo() below carries
  # the group sum the ego-wide renormalisation needs.  It is separate from
  # firstDiffJacobian() because using the one-alternative form for an ego
  # perturbation was silently wrong -- the point estimate used the ego update
  # while the Jacobian did not, giving an SE about twice the truth
  # (bootstrap/delta ratio 0.48, against 0.97 for the dyadic case).  That is
  # the regression test to re-run if either form changes.
  if (grepl("density", eff_cs_name, fixed = TRUE)) return(NULL)

  # The shift set: the perturbed effect, plus whatever a declared dependency
  # moves with it.  Built by EVALUATING THE RELATION at the perturbed state,
  # which is what lets a non-multiplicative relation work here without this
  # code changing at all.
  shifts <- .shiftSetFor(setNames(list(diff_j), eff_cs_name),
                         pa$dependencies, cs$csMat)

  A <- perturbationJacobian(shifts    = shifts,
                            direction = density,
                            csMap     = cs$changeStatsMap)
  if (length(A) == 0L) return(NULL)
  delta_u <- utilityShift(A, theta_use, length(density))

  # ── mlogit update Jacobian ────────────────────────────────────────────────
  d_scale <- if (pa$type == "tieProb") density else 1
  ## Same core A, different harness: the ego-wide update renormalises over the
  ## ego's choice set, so its derivative carries a group sum.
  J_diff  <- if (identical(pa$perturbType, "ego"))
      firstDiffJacobianEgo(Jp, changeProb, delta_u, A, cc$group_id)
  else
      firstDiffJacobian(Jp, changeProb, delta_u, A)
  (sign_j * d_scale) * J_diff   # NA sign_j propagates
}


## --------------------------------------------------------------------------
## predictSecondDiffJac — analytical Jacobian paired with predictSecondDiff
##
## A second difference is a weighted sum over four counterfactual cells,
##
##   SD = p_base - p_A - p_B + p_AB          (weights +1, -1, -1, +1)
##
## so its derivative is the same weighted sum of the cells' derivatives.  Two
## facts make that cheap.  First, each cell's probability comes from ONE update
## of the baseline with that cell's accumulated utility shift -- the mlogit
## update composes, so applying effect2 then effect1 equals applying their sum.
## Second, the weights sum to zero, so the baseline term dp_base/dtheta cancels
## and only the SHIFT part of each cell's derivative survives.  That is exactly
## what firstDiffJacobian() returns, so no new derivative is needed:
##
##   J_SD = -J_A - J_B + J_AB,   J_c = firstDiffJacobian(Jp, p, Du_c, A_c)
##
## and the base cell contributes nothing, since a zero shift gives a zero
## Jacobian.
##
## Declines the same cases predictFirstDiffJac does, for the same reasons, on
## either side.  The pairing between an update form and its derivative is NOT
## enforced anywhere -- .resolveBuiltinJac matches a Jacobian to a predictFun by
## identity, not by update form -- and an unnoticed mismatch there is what made
## the ego SE twice the truth.  So this is verified against bootstrap rather
## than assumed correct.
## --------------------------------------------------------------------------
predictSecondDiffJac <- function(cc, theta, changeProb, density, pa, cs, ...) {
  if (grepl("density", pa$effectName1, fixed = TRUE) ||
      grepl("density", pa$effectName2, fixed = TRUE)) return(NULL)

  n         <- length(density)
  theta_use <- theta[cc$effectNames]
  Jp        <- softmax_jac_rcpp(changeProb, cc$contribMat, cc$group_id)

  ## Both sides resolved by the SAME helper the first-difference Jacobian
  ## uses, so a contrast means one thing in one place.
  s1 <- .perturbationShift(pa$effectName1, pa$diff1, pa$contrast1, cs$csMat, n)
  s2 <- .perturbationShift(pa$effectName2, pa$diff2, pa$contrast2, cs$csMat, n)
  if (is.null(s1) || is.null(s2)) return(NULL)

  ## Each cell's shift set, built by the SAME helper the point estimate uses.
  ## The AB cell perturbs both effects at once, so a relation evaluated there
  ## produces the dB*dC cross term from the multiplication -- no separate
  ## rule, and no moderator to keep in step.
  cellShifts <- function(perturbations)
    .shiftSetFor(perturbations, pa$dependencies, cs$csMat)

  A_side <- cellShifts(setNames(list(s1$shift), pa$effectName1))
  B_side <- cellShifts(setNames(list(s2$shift), pa$effectName2))

  ## The AB cell perturbs both effects at once, so evaluating the relation
  ## there yields every term including dB*dC.  Nothing to assemble by hand --
  ## with no relation the two effects simply do not move each other, and the
  ## union of the two sides is the whole story.
  AB <- cellShifts(setNames(list(s1$shift, s2$shift),
                            c(pa$effectName1, pa$effectName2)))

  ## SD = p_base - p_A - p_B + p_AB.  firstDiffJacobian() already subtracts
  ## the baseline, and the weights sum to zero, so the base cell contributes
  ## nothing and only the shift parts survive.
  cells   <- list(A = A_side, B = B_side, AB = AB)
  weights <- c(A = -1, B = -1, AB = 1)
  ## MIXED perturbation types are declined, not approximated.
  ##
  ## A second difference composes two updates: effect2 moves, then effect1
  ## moves from there.  When both are the same kind the composition is again
  ## that kind -- two ego-wide renormalisations over one choice set compose
  ## into one, two one-alternative shifts likewise -- so a single harness is
  ## the right derivative for every cell.
  ##
  ## When they differ it is not.  The B cell moves only the alter-type effect
  ## and needs the one-alternative derivative; the AB cell needs the chain
  ## rule through an alter update FOLLOWED BY an ego renormalisation, which is
  ## neither harness on its own.  Applying the ego harness throughout gave an
  ## SE that did not describe the estimate: against 400 bootstrap draws on the
  ## Glasgow egoX x altX second difference it came out at 0.31-0.70 of the
  ## bootstrap SE, where the finite-difference Jacobian sits at 0.78-1.03.
  ##
  ## Finite differences are correct here, so decline until the composed
  ## derivative exists -- the same choice predictFirstDiffJac() makes for
  ## density.
  if (!identical(pa$perturbType1, pa$perturbType2)) return(NULL)
  egoWide <- identical(pa$perturbType1, "ego")

  J <- NULL
  for (nm in names(cells)) {
    A  <- perturbationJacobian(cells[[nm]], density, cs$changeStatsMap)
    if (length(A) == 0L) return(NULL)
    du <- utilityShift(A, theta_use, n)
    Jc <- if (egoWide) firstDiffJacobianEgo(Jp, changeProb, du, A, cc$group_id)
          else         firstDiffJacobian(Jp, changeProb, du, A)
    J  <- if (is.null(J)) weights[[nm]] * Jc else J + weights[[nm]] * Jc
  }

  d_scale <- if (pa$type == "tieProb") density else 1
  (s1$sign * s2$sign * d_scale) * J
}


predictSecondDiff <- function(changeContributions, theta, type,
    dependencies = NULL,
    effectName1, diff1, contrast1,
    effectName2, diff2, contrast2,
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
    thetaRaw            = theta_use,
    dependencies        = dependencies,
    csMat               = csMat,
    changeStatsMap      = cs$changeStatsMap,
    effectName2         = effectName2,
    effectContribution2 = if (!is.null(effectName2)) csMat[, effectName2] else NULL,
    diff2               = diff2,
    contrast2           = contrast2,
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
                               ## The chain Ddelta -> Du needs the RAW theta
                               ## (one entry per estimated parameter) and the
                               ## declared relations.  `theta` below is
                               ## thetaEff, folded per direction by base
                               ## effect, which is a different indexing.
                               thetaRaw,
                               dependencies = NULL,
                               csMat,
                               changeStatsMap,
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
                                      densityValue = densityValue,
                                      thetaRaw = thetaRaw,
                                      dependencies = dependencies,
                                      csMat = csMat,
                                      changeStatsMap = changeStatsMap)
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
                                effectName2, 
                                effectContribution2,
                                diff2 = NULL,
                                contrast2 = NULL,
                                type = "changeProb",
                                tieProb = NULL,
                                details = FALSE,
                                mainEffect = "riskDifference",
                                calcRiskRatio = FALSE,
                                perturbType1 = "alter",
                                perturbType2 = "alter",
                                ## Inputs for the Ddelta -> Du chain.  Both
                                ## calculateFirstDiff() calls get them, but with
                                ## DIFFERENT states: the first works from the
                                ## base, the second from the state after
                                ## effect2 has moved.
                                thetaRaw = NULL,
                                dependencies = NULL,
                                csMat = NULL,
                                changeStatsMap = NULL,
                                group_id = NULL){
  ## fd1: effect1's difference at the BASE state, so the relation is evaluated
  ## against the unshifted change statistics.
  firstDiff <- calculateFirstDiff(
    densityValue = densityValue,
    changeProb = changeProb,
    changeUtil = changeUtil,
    effectName = effectName1,
    effectContribution = effectContribution1,
    diff = diff1,
    contrast = contrast1,
    thetaRaw = thetaRaw, dependencies = dependencies,
    csMat = csMat, changeStatsMap = changeStatsMap,
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

  ## Step 1 of the second difference: effect2 moves, from the BASE state.
  ##
  ##   step 1 here   delta_A * Dd_B            (effect2 moves, effect1 at base)
  ##   step 2 fd2    Dd_A * (delta_B + Dd_B)   (effect1 moves, effect2 shifted)
  ##   -----------------------------------------------------------------
  ##   base -> AB    delta_A*Dd_B + Dd_A*delta_B + Dd_A*Dd_B
  ##
  ## The cross term falls out of the composition; neither step supplies it
  ## alone, and neither counts it twice.  Both steps go through the same
  ## calculateUtilityDiff(), which is what makes that true.
  utilDiff21 <- calculateUtilityDiff(effectName = effectName2,
                                      diff = diff2,
                                      densityValue = densityValue,
                                      thetaRaw = thetaRaw,
                                      dependencies = dependencies,
                                      csMat = csMat,
                                      changeStatsMap = changeStatsMap)

  changeProb_cf21 <- mlogit_update_r(changeProb, utilDiff21, group_id, perturbType2)
  if (type == "tieProb") {
    tieProb_cf21 <- changeProb_cf21
    idx <- which(!is.na(densityValue) & densityValue == -1)
    if (length(idx) > 0) tieProb_cf21[idx] <- 1 - changeProb_cf21[idx]
  }
  changeUtil21 <- changeUtil + utilDiff21
  rm(utilDiff21)

  ## The state after effect2 has moved.  fd2 below evaluates the relation
  ## HERE rather than at the base, so the dB*dC cross term comes out of the
  ## arithmetic.  NA marks rows outside a contrast's range -- not perturbed,
  ## so not shifted either.
  csMatShifted <- csMat
  if (!is.null(csMat) && !is.null(effectName2) &&
      effectName2 %in% colnames(csMat)) {
    .d2 <- ifelse(is.na(diff2), 0, diff2)
    csMatShifted[, effectName2] <- csMat[, effectName2] + .d2
  }
  rm(diff2)

  # Free intermediates before the second calculateFirstDiff so that
  # the second call's temporaries do not push peak memory over the limit.
  # Cost of gc() is proportional to live objects (~2 sec per 15 GB heap);
  # only worthwhile when live vectors are large enough to matter.
  if (length(densityValue) > 1e7) invisible(gc())

  # Calculate first difference of changing effect1 if effect2 was already changed
  ## fd2: effect1's difference from the state where effect2 has ALREADY moved.
  ## Hence csMatShifted -- the relation must be evaluated at that state, not at
  ## the base, or effect1's dependency uses a stale moderator.
  firstDiff2 <- calculateFirstDiff(
    densityValue = densityValue, # careful if one was density!
    changeProb = changeProb_cf21,
    changeUtil = changeUtil21,
    effectName = effectName1,
    diff = diff1,
    contrast = contrast1,
    effectContribution = effectContribution1,
    type = type,
    tieProb = tieProb_cf21,
    thetaRaw = thetaRaw, dependencies = dependencies,
    csMat = csMatShifted, changeStatsMap = changeStatsMap,
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

# --------------------------------------------------------------------------
# calculateUtilityDiff — one perturbation, as a utility shift
#
# "effectName moves by `diff`; what happens to u?"  Every step of every
# difference asks exactly this, and they must all be answered here: a second
# difference COMPOSES two such steps, so if one step read a declared relation
# and the other did not, a dependent effect would move on one and not the
# other.  That is not hypothetical -- it is what happened while the relation
# lived at one call site and intEffectNames/modEffectNames at the other, and
# it silently dropped half of every covariate-qualified second difference.
#
# The chain, in three named pieces:
#   .shiftSetFor()          which change statistics move, and by how much
#   perturbationJacobian()  that shift as A, in parameter space
#   utilityShift()          Du = A . theta
#
# The SAME A is what the Jacobian consumes.  That is the invariant: point
# estimate and derivative must describe one counterfactual, or the standard
# error belongs to a different quantity than the estimate.
#
# `diff` is in change-statistic space (absolute): +1 = one more unit of the
# feature.  NULL means a unit perturbation.
# --------------------------------------------------------------------------
calculateUtilityDiff <- function(effectName, diff, densityValue,
                                 thetaRaw, dependencies = NULL,
                                 csMat, changeStatsMap){
  if (is.null(diff)) diff <- 1  # NULL means "+1 unit" perturbation
  shifts <- .shiftSetFor(setNames(list(diff), effectName),
                         dependencies, csMat)
  utilityShift(perturbationJacobian(shifts, densityValue, changeStatsMap),
               thetaRaw, length(densityValue))
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

# --------------------------------------------------------------------------
# .perturbationShift — a perturbation spec, as a per-row shift and a sign
#
# `diff` shifts every row by the same amount.  A `contrast` is a SWAP: rows at
# contrast[1] move to contrast[2] and vice versa, rows at neither are NA and
# drop out.  Because a swap moves different rows in opposite directions, the
# reported quantity is put back on one direction by negating the rows that
# moved the other way -- that is the `sign`.
#
# One definition, used by the point estimate's Jacobians for both first and
# second differences.  It previously existed four times over
# (calculateFirstDiff, calculateSecondDiff, and once inside each Jacobian),
# with a comment in one of them reading "mirrors calculateFirstDiff" -- which
# is duplication describing itself.
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# .shiftSetFor — which change statistics a perturbation moves, and by how much
#
# Returns a named list: effect -> per-row shift.  The perturbed effects move by
# their own amount; a declared dependency moves its target by whatever the
# RELATION implies at that state.
#
# Evaluating the relation is the point.  A relation that is not a product --
# sameX ~ egoX == altX, gw(transTrip) -- has no product form to fall back to,
# so reading it here is what removes that ceiling without callers changing.
# --------------------------------------------------------------------------
.shiftSetFor <- function(perturbations, dependencies, csMat) {
  shifts <- perturbations
  if (!length(dependencies)) return(shifts)

  base <- setNames(lapply(colnames(csMat), function(k) csMat[, k]),
                   colnames(csMat))
  ## Names arrive already RESOLVED: the spec builder runs relation terms
  ## through resolveEffectName(), so there is one notion of what an effect
  ## name means.
  for (dep in dependencies) {
    ## A declared relation that cannot be located is an ERROR, not a reason to
    ## compute something else.  This used to return NULL and let the caller
    ## quietly fall back, which is how a dependency silently failing to apply
    ## went unnoticed through a whole meta-analysis: the numbers were wrong,
    ## not absent.
    missing <- setdiff(c(dep$terms, dep$target), names(base))
    if (length(missing))
      stop("Declared dependency '", dep$target, " ~ ",
           paste(dep$terms, collapse = ":"), "' names ",
           if (length(missing) > 1L) "effects that are " else "an effect that is ",
           "not in the model: ", paste(missing, collapse = ", "),
           ".\nDeclare only effects the model estimates, or drop the ",
           "dependency.", call. = FALSE)
    d <- .dependencyShift(dep$terms, base, shifts)
    if (all(d == 0 | is.na(d))) next
    shifts[[dep$target]] <-
      (if (is.null(shifts[[dep$target]])) 0 else shifts[[dep$target]]) + d
  }
  shifts
}


.perturbationShift <- function(effectName, diff, contrast, csMat, n) {
  if (is.null(contrast))
    return(list(shift = if (!is.null(diff)) diff else rep(1, n),
                sign  = rep(1L, n)))

  col <- match(effectName, colnames(csMat))
  if (is.na(col)) return(NULL)
  old <- csMat[, col]
  new <- rep(NA_real_, length(old))
  new[old == contrast[1L]] <- contrast[2L]
  new[old == contrast[2L]] <- contrast[1L]
  sg <- ifelse(!is.na(new) & new == contrast[1L], -1L, 1L)
  sg[is.na(new)] <- NA_integer_
  list(shift = new - old, sign = sg)
}


# Jacobian of the EGO-WIDE update, with respect to theta columns.
#
# The split is between a CORE and a HARNESS.  The core is A = d(Du)/d(theta):
# which parameters a perturbation reaches and on which rows.  That is
# effect-dependent and has nothing to do with ego vs alter.  The harness is
# how A is chained through the update's derivative -- and only THAT differs
# between the two forms.
#
# For the ego-wide update p'_j = w_j / S_i, with w = p exp(Du) and S the ego's
# total,
#
#   dp'_j/dtheta = [ e_j (Jp_j + p_j A_j)
#                    - p'_j sum_h e_h (Jp_h + p_h A_h) ] / S_i
#
# i.e. the same numerator terms as the one-alternative form, summed over the
# ego's alternatives because every one of them sits in the denominator.  The
# first difference subtracts the baseline Jacobian.
firstDiffJacobianEgo <- function(Jp, p, delta_u, A, group_id) {
  e <- exp(ifelse(is.na(delta_u), 0, delta_u))
  w <- p * e
  S <- .contigGroupSum(w, group_id)
  pnew <- w / S

  M <- e * Jp
  for (a in A) {
    v <- ifelse(is.na(a$val), 0, a$val)
    M[, a$col] <- M[, a$col] + e * p * v
  }
  GM <- apply(M, 2L, .contigGroupSum, group_id = group_id)
  if (is.null(dim(GM))) GM <- matrix(GM, nrow = length(p))

  (M - pnew * GM) / S - Jp
}

# Per-row group total for CONTIGUOUS groups, as the C++ kernels assume.
.contigGroupSum <- function(x, group_id) {
  if (is.null(group_id)) return(rep_len(sum(x), length(x)))
  g <- as.integer(group_id)
  ends <- c(which(diff(g) != 0L), length(x))
  rep(diff(c(0, cumsum(x)[ends])), times = diff(c(0L, ends)))
}

# Jacobian of the FIRST DIFFERENCE (p' - p) under the one-alternative update.
#
# Named for what it RETURNS, not for the update it goes through: the trailing
# -1 in `w` subtracts the baseline Jacobian, so this is d(p' - p)/dtheta and
# not dp'/dtheta.  Called mlogitUpdateJacobian() until that ambiguity caused a
# real error while deriving the second-difference version.
#
# `A` is the perturbation in parameter space: a list of list(col, val), one
# entry per parameter the perturbation reaches.  It used to be a single
# (A_col, eff_col) pair, which could not express a perturbation that moves
# more than one parameter -- eval plus creation, say -- nor one whose reach
# varies by row.  Accepting a list is also what lets a declared dependency
# contribute its own columns without this function needing to know that
# dependencies exist.
firstDiffJacobian <- function(Jp, cp, delta_u, A) {
  e_j <- exp(delta_u)
  D_j <- 1 - cp + cp * e_j
  w   <- e_j / D_j^2 - 1
  v   <- cp * e_j * (1 - cp) / D_j^2

  J <- w * Jp
  for (a in A) J[, a$col] <- J[, a$col] + v * a$val
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