##@sienaAMEDynamic
sienaAMEDynamic <- function(
    ans,
    data,
    effectName1,
    diff1 = NULL,
    contrast1 = NULL,
    interaction1 = FALSE,
    int_effectNames1 = NULL,
    mod_effectNames1 = NULL,
    effectName2 = NULL,
    diff2 = NULL,
    contrast2 = NULL,
    interaction2 = FALSE,
    int_effectNames2 = NULL,
    mod_effectNames2 = NULL,
    effects,
    algorithm,
    type = c("changeProb", "tieProb"),
    depvar = NULL,
    second = FALSE,
    level = "none",
    condition = NULL,
    sum.fun = mean,
    na.rm = TRUE,
    n3 = 500,
    useChangeContributions = FALSE, # should only be used for static case
    uncertainty = TRUE,
    uncertainty_mode = c("batch", "stream"),
    nsim = 100,
    uncertainty_sd = TRUE,
    uncertainty_ci = TRUE,
    uncertainty_probs = c(0.025, 0.5, 0.975),
    uncertainty_mcse = FALSE,
    uncertainty_mcse_batches = NULL,
    useCluster = FALSE,
    nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batch_dir = "temp",
    prefix = "simBatch_b",
    batch_size = NULL,
    combine_batch = TRUE,
    keep_batch = FALSE,
    verbose = TRUE,
    mainEffect = "riskDifference",
    details = FALSE,
    memory_scale = NULL,
    batch_unit_budget = 5e6,
    dynamic_ministep_factor = 10
){
    uncertainty_mode <- match.arg(uncertainty_mode)
    type <- match.arg(type)
    if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]

    if (is.null(memory_scale)) {
        memory_scale <- compute_memory_scale(
            data = data,
            depvar = depvar,
            dynamic = TRUE,
            n3 = n3
        )
    }
        if (is.null(batch_size)) {
            auto_batch <- auto_batch_from_budget(
                data = data,
                depvar = depvar,
                nsim = nsim,
                nbrNodes = nbrNodes,
                useCluster = useCluster,
                dynamic = TRUE,
                n3 = n3,
                unit_budget = batch_unit_budget,
                dynamic_ministep_factor = dynamic_ministep_factor
            )
            batch_size <- auto_batch$batch_size
        }
    if(!second){
        diffName <- ifelse(mainEffect == "riskDifference", 
            "firstDiff", 
            "firstRiskRatio")
        diffFun <- predictFirstDiffDynamic
        predictArgs <- list(
            ans = ans,
            data = data,
            effects = effects,
            algorithm = algorithm,
            effectName = effectName1,
            diff = diff1,
            contrast = contrast1,
            interaction = interaction1,
            int_effectNames = int_effectNames1,
            mod_effectNames = mod_effectNames1,
            type = type,
            depvar = depvar,
            n3 = n3,
            useChangeContributions = useChangeContributions,
            mainEffect = mainEffect,
            details = details
        )
    }else{
        diffName <- ifelse(mainEffect == "riskDifference", 
            "secondDiff", 
            "secondRiskRatio")
        diffFun <- predictSecondDiffDynamic
        predictArgs <- list(
            ans = ans,
            data = data,
            effects = effects,
            algorithm = algorithm,
            effectName1 = effectName1,
            diff1 = diff1,
            contrast1 = contrast1,
            interaction1 = interaction1,
            int_effectNames1 = int_effectNames1,
            mod_effectNames1 = mod_effectNames1,
            effectName2 = effectName2,
            diff2 = diff2,
            contrast2 = contrast2,
            interaction2 = interaction2,
            int_effectNames2 = int_effectNames2,
            mod_effectNames2 = mod_effectNames2,
            type = type,
            depvar = depvar,
            n3 = n3,
            useChangeContributions = useChangeContributions,
            mainEffect = mainEffect,
            details = details
        )
        ME <- "secondDiff"
    }
    
    sienaPostestimate(
        predictFun = diffFun,
        predictArgs = predictArgs,
        outcome = diffName,
        level = level,
        condition = condition,
        sum.fun = sum.fun,
        na.rm = na.rm,
        theta_hat = ans$theta,
        cov_theta = ans$covtheta,
        uncertainty = uncertainty,
        uncertainty_mode = uncertainty_mode,
        nsim = nsim,
        uncertainty_sd = uncertainty_sd,
        uncertainty_ci = uncertainty_ci,
        uncertainty_probs = uncertainty_probs,
        uncertainty_mcse = uncertainty_mcse,
        uncertainty_mcse_batches = uncertainty_mcse_batches,
        useCluster = useCluster,
        nbrNodes = nbrNodes,
        clusterType = clusterType,
        cluster = cluster,
        batch_dir = batch_dir,
        prefix = prefix,
        combine_batch = combine_batch,
        batch_size = batch_size,
        keep_batch = keep_batch,
        verbose = verbose,
        useChangeContributions = useChangeContributions,
        memory_scale = memory_scale
    )
}

predictFirstDiffDynamic <- function(ans, data, theta, effects, algorithm,
    type = "changeProb", depvar = NULL,
    effectName, diff = NULL, contrast = NULL,
    interaction = FALSE,
    int_effectNames = NULL,
    mod_effectNames = NULL,
    n3 = NULL,
    useChangeContributions = FALSE,
    details = FALSE,
    calcRiskRatio = FALSE,
    mainEffect = "riskDifference"
){
    if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
    include <- effects[["include"]]
    includedEffects <- effects[include, ]
    noRateIncluded <- includedEffects[["type"]] != "rate"
    effectNames  <- includedEffects[["shortName"]][noRateIncluded]

    # Align theta by name
    if (is.null(names(theta))) {
        # Get all effect names (including rate effects)
        allEffectNames <- ans$effects$shortName
        names(theta) <- allEffectNames
    }
    thetaNoRate <- theta[effectNames]

    df <- getDynamicChangeContributions(
        ans = ans,
        data = data,
        theta = theta,
        algorithm = algorithm,
        effects = effects,
        depvar = depvar,
        n3 = n3,
        useChangeContributions = useChangeContributions,
        returnDataFrame = TRUE
    )

    df <- widenDynamicContribution(df)
    df <- addUtilityColumn(df, effectNames, thetaNoRate)
    df <- addProbabilityColumn(df, 
        group_vars = c("chain", "period", "ministep"), 
        type = type
    )

    if (requireNamespace("data.table", quietly = TRUE) && 
        data.table::is.data.table(df)) {
        df <- df[df[["density"]] != 0]
    } else {
        df <- df[df[["density"]] != 0, , drop = FALSE]
    }

    fd <- calculateFirstDiff(
        densityValue = df[["density"]],
        changeProb = df[["changeProb"]],
        changeUtil = df[["changeUtil"]],
        effectName = effectName,
        effectContribution = df[[effectName]],
        diff = diff, contrast = contrast,
        interaction = interaction,
        int_effectNames = int_effectNames,
        mod_effectNames = mod_effectNames,
        modContribution = df[[mod_effectNames]],
        effectNames = effectNames,
        theta = thetaNoRate,
        type = type,
        tieProb = df[["tieProb"]],
        details = details,
        calcRiskRatio = calcRiskRatio,
        mainEffect = mainEffect
    )
    if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
        df[, (names(fd)) := fd]
    } else {
        df[names(fd)] <- fd
    }
    rm(fd)

    df <- conditionalReplace(df, df[["density"]] == -1, setdiff(effectNames, "density"), function(x) x * -1)
    df
}

predictSecondDiffDynamic <- function(ans, data, theta, effects, algorithm,
    type = "changeProb", depvar = NULL,
    effectName1, 
    diff1 = NULL, contrast1 = NULL,
    interaction1 = FALSE,
    int_effectNames1 = NULL,
    mod_effectNames1 = NULL,
    effectName2,
    diff2 = NULL, contrast2 = NULL,
    interaction2 = FALSE,
    int_effectNames2 = NULL,
    mod_effectNames2 = NULL,
    n3 = NULL,
    useChangeContributions = FALSE,
    calcRiskRatio = FALSE,
    mainEffect = "riskDifference",
    details = FALSE
){
    if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
    include <- effects[["include"]]
    includedEffects <- effects[include, ]
    noRateIncluded <- includedEffects[["type"]] != "rate"
    effectNames  <- includedEffects[["shortName"]][noRateIncluded]

    if (is.null(names(theta))) {
        # Get all effect names (including rate effects)
        allEffectNames <- ans$effects$shortName
        names(theta) <- allEffectNames
    }
    thetaNoRate <- theta[effectNames]

    df <- getDynamicChangeContributions(
        ans = ans,
        data = data,
        theta = theta,
        algorithm = algorithm,
        effects = effects,
        depvar = depvar,
        n3 = n3,
        useChangeContributions = useChangeContributions,
        returnDataFrame = TRUE
    )  
    df <- widenDynamicContribution(df)
    df <- addUtilityColumn(df, effectNames, thetaNoRate)
    df <- addProbabilityColumn(df, 
        group_vars = c("chain", "period", "ministep"), type = type)
    if (requireNamespace("data.table", quietly = TRUE) && 
        data.table::is.data.table(df)) {
        df <- df[df[["density"]] != 0]
    } else {
        df <- df[df[["density"]] != 0, , drop = FALSE]
    }
    sd <- calculateSecondDiff(
        densityValue = df[["density"]],
        changeProb = df[["changeProb"]],
        changeUtil = df[["changeUtil"]],
        effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
        effectContribution1 = df[[effectName1]],
        interaction1 = interaction1,
        int_effectNames1 = int_effectNames1,
        mod_effectNames1 = mod_effectNames1,
        modContribution1 = df[[mod_effectNames1]],
        effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
        effectContribution2 = df[[effectName2]],
        interaction2 = interaction2,
        int_effectNames2 = int_effectNames2,
        mod_effectNames2 = mod_effectNames2,
        modContribution2 = df[[mod_effectNames2]],
        effectNames = effectNames,
        theta = thetaNoRate,
        type = type,
        tieProb = df[["tieProb"]],
        details = details,
        calcRiskRatio = calcRiskRatio,
        mainEffect = mainEffect
    )

    if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
        df[, (names(sd)) := sd]
    } else {
        df[names(sd)] <- sd
    }
    rm(sd)

    df <- conditionalReplace(df, df[["density"]] == -1, 
        setdiff(effectNames, "density"), function(x) x * -1)
    df
}
