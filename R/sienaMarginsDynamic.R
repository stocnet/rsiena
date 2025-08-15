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
    useTieProb = TRUE,
    depvar = NULL,
    second = FALSE,
    level = "none",
    condition = NULL,
    sum.fun = mean,
    na.rm = TRUE,
    n3 = 500,
    useChangeContributions = TRUE,
    uncertainty = TRUE,
    nsim = 100,
    useCluster = FALSE,
    nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batch_dir = "temp",
    prefix = "simBatch_b",
    batch_size = 1,
    combine_batch = TRUE,
    keep_batch = FALSE,
    verbose = TRUE
){
    if(!second){
        diffName <- "firstDiff"
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
            useTieProb = useTieProb,
            depvar = depvar,
            useChangeContributions = FALSE
        )
    }else{
        diffName <- "secondDiff"
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
            useTieProb = useTieProb,
            depvar = depvar,
            useChangeContributions = FALSE
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
        theta_hat = ans[["theta"]],
        cov_theta = ans[["covtheta"]],
        uncertainty = uncertainty,
        nsim = nsim,
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
        useChangeContributions = useChangeContributions
    )
}

predictFirstDiffDynamic <- function(ans, data, theta, effects, algorithm,
    useTieProb = TRUE, depvar = NULL,
    effectName, diff = NULL, contrast = NULL,
    interaction = FALSE,
    int_effectNames = NULL,
    mod_effectNames = NULL,
    n3 = NULL,
    useChangeContributions = FALSE
){
    if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
    include <- effects[["include"]]
    includedEffects <- effects[include, ]
    if(length(theta) > length(includedEffects)) {
        noRateIncluded <- includedEffects[["type"]] != "rate"
        thetaNoRate <- theta[noRateIncluded]
        effectNames  <- includedEffects[["shortName"]][noRateIncluded]
    } else {
        thetaNoRate <- theta
        effectNames  <- includedEffects[["shortName"]]
    }

    df <- getChangeContributionsDynamic(
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

    df <- widenContribution(df)
    df <- addUtilityColumn(df, effectNames, thetaNoRate)
    df <- addProbabilityColumn(df, group_vars=c("chain", "period", "ministep"), useTieProb = useTieProb)

    df[["firstDiff"]] <- calculateFirstDiff(
        densityValue = df[["density"]],
        changeProb = df[["changeProb"]],
        effectName = effectName,
        effectContribution = df[[effectName]],
        diff = diff, contrast = contrast,
        interaction = interaction,
        int_effectNames = int_effectNames,
        mod_effectNames = mod_effectNames,
        modContribution = df[[mod_effectNames]],
        effectNames = effectNames,
        theta = thetaNoRate,
        useTieProb = useTieProb,
        tieProb = df[["tieProb"]]
    )
    df <- subset(df, df[["density"]] != 0)
    df <- conditionalReplace(df, df[["density"]] == -1, setdiff(effectNames, "density"), function(x) x * -1)
    df
}

predictSecondDiffDynamic <- function(ans, data, theta, effects, algorithm,
    useTieProb = TRUE, depvar = NULL,
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
    useChangeContributions = FALSE
){
    if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
    include <- effects[["include"]]
    includedEffects <- effects[include, ]
    if(length(theta) > length(includedEffects)) {
        noRateIncluded <- includedEffects[["type"]] != "rate"
        thetaNoRate <- theta[noRateIncluded]
        effectNames  <- includedEffects[["shortName"]][noRateIncluded]
    } else {
        thetaNoRate <- theta
        effectNames  <- includedEffects[["shortName"]]
    }

    df <- getChangeContributionsDynamic(
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
    df <- widenContribution(df)
    df <- addUtilityColumn(df, effectNames, thetaNoRate)
    df <- addProbabilityColumn(df, group_vars = c("chain", "period", "ministep"), useTieProb = useTieProb)

    df <- subset(df, df[["density"]] != 0)

    df[["secondDiff"]] <- calculateSecondDiff(
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
        useTieProb = useTieProb,
        tieProb = df[["tieProb"]]# if not tieProb NULL?
    )
  df <- conditionalReplace(df, df[["density"]] == -1, setdiff(effectNames, "density"), function(x) x * -1)
  df
}
