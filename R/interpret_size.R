##@interpret_size


interpret_size <- function(object, ...) UseMethod("interpret_size", object)

##@interpret_size  method for sienaEffects objects
# Wrapper: delegates RI computation to relativeImportance.sienaEffects.
# interpret_size may later attach additional non-RI measures (entropy, etc.).
interpret_size.sienaEffects <- function(object, data,
		theta, getChangeStats = FALSE, ...)
{
    relativeImportance(object, data = data, theta = theta,
        getChangeStats = getChangeStats, ...)
}

##@interpret_size  method for sienaFit objects
# Orchestrator for effect-size measures at observation moments.
# Calls relativeImportance.sienaFit for the RI computation
# and optionally attaches change statistics.
interpret_size.sienaFit <- function(object, data, getChangeStats=FALSE,
    depvar = NULL,
    effects = NULL,
    dynamic = FALSE,
    distFun = c("L1D", "KLD"),
    algorithm = NULL,
    n3 = 500,
    useChangeContributions = FALSE,
    uncertainty = FALSE,
    level = "ego", # should be egoChain for dynamic RI? check if we can auto-translate ego -> egoChain when dynamic=TRUE
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
    ...) {

    # what happens for multiple depvars?  Should we return a list of sienaRI 
    # objects, or a single sienaRI object with an added "dependentVariable" column?  
    # Delegate to relativeImportance — which does validation, contribution
    # extraction, RI computation, and optional uncertainty
    result <- relativeImportance(object, data = data,
        depvar           = depvar,
        effects          = effects,
        dynamic          = dynamic,
        distFun          = distFun,
        algorithm        = algorithm,
        n3               = n3,
        useChangeContributions = useChangeContributions,
        uncertainty      = uncertainty,
        level            = level,
        sum_fun          = sum_fun,
        na.rm            = na.rm,
        nsim             = nsim,
        uncertaintyMode  = uncertaintyMode,
        uncertaintySd    = uncertaintySd,
        uncertaintyCi    = uncertaintyCi,
        uncertaintyProbs = uncertaintyProbs,
        useCluster       = useCluster,
        nbrNodes         = nbrNodes,
        clusterType      = clusterType,
        batchDir         = batchDir,
        batchSize        = batchSize,
        keepBatch        = keepBatch,
        verbose          = verbose,
        memoryScale      = memoryScale,
        batchUnitBudget  = batchUnitBudget,
        dynamicMinistepFactor = dynamicMinistepFactor,
        getChangeStats   = getChangeStats,
        ...)

    result
}

# ============================================================================
# Dynamic RI (formerly in sienaRIDynamics.r)
# ============================================================================

##@interpret_size_dynamics  Deprecated — use interpret_size(dynamic=TRUE)
interpret_size_dynamics <- function(data, ans=NULL, theta=NULL,
				algorithm=NULL, effects=NULL, depvar=NULL, intervalsPerPeriod=NULL,
				n3 = NULL, useChangeContributions = FALSE,
				silent = TRUE, seed = NULL,
				uncertainty = FALSE, ...)
{
	.Deprecated("interpret_size",
		msg = paste("interpret_size_dynamics() is deprecated.",
			"Use interpret_size(object, data, dynamic=TRUE, algorithm=...) instead."))
	if (is.null(ans))
		stop("interpret_size_dynamics() now requires a sienaFit object 'ans'.")
	if (!inherits(ans, "sienaFit"))
		stop("'ans' is not a legitimate Siena fit object")
	if (is.null(algorithm)) algorithm <- ans$x
	interpret_size(ans, data = data, depvar = depvar,
		effects = effects, dynamic = TRUE, algorithm = algorithm,
		n3 = n3, useChangeContributions = useChangeContributions,
		uncertainty = uncertainty, ...)
}
