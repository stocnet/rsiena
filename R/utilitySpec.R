## Spec-based utility framework for structured marginal effects.
##
## Works with **change statistics** in attachContributions(flip=TRUE) style:
## the density column is ±1 (the sign), all other columns are unsigned
## change statistics δ_k = |s_k(x+ij) - s_k(x-ij)|.
##
## The utility for dyad ij is:
##
##   U_ij = δ_density × (θ_density + Σ_{k≠density} δ_k × θ_k)
##
## or equivalently: U_ij = contribMat %*% θ  (since contrib_k = δ_density × δ_k).
##
## Density is implicit — never declared in the spec, always interacts with
## everything.  Perturbation of density: utilDiff = −2 × U_ij.
## Perturbation of effect k by Δ:
##   utilDiff = δ_density × Δ × (θ_k + Σ_int Π_{δ_other} × θ_int)
##
## State updates are trivial: δ[k] += Δ.  No density multiplication needed,
## eliminating the moderator-update bug class.
##
## A utilitySpec declares interaction decompositions using R formula syntax:
##
##   spec <- utilitySpec(transRecTrip ~ factor(recip) : transTrip)
##
## This maps δ_recip × δ_transTrip → the SIENA effect transRecTrip.
## factor() marks a binary factor (default contrast c(0,1)).
## Bare names are continuous (default perturbation diff = +1).

# ---- Constructor -------------------------------------------------------

##@utilitySpec Constructor
utilitySpec <- function(..., contrasts = list(), diffs = list()) {
  formulas <- list(...)
  if (length(formulas) == 1L && is.list(formulas[[1L]]) &&
      !inherits(formulas[[1L]], "formula"))
    formulas <- formulas[[1L]]

  interactions   <- list()
  factorEffects  <- list()

  for (f in formulas) {
    stopifnot(inherits(f, "formula"), length(f) == 3L)
    intName  <- deparse(f[[2]])
    rhsLabel <- attr(terms(f), "term.labels")
    if (length(rhsLabel) != 1L)
      stop("RHS should be a single interaction term (use ':' to separate factors).")
    parsed <- .parseInteractionLabel(rhsLabel)
    # density is always implicit — warn if user tried to include it
    denIdx <- grep("density", parsed$names, ignore.case = TRUE)
    if (length(denIdx)) {
      message("Note: density is implicit and cannot appear in the spec. Ignoring.")
      parsed$names    <- parsed$names[-denIdx]
      parsed$isFactor <- parsed$isFactor[-denIdx]
    }
    interactions[[intName]] <- list(
      intEffectName = intName,
      factors       = parsed$names,
      isFactor      = parsed$isFactor
    )
    for (i in seq_along(parsed$names)) {
      nm <- parsed$names[i]
      if (parsed$isFactor[i] && !(nm %in% names(factorEffects)))
        factorEffects[[nm]] <- c(0L, 1L)
    }
  }

  for (nm in names(contrasts)) {
    if (grepl("density", nm, ignore.case = TRUE)) {
      message("Note: density contrast is always c(-1, 1). Ignoring user override.")
      next
    }
    factorEffects[[nm]] <- contrasts[[nm]]
  }

  perturb <- list()
  for (nm in names(factorEffects))
    perturb[[nm]] <- list(type = "contrast", levels = factorEffects[[nm]])
  for (nm in names(diffs))
    perturb[[nm]] <- list(type = "diff", value = diffs[[nm]])

  allFactors <- unique(unlist(lapply(interactions, `[[`, "factors")))
  factorIndex <- setNames(lapply(allFactors, function(f) {
    names(Filter(function(t) f %in% t$factors, interactions))
  }), allFactors)

  structure(
    list(interactions = interactions,
         perturb      = perturb,
         factorIndex  = factorIndex),
    class = "utilitySpec"
  )
}

# ---- Internal: parse interaction label ---------------------------------

.parseInteractionLabel <- function(label) {
  parts <- strsplit(label, ":", fixed = TRUE)[[1]]
  n     <- length(parts)
  nms   <- character(n)
  isF   <- logical(n)
  for (i in seq_len(n)) {
    p <- trimws(parts[i])
    m <- regexec("^factor\\(\\s*([^,)]+?)\\s*\\)$", p)[[1]]
    if (m[1] > 0) {
      nms[i] <- trimws(regmatches(p, list(m))[[1]][2])
      isF[i] <- TRUE
    } else {
      nms[i] <- p
      isF[i] <- FALSE
    }
  }
  list(names = nms, isFactor = isF)
}

# ---- Name resolution ---------------------------------------------------

##@resolveSpec Resolve short effect names to change-stat column names
resolveSpec <- function(spec, effectNames) {
  resolve1 <- function(nm) resolveEffectName(nm, effectNames)

  ints <- lapply(spec$interactions, function(t) {
    list(intEffectName = resolve1(t$intEffectName),
         factors       = vapply(t$factors, resolve1, character(1),
                                USE.NAMES = FALSE),
         isFactor      = t$isFactor)
  })
  names(ints) <- vapply(ints, `[[`, character(1), "intEffectName")

  allFactors <- unique(unlist(lapply(ints, `[[`, "factors")))
  fIdx <- setNames(lapply(allFactors, function(f) {
    names(Filter(function(t) f %in% t$factors, ints))
  }), allFactors)

  prb <- list()
  for (nm in names(spec$perturb))
    prb[[resolve1(nm)]] <- spec$perturb[[nm]]

  structure(
    list(interactions = ints, perturb = prb, factorIndex = fIdx),
    class = "resolvedSpec"
  )
}

# ---- Conversion helpers ------------------------------------------------

# contribToChangeStats is now defined in postestimate.R

# ---- Core: utility difference ------------------------------------------

## Utility difference when perturbing effect k by Δ (non-density only).
## Density perturbation is -2*U, handled directly in calculateFirstDiffSpec.
##
##   utilDiff = δ_density × Δ × (θ_k + Σ_int Π_{δ_other} × θ_int)
##
## where δ_other are change-statistic columns (unsigned).
##
##@calculateUtilityDiffSpec Utility difference (spec-based)
calculateUtilityDiffSpec <- function(rspec, effectName, delta,
                                     theta, density, changeStats,
                                     effectNames) {
  effectNum <- which(effectNames == effectName)
  utilDiff  <- density * delta * theta[effectNum]
  intNames  <- rspec$factorIndex[[effectName]]
  for (intName in intNames) {
    term   <- rspec$interactions[[intName]]
    intNum <- which(effectNames == term$intEffectName)
    others <- setdiff(term$factors, effectName)
    modProduct <- if (length(others) == 1L) {
      changeStats[, others]
    } else {
      Reduce(`*`, lapply(others, function(f) changeStats[, f]))
    }
    utilDiff <- utilDiff + density * modProduct * delta * theta[intNum]
  }
  utilDiff
}

# ---- Perturbation resolution -------------------------------------------

## Determine delta for effectName from spec + optional overrides.
## Returns list(delta, contrast, oldCS, newCS).
##@.resolvePerturbation Internal perturbation resolver
.resolvePerturbation <- function(rspec, effectName,
                                 diff, contrast, changeStats) {
  if (!is.null(contrast))
    return(.computeContrastDiff(changeStats[, effectName], contrast))
  if (!is.null(diff))
    return(list(delta = diff, contrast = NULL, oldCS = NULL, newCS = NULL))

  prb <- rspec$perturb[[effectName]]
  if (!is.null(prb) && prb$type == "contrast")
    return(.computeContrastDiff(changeStats[, effectName], prb$levels))

  d <- if (!is.null(prb) && prb$type == "diff") prb$value else 1
  list(delta = d, contrast = NULL, oldCS = NULL, newCS = NULL)
}

## Convert contrast levels to a delta vector.
##@.computeContrastDiff Internal contrast helper
.computeContrastDiff <- function(changeStatistic, levels) {
  old <- changeStatistic
  new <- rep(NA_real_, length(old))
  new[old == levels[1]] <- levels[2]
  new[old == levels[2]] <- levels[1]
  list(delta    = new - old,
       contrast = levels,
       oldCS    = old,
       newCS    = new)
}

# ---- First difference --------------------------------------------------

##@calculateFirstDiffSpec First difference (spec-based)
calculateFirstDiffSpec <- function(rspec, effectName,
                                   density, changeProb, changeUtil,
                                   changeStats, theta, effectNames,
                                   diff = NULL, contrast = NULL,
                                   type = "changeProb", tieProb = NULL,
                                   mainEffect = "riskDifference",
                                   details = FALSE, calcRiskRatio = FALSE,
                                   perturbType = "alter", group_id = NULL) {

  isDensity <- grepl("density", effectName, fixed = TRUE)

  if (isDensity) {
    # Density perturbation: flip ±1, utilDiff = -2U
    ct <- if (!is.null(contrast)) contrast else c(-1L, 1L)
    newDensity <- ifelse(density == ct[1], ct[2],
                         ifelse(density == ct[2], ct[1], NA_real_))
    utilDiff <- ifelse(is.na(newDensity), NA_real_, -2 * changeUtil)
    contrastInfo <- list(contrast = ct, newCS = newDensity)
  } else {
    # Regular perturbation via spec
    pert <- .resolvePerturbation(rspec, effectName, diff, contrast,
                                 changeStats)
    utilDiff <- calculateUtilityDiffSpec(rspec, effectName, pert$delta,
                                         theta, density, changeStats,
                                         effectNames)
    contrastInfo <- if (!is.null(pert$contrast))
      list(contrast = pert$contrast, newCS = pert$newCS) else NULL
  }

  # Counterfactual probability
  changeProbCf <- mlogit_update_r(changeProb, utilDiff, group_id, perturbType)
  changeProbCf[density == 0] <- NA

  if (type == "tieProb") {
    tieProbCf <- changeProbCf
    idx <- which(!is.na(density) & density == -1)
    if (length(idx)) tieProbCf[idx] <- 1 - changeProbCf[idx]
    firstDiff <- tieProbCf - tieProb
  } else {
    firstDiff <- changeProbCf - changeProb
  }

  # Contrast flip
  if (!is.null(contrastInfo)) {
    flipIdx <- which(contrastInfo$newCS == min(contrastInfo$contrast))
    firstDiff[flipIdx] <- -firstDiff[flipIdx]
  }

  # Risk ratio
  if (calcRiskRatio || mainEffect == "riskRatio") {
    firstRiskRatio <- if (type == "tieProb") tieProbCf / tieProb
                      else changeProbCf / changeProb
    if (!is.null(contrastInfo)) {
      flipIdx <- which(contrastInfo$newCS == min(contrastInfo$contrast))
      firstRiskRatio[flipIdx] <- 1 / firstRiskRatio[flipIdx]
    }
  }

  # Output
  if (details) {
    out <- data.frame(firstDiff = firstDiff, utilDiff = utilDiff,
                      newChangeProb = changeProbCf, oldChangeProb = changeProb)
    if (type == "tieProb") {
      out$newTieProb <- tieProbCf
      out$oldTieProb <- tieProb
    }
    if (calcRiskRatio || mainEffect == "riskRatio")
      out$firstRiskRatio <- firstRiskRatio
    out
  } else if (mainEffect == "riskRatio") {
    list(firstRiskRatio = firstRiskRatio)
  } else {
    list(firstDiff = firstDiff)
  }
}

# ---- Second difference -------------------------------------------------

##@calculateSecondDiffSpec Second difference (spec-based)
calculateSecondDiffSpec <- function(rspec, effectName1, effectName2,
                                    density, changeProb, changeUtil,
                                    changeStats, theta, effectNames,
                                    diff1 = NULL, contrast1 = NULL,
                                    diff2 = NULL, contrast2 = NULL,
                                    type = "changeProb", tieProb = NULL,
                                    mainEffect = "riskDifference",
                                    details = FALSE, calcRiskRatio = FALSE,
                                    perturbType1 = "alter",
                                    perturbType2 = "alter",
                                    group_id = NULL) {

  if (grepl("density", effectName2, fixed = TRUE))
    stop("density as the moderator (effectName2) is not supported")

  # Step 1: firstDiff of effectName1 at baseline
  fd1 <- calculateFirstDiffSpec(rspec, effectName1,
      density, changeProb, changeUtil, changeStats, theta, effectNames,
      diff = diff1, contrast = contrast1,
      type = type, tieProb = tieProb,
      mainEffect = mainEffect, details = details,
      calcRiskRatio = calcRiskRatio,
      perturbType = perturbType1, group_id = group_id)

  # Step 2: perturb effectName2
  pert2 <- .resolvePerturbation(rspec, effectName2, diff2, contrast2,
                                changeStats)
  contrastInfo2 <- if (!is.null(pert2$contrast))
    list(contrast = pert2$contrast, newCS = pert2$newCS) else NULL
  utilDiff2 <- calculateUtilityDiffSpec(rspec, effectName2, pert2$delta,
                                         theta, density, changeStats,
                                         effectNames)
  changeProbCf2 <- mlogit_update_r(changeProb, utilDiff2, group_id,
                                    perturbType2)
  changeUtilCf2 <- changeUtil + utilDiff2

  tieProbCf2 <- NULL
  if (type == "tieProb") {
    tieProbCf2 <- changeProbCf2
    idx <- which(!is.na(density) & density == -1)
    if (length(idx)) tieProbCf2[idx] <- 1 - changeProbCf2[idx]
  }

  # Step 3: update state — only the perturbed column
  delta2 <- pert2$delta
  delta2[is.na(delta2)] <- 0
  changeStats2 <- changeStats
  changeStats2[, effectName2] <- changeStats[, effectName2] + delta2

  # Step 4: firstDiff of effectName1 at shifted state
  fd2 <- calculateFirstDiffSpec(rspec, effectName1,
      density, changeProbCf2, changeUtilCf2, changeStats2, theta, effectNames,
      diff = diff1, contrast = contrast1,
      type = type, tieProb = tieProbCf2,
      mainEffect = mainEffect, details = details,
      calcRiskRatio = calcRiskRatio,
      perturbType = perturbType1, group_id = group_id)

  # Step 5: second difference
  secondDiff <- fd2[["firstDiff"]] - fd1[["firstDiff"]]
  if (!is.null(contrastInfo2)) {
    flipIdx <- which(contrastInfo2$newCS == min(contrastInfo2$contrast))
    secondDiff[flipIdx] <- -secondDiff[flipIdx]
  }

  if (mainEffect == "riskRatio") {
    secondRiskRatio <- fd2[["firstRiskRatio"]] / fd1[["firstRiskRatio"]]
    if (!is.null(contrastInfo2)) {
      flipIdx <- which(contrastInfo2$newCS == min(contrastInfo2$contrast))
      secondRiskRatio[flipIdx] <- 1 / secondRiskRatio[flipIdx]
    }
  }

  if (details) {
    out <- data.frame(
      changeProb_base = changeProb,
      changeProb_main = fd1[["newChangeProb"]],
      changeProb_mod  = changeProbCf2,
      changeProb_both = fd2[["newChangeProb"]],
      firstDiff1      = fd1[["firstDiff"]],
      firstDiff2      = fd2[["firstDiff"]],
      secondDiff      = secondDiff
    )
    if (type == "tieProb") {
      out$tieProb_base <- tieProb
      out$tieProb_main <- fd1[["newTieProb"]]
      out$tieProb_mod  <- tieProbCf2
      out$tieProb_both <- fd2[["newTieProb"]]
    }
    if (mainEffect == "riskRatio" || calcRiskRatio) {
      out$firstRiskRatio1 <- fd1[["firstRiskRatio"]]
      out$firstRiskRatio2 <- fd2[["firstRiskRatio"]]
      out$secondRiskRatio <- secondRiskRatio
    }
    out
  } else if (mainEffect == "riskRatio") {
    list(secondRiskRatio = secondRiskRatio)
  } else {
    list(secondDiff = secondDiff)
  }
}

# ---- Predict wrappers --------------------------------------------------

##@predictFirstDiffStaticSpec Thin wrapper for spec-based first difference
predictFirstDiffStaticSpec <- function(theta, staticContributions,
    type = "changeProb", rspec,
    effectName, diff = NULL, contrast = NULL,
    details = FALSE, calcRiskRatio = FALSE,
    mainEffect = "riskDifference", perturbType = "alter",
    massContrasts = FALSE)
{
  effectNames <- staticContributions$effectNames
  thetaUse    <- theta[effectNames]
  contribMat  <- staticContributions$contribMat
  densityName <- resolveEffectName("density", effectNames)
  density     <- contribMat[, densityName]

  changeStats <- contribToChangeStats(contribMat, effectNames)
  utility     <- calculateUtility(contribMat, thetaUse)
  changeProb  <- calculateChangeProb(utility, staticContributions$group_id)
  tieProb     <- if (type == "tieProb")
    calculateTieProb(changeProb, density) else NULL

  fd <- calculateFirstDiffSpec(rspec, effectName,
      density, changeProb, utility, changeStats, thetaUse, effectNames,
      diff = diff, contrast = contrast,
      type = type, tieProb = tieProb,
      mainEffect = mainEffect, details = details,
      calcRiskRatio = calcRiskRatio, perturbType = perturbType,
      group_id = staticContributions$group_id)

  keep <- density != 0L
  out  <- data.frame(groupColsList(staticContributions, keep),
                     stringsAsFactors = FALSE)
  out  <- attachContributions(out, effectNames,
                              contribMat[keep, , drop = FALSE], flip = TRUE)
  out[names(fd)] <- lapply(fd, `[`, keep)

  if (details) {
    out$changeUtil <- utility[keep]
    out$changeProb <- changeProb[keep]
    if (!is.null(tieProb)) out$tieProb <- tieProb[keep]
  }

  if (massContrasts) {
    diffColName <- intersect(c("firstDiff", "firstRiskRatio"), names(fd))[1L]
    if (!is.na(diffColName) && diffColName == "firstDiff") {
      mc <- computeMassContrasts(
        firstDiff = fd[["firstDiff"]][keep],
        density   = density[keep],
        ego       = staticContributions$ego[keep],
        period    = staticContributions$period[keep],
        group     = if (!is.null(staticContributions$group))
                      staticContributions$group[keep] else rep(1L, sum(keep)),
        type      = type)
      out$massCreation    <- mc$massCreation
      out$massDissolution <- mc$massDissolution
    }
  }
  out
}

##@predictSecondDiffStaticSpec Thin wrapper for spec-based second difference
predictSecondDiffStaticSpec <- function(theta, staticContributions,
    type = "changeProb", rspec,
    effectName1, diff1 = NULL, contrast1 = NULL,
    effectName2,  diff2 = NULL, contrast2 = NULL,
    mainEffect = "riskDifference", details = FALSE,
    perturbType1 = "alter", perturbType2 = "alter",
    massContrasts = FALSE)
{
  effectNames <- staticContributions$effectNames
  thetaUse    <- theta[effectNames]
  contribMat  <- staticContributions$contribMat
  densityName <- resolveEffectName("density", effectNames)
  density     <- contribMat[, densityName]

  changeStats <- contribToChangeStats(contribMat, effectNames)
  utility     <- calculateUtility(contribMat, thetaUse)
  changeProb  <- calculateChangeProb(utility, staticContributions$group_id)
  tieProb     <- if (type == "tieProb")
    calculateTieProb(changeProb, density) else NULL

  sd <- calculateSecondDiffSpec(rspec, effectName1, effectName2,
      density, changeProb, utility, changeStats, thetaUse, effectNames,
      diff1 = diff1, contrast1 = contrast1,
      diff2 = diff2, contrast2 = contrast2,
      type = type, tieProb = tieProb,
      mainEffect = mainEffect, details = details,
      perturbType1 = perturbType1, perturbType2 = perturbType2,
      group_id = staticContributions$group_id)

  keep <- density != 0L
  out  <- data.frame(groupColsList(staticContributions, keep),
                     stringsAsFactors = FALSE)
  out  <- attachContributions(out, effectNames,
                              contribMat[keep, , drop = FALSE], flip = TRUE)
  out[names(sd)] <- lapply(sd, `[`, keep)

  if (details) {
    out$changeUtil <- utility[keep]
    out$changeProb <- changeProb[keep]
    if (!is.null(tieProb)) out$tieProb <- tieProb[keep]
  }

  if (massContrasts) {
    diffColName <- intersect(c("secondDiff", "secondRiskRatio"), names(sd))[1L]
    if (!is.na(diffColName) && diffColName == "secondDiff") {
      mc <- computeMassContrasts(
        firstDiff = sd[["secondDiff"]][keep],
        density   = density[keep],
        ego       = staticContributions$ego[keep],
        period    = staticContributions$period[keep],
        group     = if (!is.null(staticContributions$group))
                      staticContributions$group[keep] else rep(1L, sum(keep)),
        type      = type)
      out$massCreation    <- mc$massCreation
      out$massDissolution <- mc$massDissolution
    }
  }
  out
}

# ---- Print method ------------------------------------------------------

##@print.utilitySpec Print
print.utilitySpec <- function(x, ...) {
  cat("utilitySpec\n")
  cat("  Density: implicit factor c(-1, 1), interacts with all effects\n")
  if (length(x$interactions) == 0L) {
    cat("  (all effects atomic)\n")
  } else {
    cat("  Interactions:\n")
    for (nm in names(x$interactions)) {
      t <- x$interactions[[nm]]
      fstr <- ifelse(t$isFactor,
                     paste0("factor(", t$factors, ")"),
                     t$factors)
      cat(sprintf("    %s ~ %s\n", nm, paste(fstr, collapse = " : ")))
    }
  }
  if (length(x$perturb) > 0L) {
    cat("  Perturbation defaults:\n")
    for (nm in names(x$perturb)) {
      p <- x$perturb[[nm]]
      if (p$type == "contrast")
        cat(sprintf("    %s: contrast c(%s)\n", nm,
                    paste(p$levels, collapse = ", ")))
      else
        cat(sprintf("    %s: diff = %s\n", nm, p$value))
    }
  }
  invisible(x)
}
