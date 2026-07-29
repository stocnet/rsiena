## postestTargets.R — target specification for postestimation marginal effects
##
## A "target" is one perturbation to evaluate: an effect, how it is shifted
## (a step or a contrast), and optionally a second effect for a cross-partial.
##
## The object is a data.frame with one row per target, built and modified the
## same way an effects specification is:
##
##   mytg <- make_targets(ans)                       # all non-rate effects
##   mytg <- set_target(mytg, transTrip, diff = 2)   # include and tune
##   mytg <- set_second_diff(mytg, c(density, recip))
##
## Overridable settings
## --------------------
## level, condition, accumulated, rateWeight and massContrasts may be set once
## for the whole model and overridden per target.  Distinguishing "not set for
## this target" from "set to NULL for this target" matters: an explicit NULL
## condition on a target must switch conditioning OFF even when the model sets
## one.  A data.frame cell cannot express that difference, so explicitly-set
## values live in a `.overrides` list-column holding a named list per row --
## a name is present iff the user set that field for that target.
##
## This mirrors the `"level" %in% names(spec)` test the estimation path already
## uses, and lowers directly onto it.


## Effects whose perturbation default is a contrast rather than a step,
## because their change statistic is not on a scale a unit step means
## anything on: density's only ever takes +/-1, and reciprocity's is whether
## the reverse tie is there, so 0 or 1.
.targetDefaultContrast <- list(density = c(-1, 1),
                               recip   = c(0, 1))

##@make_postest_targets Postestimation
make_postest_targets <- function(x, ...) UseMethod("make_postest_targets", x)

##@make_postest_targets.sienaFit Postestimation
##
## Holds three things: the target rows, the model-level defaults those rows
## inherit, and the effects/depvar the targets were enumerated from.  Inputs
## enter where they are used -- `effects` here because targets are derived from
## it; `algorithm` is NOT here because it is only needed for dynamic
## simulation, and so lives in set_postest_algo_saom().
##
## `effects` is required rather than defaulted from the fit: the fit does carry
## it, but edge cases (conditional estimation, stripped effect objects) make
## the implicit route unreliable, so the caller states it.
make_postest_targets.sienaFit <- function(x, data = NULL, effects = NULL,
                                  depvar = NULL,
                                  dynamic = FALSE,
                                  type = c("changeProb", "tieProb"),
                                  mainEffect = c("riskDifference", "riskRatio"),
                                  level = "period",
                                  condition = NULL,
                                  egoNormalize = TRUE,
                                  accumulated = FALSE,
                                  rateWeight = FALSE,
                                  massContrasts = NULL,
                                  dependencies = NULL,
                                  includeDefaults = TRUE, ...) {
    if (!is.logical(includeDefaults) || length(includeDefaults) != 1L ||
        is.na(includeDefaults))
        stop("'includeDefaults' must be a single non-NA logical.", call. = FALSE)
    type       <- match.arg(type)
    mainEffect <- match.arg(mainEffect)
    if (!is.logical(dynamic) || length(dynamic) != 1L || is.na(dynamic))
        stop("'dynamic' must be a single non-NA logical.", call. = FALSE)

    if (isTRUE(accumulated) && !dynamic)
        stop("accumulated = TRUE requires dynamic = TRUE: an accumulated ",
             "marginal effect sums over the ministeps of a simulated chain, ",
             "which a static evaluation does not have.", call. = FALSE)
    if (isTRUE(rateWeight) && dynamic)
        message("Note: rate-weighting is absorbed by the simulation when ",
                "dynamic = TRUE (actors are drawn proportional to lambda), ",
                "so 'rateWeight' has no additional effect here.")

    if (is.null(effects))
        stop("'effects' must be supplied -- the effects object the model was ",
             "fitted with, e.g. effects = mymodel.", call. = FALSE)
    if (!inherits(effects, "sienaEffects"))
        stop("'effects' must be a sienaEffects object.", call. = FALSE)
    if (is.null(depvar)) {
        nm <- unique(effects[["name"]])
        depvar <- nm[1L]
    }
    hasRate <- !("basicRate" %in% names(effects)) || any(effects$basicRate)

    reg <- buildEffectNameRegistry(effects, depvar = depvar,
                                   includeOnly = TRUE, append_parm = FALSE)
    if (nrow(reg) == 0L)
        stop("No included effects for dependent variable '", depvar, "'.",
             call. = FALSE)

    reg <- reg[!duplicated(reg$base_name), , drop = FALSE]

    ## Rate effects are dropped rather than listed-but-unselected: they have no
    ## meaningful first difference in this framework, so offering them as
    ## selectable targets would be misleading.  They also share the short name
    ## "Rate" across periods, which would put duplicate identifiers in the
    ## table.
    reg <- reg[reg$effect_type != "rate", , drop = FALSE]
    if (nrow(reg) == 0L)
        stop("No non-rate effects for dependent variable '", depvar,
             "'; there is nothing to perturb.", call. = FALSE)

    isRate <- rep(FALSE, nrow(reg))
    n      <- nrow(reg)

    ## Declared dependencies shape the DEFAULT selection.  An effect that is
    ## declared as the consequence of others is not an independent target:
    ## perturbing it on its own contradicts the relationship just stated.  So
    ## it is present in the table but not selected by default; the user can
    ## still select it explicitly if that is genuinely what they want.
    ##
    ## Parsing here also validates the declaration against the model at
    ## construction time, rather than several steps later at lowering.
    derived <- character(0)
    if (length(dependencies) > 0L) {
        for (d in dependencies) {
            pd  <- .parseDependency(d)
            unk <- setdiff(c(pd$target, pd$terms), reg$short_raw)
            if (length(unk))
                stop("Dependency ", deparse(d), " names effect(s) not in the ",
                     "model: ", paste(unk, collapse = ", "), ".\n",
                     "Available: ", paste(reg$short_raw, collapse = ", "), ".",
                     call. = FALSE)
            derived <- c(derived, pd$target)
        }
        derived <- unique(derived)
    }

    diff1 <- vector("list", n)
    cont1 <- vector("list", n)
    for (i in seq_len(n)) {
        if (isRate[i]) next                      # rate: no perturbation
        dc <- .targetDefaultContrast[[reg$short_raw[i]]]
        if (!is.null(dc)) cont1[[i]] <- dc       # contrast-defaulted effect
        else              diff1[[i]] <- 1        # unit step
    }

    tg <- list(
        name            = paste0(.identNames(reg$short_raw,
                                             reg$short_with_covar), "_fd"),
        effectName1     = reg$short_raw,
        effectType      = reg$effect_type,
        .covar1           = reg$covar1,
        .covar2           = reg$covar2,
        .qual1          = reg$short_with_covar,
        .qual2          = rep(NA_character_, n),
        include         = includeDefaults & !isRate &
                          !(reg$short_raw %in% derived),
        second          = rep(FALSE, n),
        effectName2     = rep(NA_character_, n),
        diff1           = diff1,
        contrast1       = cont1,
        diff2           = vector("list", n),
        contrast2       = vector("list", n),
        interaction1    = rep(FALSE, n),
        intEffectNames1 = vector("list", n),
        modEffectNames1 = vector("list", n),
        perturbType1    = vector("list", n),
        perturbType2    = vector("list", n),
        interaction2    = rep(FALSE, n),
        intEffectNames2 = vector("list", n),
        modEffectNames2 = vector("list", n),
        returnDecisionDetails = rep(FALSE, n),
        .seq            = as.numeric(seq_len(n)),
        .overrides      = replicate(n, list(), simplify = FALSE)
    )
    attr(tg, "row.names") <- .set_row_names(n)
    attr(tg, "depvar")       <- depvar
    attr(tg, "effects")      <- effects
    attr(tg, "hasRate")      <- hasRate
    attr(tg, "dependencies") <- dependencies
    ## Model-level defaults; per-target entries in .overrides win over these.
    attr(tg, "defaults") <- list(
        dynamic = dynamic,
        type = type, mainEffect = mainEffect, level = level,
        condition = condition, egoNormalize = egoNormalize,
        accumulated = accumulated, rateWeight = rateWeight,
        massContrasts = massContrasts
    )
    attr(tg, "version") <- utils::packageDescription("RSiena", fields = "Version")
    class(tg) <- c("sienaPostestTargets", "data.frame")
    tg
}


## Resolve bare symbols or characters to effect names, as set_effect() does.
.targetNames <- function(subst, env) {
    val <- tryCatch(eval(subst, envir = env), error = function(e) NULL)
    if (is.character(val)) return(val)
    if (is.call(subst) && identical(as.character(subst[[1L]]), "c"))
        return(vapply(as.list(subst)[-1L], deparse, character(1L)))
    deparse(subst)
}

## Rows addressable by an effect name: the ones enumerated from the model.
## Second differences are added later and are addressed by their target name.
.baseRows <- function(x) which(!x$second)

## short name where that is unique, qualified short name where it is not.
.identNames <- function(short, qual) {
    dup <- short %in% short[duplicated(short)]
    ifelse(dup, qual, short)
}

## The identifying name of every row, second-difference rows included (those
## already store idents, resolved when they were added).
.rowIdent <- function(x) {
    ident <- x$effectName1
    b <- .baseRows(x)
    ident[b] <- .identNames(x$effectName1[b], x$.qual1[b])
    ident
}

.targetRow <- function(x, nm, covar1 = NULL, covar2 = NULL) {
    b   <- .baseRows(x)
    hit <- b[x$effectName1[b] == nm]
    ## Also accept the qualified form, so that whatever print() showed can be
    ## typed straight back in.
    if (length(hit) == 0L)
        hit <- b[!is.na(x$.qual1[b]) & x$.qual1[b] == nm]
    if (length(hit) == 0L)
        stop("Effect '", nm, "' is not among the targets: ",
             paste(unique(.rowIdent(x)[b]), collapse = ", "), ".",
             call. = FALSE)

    if (!is.null(covar1)) hit <- hit[x$.covar1[hit] == covar1]
    if (!is.null(covar2)) hit <- hit[x$.covar2[hit] == covar2]
    if (length(hit) == 0L)
        stop("No target for effect '", nm, "' with ",
             if (!is.null(covar1)) paste0("covar1 = \"", covar1, "\"") else "",
             if (!is.null(covar1) && !is.null(covar2)) ", " else "",
             if (!is.null(covar2)) paste0("covar2 = \"", covar2, "\"") else "",
             ". Available for '", nm, "': ",
             .describeCandidates(x, b[x$effectName1[b] == nm]), ".",
             call. = FALSE)
    if (length(hit) > 1L)
        stop("Effect '", nm, "' matches ", length(hit), " targets, so it does ",
             "not identify one: ", .describeCandidates(x, hit), ".\n",
             "Name the covariate, as set_effect() does -- e.g. ",
             if (nzchar(x$.covar1[hit[1L]]))
                 paste0("covar1 = \"", x$.covar1[hit[1L]], "\"")
             else paste0("use the qualified name '", x$.qual1[hit[1L]], "'"),
             ".", call. = FALSE)
    hit
}

## How the candidate targets differ, in the terms used to tell them apart.
.describeCandidates <- function(x, rows) {
    paste(vapply(rows, function(i) {
        cv <- .covarLabel(x$.covar1[i], x$.covar2[i])
        if (nzchar(cv)) paste0(x$effectName1[i], " on ", cv)
        else x$.qual1[i]
    }, character(1L)), collapse = ", ")
}


##@set_target Postestimation
set_target <- function(x, ...) UseMethod("set_target", x)

##@set_target.sienaPostestTargets Postestimation
set_target.sienaPostestTargets <- function(x, shortNames,
                                           diff = NULL, contrast = NULL,
                                           covar1 = NULL, covar2 = NULL,
                                           perturbType = NULL,
                                           interaction = NULL,
                                           intEffectNames = NULL,
                                           modEffectNames = NULL,
                                           returnDecisionDetails = NULL,
                                           name = NULL,
                                           level, condition, accumulated,
                                           rateWeight, massContrasts,
                                           include = TRUE, verbose = TRUE,
                                           ...) {
    if (!hasArg(shortNames))
        stop("set_target needs one or more shortNames.", call. = FALSE)
    nms <- .targetNames(substitute(shortNames), parent.frame())

    if (!is.null(diff) && !is.null(contrast))
        stop("Supply either 'diff' or 'contrast', not both.", call. = FALSE)
    ## covar1/covar2 pick out WHICH target a short name refers to, so like `name`
    ## they only mean anything when one target is being addressed.
    if ((!is.null(covar1) || !is.null(covar2)) && length(nms) != 1L)
        stop("'covar1'/'covar2' identify a single target, so they can only be ",
             "given with a single short name.", call. = FALSE)
    ## `name` labels the output for this target, as an effectList entry's list
    ## name does.  Only meaningful for a single target at a time.
    if (!is.null(name)) {
        if (length(nms) != 1L)
            stop("'name' can only be given when setting a single target.",
                 call. = FALSE)
        if (!is.character(name) || length(name) != 1L || is.na(name))
            stop("'name' must be a single non-NA string.", call. = FALSE)
    }

    ## Only fields the caller actually named become overrides; an unset field
    ## must inherit from the model, and an explicit NULL must override it.
    ## `ov[[f]] <- NULL` would DELETE the entry, making an explicit
    ## `condition = NULL` indistinguishable from not setting it at all --
    ## exactly the distinction this mechanism exists to preserve.
    ## `ov[f] <- list(value)` stores NULL as a present entry.
    ov <- list()
    for (f in c("level", "condition", "accumulated", "rateWeight",
                "massContrasts"))
        if (!missing_arg_named(f, environment())) ov[f] <- list(get(f))

    for (nm in nms) {
        i <- .targetRow(x, nm, covar1, covar2)
        ## Re-tuning an already-selected target must not reorder it; only a
        ## newly selected one moves to the end.
        if (include && !isTRUE(x$include[i]))
            x$.seq[i] <- max(c(x$.seq, 0), na.rm = TRUE) + 1
        x$include[i] <- include
        ## `x$col[[i]] <- NULL` would REMOVE the element and shift the column;
        ## `x$col[i] <- list(NULL)` sets it to NULL in place.
        if (!is.null(diff))     { x$diff1[[i]] <- diff;     x$contrast1[i] <- list(NULL) }
        if (!is.null(contrast)) { x$contrast1[[i]] <- contrast; x$diff1[i] <- list(NULL) }
        if (!is.null(perturbType)) x$perturbType1[[i]] <- perturbType
        ## Compound-effect perturbation.  This is the pre-dependency-layer
        ## mechanism: the moderating effects are named explicitly rather than
        ## derived from a declared relationship.
        if (!is.null(interaction))    x$interaction1[i]       <- interaction
        if (!is.null(intEffectNames)) x$intEffectNames1[[i]]  <- intEffectNames
        if (!is.null(modEffectNames)) x$modEffectNames1[[i]]  <- modEffectNames
        if (!is.null(returnDecisionDetails))
            x$returnDecisionDetails[i] <- returnDecisionDetails
        if (!is.null(name)) {
            if (name %in% x$name[-i])
                stop("A target named '", name, "' already exists.",
                     call. = FALSE)
            x$name[i] <- name
        }
        if (length(ov))
            x$.overrides[[i]][names(ov)] <- ov
        if (verbose)
            message("set_target: ", if (include) "included" else "excluded",
                    " '", nm, "'",
                    if (length(ov))
                        paste0(" (", paste(names(ov), collapse = ", "),
                               " overridden)") else "")
    }
    x
}


## missing() cannot be called on a variable in a loop, so test the call frame.
missing_arg_named <- function(nm, env) {
    do.call(missing, list(as.name(nm)), envir = env)
}


## --------------------------------------------------------------------------
## .parsePerturbList — the per-effect form of a multi-way difference
##
##   list(transTrip = list(diff = 1),
##        recip     = list(contrast = c(0, 1)))
##
## Each element says what is done to ONE effect's change statistic, keyed by
## that effect's short name and using the same vocabulary as set_target().
## The numbered alternative (diff1 =, contrast2 = ...) makes the reader carry
## the mapping from suffix to position in their head, and does not extend:
## a third-order difference would need diff3/contrast3/perturbType3 formals,
## whereas a list simply gains a third element.
##
## Returns a list of per-effect specs, in the order given, names = effects.
## --------------------------------------------------------------------------
.perturbKeys <- c("diff", "contrast", "covar1", "covar2", "perturbType",
                  "interaction", "intEffectNames", "modEffectNames")

.parsePerturbList <- function(spec) {
    nms <- names(spec)
    if (is.null(nms) || any(!nzchar(nms)))
        stop("Every element of the perturbation list must be named with an ",
             "effect short name, e.g. list(transTrip = list(diff = 1), ",
             "recip = list(contrast = c(0, 1))).", call. = FALSE)
    ## Repeated names are allowed: crossing an effect with itself is a real
    ## quantity (the same base effect entering through two different compound
    ## interactions, say), and the numbered form has always permitted it.

    lapply(seq_along(spec), function(i) {
        e  <- spec[[i]]
        nm <- nms[i]
        ## NULL means "whatever this effect's default perturbation is",
        ## matching an unset diff/contrast in the numbered form.
        if (is.null(e)) return(list())
        if (!is.list(e))
            stop("Perturbation for '", nm, "' must be a list, e.g. ",
                 nm, " = list(diff = 1) or ", nm,
                 " = list(contrast = c(0, 1)); got ", class(e)[1L], ".",
                 call. = FALSE)
        bad <- setdiff(names(e), .perturbKeys)
        if (length(bad))
            stop("Unknown perturbation setting(s) for '", nm, "': ",
                 paste(bad, collapse = ", "), ". Available: ",
                 paste(.perturbKeys, collapse = ", "), ".", call. = FALSE)
        if (!is.null(e$diff) && !is.null(e$contrast))
            stop("Supply either 'diff' or 'contrast' for '", nm,
                 "', not both.", call. = FALSE)
        e
    })
}


##@set_second_diff Postestimation
set_second_diff <- function(x, ...) UseMethod("set_second_diff", x)

##@set_second_diff.sienaPostestTargets Postestimation
##
## `perturbations` takes either the per-effect list form (preferred) or the
## older pair-of-names form with numbered arguments; see .parsePerturbList.
set_second_diff.sienaPostestTargets <- function(x, perturbations,
                                                diff1 = NULL, diff2 = NULL,
                                                contrast1 = NULL,
                                                contrast2 = NULL,
                                                interaction1 = NULL,
                                                intEffectNames1 = NULL,
                                                modEffectNames1 = NULL,
                                                interaction2 = NULL,
                                                intEffectNames2 = NULL,
                                                modEffectNames2 = NULL,
                                                perturbType1 = NULL,
                                                perturbType2 = NULL,
                                                returnDecisionDetails = FALSE,
                                                name = NULL,
                                                include = TRUE, verbose = TRUE,
                                                ...) {
    if (!hasArg(perturbations))
        stop("set_second_diff needs two effects to cross, e.g.\n",
             "  set_second_diff(tg, list(transTrip = list(diff = 1),\n",
             "                           recip = list(contrast = c(0, 1))))",
             call. = FALSE)

    ## Which form was used?  A bare c(a, b) of effect names does not evaluate
    ## (the symbols are effect names, not objects), so failing to evaluate is
    ## itself the signal for the older form; so is a character vector.
    .spec <- substitute(perturbations)
    .val  <- tryCatch(eval(.spec, parent.frame()), error = function(e) NULL)

    if (is.list(.val)) {
        if (any(vapply(list(diff1, diff2, contrast1, contrast2, perturbType1,
                            perturbType2, interaction1, interaction2,
                            intEffectNames1, intEffectNames2,
                            modEffectNames1, modEffectNames2),
                       Negate(is.null), logical(1L))))
            stop("The perturbation list already says what is done to each ",
                 "effect; the numbered arguments (diff1, contrast2, ...) are ",
                 "the older way of saying the same thing. Use one or the ",
                 "other.", call. = FALSE)
        if (length(.val) != 2L)
            stop("A second difference crosses exactly two effects; the ",
                 "perturbation list has ", length(.val), ".",
                 if (length(.val) > 2L)
                     " Higher-order differences are not supported yet."
                 else "", call. = FALSE)
        parsed <- .parsePerturbList(.val)
        nms    <- names(.val)
        p1 <- parsed[[1L]]; p2 <- parsed[[2L]]
        covars1 <- list(p1$covar1, p2$covar1)
        covars2 <- list(p1$covar2, p2$covar2)
        diff1 <- p1$diff; contrast1 <- p1$contrast; perturbType1 <- p1$perturbType
        diff2 <- p2$diff; contrast2 <- p2$contrast; perturbType2 <- p2$perturbType
        interaction1 <- p1$interaction; interaction2 <- p2$interaction
        intEffectNames1 <- p1$intEffectNames; intEffectNames2 <- p2$intEffectNames
        modEffectNames1 <- p1$modEffectNames; modEffectNames2 <- p2$modEffectNames
    } else {
        nms <- .targetNames(.spec, parent.frame())
        if (length(nms) != 2L)
            stop("set_second_diff needs exactly two effects, got ",
                 length(nms), ".", call. = FALSE)
        if (!is.null(diff1) && !is.null(contrast1))
            stop("Supply either 'diff1' or 'contrast1', not both.",
                 call. = FALSE)
        if (!is.null(diff2) && !is.null(contrast2))
            stop("Supply either 'diff2' or 'contrast2', not both.",
                 call. = FALSE)
        covars1 <- list(NULL, NULL)
        covars2 <- list(NULL, NULL)
    }

    ## Resolve each component to a row now, so an ambiguous short name is
    ## caught here rather than becoming a silently-wrong perturbation.  The
    ## row supplies both spellings: the readable one for naming and printing,
    ## the qualified one for the engine.
    rows  <- vapply(seq_along(nms),
                    function(k) .targetRow(x, nms[k], covars1[[k]], covars2[[k]]),
                    integer(1L))
    ident <- .rowIdent(x)[rows]
    quals <- x$.qual1[rows]
    nms   <- ident

    if (is.null(name)) name <- paste0(paste(nms, collapse = "_"), "_sd")
    if (name %in% x$name)
        stop("A target named '", name, "' already exists.", call. = FALSE)

    ## Default the first component the same way make_targets() would.
    if (is.null(diff1) && is.null(contrast1)) {
        dc <- .targetDefaultContrast[[nms[1L]]]
        if (!is.null(dc)) contrast1 <- dc else diff1 <- 1
    }
    if (is.null(diff2) && is.null(contrast2)) {
        dc <- .targetDefaultContrast[[nms[2L]]]
        if (!is.null(dc)) contrast2 <- dc else diff2 <- 1
    }

    ## Append by extending each column directly.  rbind() on a data.frame with
    ## list columns recycles them elementwise and warns; building the row as a
    ## list and concatenating per column is exact.
    atomic_new <- list(
        name = name, effectName1 = nms[1L], effectType = NA_character_,
        ## The qualified names of the two components, for lowering.  The
        ## covariate columns stay empty: a second difference has no covariate
        ## of its own, and its components' are read off their own rows.
        .covar1 = "", .covar2 = "",
        .qual1 = quals[1L], .qual2 = quals[2L],
        include = include, second = TRUE, effectName2 = nms[2L],
        interaction1 = isTRUE(interaction1),
        interaction2 = isTRUE(interaction2),
        returnDecisionDetails = isTRUE(returnDecisionDetails),
        .seq = max(c(x$.seq, 0), na.rm = TRUE) + 1
    )
    list_new <- list(
        diff1 = diff1, contrast1 = contrast1,
        diff2 = diff2, contrast2 = contrast2,
        intEffectNames1 = intEffectNames1, modEffectNames1 = modEffectNames1,
        intEffectNames2 = intEffectNames2, modEffectNames2 = modEffectNames2,
        perturbType1 = perturbType1, perturbType2 = perturbType2,
        .overrides = list()
    )

    out <- unclass(x)
    for (cn in names(out)) {
        if (cn %in% names(atomic_new))
            out[[cn]] <- c(out[[cn]], atomic_new[[cn]])
        else
            out[[cn]] <- c(out[[cn]], list_new[cn])   # single-bracket keeps NULL
    }
    attr(out, "row.names") <- .set_row_names(length(out$name))
    attr(out, "depvar")    <- attr(x, "depvar")
    attr(out, "version")   <- attr(x, "version")
    class(out) <- c("sienaPostestTargets", "data.frame")

    if (verbose)
        message("set_second_diff: added '", name, "' (",
                nms[1L], " x ", nms[2L], ")")
    out
}



## "gender", or "gender x age" when an effect carries two covariates.
.covarLabel <- function(c1, c2) {
    c1 <- ifelse(is.na(c1), "", c1)
    c2 <- ifelse(is.na(c2), "", c2)
    ifelse(c1 == "", "",
           ifelse(c2 == "", c1, paste(c1, c2, sep = " x ")))
}

## What to call an effect on screen: the readable short name, plus the
## covariate that says WHICH one when the model carries it more than once.
## Where there is no covariate to name (two unspInt terms, say), the
## qualified name is the only thing that tells them apart, so use that.
.rowLabel <- function(short, ident, covar1, covar2) {
    cv <- .covarLabel(covar1, covar2)
    if (nzchar(cv)) paste0(short, " (", cv, ")") else ident
}

## Same, for a second difference's components, which are stored by their
## identifying name and so resolve to exactly one row.
.targetEffectLabel <- function(x, nm) {
    b     <- .baseRows(x)
    ident <- .rowIdent(x)
    hit   <- b[ident[b] == nm]
    if (length(hit) != 1L) return(nm)
    .rowLabel(x$effectName1[hit], ident[hit], x$.covar1[hit], x$.covar2[hit])
}

## What is done to the change statistic.  A contrast names two levels the
## statistic is held at; a difference steps it from wherever it is.
.perturbLabel <- function(diff, contrast) {
    if (!is.null(contrast)) {
        if (length(contrast) == 2L)
            paste0(contrast[1L], " -> ", contrast[2L])
        else paste0("contrast (", paste(contrast, collapse = ", "), ")")
    } else if (!is.null(diff)) {
        paste0(if (diff >= 0) "+" else "", diff)
    } else "-"
}

##@print.sienaPostestTargets Postestimation
print.sienaPostestTargets <- function(x, ...) {
    ## Work with row indices into x: a row's readable name depends on whether
    ## its short name is unique in the WHOLE table, which a subset cannot say.
    ## Ordered as .targetsToEffectList() lowers, so the printed order is the
    ## order the results come back in.
    idx <- which(x$include)
    idx <- idx[order(x$.seq[idx])]
    inc <- x[idx, , drop = FALSE]
    identAll <- .rowIdent(x)
    d <- attr(x, "defaults")
    cat("RSiena postestimation targets  (",
        nrow(inc), " of ", nrow(x), " selected",
        if (!is.null(attr(x, "depvar")))
            paste0("; depvar '", attr(x, "depvar"), "'") else "",
        ")\n", sep = "")
    if (!is.null(d)) {
        cat("  estimand : ",
            if (identical(d$mainEffect, "riskRatio")) "risk ratio in "
            else "risk difference in ",
            if (identical(d$type, "tieProb")) "tie probability"
            else "change probability",
            if (isTRUE(d$dynamic)) ", dynamic (simulated chains)"
            else ", static (observed waves)",
            "\n", sep = "")
        extra <- c(
            if (isTRUE(d$accumulated)) "accumulated over ministeps",
            if (isTRUE(d$rateWeight))  "rate-weighted",
            if (isFALSE(d$egoNormalize)) "not ego-normalised",
            if (!is.null(d$massContrasts)) "with mass contrasts")
        cat("  reported : by ", d$level,
            if (!is.null(d$condition))
                paste0(", conditional on ", paste(d$condition, collapse = " + "))
            else "",
            if (length(extra)) paste0(", ", paste(extra, collapse = ", ")) else "",
            "\n", sep = "")
    }

    dep <- attr(x, "dependencies")
    if (length(dep)) {
        cat("  depends  : ", deparse(dep[[1L]]), "\n", sep = "")
        for (dd in dep[-1L])
            cat("             ", deparse(dd), "\n", sep = "")
    }

    if (nrow(inc) == 0L) {
        cat("  (none selected)\n")
        return(invisible(x))
    }

    anyOv <- any(vapply(inc$.overrides, length, integer(1L)) > 0L)
    if (anyOv && !is.null(d))
        cat("  (the above holds for every target except where [bracketed])\n")

    cat("\n")
    for (i in seq_len(nrow(inc))) {
        p1 <- .perturbLabel(inc$diff1[[i]], inc$contrast1[[i]])
        if (isTRUE(inc$second[i])) {
            what <- paste0(.targetEffectLabel(x, inc$effectName1[i]), " ", p1,
                           "  x  ",
                           .targetEffectLabel(x, inc$effectName2[i]), " ",
                           .perturbLabel(inc$diff2[[i]], inc$contrast2[[i]]))
            kind <- "second difference"
        } else {
            what <- paste0(.rowLabel(inc$effectName1[i], identAll[idx[i]],
                                     inc$.covar1[i], inc$.covar2[i]), " ", p1)
            kind <- "first difference"
        }
        ov  <- inc$.overrides[[i]]
        ovs <- if (length(ov))
            paste0("  [", paste(vapply(names(ov), function(f) {
                v <- ov[[f]]
                paste0(f, " = ", if (is.null(v)) "NULL"
                                 else paste(v, collapse = " + "))
            }, character(1L)), collapse = ", "), "]") else ""
        cat(sprintf("  %-22s %-34s %s%s\n",
                    inc$name[i], what, kind, ovs))
    }
    invisible(x)
}


## --------------------------------------------------------------------------
## .targetsToEffectList — lower a targets object onto the internal spec format
##
## Returns list(effectList, defaults):
##   effectList — named list, one entry per INCLUDED target, in exactly the
##                shape marginalEffects() already consumes
##   defaults   — the model-level settings, to be assigned as the call-level
##                variables the per-spec overrides fall back to
##
## Override semantics are preserved exactly as the estimation path applies
## them (R/sienaMargins.r, "Per-spec level/condition overrides"):
##
##   level, condition   -> present-in-spec wins, INCLUDING an explicit NULL
##   accumulated, rateWeight -> OR'd with the call-level value
##   massContrasts      -> non-NULL in spec wins
## --------------------------------------------------------------------------
.targetsToEffectList <- function(tg) {
    if (!inherits(tg, "sienaPostestTargets"))
        stop("'targets' must be a sienaPostestTargets object, as returned by ",
             "make_postest_targets().", call. = FALSE)
    keep <- which(tg$include)
    keep <- keep[order(tg$.seq[keep])]
    if (length(keep) == 0L)
        stop("No targets selected: every row has include = FALSE.",
             call. = FALSE)

    out <- vector("list", length(keep))
    for (k in seq_along(keep)) {
        i <- keep[k]
        e <- list(
            effectName1     = if (!is.na(tg$.qual1[i])) tg$.qual1[i]
                              else tg$effectName1[i],
            diff1           = tg$diff1[[i]],
            contrast1       = tg$contrast1[[i]],
            interaction1    = isTRUE(tg$interaction1[i]),
            intEffectNames1 = tg$intEffectNames1[[i]],
            modEffectNames1 = tg$modEffectNames1[[i]],
            second          = isTRUE(tg$second[i]),
            perturbType1    = tg$perturbType1[[i]],
            perturbType2    = tg$perturbType2[[i]],
            returnDecisionDetails = isTRUE(tg$returnDecisionDetails[i])
        )
        if (isTRUE(tg$second[i])) {
            e$effectName2     <- if (!is.na(tg$.qual2[i])) tg$.qual2[i]
                                 else tg$effectName2[i]
            e$diff2           <- tg$diff2[[i]]
            e$contrast2       <- tg$contrast2[[i]]
            e$interaction2    <- isTRUE(tg$interaction2[i])
            e$intEffectNames2 <- tg$intEffectNames2[[i]]
            e$modEffectNames2 <- tg$modEffectNames2[[i]]
        }
        ov <- tg$.overrides[[i]]
        for (f in names(ov)) e[f] <- list(ov[[f]])

        if (!("massContrasts" %in% names(ov)) &&
            !is.null(d0 <- attr(tg, "defaults")$massContrasts))
            e$massContrasts <- d0

        e <- .applyDependencies(e, attr(tg, "dependencies"))
        out[[k]] <- e
    }
    names(out) <- tg$name[keep]
    defs <- attr(tg, "defaults")
    defs$massContrasts <- NULL
    list(effectList = out, defaults = defs)
}


##@set_dependency Postestimation
##
## Declare how effects relate, so a perturbation propagates to the change
## statistics that move with it.  Model-level: a relationship holds regardless
## of which effect is being perturbed, so it is stated once.
##
## A declared dependency is applied when the targets are lowered: perturbing an
## effect that participates in one also moves the dependent effect's change
## statistic.  Only multiplicative relationships are supported; other forms are
## rejected rather than approximated.
set_dependency <- function(x, ...) UseMethod("set_dependency", x)

##@set_dependency.sienaPostestTargets Postestimation
set_dependency.sienaPostestTargets <- function(x, ..., verbose = TRUE) {
    deps <- list(...)
    if (length(deps) == 0L)
        stop("set_dependency needs at least one dependency.", call. = FALSE)
    ok <- vapply(deps, function(d) inherits(d, "formula"), logical(1L))
    if (!all(ok))
        stop("Dependencies must be one-sided formulas, e.g. ",
             "egoXaltX ~ egoX:altX.", call. = FALSE)
    for (d in deps)
        if (length(d) != 3L)
            stop("A dependency needs a left-hand side naming the dependent ",
                 "effect: egoXaltX ~ egoX:altX.", call. = FALSE)

    for (d in deps) {
        lhs <- tryCatch(.parseDependency(d)$target, error = function(e) NULL)
        if (!is.null(lhs) && lhs %in% x$effectName1[x$include])
            warning("'", lhs, "' is a selected target and is now declared as ",
                    "derived from other effects. Perturbing it directly holds ",
                    "its components fixed, which the declaration says is not ",
                    "achievable. Drop it with set_target(x, ", lhs,
                    ", include = FALSE), or keep it if a direct perturbation ",
                    "of its change statistic is what you intend.",
                    call. = FALSE)
    }

    cur <- attr(x, "dependencies")
    attr(x, "dependencies") <- c(cur, deps)
    if (verbose)
        message("set_dependency: ", length(deps), " added (",
                length(cur) + length(deps), " total).")
    x
}


## --------------------------------------------------------------------------
## .parseDependency — one dependency formula -> list(target, op, terms)
##
## Iteration 1 accepts pure interactions only:
##
##     egoXaltX ~ egoX:altX
##
## meaning "egoXaltX's change statistic is the product of egoX's and altX's".
## --------------------------------------------------------------------------
.parseDependency <- function(f) {
    if (!inherits(f, "formula") || length(f) != 3L)
        stop("A dependency must have the form target ~ a * b.", call. = FALSE)
    lhs <- all.vars(f[[2L]])
    if (length(lhs) != 1L)
        stop("The left-hand side of a dependency must name exactly one ",
             "effect.", call. = FALSE)
    rhs <- f[[3L]]
    op  <- if (is.call(rhs)) as.character(rhs[[1L]]) else ""
    if (identical(op, "*"))
        stop("Use ':' rather than '*' for a dependency: ", deparse(f),
             ". In R's formula algebra 'a * b' means 'a + b + a:b' -- the main ",
             "effects as well as their interaction -- whereas a dependency ",
             "states only the product. Write: ",
             deparse(f[[2L]]), " ~ ",
             paste(vapply(as.list(rhs)[-1L], deparse, character(1L)),
                   collapse = ":"), call. = FALSE)
    if (!identical(op, ":"))
        stop("Only pure-interaction dependencies are supported so far: ",
             "target ~ a:b. Got: ", deparse(f), call. = FALSE)
    terms <- vapply(as.list(rhs)[-1L], function(z) {
        v <- all.vars(z)
        if (length(v) != 1L)
            stop("Each side of ':' must be a single effect name. Got: ",
                 deparse(f), call. = FALSE)
        v
    }, character(1L))
    list(target = lhs, op = "*", terms = terms)
}


## --------------------------------------------------------------------------
## .applyDependencies — turn declarations into the compound-effect arguments
##
## Given a target perturbing effect X and a dependency `A ~ X * C`, perturbing
## X also moves A's change statistic by delta * s_C.  That is exactly what the
## existing compound-effect path computes:
##
##     Delta u = d * delta * (theta_X + sum_k s_{mod_k} * theta_{int_k})
##
## so the declaration is resolved INTO intEffectNames / modEffectNames rather
## than given a separate numeric path.  Same arithmetic, same tested code --
## the dependency layer only removes the need to work the arguments out by hand.
##
## The relation is symmetric in its two right-hand terms: whichever one is
## being perturbed, the other is the moderator.
## --------------------------------------------------------------------------
.applyDependencies <- function(e, deps) {
    if (length(deps) == 0L) return(e)

    ## Resolve for ONE side of the target.  A second difference perturbs two
    ## effects, and a dependency can cover either -- resolving only the first
    ## would silently drop the declaration for the second half.
    resolve1 <- function(x) {
        ints <- character(0); mods <- character(0)
        for (d in deps) {
            p <- .parseDependency(d)
            hit <- which(p$terms == x)
            if (length(hit) == 0L) next
            if (length(hit) == 2L)
                stop("Dependency '", deparse(d), "' uses '", x,
                     "' on both sides; a self-product is not supported.",
                     call. = FALSE)
            ints <- c(ints, p$target)
            mods <- c(mods, p$terms[-hit])
        }
        if (length(ints) == 0L) NULL else list(ints = ints, mods = mods)
    }

    r1 <- resolve1(e$effectName1)
    if (!is.null(r1)) {
        if (isTRUE(e$interaction1))
            stop("Target '", e$effectName1, "' sets interaction arguments ",
                 "explicitly and is also covered by a declared dependency. ",
                 "Use one or the other.", call. = FALSE)
        e$interaction1    <- TRUE
        e$intEffectNames1 <- r1$ints
        e$modEffectNames1 <- r1$mods
    }

    if (isTRUE(e$second) && !is.null(e$effectName2)) {
        r2 <- resolve1(e$effectName2)
        if (!is.null(r2)) {
            if (isTRUE(e$interaction2))
                stop("The second effect of target pair '", e$effectName1,
                     " x ", e$effectName2, "' sets interaction arguments ",
                     "explicitly and is also covered by a declared ",
                     "dependency. Use one or the other.", call. = FALSE)
            e$interaction2    <- TRUE
            e$intEffectNames2 <- r2$ints
            e$modEffectNames2 <- r2$mods
        }
    }
    e
}
