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


## Effects whose perturbation default is a contrast rather than a step.
## density's change statistic only takes +/-1, so a unit step is meaningless.
.targetDefaultContrast <- list(density = c(-1, 1))


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
                                  type = c("changeProb", "tieProb"),
                                  mainEffect = c("riskDifference", "riskRatio"),
                                  level = "period",
                                  condition = NULL,
                                  egoNormalize = TRUE,
                                  accumulated = FALSE,
                                  rateWeight = FALSE,
                                  massContrasts = NULL,
                                  combineSameLevel = FALSE,
                                  dependencies = NULL,
                                  includeDefaults = TRUE, ...) {
    if (!is.logical(includeDefaults) || length(includeDefaults) != 1L ||
        is.na(includeDefaults))
        stop("'includeDefaults' must be a single non-NA logical.", call. = FALSE)
    type       <- match.arg(type)
    mainEffect <- match.arg(mainEffect)

    if (is.null(effects))
        stop("'effects' must be supplied -- the effects object the model was ",
             "fitted with, e.g. effects = mymodel.", call. = FALSE)
    if (!inherits(effects, "sienaEffects"))
        stop("'effects' must be a sienaEffects object.", call. = FALSE)
    if (is.null(depvar)) {
        nm <- unique(effects[["name"]])
        depvar <- nm[1L]
    }

    ## Checked here rather than deep inside marginalEffects(): a stripped
    ## effects object cannot drive the dynamic path, and the caller is looking
    ## at the argument right now.
    hasRate <- !("basicRate" %in% names(effects)) || any(effects$basicRate)

    reg <- buildEffectNameRegistry(effects, depvar = depvar,
                                   includeOnly = TRUE, append_parm = FALSE)
    if (nrow(reg) == 0L)
        stop("No included effects for dependent variable '", depvar, "'.",
             call. = FALSE)

    ## One row per underlying effect: creation/endow/eval variants of the same
    ## effect are a single perturbation target, not three.
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

    ## vector("list", n) is already all-NULL, so only NON-NULL defaults are
    ## assigned.  Writing `lst[[i]] <- NULL` would DELETE element i and shift
    ## every later element up, silently misaligning the column.
    diff1 <- vector("list", n)
    cont1 <- vector("list", n)
    for (i in seq_len(n)) {
        if (isRate[i]) next                      # rate: no perturbation
        dc <- .targetDefaultContrast[[reg$short_raw[i]]]
        if (!is.null(dc)) cont1[[i]] <- dc       # contrast-defaulted effect
        else              diff1[[i]] <- 1        # unit step
    }

    ## Assembled as a plain list and stamped afterwards: `$<-` on a data.frame
    ## silently DROPS NULL entries from a list column, which would shorten
    ## every column whose default is NULL.
    ## Default names encode the KIND of quantity, not just the effect:
    ##   <effect>_fd            first difference
    ##   <effect1>_<effect2>_sd second difference (set_second_diff)
    ## so a result list is self-describing.  `name =` overrides per target.
    tg <- list(
        name            = paste0(reg$short_raw, "_fd"),
        effectName1     = reg$short_raw,
        effectType      = reg$effect_type,
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
        ## Selection order.  Starting from the full set, targets keep the
        ## model's effect order; when building up from includeDefaults = FALSE,
        ## each newly included target is appended, so results appear in the
        ## order the user selected them rather than in table order.
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
        type = type, mainEffect = mainEffect, level = level,
        condition = condition, egoNormalize = egoNormalize,
        accumulated = accumulated, rateWeight = rateWeight,
        massContrasts = massContrasts, combineSameLevel = combineSameLevel
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

.targetRow <- function(x, nm) {
    i <- match(nm, x$effectName1)
    if (is.na(i))
        stop("Effect '", nm, "' is not among the targets: ",
             paste(unique(x$effectName1), collapse = ", "), ".", call. = FALSE)
    i
}


##@set_target Postestimation
set_target <- function(x, ...) UseMethod("set_target", x)

##@set_target.sienaPostestTargets Postestimation
set_target.sienaPostestTargets <- function(x, shortNames,
                                           diff = NULL, contrast = NULL,
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
        i <- .targetRow(x, nm)
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


##@set_second_diff Postestimation
set_second_diff <- function(x, ...) UseMethod("set_second_diff", x)

##@set_second_diff.sienaPostestTargets Postestimation
set_second_diff.sienaPostestTargets <- function(x, shortNames,
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
    if (!hasArg(shortNames))
        stop("set_second_diff needs a pair of shortNames.", call. = FALSE)
    nms <- .targetNames(substitute(shortNames), parent.frame())
    if (length(nms) != 2L)
        stop("set_second_diff needs exactly two shortNames, got ",
             length(nms), ".", call. = FALSE)
    if (!is.null(diff1) && !is.null(contrast1))
        stop("Supply either 'diff1' or 'contrast1', not both.", call. = FALSE)
    if (!is.null(diff2) && !is.null(contrast2))
        stop("Supply either 'diff2' or 'contrast2', not both.", call. = FALSE)

    for (nm in nms) .targetRow(x, nm)   # validate both exist

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


##@print.sienaPostestTargets Postestimation
print.sienaPostestTargets <- function(x, ...) {
    inc <- x[x$include, , drop = FALSE]
    d <- attr(x, "defaults")
    cat("RSiena postestimation targets  (",
        nrow(inc), " of ", nrow(x), " selected",
        if (!is.null(attr(x, "depvar")))
            paste0("; depvar '", attr(x, "depvar"), "'") else "",
        ")\n", sep = "")
    if (!is.null(d))
        cat("  defaults : type = ", d$type, ", level = ", d$level,
            if (!is.null(d$condition))
                paste0(", condition = ", paste(d$condition, collapse = "+"))
            else "",
            ", mainEffect = ", d$mainEffect, "\n", sep = "")
    dep <- attr(x, "dependencies")
    if (length(dep)) {
        ## Show the declarations themselves -- what was declared is the whole
        ## point, and a count tells the reader nothing they can check.
        cat("  depends  : ", deparse(dep[[1L]]), "\n", sep = "")
        for (d in dep[-1L])
            cat("             ", deparse(d), "\n", sep = "")
    }
    if (nrow(inc) == 0L) {
        cat("  (none selected)\n")
        return(invisible(x))
    }
    for (i in seq_len(nrow(inc))) {
        shift <- if (!is.null(inc$contrast1[[i]]))
            paste0("contrast = (", paste(inc$contrast1[[i]], collapse = ", "), ")")
        else if (!is.null(inc$diff1[[i]]))
            paste0("diff = ", inc$diff1[[i]])
        else "-"
        lbl <- if (isTRUE(inc$second[i]))
            paste0(inc$effectName1[i], " x ", inc$effectName2[i])
        else inc$effectName1[i]
        ovn <- names(inc$.overrides[[i]])
        cat(sprintf("  %-24s %-24s%s\n", lbl, shift,
                    if (length(ovn))
                        paste0("[", paste(ovn, collapse = ", "), "]") else ""))
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
##
## The level/condition vs accumulated/rateWeight asymmetry is almost certainly
## accidental, but it is current behaviour and equivalence testing depends on
## reproducing it.  Change it deliberately later, not here.
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
            effectName1     = tg$effectName1[i],
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
            e$effectName2     <- tg$effectName2[i]
            e$diff2           <- tg$diff2[[i]]
            e$contrast2       <- tg$contrast2[[i]]
            e$interaction2    <- isTRUE(tg$interaction2[i])
            e$intEffectNames2 <- tg$intEffectNames2[[i]]
            e$modEffectNames2 <- tg$modEffectNames2[[i]]
        }
        ## Overrides are spliced in by name, so a field is "present in the
        ## spec" iff the user set it -- which is what the consumer tests.
        ov <- tg$.overrides[[i]]
        for (f in names(ov)) e[f] <- list(ov[[f]])

        ## massContrasts is resolved PER SPEC, not from a call-level variable:
        ##   eff_massC <- if (!is.null(spec$massContrasts)) spec$massContrasts
        ##                else (perturbType == "ego") || ...
        ## so a model-level setting has to be written into each entry, which is
        ## what the flat single-effect wrapper effectively does.  A per-target
        ## override (already spliced above) wins.
        if (!("massContrasts" %in% names(ov)) &&
            !is.null(d0 <- attr(tg, "defaults")$massContrasts))
            e$massContrasts <- d0

        e <- .applyDependencies(e, attr(tg, "dependencies"))
        out[[k]] <- e
    }
    names(out) <- tg$name[keep]
    ## massContrasts is carried on the entries (above), not as a call-level
    ## variable -- assigning both would double-apply it.
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

    ## A declared-derived effect that is CURRENTLY a selected target is
    ## flagged, not deselected.  make_postest_targets() leaves such an effect
    ## out of the DEFAULT selection, but once the object exists a selection may
    ## be deliberate, and silently dropping it would override an explicit
    ## choice.  So: say what the consequence is and leave the decision.
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
##
## `:` and not `*`, following R's formula algebra: `a * b` expands to
## `a + b + a:b`, i.e. the main effects AND their interaction, which is not
## what a dependency states.  `a:b` is the pure interaction.  Reusing formula
## syntax obliges us to reuse its meaning, so `*` is rejected rather than
## silently read as `:`.
##
## Anything else is rejected rather than approximated -- silently treating an
## unsupported relationship as multiplicative is the failure this layer exists
## to remove.
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
