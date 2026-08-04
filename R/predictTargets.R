## predictTargets.R — target specification for predicted probabilities
##
## The same idea as postestTargets.R, and much smaller.  A marginal-effect
## target says WHICH EFFECT to perturb and HOW; a prediction has nothing to
## perturb, so a prediction target says only WHAT TO CONDITION ON and AT WHAT
## LEVEL.  One row per requested stratification:
##
##   ptg <- make_predict_targets(ans, effects = myeff, type = "tieProb")
##   ptg <- set_condition(ptg, recip)                       # by reciprocation
##   ptg <- set_condition(ptg, density, level = "egoChoice")# ego tie profile
##
## Everything about what the estimand IS -- type, dynamic, level, condition,
## egoNormalize, accumulated, rateWeight -- lives here, exactly as it does on
## the marginal-effects object, and per-row overrides use the same `.overrides`
## convention: a name is present in that list iff the user set it for that row,
## which is what distinguishes "not set" from "set to NULL".
##
## The estimand rules (accumulated needs dynamic, rateWeight is absorbed by a
## dynamic run) are enforced by .checkEstimandCombo(), shared with
## make_marginal_targets() so the two entry points cannot drift apart.


##@make_predict_targets Postestimation
make_predict_targets <- function(x, ...) UseMethod("make_predict_targets", x)

##@make_predict_targets.sienaFit Postestimation
make_predict_targets.sienaFit <- function(x, data = NULL, effects = NULL,
                                          depvar = NULL,
                                          type = c("changeProb", "tieProb"),
                                          dynamic = FALSE,
                                          level = "period",
                                          condition = NULL,
                                          egoNormalize = TRUE,
                                          accumulated = FALSE,
                                          rateWeight = FALSE,
                                          na.rm = TRUE, ...) {
    type <- match.arg(type)
    .checkEstimandCombo(dynamic, accumulated, rateWeight, "prediction")
    .checkLevel(level)

    if (is.null(effects))
        stop("'effects' must be supplied -- the effects object the model was ",
             "fitted with, e.g. effects = mymodel.", call. = FALSE)
    if (!inherits(effects, "sienaEffects"))
        stop("'effects' must be a sienaEffects object.", call. = FALSE)
    if (is.null(depvar)) depvar <- unique(effects[["name"]])[1L]

    ## One row, carrying the model-level request.  set_condition() appends.
    tg <- list(name       = "prediction",
               condition  = I(list(condition)),
               .seq       = 1,
               .overrides = I(list(list())))
    attr(tg, "row.names") <- .set_row_names(1L)
    attr(tg, "depvar")    <- depvar
    attr(tg, "effects")   <- effects
    attr(tg, "defaults")  <- list(
        type = type, dynamic = dynamic, level = level,
        condition = condition, egoNormalize = egoNormalize,
        accumulated = accumulated, rateWeight = rateWeight, na.rm = na.rm)
    attr(tg, "version") <- utils::packageDescription("RSiena", fields = "Version")
    class(tg) <- c("sienaPredictTargets", "data.frame")
    tg
}


##@set_condition Postestimation
set_condition <- function(x, ...) UseMethod("set_condition", x)

##@set_condition.sienaPredictTargets Postestimation
##
## Appends a row rather than retuning one: several stratifications of the same
## prediction are the normal case, and are what the old interface could only
## get by calling predict() once per condition.
##
## `condition` takes bare effect short names, as set_target() does.  An
## explicit NULL asks for the unconditional prediction.
set_condition.sienaPredictTargets <- function(x, condition, level,
                                              egoNormalize, accumulated,
                                              rateWeight,
                                              name = NULL, verbose = TRUE,
                                              ...) {
    cond <- if (missing(condition)) NULL
            else .targetNames(substitute(condition), parent.frame())
    ## `.targetNames` deparses a bare NULL to the string "NULL"; that is the
    ## request for an unconditional row, not a condition on an effect called
    ## NULL.
    if (identical(cond, "NULL")) cond <- NULL

    ov <- list()
    for (f in c("level", "egoNormalize", "accumulated", "rateWeight"))
        if (!missing_arg_named(f, environment())) ov[f] <- list(get(f))
    if (!is.null(ov$level)) .checkLevel(ov$level)

    if (is.null(name))
        name <- if (is.null(cond)) "unconditional"
                else paste0("by_", paste(cond, collapse = "_"))
    if (name %in% x$name)
        stop("A prediction target named '", name, "' already exists.",
             call. = FALSE)

    row <- list(name = name, condition = I(list(cond)),
                .seq = max(c(x$.seq, 0)) + 1, .overrides = I(list(ov)))
    out <- rbind_predict_rows(x, row)
    if (verbose)
        message("set_condition: added '", name, "'",
                if (length(ov))
                    paste0(" (", paste(names(ov), collapse = ", "),
                           " overridden)") else "")
    out
}

## rbind() on a data.frame with list-columns drops the I() wrapper and can
## unlist a single-element list, so the rows are stitched by hand.
rbind_predict_rows <- function(x, row) {
    out <- list(name      = c(x$name, row$name),
                condition = c(unclass(x$condition), unclass(row$condition)),
                .seq      = c(x$.seq, row$.seq),
                .overrides = c(unclass(x$.overrides), unclass(row$.overrides)))
    out$condition  <- I(out$condition)
    out$.overrides <- I(out$.overrides)
    attributes(out) <- c(attributes(out),
                         attributes(x)[setdiff(names(attributes(x)),
                                               c("names", "row.names", "class"))])
    attr(out, "row.names") <- .set_row_names(length(out$name))
    class(out) <- c("sienaPredictTargets", "data.frame")
    out
}


##@print.sienaPredictTargets Postestimation
print.sienaPredictTargets <- function(x, ...) {
    d <- attr(x, "defaults")
    cat("RSiena prediction targets  (", nrow(x), " requested; depvar '",
        attr(x, "depvar"), "')\n", sep = "")
    cat("  quantity : ",
        if (identical(d$type, "tieProb"))
            "probability the tie EXISTS after the ministep"
        else "probability ego CHANGES the tie",
        if (isTRUE(d$dynamic)) ", dynamic (simulated chains)"
        else ", static (observed waves)", "\n", sep = "")
    cat("  reported : by ", d$level,
        if (isTRUE(d$accumulated)) ", accumulated over ministeps" else "",
        "\n\n", sep = "")

    for (i in order(x$.seq)) {
        ov   <- x$.overrides[[i]]
        cond <- x$condition[[i]]
        desc <- if (is.null(cond)) "unconditional"
                else paste("conditional on", paste(cond, collapse = " + "))
        extra <- if (length(ov))
            paste0("  [", paste(names(ov), "=",
                                vapply(ov, function(v)
                                    if (is.null(v)) "NULL"
                                    else paste(v, collapse = ","),
                                    character(1L)),
                                collapse = ", "), "]") else ""
        cat(sprintf("  %-22s %s%s\n", x$name[i], desc, extra))
    }
    invisible(x)
}


## Lower to the per-row settings the prediction path consumes: model defaults
## with that row's overrides applied.  Same resolution order as the
## marginal-effects object.
.predictRowSettings <- function(x, i) {
    d  <- attr(x, "defaults")
    ov <- x$.overrides[[i]]
    d$condition <- x$condition[[i]]
    for (nm in names(ov)) d[[nm]] <- ov[[nm]]
    d
}
