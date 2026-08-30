## chainAccumulate.R — period-level quantities from simulated ministep chains
##
## Everything else in the postestimation path is ROW-WISE: a flattened chain
## frame goes in, one number per (ministep, alternative) comes out, and the
## aggregation layer averages it.  The quantities here are not row-wise.  They
## are defined per DYAD over a whole period, and getting them requires walking
## the chain in ministep order while carrying state for each dyad.  That is the
## one structural thing this file adds; the utilities, probabilities and
## aggregation it leans on all already exist.
##
## The quantity
## ------------
## Fix a dyad (i,j) and a period.  Let sigma_1 < ... < sigma_K be the ministeps
## at which i has an opportunity and j is a permitted alternative.  Write
##
##   c_k = P(i chooses j | x_ij = 0)      d_k = P(i chooses j | x_ij = 1)
##
## and propagate
##
##   pi_0 = x_ij(t_0),   pi_k = c_k + pi_{k-1} * rho_k,   rho_k = 1 - c_k - d_k
##
## so pi_K is the model-implied probability the tie exists at the end of the
## period.  See docs/design/period_level_tie_probability.md.
##
## Why BOTH branches, every step
## -----------------------------
## c_k and d_k are not two counterfactuals.  They are the two rows of the
## dyad's transition kernel, and both are needed at every step because after
## step 1 the dyad's state is no longer known -- it is pi_{k-1}.  The recursion
## is exactly the mixture:
##
##   pi_k = (1 - pi_{k-1}) c_k + pi_{k-1} (1 - d_k)
##
## Only one of the two is the probability the simulation stored; the other is
## the same decision re-evaluated with the focal tie in the opposite state.
## This is why the one-step readout r_K (which conditions on the REALISED
## state, and so needs only the stored branch) is unbiased while pi_K is not:
## see .accumulatorRegistry below, and section 3 of the design note.
##
## A counterfactual, when there is one, is a separate axis: a statistic pinned
## at a stipulated LEVEL and held there for the whole period.  Each level gets
## its own recursion, and the reported contrast is the difference of endpoints.
## Levels are "branches" in the code below; the two kernel arms are not.


## --------------------------------------------------------------------------
## periodStartNetworks — pi_0 for every period
##
## The chain frame carries no network state, and it does not need to: the
## recursion only ever needs the value at t_0, which is an observed wave.
## RSiena codes structurally determined ties as 10/11; those are still tie
## states for this purpose, so they are folded back to 0/1.  Missing entries
## stay NA and, because NA propagates through the recursion untouched, drop
## out pairwise at the end without needing a mask.
## --------------------------------------------------------------------------
periodStartNetworks <- function(data, depvar) {
    dv <- data[["depvars"]][[depvar]]
    if (is.null(dv))
        stop("periodStartNetworks: no dependent variable '", depvar, "'.",
             call. = FALSE)
    nPeriods <- dim(dv)[3L] - 1L
    lapply(seq_len(nPeriods), function(p) {
        m <- dv[, , p]
        storage.mode(m) <- "double"
        m[!is.na(m) & m >= 10] <- m[!is.na(m) & m >= 10] - 10
        m
    })
}


## --------------------------------------------------------------------------
## periodDirectionDelta — the utility change from conditioning on a tie state
##
## Returns, per arm, the per-row shift that moves a decision to
##
##   arm 0 = tie absent  (the toggle would CREATE)
##   arm 1 = tie present (the toggle would DISSOLVE)
##
## relative to the utility the simulation stored.  This is the only genuinely
## period-specific piece of the perturbation: the ministep marginal effect
## never varies the tie state, it varies statistics.
##
## Direction-dependent theta -- creation vs endowment/maintenance -- is NOT
## handled here.  calculateUtility() already owns that split, so setting the
## density column to the arm's direction and calling it gives the right answer
## for evaluation, creation and endowment specifications alike, with one
## definition of what those types mean.  An earlier version inlined a closed
## form for speed; it was ~3x faster on this step and duplicated exactly the
## logic that would have to be fixed twice.
##
## Computed ONCE per arm and reused across every cell of a contrast, because it
## does not depend on the counterfactual at all -- which more than pays back
## the matrix rebuild that calling calculateUtility() costs.
## --------------------------------------------------------------------------
periodDirectionDelta <- function(cs, thetaEff, u, permitted = NULL) {
    dI <- cs$densityIdx
    lapply(c(1, -1), function(sigma) {
        m <- cs$csMat
        m[, dI] <- sigma
        calculateUtility(m, thetaEff, permitted, dI) - u
    })
}


## --------------------------------------------------------------------------
## branchProbabilities — c_k and d_k for one cell of a counterfactual
##
## `shifts` is a named list of effect -> the shift applied to that effect's
## change statistic, exactly as the per-ministep marginal effect expresses a
## perturbation.  A LEVEL is a shift too (level - factual value), so levels and
## diffs are one mechanism rather than two.  An empty `shifts` is the factual
## condition, where c/d reduce to the stored probability on whichever arm
## matches the row's factual direction -- an identity the tests assert.
##
## Declared `dependencies` are expanded by .shiftSetFor(), the SAME function
## the ministep path uses, so there is one definition of what a perturbation
## does and one place to fix when it is wrong.
##
## The softmax is patched, not rebuilt.  The shift over the ego's choice set
## splits into a part common to every alternative and a spike on the focal one:
##
##   delta_h = a_h                  h != focal   ego-wide: an ego statistic
##   delta_j = a_j + f_j                         plus this arm's tie state on j
##
## and with w_h = p_h e^{a_h}, S = sum_h w_h,
##
##   p'_j = w_j e^{f_j} / (S - w_j + w_j e^{f_j})
##
## which is the ego kernel applied with `a` then the alter kernel with `f` --
## i.e. applyUtilityShift().  perturbType only decides whether `a` is nonzero:
## a dyadic statistic touches the focal alternative alone (a = 0, one alter
## update), an ego statistic touches them all.
## --------------------------------------------------------------------------
branchProbabilities <- function(cc, theta, shifts = list(),
                                dependencies = list(), perturbTypes = NULL,
                                dirDelta = NULL) {
    cs <- cc$changeStats
    if (is.null(cs))
        cs <- contribToChangeStats(cc$contribMat, cc$effectNames,
                                   direction = cc$direction)
    theta_use <- theta[cc$effectNames]
    names(theta_use) <- cc$effectNames
    thetaEff  <- buildThetaEff(theta_use, cs$changeStatsMap)

    dI <- cs$densityIdx
    if (is.null(dI) || is.na(dI))
        stop("branchProbabilities: the period-level recursion needs a density ",
             "column to read the tie direction from.", call. = FALSE)

    u    <- if (!is.null(cc$changeUtility) && !all(is.na(cc$changeUtility)))
                as.numeric(cc$changeUtility)
            else calculateUtility(cs$csMat, thetaEff, cc$permitted, dI)
    p    <- if (!is.null(cc$changeProbability) &&
                !all(is.na(cc$changeProbability)))
                as.numeric(cc$changeProbability)
            else calculateChangeProb(u, cc$group_id)
    dens <- as.numeric(cs$density)
    n    <- length(u)

    ## Independent of the counterfactual, so a caller evaluating several cells
    ## computes it once and passes it in.
    if (is.null(dirDelta))
        dirDelta <- periodDirectionDelta(cs, thetaEff, u, cc$permitted)

    ## ---- the counterfactual, via the SHARED builder ----------------------
    ## Names are resolved once, then handed to .shiftSetFor() -- the same
    ## function the per-ministep marginal effect uses.  Declared dependencies
    ## are expanded THERE and nowhere else, so a dependency fixed for one
    ## quantity is fixed for both, and there is one definition of what a
    ## perturbation does.
    shiftSet <- list()
    if (length(shifts)) {
        for (nm in names(shifts)) {
            resolved <- resolveEffectName(nm, cs$csNames)
            if (length(resolved) != 1L)
                stop("branchProbabilities: '", nm, "' resolves to ",
                     length(resolved), " change-statistic columns; exactly ",
                     "one is required.", call. = FALSE)
            if (!is.null(shiftSet[[resolved]]))
                stop("branchProbabilities: '", resolved, "' is shifted twice; ",
                     "a second difference needs two distinct effects.",
                     call. = FALSE)
            shiftSet[[resolved]] <- shifts[[nm]]
        }
        shiftSet <- .shiftSetFor(shiftSet, dependencies, cs$csMat)
    }

    ## ---- the ego-wide part ------------------------------------------------
    ## A shift on an EGO-specific statistic moves every alternative in the
    ## ego's choice set, not just the focal one, so it enters the block sum.
    ## Evaluated at each row's OWN direction, since only the focal row's tie
    ## state is being varied.  `dens` carries the factor, so this is
    ## identically zero on the no-change alternative -- which matters, because
    ## in the ego kernel every row enters the shared denominator.
    shiftAll <- 0
    if (length(shiftSet)) {
        egoNames <- Filter(function(nm)
            identical(resolvePerturbType(.pinCompositeName(cc, nm),
                                         cc$effectInteractionTypes,
                                         if (is.null(perturbTypes)) NULL
                                         else perturbTypes[[nm]]), "ego"),
            names(shiftSet))
        if (length(egoNames))
            shiftAll <- utilityShift(
                perturbationJacobian(shiftSet[egoNames], dens,
                                     cs$changeStatsMap), theta_use, n)
    }

    onArm <- function(a) {
        sigma <- if (a == 0L) 1 else -1
        ## Two independent pieces: the tie state this arm conditions on
        ## (period-specific), and the counterfactual shift (shared).  The shift
        ## is evaluated at THIS arm's direction, so perturbationJacobian's
        ## creation/endowment masking is applied per arm rather than assumed
        ## to negate.
        du <- dirDelta[[if (a == 0L) 1L else 2L]]
        if (length(shiftSet))
            du <- du + utilityShift(
                perturbationJacobian(shiftSet, rep.int(sigma, n),
                                     cs$changeStatsMap), theta_use, n)
        du[!is.finite(du)] <- NA_real_
        ## shiftAll = 0 makes applyUtilityShift() skip its ego kernel, so a
        ## dyadic shift costs exactly one alter update.
        applyUtilityShift(p, shiftAll = shiftAll, shiftFocal = du - shiftAll,
                          group_id = cc$group_id)
    }
    ## c and d are CHANGE probabilities: P(ego picks j | tie absent) and
    ## P(ego picks j | tie present).  Converting either of them to a tie
    ## probability here would double-count -- the recursion does that itself,
    ## mixing them as (1 - pi)*c + pi*(1 - d).
    list(c = onArm(0L), d = onArm(1L))
}


## --------------------------------------------------------------------------
## .chainKeep — the focal rows
##
## A row is a focal row for dyad (i,j) iff it is a real alternative (not the
## no-change option), the alternative is permitted, and the direction is known.
## Everything downstream assumes this mask has been applied.
## --------------------------------------------------------------------------
.chainKeep <- function(cc) {
    dens <- as.numeric(cc$changeStats$density)
    perm <- if (is.null(cc$permitted)) TRUE else cc$permitted
    !is.na(dens) & dens != 0 & perm & cc$choice != cc$ego
}


## --------------------------------------------------------------------------
## Accumulator registry
##
## An accumulator declares the extra per-dyad state matrices it needs, how to
## update them from one ministep block, and what to emit when the chain ends.
## The recursion state `pi` is always maintained by the core, since everything
## registered here is defined in terms of it.
##
## Adding a quantity that falls out of the same walk -- churn, timing profiles,
## first-passage -- means adding an entry here, not touching the core.
##
##   fields   extra n x n matrices, initialised to 0
##   update   (blk) -> named list of vectors, one per field.  blk carries the
##            block's c, d, dens and the dyad's pi before and after this
##            ministep.  PURE: it never sees the state.
##   combine  "add" (accumulate over ministeps) or "set" (keep the last value)
##   emit     (S, idx) -> named list of vectors, idx a two-column (i,j) matrix
##
## `update` is deliberately a pure function of the block rather than a state
## transformer.  Handing an n x n matrix to a function bumps its reference
## count, so every write inside that function copies the whole matrix; on
## Glasgow-sized chains that one detail was ~90% of the total runtime (measured,
## 10x on a microbenchmark).  Returning values and letting the core write them
## in its own frame keeps the writes in place and the registry abstraction
## intact.
## --------------------------------------------------------------------------
.accumulatorRegistry <- function() list(

    ## The endpoint probability itself.
    piK = list(
        fields = character(0),
        update = function(S, i, j, blk) S,
        emit   = function(S, idx) list(piK = S$pi[idx])
    ),

    ## The one-step readout r_K: the tie probability implied by the ministep
    ## the dyad last faced, conditioned on the state that decision ACTUALLY
    ## faced rather than on a carried probability.  Because
    ## r_k = E[x_ij(sigma_k) | F_{sigma_k - 1}] exactly, E[r_K] is the period-end
    ## tie probability with no approximation -- which makes it the validation
    ## control: it must reproduce fitted(), and it fails loudly if focal-row
    ## selection, direction handling, the K = 0 path or chain alignment are
    ## wrong.  On the factual condition it equals calculateTieProb(p, density).
    rK = list(
        fields  = "r",
        combine = "set",
        update  = function(blk) list(r = ifelse(blk$dens > 0, blk$c, 1 - blk$d)),
        emit    = function(S, idx) list(rK = S$r[idx])
    ),

    ## Expected creations and dissolutions over the period.  These satisfy
    ## nPlus - nMinus = piK - pi_0 by construction, which the tests check --
    ## it is a free consistency check on the recursion.
    counts = list(
        fields  = c("nPlus", "nMinus"),
        combine = "add",
        update  = function(blk) list(nPlus  = (1 - blk$piPrev) * blk$c,
                                     nMinus = blk$piPrev * blk$d),
        emit    = function(S, idx) list(nCreations    = S$nPlus[idx],
                                        nDissolutions = S$nMinus[idx])
    ),

    ## Expected time-on, in ministeps: sum_k pi_k.
    occupancy = list(
        fields  = "occ",
        combine = "add",
        update  = function(blk) list(occ = blk$piNew),
        emit    = function(S, idx) list(occupancy = S$occ[idx])
    )
)


## --------------------------------------------------------------------------
## chainAccumulate — the streaming core
##
## Walks the chain in ministep order holding one n x n state matrix per branch
## per field.  A ministep updates exactly one row of those matrices -- the
## acting ego's -- so the inner operation is a vectorised write of length
## (choice set size), not a per-dyad loop.  Nothing crosses chains, so this
## survives batching unchanged and a chain's state is discarded as soon as it
## ends.
##
## The rows of a flattened chain frame arrive in ministep order already (they
## are emitted that way by the C++ flattener), so no sort is required.  That is
## checked once, cheaply, rather than assumed.
##
## Dyads that never came up (K = 0) are absent from the chain entirely.  By
## default they are not emitted -- their pi_K is x_ij(t_0) by definition, so
## nothing is estimated for them.  includeUnvisited = TRUE seeds every dyad
## from x0 instead, at the cost of an n x n row block per chain.
##
## `branches` is a named list of list(c = , d = ) full-length vectors, one per
## counterfactual condition.  One branch is a prediction; two are a contrast.
## --------------------------------------------------------------------------
chainAccumulate <- function(cc, branches, x0,
                            accumulators = "piK",
                            weights = NULL,
                            keep = NULL,
                            includeUnvisited = FALSE) {

    if (is.null(cc$chain))
        stop("chainAccumulate: the chain frame has no 'chain' column; the ",
             "period-level recursion needs simulated chains (dynamic = TRUE).",
             call. = FALSE)
    if (length(branches) == 0L)
        stop("chainAccumulate: no branches supplied.", call. = FALSE)

    reg <- .accumulatorRegistry()
    unknown <- setdiff(accumulators, names(reg))
    if (length(unknown))
        stop("chainAccumulate: unknown accumulator(s): ",
             paste(unknown, collapse = ", "), ". Available: ",
             paste(names(reg), collapse = ", "), ".", call. = FALSE)
    accs <- reg[accumulators]
    ## Only accumulators that carry state need touching per ministep; piK
    ## declares no fields because the recursion state IS its output, so it
    ## costs nothing beyond its emit().
    accsU  <- accs[vapply(accs, function(a) length(a$fields) > 0L, logical(1L))]
    nAcc   <- length(accsU)
    accAdd <- vapply(accsU, function(a) !identical(a$combine, "set"), logical(1L))
    ## update() returns its fields in declared order, so the inner loop indexes
    ## positionally and never calls names() per ministep.
    accFields <- lapply(accsU, `[[`, "fields")

    if (is.null(keep)) keep <- .chainKeep(cc)
    idx <- which(keep)
    if (length(idx) == 0L)
        stop("chainAccumulate: no focal rows survive the keep mask.",
             call. = FALSE)

    ego    <- cc$ego[idx]
    alt    <- cc$choice[idx]
    chain  <- cc$chain[idx]
    period <- cc$period[idx]
    group  <- if (is.null(cc$group)) rep(1L, length(idx)) else cc$group[idx]
    gid    <- cc$group_id[idx]
    dens   <- as.numeric(cc$changeStats$density)[idx]

    ## Ministep order is what makes the walk correct; assert it rather than
    ## trusting the producer.
    if (is.unsorted(gid))
        stop("chainAccumulate: rows are not in ministep order.", call. = FALSE)

    bn <- names(branches)
    if (is.null(bn) || any(bn == ""))
        stop("chainAccumulate: 'branches' must be named.", call. = FALSE)
    if (!is.null(weights)) {
        if (is.null(names(weights)) || !setequal(names(weights), bn))
            stop("chainAccumulate: 'weights' must be named for exactly the ",
                 "branches supplied.", call. = FALSE)
        weights <- weights[bn]
    }
    bc <- lapply(branches, function(b) as.numeric(b$c)[idx])
    bd <- lapply(branches, function(b) as.numeric(b$d)[idx])

    nA <- nrow(x0[[1L]])
    fields <- unique(unlist(lapply(accs, `[[`, "fields")))

    ## Block boundaries: one block per ministep.  rle() over the already-sorted
    ## group_id is cheaper than split() and keeps the blocks contiguous, which
    ## is what lets the state matrices be written a row at a time.
    r      <- rle(gid)
    ends   <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    nBlk   <- length(starts)

    out       <- vector("list", 0L)
    curChain  <- NA_integer_
    curPeriod <- NA_integer_
    curGroup  <- NA_integer_
    S  <- NULL
    KM <- NULL

    ## Seeding from x0 would otherwise carry the diagonal in as a set of
    ## self-"dyads" with K = 0 and pi_K = 0.  There is no such dyad; the
    ## no-change alternative is a choice, not a tie.
    base0 <- if (includeUnvisited) {
        b0 <- x0
        for (p in seq_along(b0)) diag(b0[[p]]) <- NA_real_
        b0
    } else NULL

    newState <- function(pr) {
        st <- setNames(vector("list", length(bn)), bn)
        for (b in bn) {
            s <- list(pi = if (includeUnvisited) base0[[pr]]
                           else matrix(NA_real_, nA, nA))
            for (f in fields) s[[f]] <- matrix(0, nA, nA)
            st[[b]] <- s
        }
        st
    }

    emitChain <- function(S, KM, ch, pr, gp) {
        ## A dyad is present iff the recursion has a value for it: touched at
        ## least once, or seeded when includeUnvisited.  Structural NA in x0
        ## propagates and drops the dyad here, which is the pairwise deletion
        ## the design note asks for.
        seen <- !is.na(S[[bn[1L]]]$pi)
        if (!any(seen)) return(NULL)
        ij <- which(seen, arr.ind = TRUE, useNames = FALSE)
        ## pi_0 rides along: the change-vs-tie readout is keyed on the dyad's
        ## initial state, and every consumer would otherwise have to look it up
        ## against x0 again.
        res <- list(chain  = rep.int(ch, nrow(ij)),
                    group  = rep.int(gp, nrow(ij)),
                    period = rep.int(pr, nrow(ij)),
                    ego    = ij[, 1L], choice = ij[, 2L],
                    K      = KM[ij],
                    pi0    = x0[[pr]][ij])
        ## With weights, the cells are reduced here: every accumulated quantity
        ## becomes sum_cells weight * value, so a first difference (-1, +1) and
        ## a second difference (+1, -1, -1, +1) differ only in the weights.
        ## Without weights each branch keeps its own column.
        if (is.null(weights)) {
            for (b in bn) for (a in accs) {
                v <- a$emit(S[[b]], ij)
                for (nm in names(v))
                    res[[if (length(bn) > 1L) paste(nm, b, sep = "..")
                         else nm]] <- v[[nm]]
            }
        } else {
            for (a in accs) {
                acc <- NULL
                for (b in bn) {
                    v <- a$emit(S[[b]], ij)
                    if (is.null(acc))
                        acc <- lapply(v, function(x) x * weights[[b]])
                    else
                        for (nm in names(v))
                            acc[[nm]] <- acc[[nm]] + v[[nm]] * weights[[b]]
                }
                for (nm in names(acc)) res[[nm]] <- acc[[nm]]
            }
        }
        ## A plain list, not a data.frame: the caller concatenates column-wise
        ## once at the end.  rbind()-ing one frame per chain is quadratic in
        ## the number of chains and was the second-largest cost at n3 = 1000.
        res
    }

    for (blk in seq_len(nBlk)) {
        s <- starts[blk]; e <- ends[blk]
        rows <- s:e
        ch <- chain[s]; pr <- period[s]

        if (!identical(ch, curChain) || !identical(pr, curPeriod)) {
            if (!is.null(S)) {
                em <- emitChain(S, KM, curChain, curPeriod, curGroup)
                if (!is.null(em)) out[[length(out) + 1L]] <- em
            }
            S  <- newState(pr)
            KM <- matrix(0L, nA, nA)
            x0p <- x0[[pr]]
            curChain <- ch; curPeriod <- pr; curGroup <- group[s]
        }

        ## A ministep is one ego's whole choice set, so every write below is a
        ## slice of ONE row.  M[i, j] with scalar i takes R's fast row path;
        ## the equivalent M[cbind(rep(i, n), j)] goes through generic
        ## matrix-indexing and dominated the runtime on Glasgow-sized chains.
        i  <- ego[s]
        j  <- alt[rows]
        seed <- x0p[i, j]
        KM[i, j] <- KM[i, j] + 1L

        for (b in bn) {
            ck <- bc[[b]][rows]
            dk <- bd[[b]][rows]
            prev <- S[[b]]$pi[i, j]
            na <- is.na(prev)
            if (any(na)) prev[na] <- seed[na]
            new  <- ck + prev * (1 - ck - dk)
            S[[b]]$pi[i, j] <- new
            if (nAcc > 0L) {
                blkCtx <- list(c = ck, d = dk, dens = dens[rows],
                               piPrev = prev, piNew = new)
                for (ai in seq_len(nAcc)) {
                    dv <- accsU[[ai]]$update(blkCtx)
                    fns <- accFields[[ai]]
                    if (accAdd[ai]) {
                        for (fi in seq_along(fns))
                            S[[b]][[fns[fi]]][i, j] <-
                                S[[b]][[fns[fi]]][i, j] + dv[[fi]]
                    } else {
                        for (fi in seq_along(fns))
                            S[[b]][[fns[fi]]][i, j] <- dv[[fi]]
                    }
                }
            }
        }
    }
    if (!is.null(S)) {
        em <- emitChain(S, KM, curChain, curPeriod, curGroup)
        if (!is.null(em)) out[[length(out) + 1L]] <- em
    }

    if (length(out) == 0L) return(NULL)
    nms <- names(out[[1L]])
    res <- setNames(lapply(nms, function(nm)
                           unlist(lapply(out, .subset2, nm),
                                  use.names = FALSE, recursive = FALSE)),
                    nms)
    attr(res, "row.names") <- .set_row_names(length(res[[1L]]))
    class(res) <- "data.frame"
    res
}


## --------------------------------------------------------------------------
## The two predictFun wrappers
##
## These are the shape .evalBatch() expects, and they exploit the branch it
## already has for prediction functions that return a data.frame: the frame is
## handed straight to aggAccum() at the spec's level and condition.  That is
## what lets a per-DYAD quantity ride the per-MINISTEP aggregation machinery
## without changing it -- "egoChoice" already means (period, ego, choice).
##
## Note they ignore `outcomesOnly` and `baseline`.  Both exist to let row-wise
## predictions share a precomputed softmax; there is nothing to share here,
## because the recursion needs the opposite branch that the baseline does not
## carry.
## --------------------------------------------------------------------------

## csName -> the composite effect name resolvePerturbType() keys on.
## changeStats drops the _eval/_endow/_creation suffix that effectNames carry,
## so the two vocabularies have to be joined before the interactionType lookup.
## Falls through to the csName itself, which makes resolvePerturbType() default
## to "alter" -- the safe choice for anything unrecognised.
.pinCompositeName <- function(cc, csName) {
    nms <- names(cc$effectInteractionTypes)
    if (is.null(nms)) return(csName)
    hit <- nms[sub("_(eval|endow|creation)$", "", nms) == csName]
    if (length(hit) == 0L) csName else hit[1L]
}


## --------------------------------------------------------------------------
## .periodCells — a contrast as cells and weights
##
## A sustained counterfactual is a weighted sum over cells, where a cell is one
## complete set of stipulated change statistics held for the whole period.  The
## weights are the ones RSiena already uses for perturbation designs
## (.cellWeights in perturbation.R):
##
##   one effect   base -1, A +1                    the first difference
##   two effects  base +1, A -1, B -1, AB +1       the cross-partial
##
## so the second difference needs no separate machinery.
##
## Each effect is given EITHER a contrast (two stipulated levels) or a diff (a
## sustained shift from whatever the statistic factually is).  Both are emitted
## as SHIFTS, because that is what the shared builder consumes:
##
##   contrast c(a, b)   base: a - s_row     A: b - s_row
##   diff d             base: 0             A: d
##
## Which one is meaningful depends on the statistic.  A level is only defined
## for something that can be HELD -- x_ji for recip, a covariate for altX.  For
## a derived count such as transTrip a level is not a state anyone could
## arrange, while a diff is: "one more two-path than actually available, at
## every decision", which is an inert extra node.  The caller chooses; this
## function only turns the choice into cells.
## --------------------------------------------------------------------------
.periodCells <- function(effectNames, contrasts = NULL, diffs = NULL, csMat) {
    n <- length(effectNames)
    if (n < 1L || n > 2L)
        stop("A period-level contrast takes one or two effects; got ", n, ".",
             call. = FALSE)
    if (is.null(contrasts)) contrasts <- vector("list", n)
    if (is.null(diffs))     diffs     <- vector("list", n)
    if (!is.list(contrasts)) contrasts <- list(contrasts)
    if (!is.list(diffs))     diffs     <- list(diffs)
    length(contrasts) <- n; length(diffs) <- n

    for (k in seq_len(n)) {
        hasC <- !is.null(contrasts[[k]]); hasD <- !is.null(diffs[[k]])
        if (hasC && hasD)
            stop("Effect '", effectNames[k], "': supply either a contrast or ",
                 "a diff, not both.", call. = FALSE)
        if (!hasC && !hasD)
            stop("Effect '", effectNames[k], "' needs either contrast = ",
                 "c(<reference>, <comparison>) -- two levels held for the ",
                 "period -- or diff = <d>, a sustained shift from the factual ",
                 "value.", call. = FALSE)
        if (hasC && length(contrasts[[k]]) != 2L)
            stop("Effect '", effectNames[k], "': contrast must be two levels.",
                 call. = FALSE)
        if (hasD && length(diffs[[k]]) != 1L)
            stop("Effect '", effectNames[k], "': diff must be a single shift.",
                 call. = FALSE)
    }

    w <- .cellWeights(n)
    ## Which effects sit at their COMPARISON state in each cell; the rest sit
    ## at their reference state.  These names are .cellWeights()'s own.
    upIn <- if (n == 1L) list(base = integer(0), A = 1L)
            else list(base = integer(0), A = 1L, B = 2L, AB = c(1L, 2L))

    lapply(names(w), function(cn) {
        up <- upIn[[cn]]
        sh <- lapply(seq_len(n), function(k) {
            if (!is.null(diffs[[k]]))
                if (k %in% up) diffs[[k]] else 0
            else {
                lev <- contrasts[[k]][if (k %in% up) 2L else 1L]
                lev - csMat[, effectNames[k]]
            }
        })
        list(name = cn, weight = unname(w[cn]),
             shifts = setNames(sh, effectNames))
    })
}


## --------------------------------------------------------------------------
## The period-level readout
##
## "tie vs change" and "period vs ministep" are two INDEPENDENT decisions, and
## conflating them was a mistake worth not repeating.  The recursion always
## carries pi_k, the probability the tie EXISTS -- its (1 - d) term is the
## state mixture and is not optional.  What IS a free choice is what gets
## reported at the end:
##
##   tieProb     pi_K                    the tie exists at t1
##   changeProb  pi_K      if x_ij(t0)=0 the tie's state CHANGED over the period
##               1 - pi_K  if x_ij(t0)=1
##
## which is the same transformation calculateTieProb() applies per ministep,
## keyed on the dyad's initial state instead of on the ministep's direction.
##
## For a CONTRAST the weights sum to zero, so 1 - pi maps to -(contrast): the
## change reading is the tie reading with its sign flipped on dyads that start
## with the tie present.  Not a relabelling -- the maintenance stratum enters
## the aggregate with the opposite sign.
## --------------------------------------------------------------------------
.periodReadout <- function(value, pi0, quantity, isContrast) {
    if (identical(quantity, "tieProb")) return(value)
    if (!identical(quantity, "changeProb"))
        stop("Unknown quantity '", quantity, "'; expected tieProb or ",
             "changeProb.", call. = FALSE)
    flip <- !is.na(pi0) & pi0 == 1
    if (isContrast) ifelse(flip, -value, value)
    else            ifelse(flip, 1 - value, value)
}


##@predictPeriodTieProb Postestimation
##
## The period-level LEVEL: no counterfactual, just the recursion on the factual
## condition.  `quantity` chooses the readout -- see .periodReadout().
predictPeriodTieProb <- function(changeContributions, theta, x0,
                                 quantity         = c("tieProb", "changeProb"),
                                 accumulators     = "piK",
                                 includeUnvisited = FALSE,
                                 baseline = NULL, outcomesOnly = FALSE, ...) {
    quantity <- match.arg(quantity)
    br  <- branchProbabilities(changeContributions, theta)
    res <- chainAccumulate(changeContributions, list(factual = br), x0,
                           accumulators     = union("piK", accumulators),
                           includeUnvisited = includeUnvisited)
    res$piK <- .periodReadout(res$piK, res$pi0, quantity, isContrast = FALSE)
    names(res)[names(res) == "piK"] <-
        if (identical(quantity, "tieProb")) "periodTieProb" else "periodChangeProb"
    res
}

##@predictPeriodDiff Postestimation
##
## A SUSTAINED counterfactual: each effect is held at a stipulated level
## (`contrasts`) or shifted from its factual value (`diffs`) at every one of
## ego's opportunities.  One effect gives the first difference, two the
## cross-partial; .periodCells() turns either into cells and weights, so second
## differences need no separate code path.
##
## Every cell is carried through ONE walk of the chain -- one n x n state
## matrix each -- so a second difference costs four matrices, not four passes.
## Every accumulated quantity is contrasted, not just piK, so the contrast in
## expected creations or in occupancy stays consistent with the endpoint.
predictPeriodDiff <- function(changeContributions, theta, x0,
                              effectNames, contrasts = NULL, diffs = NULL,
                              quantity         = c("tieProb", "changeProb"),
                              dependencies     = list(),
                              perturbTypes     = NULL,
                              accumulators     = "piK",
                              includeUnvisited = FALSE,
                              baseline = NULL, outcomesOnly = FALSE, ...) {
    quantity <- match.arg(quantity)
    cs <- changeContributions$changeStats
    if (is.null(cs))
        cs <- contribToChangeStats(changeContributions$contribMat,
                                   changeContributions$effectNames,
                                   direction = changeContributions$direction)
    resolved <- vapply(effectNames, function(nm) {
        r <- resolveEffectName(nm, cs$csNames)
        if (length(r) != 1L)
            stop("'", nm, "' resolves to ", length(r), " change-statistic ",
                 "columns; exactly one is required.", call. = FALSE)
        r
    }, character(1L), USE.NAMES = FALSE)
    cells <- .periodCells(resolved, contrasts, diffs, cs$csMat)

    ## The tie-state conditioning does not depend on the counterfactual, so it
    ## is computed once here rather than once per cell.
    theta_use <- theta[changeContributions$effectNames]
    names(theta_use) <- changeContributions$effectNames
    thetaEff  <- buildThetaEff(theta_use, cs$changeStatsMap)
    u <- if (!is.null(changeContributions$changeUtility) &&
             !all(is.na(changeContributions$changeUtility)))
             as.numeric(changeContributions$changeUtility)
         else calculateUtility(cs$csMat, thetaEff,
                               changeContributions$permitted, cs$densityIdx)
    dd <- periodDirectionDelta(cs, thetaEff, u, changeContributions$permitted)

    branches <- setNames(
        lapply(cells, function(cl)
               branchProbabilities(changeContributions, theta, cl$shifts,
                                   dependencies = dependencies,
                                   perturbTypes = perturbTypes,
                                   dirDelta     = dd)),
        vapply(cells, `[[`, character(1L), "name"))
    weights <- setNames(vapply(cells, `[[`, numeric(1L), "weight"),
                        names(branches))

    res <- chainAccumulate(changeContributions, branches, x0,
                           accumulators     = union("piK", accumulators),
                           weights          = weights,
                           includeUnvisited = includeUnvisited)
    res$piK <- .periodReadout(res$piK, res$pi0, quantity, isContrast = TRUE)
    outName <- if (length(effectNames) == 1L) "periodFirstDiff"
               else "periodSecondDiff"
    names(res)[names(res) == "piK"] <- outName
    for (nm in setdiff(names(res), c("chain", "group", "period", "ego",
                                     "choice", "K", "pi0", outName)))
        names(res)[names(res) == nm] <- paste0(nm, "Diff")
    res
}
