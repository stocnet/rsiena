## perturbation.R — a perturbation expressed where it acts
##
## A marginal effect perturbs change statistics and asks what that does to the
## choice probability.  Two spaces, and the mathematics is very different in
## each -- conflating them is easy and was the source of several errors here.
##
## ---------------------------------------------------------------------------
## What the counterfactual IS -- and what it deliberately leaves out
## ---------------------------------------------------------------------------
## The change statistics are not independent quantities.  In the data they are
## all functions of the same network x, so they co-vary whether or not anyone
## perturbs anything.  A counterfactual that changed x would move many of them
## at once, in ways specific to the network.
##
## What is computed here is a SIMPLIFICATION of that, and it is a modelling
## choice rather than a derivation:
##
##   perturb in delta-space, not in x-space
##   DO propagate DIRECT functional dependence -- where one delta is a
##      mathematical function of others (a derived statistic: unspInt from its
##      components, and later gw(), sameX, ...)
##   DO NOT propagate indirect co-variation -- changing reciprocity in a real
##      network would also move transitivity; that is not followed
##
## So the delta vector a cell describes need not correspond to any network that
## could exist.  That is intended, and it is what makes the quantity a
## ceteris-paribus marginal effect rather than a simulation.
##
## Once that choice is made, the arithmetic is exact rather than approximate:
## u = sum_k theta_k delta_k is linear in the delta VECTOR, so for whatever
## Ddelta a cell specifies,
##
##     Du = Ddelta . theta            (1)   exact, no truncation
##     d(Du)/d(theta) = Ddelta        (2)   since (1) is linear in theta
##
## Ddelta is NOT a derivative -- it is a coefficient matrix that happens to be
## its own Jacobian.  Note what (1) and (2) do NOT claim: Ddelta is generally
## non-linear in the size of the perturbation, because direct dependencies make
## it so (for A = B:C perturbed on both sides, Ddelta_A carries a dB*dC term).
## The non-linearity is upstream, in constructing Ddelta, and never reaches (1)
## or (2).  That is why `.dependencyShift()` EVALUATES a relation at the
## perturbed state instead of differentiating it: evaluation is exact for any
## relation however non-linear.
##
## ---------------------------------------------------------------------------
## THE INVARIANT: one Ddelta, shared
## ---------------------------------------------------------------------------
## Because the counterfactual is a choice, "correct" means CONSISTENT: the
## point estimate, the Jacobian and the delta-method SE must describe the same
## counterfactual.  If the Jacobian propagates a dependency that the point
## estimate does not, the SE is the standard error of a different quantity --
## and it will look entirely reasonable.
##
## That is not hypothetical.  It is exactly how the ego-perturbation SE came to
## be twice its true value: the point estimate used the ego-wide update and the
## Jacobian used the one-alternative derivative, and nothing checked they
## agreed.
##
## Hence: Ddelta is built ONCE, by `.shiftSetFor()`, and every consumer takes
## it from there.  Any new consumer must do the same rather than reconstruct
## it.
##
## ---------------------------------------------------------------------------
## In PROBABILITY space: non-linear, so they do NOT coincide
## ---------------------------------------------------------------------------
## p = softmax(u) is not linear in u, so a finite utility shift gives a genuine
## discrete difference Dp, while d(Dp)/d(theta) needs the chain rule through
## the update.  Nothing is shared between the two there, and that is where the
## real derivation work sits -- `firstDiffJacobian()` and its ego-wide
## counterpart, the HARNESS.
##
## The split to keep in mind:
##
##     CORE     A = Ds        which statistics move, and by how much.
##                            Effect-dependent.  Identical for every
##                            perturbation mode, because it is about the
##                            counterfactual, not about the probability model.
##
##     HARNESS  chain A through the update's derivative.  This is where
##                            one-alternative and ego-wide differ, and the only
##                            place they differ.
##
## ---------------------------------------------------------------------------
## Why A needs more columns than Ds has statistics
## ---------------------------------------------------------------------------
## csMat holds one column per BASE effect, in eval-space; creation and
## endowment live in the theta mapping rather than in extra columns
## (`buildThetaEff`):
##
##     thetaEff[k, "creation"]    = theta_eval + theta_creation
##     thetaEff[k, "dissolution"] = theta_eval + theta_endow
##
## and `calculateUtility` picks the column per row via the direction d.  So in
## (2) the derivative is taken with respect to the ESTIMATED parameters, and
## perturbing one base effect moves several of them -- which ones depending on
## the row: a creating ministep responds to theta_creation, a dissolving one to
## theta_endow.  A therefore has one column per estimated parameter with a
## row-dependent mask, which a single column index and a scalar multiplier
## cannot express.


## --------------------------------------------------------------------------
## .parameterMask — which rows an estimated parameter acts on
##
## eval acts on every row that changes anything; creation only where a tie is
## created, endowment only where one is dropped.  Rows with d = 0 are the
## no-change alternative and contribute to nothing.
##
## The mask comes from the model's semantics via the direction vector, NOT
## from inspecting which change statistics happen to be non-zero: a change
## statistic may legitimately be zero, so reading the mask off the data would
## agree with this definition on most fixtures and diverge silently on others.
## --------------------------------------------------------------------------
.parameterMask <- function(type, direction) {
    switch(type,
           eval     = direction != 0L,
           creation = direction ==  1L,
           endow    = direction == -1L,
           stop("Unknown effect type '", type, "'; expected eval, creation ",
                "or endow.", call. = FALSE))
}


## --------------------------------------------------------------------------
## .effectColumns — A-columns for shifting ONE effect's change statistic
##
## Split out so that perturbationJacobian() below is purely about combining
## shifts: this is the only place that knows how an estimated parameter maps
## onto rows.
## --------------------------------------------------------------------------
.effectColumns <- function(effectName, shift, direction, csMap, knownDirection) {
    members <- which(csMap$bases == effectName)
    if (length(members) == 0L) return(list())
    ## NA is NOT zeroed here.  It marks a row the perturbation does not reach
    ## -- a contrast row outside the contrast's range -- and must propagate
    ## through Du so the row drops out of the aggregation.  Zeroing it would
    ## return the row unchanged: a valid-looking number that stays in the
    ## average and pulls it toward zero.

    if (!knownDirection) {
        ## No direction to condition on: utility is a plain inner product
        ## against the creation column, so every row responds identically and
        ## only eval parameters are in play.  A creation/endow effect here
        ## would be unrepresentable, and contribToChangeStats() already raises
        ## that case rather than letting it through.
        keep <- members[csMap$types[members] == "eval"]
        return(lapply(keep, function(j) list(col = j, val = shift)))
    }

    out <- vector("list", length(members))
    n   <- 0L
    for (j in members) {
        mask <- .parameterMask(csMap$types[j], direction)
        ## isTRUE(): an NA direction makes any() return NA, and `if (!NA)` is
        ## an error rather than a decision.
        if (!isTRUE(any(mask, na.rm = TRUE))) next
        ## direction carries the sign with which the change statistic enters
        ## utility, so it multiplies the shift as well.
        n <- n + 1L
        out[[n]] <- list(col = j, val = direction * shift * as.numeric(mask))
    }
    if (n == 0L) list() else out[seq_len(n)]
}


## --------------------------------------------------------------------------
## perturbationJacobian — A for a set of change-statistic shifts
##
##   shifts     named list: effect name -> the shift applied to THAT effect's
##              change statistic (scalar or one value per row)
##   direction  the d vector: +1 creating, -1 dissolving, 0 no change; NULL or
##              all-NA when the model carries no density effect and no
##              constraint-derived direction
##   csMap      changeStatsMap from contribToChangeStats()
##
## Taking a SET of shifts rather than one effect plus interaction arguments is
## what makes this general.  A perturbation of one effect is a one-element
## set; a declared dependency simply adds an entry, and so would any future
## dependency form:
##
##   plain            list(recip = delta)
##   A = B:C          list(B = delta, A = delta * s_C)
##   A = f(B)         list(B = delta, A = f(s_B + delta) - f(s_B))
##
## All three have the same shape because SAOM is linear in the parameters:
## however non-linear a dependency is in the STATISTICS, the theta-derivative
## is still just the shift, so non-linearity is absorbed when the entry is
## computed and never reaches this structure.
##
## Returns a list of list(col = <index into theta>, val = <one per row>),
## with contributions to the same parameter summed.  Sparse by construction:
## only parameters actually reached appear.
## --------------------------------------------------------------------------
perturbationJacobian <- function(shifts, direction, csMap) {
    knownDirection <- !is.null(direction) && !all(is.na(direction))

    cols <- list()
    for (nm in names(shifts)) {
        for (a in .effectColumns(nm, shifts[[nm]], direction, csMap,
                                 knownDirection)) {
            key <- as.character(a$col)
            ## Two shifts may reach the same parameter -- an effect that is
            ## both perturbed and depended upon, say.  They add; keeping them
            ## as separate entries would double-count in the consumer.
            cols[[key]] <- if (is.null(cols[[key]])) a
                           else list(col = a$col, val = cols[[key]]$val + a$val)
        }
    }
    unname(cols)
}


## --------------------------------------------------------------------------
## utilityShift — Du implied by a perturbation, given theta
##
## The same A that serves as the Jacobian also produces the utility shift,
## which is the point: Du = A . theta.  Computing them from one object is what
## keeps them consistent by construction rather than by matching two
## derivations.
## --------------------------------------------------------------------------
utilityShift <- function(A, theta, nrow) {
    du <- numeric(nrow)
    for (a in A) du <- du + a$val * theta[a$col]
    du
}


## --------------------------------------------------------------------------
## Cells — the counterfactual states a target compares
##
## A target is a set of cells and a weight per cell.  The estimand is the
## weighted sum of the outcome across them:
##
##   first difference    {base, A}              -1, +1
##   second difference   {base, A, B, AB}       +1, -1, -1, +1
##   third order         eight cells            alternating
##
## The code already evaluates these states -- a second difference computes the
## first difference of effect1, shifts effect2, and computes it again, which is
## four states -- but it does so by patching increments onto intermediate
## values.  Naming the states is what removes that: a cell says what the change
## statistics ARE, so nothing can be stale, and the same object yields both the
## point estimate and the derivative.
##
## A cell holds SHIFTS from the base state rather than absolute values, because
## that is what both the utility shift and the Jacobian need, and because it
## keeps the base state untouched (see .dependencyShift below for why that
## matters).
## --------------------------------------------------------------------------

## The weights that turn cell outcomes into the estimand.  Sign pattern only;
## the caller decides what "outcome" means.
.cellWeights <- function(nEffects) {
    if (nEffects == 1L) return(c(base = -1, A = 1))
    if (nEffects == 2L) return(c(base = 1, A = -1, B = -1, AB = 1))
    stop("Only first and second differences are supported; got ", nEffects,
         " perturbed effects.", call. = FALSE)
}


## --------------------------------------------------------------------------
## .dependencyShift — what a declared dependency does in a given cell
##
## For `A ~ B:C`, A's change statistic follows the product of B's and C's.  In
## a cell where B and C are at s_B + dB and s_C + dC, the induced shift is
##
##   (s_B + dB)(s_C + dC) - s_B * s_C  =  dB*s_C + dC*s_B + dB*dC
##
## including the dB*dC cross term, which appears in the AB cell of a second
## difference and nowhere else.  Evaluating the relation AT THE CELL'S STATE is
## what produces it; no separate cross-term rule is needed.
##
## Expressed as a DELTA against the base product rather than by overwriting
## A's column with the product.  The declaration is the caller's assertion
## about how A moves, not a claim that A equals the product exactly in the
## data; overwriting would silently change the base state wherever the two
## disagree.  As a delta, a base that does not match is left alone and only the
## induced change is applied.
##
## Non-multiplicative forms slot in here and nowhere else: `sameX ~ egoX ==
## altX` or `gw(transTrip)` differ only in the function evaluated at the two
## states.  The Jacobian is unaffected either way, because the relation
## involves no parameters.
## --------------------------------------------------------------------------
.dependencyShift <- function(terms, baseStats, cellShifts) {
    shifted <- lapply(terms, function(tm) {
        s <- baseStats[[tm]]
        if (is.null(s))
            stop("Dependency term '", tm, "' has no change-statistic column.",
                 call. = FALSE)
        s + (if (is.null(cellShifts[[tm]])) 0 else cellShifts[[tm]])
    })
    base <- Reduce(`*`, lapply(terms, function(tm) baseStats[[tm]]))
    Reduce(`*`, shifted) - base
}


## --------------------------------------------------------------------------
## perturbationCells — the cells of a target, each as a shift set
##
##   perturbations  named list, one entry per perturbed effect: the per-row
##                  shift applied to that effect's change statistic
##   baseStats      named list of the base change statistics, needed to
##                  evaluate dependencies at each cell's state
##   dependencies   list of list(target = <effect>, terms = <effects>)
##
## Returns a list of cells: list(name, weight, shifts), where `shifts` is in
## exactly the form perturbationJacobian() consumes.
## --------------------------------------------------------------------------
perturbationCells <- function(perturbations, baseStats = list(),
                              dependencies = list()) {
    nms <- names(perturbations)
    w   <- .cellWeights(length(nms))

    ## Which perturbed effects are active in each cell.
    ## Written out rather than switch(length(nms), ...): switch() on a number
    ## selects by POSITION and ignores the names, so the branch labels would
    ## be decorative and correct only by accident of ordering.
    members <- if (length(nms) == 1L)
        list(base = character(0), A = nms[1L])
    else
        list(base = character(0), A = nms[1L], B = nms[2L], AB = nms)

    lapply(names(w), function(cellName) {
        active <- members[[cellName]]
        shifts <- perturbations[active]

        ## Dependencies are evaluated at THIS cell's state, which is what makes
        ## the cross term fall out rather than needing to be added.
        for (dep in dependencies) {
            d <- .dependencyShift(dep$terms, baseStats, shifts)
            ## `all(d == 0)` is NA as soon as any element is NA, and `if (NA)`
            ## errors.  A shift that is zero or unreachable everywhere moves
            ## nothing, so skipping it is right; NA must not turn that into a
            ## crash.
            if (all(d == 0 | is.na(d))) next
            shifts[[dep$target]] <-
                (if (is.null(shifts[[dep$target]])) 0 else shifts[[dep$target]]) + d
        }
        list(name = cellName, weight = unname(w[cellName]), shifts = shifts)
    })
}


## --------------------------------------------------------------------------
## applyUtilityShift — the incremental multinomial-logit update
##
## Implements eq. (mlogitupdate) of the paper (Ben-Akiva & Lerman 1985, p.114),
## which is valid for ANY linear-in-parameters preference shift:
##
##   p_ij(x; Df_i) = p_ij(x) exp(Df_ij) / sum_h p_ih(x) exp(Df_ih)
##
## There is ONE formula.  The two code paths that exist today are not two
## models; they are this formula specialised to the two shapes the shift takes:
##
##   alter-specific (dyadic) effect
##       The counterfactual perturbs ONE alternative, and each row carries its
##       own counterfactual -- n of them, not one.  Df is nonzero only at the
##       focal alternative, the denominator collapses to 1 - p + p*exp(Df), and
##       no group sum is needed.
##
##   ego-specific effect
##       Perturbing an ego attribute is not coherent for a single alter: the
##       change statistic is constant across alters, so the counterfactual is
##       an ego-wide shift -- all creation alternatives up by the same amount,
##       all dissolution alternatives down by it, the no-change option
##       unchanged.  One counterfactual per ego, and the denominator needs the
##       full sum.
##
## Splitting the shift into the part common to every alternative of the ego and
## the part specific to the focal alternative covers both, AND the mixed case a
## second difference can produce (an ego-wide shift plus a dyadic increment on
## one alternative), which neither existing path expresses:
##
##   w_h    = p_h * exp(shiftAll_h)          the ego-wide part
##   S_i    = sum_h w_h                      once per ego
##   p'_j   = w_j*exp(shiftFocal_j) / ( S_i + w_j*(exp(shiftFocal_j) - 1) )
##
## The denominator is the group sum with the focal alternative's own term
## replaced by its shifted value, which is what makes one pass per ego enough
## rather than one per counterfactual.
##
## It reduces to both existing paths exactly, and both reductions are tested:
##   shiftAll   = 0  ->  S_i = 1, denominator = 1 - p + p*exp(shiftFocal)
##   shiftFocal = 0  ->  p'_j = w_j / S_i
## --------------------------------------------------------------------------
applyUtilityShift <- function(p, shiftAll = 0, shiftFocal = 0,
                              group_id = NULL) {
    n  <- length(p)
    sa <- rep_len(shiftAll,   n)
    sf <- rep_len(shiftFocal, n)

    ## The unified formula is EXACTLY the two existing kernels composed, in
    ## this order:
    ##
    ##   q_j  = w_j / S_i                          (ego kernel, shiftAll)
    ##   p'_j = q_j e^sf / (1 - q_j + q_j e^sf)    (alter kernel, shiftFocal)
    ##        = w_j e^sf / (S_i + w_j (e^sf - 1))
    ##
    ## so no new kernel is needed and both compiled paths are reused as they
    ## are.  Verified against the closed form to 3e-15.
    ##
    ## The ORDER is forced, not chosen: the ego-wide part is common to every
    ## alternative and so belongs inside every denominator, which is only true
    ## if it is applied first.  Reversing it is wrong by ~2e-2, not by
    ## rounding.
    ##
    ## Each branch below costs exactly what it costs today -- one compiled
    ## call for a purely dyadic or purely ego-wide shift.  Only the mixed
    ## case, which no current path can express at all, pays for two.
    egoWide  <- any(sa != 0, na.rm = TRUE)
    focal    <- any(sf != 0, na.rm = TRUE) || !egoWide

    out <- p
    if (egoWide) out <- mlogit_update_r(out, sa, group_id, "ego")
    if (focal)   out <- mlogit_update_r(out, sf, group_id, "alter")
    out
}
