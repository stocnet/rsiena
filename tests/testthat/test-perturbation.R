# test-perturbation.R
#
# The perturbation representation (R/perturbation.R): Ds in change-statistic
# space, A in parameter space, Du = A . theta.
#
# Correctness is checked DIRECTLY rather than against another implementation:
# the Du implied by A must equal the utility difference obtained by actually
# shifting the change statistic and recomputing utility.  That check shares no
# code with the representation, so agreement is evidence rather than tautology,
# and it is the only kind of check that would catch a wrong mask.

cs_fixture <- function(fitName, dataName, effName, depvar) {
  ans <- load_fixture(fitName); dat <- load_fixture(dataName)
  eff <- load_fixture(effName)
  if (is.null(ans) || is.null(dat) || is.null(eff)) return(NULL)
  w <- RSiena:::getStaticChangeContributions(ans = ans, data = dat,
         effects = eff, depvar = depvar, returnWide = TRUE)
  cs <- RSiena:::contribToChangeStats(w$contribMat, w$effectNames,
                                      theta = coef(ans))
  ## Effects are addressed by their RESOLVED (depvar-prefixed) name here --
  ## "mynet_cm_recip", not "recip" -- because that is what csMap$bases holds
  ## and what resolveEffectName() produces in the calling path.
  resolve <- function(short)
    grep(paste0("(^|_)", short, "$"), unique(cs$changeStatsMap$bases),
         value = TRUE)[1L]
  list(fit = ans, wide = w, cs = cs, theta = coef(ans), resolve = resolve)
}


test_that("A is zero on no-change rows", {
  skip_on_cran()
  fx <- cs_fixture("ans_cm", "mydata_cm", "mymodel_cm", "mynet_cm")
  skip_if(is.null(fx), "cm fixtures unavailable")

  d <- fx$cs$density
  skip_if(!any(d == 0L), "fixture has no no-change rows")

  ## A no-change alternative contributes nothing to utility, so it must
  ## contribute nothing to the derivative either. Broadcasting the shift
  ## across all rows would pass every other test in this file.
  A <- RSiena:::perturbationJacobian(
           setNames(list(1), fx$resolve("recip")), d, fx$cs$changeStatsMap)
  expect_gt(length(A), 0L)
  for (a in A) expect_true(all(a$val[d == 0L] == 0))
})


test_that("creation and endowment parameters are reached, on their own rows only", {
  skip_on_cran()
  fx <- cs_fixture("ans_cm", "mydata_cm", "mymodel_cm", "mynet_cm")
  skip_if(is.null(fx), "cm fixtures unavailable")

  csMap <- fx$cs$changeStatsMap
  d     <- fx$cs$density

  ## recip carries eval + creation, transTrip carries eval + endow.
  A_recip <- RSiena:::perturbationJacobian(setNames(list(1), fx$resolve("recip")), d, csMap)
  types <- vapply(A_recip, function(a) csMap$types[a$col], character(1L))
  expect_setequal(types, c("eval", "creation"))

  cre <- A_recip[[which(types == "creation")]]
  ## Creation acts only where a tie is created; nowhere else.
  expect_true(all(cre$val[d != 1L] == 0))
  expect_true(any(cre$val[d == 1L] != 0))

  A_tt  <- RSiena:::perturbationJacobian(setNames(list(1), fx$resolve("transTrip")), d, csMap)
  types_tt <- vapply(A_tt, function(a) csMap$types[a$col], character(1L))
  expect_setequal(types_tt, c("eval", "endow"))
  end <- A_tt[[which(types_tt == "endow")]]
  expect_true(all(end$val[d != -1L] == 0))
  expect_true(any(end$val[d == -1L] != 0))
})


test_that("Du from A equals the utility difference from actually shifting the statistic", {
  skip_on_cran()
  for (nm in list(c("ans", "mydata", "mymodel", "mynet"),
                  c("ans_cm", "mydata_cm", "mymodel_cm", "mynet_cm"))) {
    fx <- cs_fixture(nm[1], nm[2], nm[3], nm[4])
    if (is.null(fx)) next

    csMap <- fx$cs$changeStatsMap
    d     <- fx$cs$density
    csMat <- fx$cs$csMat
    thEff <- RSiena:::buildThetaEff(fx$theta, csMap)
    dIdx  <- fx$cs$densityIdx

    ## density excluded deliberately: perturbing it changes the DIRECTION of
    ## the ministep, not just a change statistic, so it is not a shift in this
    ## representation at all. It gets its own route later in step 6.
    effs <- unique(csMap$bases)
    effs <- effs[!grepl("density", effs, fixed = TRUE)]
    for (eff in effs) {
      A <- RSiena:::perturbationJacobian(setNames(list(1), eff), d, csMap)
      if (length(A) == 0L) next
      du_A <- RSiena:::utilityShift(A, fx$theta, length(d))

      ## The independent check: shift the change statistic itself and let the
      ## existing utility function say what happened. This shares no code with
      ## perturbationJacobian(), so agreement is evidence rather than tautology.
      k <- match(eff, colnames(csMat))
      shifted <- csMat
      shifted[, k] <- shifted[, k] + 1
      du_direct <-
        RSiena:::calculateUtility(shifted, thEff, NULL, dIdx) -
        RSiena:::calculateUtility(csMat,   thEff, NULL, dIdx)

      expect_equal(du_A, du_direct, tolerance = 1e-12,
        info = paste0(nm[1], " / ", eff,
                      ": A implies a different utility shift than shifting ",
                      "the change statistic does"))
    }
  }
})


# ── Dependencies: an effect that moves because another one did ───────────────
#
# A declared dependency A = B:C means perturbing B by delta also moves A's
# change statistic, by delta * s_C. In this representation that is not a
# special case: it is simply a second entry in the shift set. These check that
# the entry lands on A's parameters, with the right per-row multiplier, and
# that the resulting Du still matches what shifting the statistics directly
# produces -- the same independent check used for the plain cases.

test_that("a dependency contributes its own parameter columns", {
  skip_on_cran()
  fx <- cs_fixture("ans2", "mydata2", "mymodel2", "mynet2")
  skip_if(is.null(fx), "ans2 fixtures unavailable")

  csMap <- fx$cs$changeStatsMap
  d     <- fx$cs$density
  csMat <- fx$cs$csMat

  b   <- fx$resolve("transTrip")     # perturbed
  dep <- fx$resolve("transRecTrip")  # moves with it
  mod <- fx$resolve("recip")         # the moderator
  skip_if(any(is.na(c(b, dep, mod))), "model lacks the interaction trio")

  delta <- 1
  plain <- RSiena:::perturbationJacobian(setNames(list(delta), b), d, csMap)
  withDep <- RSiena:::perturbationJacobian(
      setNames(list(delta, delta * csMat[, mod]), c(b, dep)), d, csMap)

  ## The dependency reaches strictly more parameters than the bare shift.
  expect_gt(length(withDep), length(plain))
  depCols <- which(csMap$bases == dep)
  expect_true(any(vapply(withDep, function(a) a$col %in% depCols, logical(1L))))

  ## and the column it adds is delta * s_mod, carrying the direction sign.
  a_dep <- Filter(function(a) a$col %in% depCols, withDep)[[1]]
  expect_equal(a_dep$val, d * delta * csMat[, mod] * as.numeric(d != 0L))
})

test_that("a dependency's Du matches shifting both statistics directly", {
  skip_on_cran()
  fx <- cs_fixture("ans2", "mydata2", "mymodel2", "mynet2")
  skip_if(is.null(fx), "ans2 fixtures unavailable")

  csMap <- fx$cs$changeStatsMap
  d     <- fx$cs$density
  csMat <- fx$cs$csMat
  thEff <- RSiena:::buildThetaEff(fx$theta, csMap)
  dIdx  <- fx$cs$densityIdx

  b   <- fx$resolve("transTrip")
  dep <- fx$resolve("transRecTrip")
  mod <- fx$resolve("recip")
  skip_if(any(is.na(c(b, dep, mod))), "model lacks the interaction trio")

  delta   <- 1
  shifts  <- setNames(list(delta, delta * csMat[, mod]), c(b, dep))
  A       <- RSiena:::perturbationJacobian(shifts, d, csMap)
  du_A    <- RSiena:::utilityShift(A, fx$theta, length(d))

  ## Independent: move BOTH change statistics and let calculateUtility say
  ## what happened. Shares no code with perturbationJacobian().
  shifted <- csMat
  for (nm in names(shifts)) shifted[, nm] <- shifted[, nm] + shifts[[nm]]
  du_direct <- RSiena:::calculateUtility(shifted, thEff, NULL, dIdx) -
               RSiena:::calculateUtility(csMat,   thEff, NULL, dIdx)

  expect_equal(du_A, du_direct, tolerance = 1e-12)
})

test_that("shifts reaching the same parameter are summed, not duplicated", {
  skip_on_cran()
  fx <- cs_fixture("ans", "mydata", "mymodel", "mynet")
  skip_if(is.null(fx), "base fixtures unavailable")

  csMap <- fx$cs$changeStatsMap
  d     <- fx$cs$density
  eff   <- fx$resolve("transTrip")

  ## An effect that is both perturbed and depended upon would otherwise
  ## appear twice and be counted twice by the consumer.
  once  <- RSiena:::perturbationJacobian(setNames(list(1), eff), d, csMap)
  twice <- RSiena:::perturbationJacobian(setNames(list(2), eff), d, csMap)
  expect_length(once, length(twice))
  expect_equal(twice[[1]]$val, 2 * once[[1]]$val)
})


# ── Cells: the counterfactual states a target compares ───────────────────────
#
# The point of naming the states is that a dependency is evaluated AT a state
# rather than patched onto an intermediate value. The decisive property is the
# cross term: for A ~ B:C with both B and C perturbed, the AB cell must carry
# dB*dC, and it must do so because the product was evaluated at the shifted
# state -- not because a cross-term rule was added.

test_that("a first difference has two cells with opposite weights", {
  cells <- RSiena:::perturbationCells(list(recip = 1))
  expect_length(cells, 2L)
  expect_equal(vapply(cells, `[[`, "", "name"), c("base", "A"))
  expect_equal(vapply(cells, `[[`, 0, "weight"), c(-1, 1))
  ## The base cell perturbs nothing; that is what makes it the base.
  expect_length(cells[[1]]$shifts, 0L)
  expect_equal(cells[[2]]$shifts$recip, 1)
})

test_that("a second difference has four cells in the difference-of-differences pattern", {
  cells <- RSiena:::perturbationCells(list(transTrip = 1, recip = 2))
  expect_equal(vapply(cells, `[[`, "", "name"), c("base", "A", "B", "AB"))
  ## +1 -1 -1 +1 is (p_AB - p_B) - (p_A - p_base).
  expect_equal(vapply(cells, `[[`, 0, "weight"), c(1, -1, -1, 1))
  expect_equal(names(cells[[2]]$shifts), "transTrip")
  expect_equal(names(cells[[3]]$shifts), "recip")
  expect_setequal(names(cells[[4]]$shifts), c("transTrip", "recip"))
})

test_that("a dependency's cross term appears in the AB cell and nowhere else", {
  s_B <- c(1, 2, 3); s_C <- c(10, 20, 30)
  base <- list(B = s_B, C = s_C, A = s_B * s_C)
  dB <- 1; dC <- 2

  cells <- RSiena:::perturbationCells(
      list(B = dB, C = dC), baseStats = base,
      dependencies = list(list(target = "A", terms = c("B", "C"))))
  byName <- setNames(cells, vapply(cells, `[[`, "", "name"))

  ## Perturbing B alone moves A by dB * s_C -- no cross term.
  expect_equal(byName$A$shifts$A, dB * s_C)
  ## Perturbing C alone moves A by dC * s_B.
  expect_equal(byName$B$shifts$A, dC * s_B)
  ## Both: the full product difference, which INCLUDES dB*dC.
  expect_equal(byName$AB$shifts$A, dB * s_C + dC * s_B + dB * dC)
  ## and that is exactly (s_B+dB)(s_C+dC) - s_B*s_C, i.e. the relation
  ## evaluated at the cell's state rather than a rule that was added.
  expect_equal(byName$AB$shifts$A, (s_B + dB) * (s_C + dC) - s_B * s_C)
  ## The base cell moves nothing at all.
  expect_length(byName$base$shifts, 0L)
})

test_that("a dependency is applied as a delta, leaving a mismatched base alone", {
  ## The declaration asserts how A MOVES, not that A equals the product in the
  ## data. Where they disagree, overwriting A with the product would silently
  ## change the base state; a delta must not.
  s_B <- c(1, 2); s_C <- c(5, 7)
  base <- list(B = s_B, C = s_C, A = c(999, 999))   # deliberately not s_B*s_C

  cells <- RSiena:::perturbationCells(
      list(B = 1), baseStats = base,
      dependencies = list(list(target = "A", terms = c("B", "C"))))
  byName <- setNames(cells, vapply(cells, `[[`, "", "name"))
  ## The shift depends only on the relation, not on A's actual base value.
  expect_equal(byName$A$shifts$A, 1 * s_C)
})

test_that("higher-order targets are refused rather than silently truncated", {
  expect_error(RSiena:::perturbationCells(list(a = 1, b = 1, c = 1)),
               "Only first and second differences")
})


# ── Cells reproduce the existing second difference ───────────────────────────
#
# Routing second differences through cells has to be a PURE REFACTOR: the same
# numbers, computed from states instead of by patching increments. This is the
# check that licenses the switch.
#
# It rests on the mlogit update composing additively -- applying shift a then
# shift b equals applying a+b -- so that the current sequential computation
# (shift effect2, then take effect1's difference there) equals a single update
# per cell. That property is verified directly below, for both perturbation
# types, because the ego update renormalises and had no obvious reason to
# compose.

test_that("the mlogit update composes additively, for both perturbation types", {
  set.seed(7)
  n <- 12; g <- rep(1:3, each = 4)
  p <- as.vector(unlist(tapply(runif(n), g, function(z) z / sum(z))))
  a <- rnorm(n, 0, 0.4); b <- rnorm(n, 0, 0.4)
  for (pt in c("alter", "ego")) {
    seq2 <- RSiena:::mlogit_update_r(RSiena:::mlogit_update_r(p, a, g, pt), b, g, pt)
    one  <- RSiena:::mlogit_update_r(p, a + b, g, pt)
    expect_equal(seq2, one, tolerance = 1e-12,
      info = paste("perturbType", pt, "must compose for cells to be exact"))
  }
})

test_that("cells reproduce marginalEffects' second difference exactly", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mydata2 <- load_fixture("mydata2")
  mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mydata2) || is.null(mymodel2),
          "ans2 fixtures unavailable")

  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             ## egoChoice is the row-level aggregation: one
                             ## row per (period, ego, choice). "none" would
                             ## collapse to a single number.
                             type = "changeProb", level = "egoChoice",
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_second_diff(tg,
          list(transTrip = list(diff = 1), recip = list(diff = 1)), name = "sd"))
  ref <- suppressMessages(marginalEffects(ans2, mydata2, targets = tg,
           control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
           control_algo = set_postest_algo_saom(verbose = FALSE)))

  ## Independent cell computation, sharing no code with calculateSecondDiff().
  w  <- RSiena:::getStaticChangeContributions(ans = ans2, data = mydata2,
          effects = mymodel2, depvar = "mynet2", returnWide = TRUE)
  th <- coef(ans2)
  cs <- RSiena:::contribToChangeStats(w$contribMat, w$effectNames, theta = th)
  thEff <- RSiena:::buildThetaEff(th, cs$changeStatsMap)
  u  <- RSiena:::calculateUtility(cs$csMat, thEff, NULL, cs$densityIdx)
  p0 <- RSiena:::calculateChangeProb(u, w$group_id)

  nm <- function(short)
    grep(paste0("(^|_)", short, "$"), colnames(cs$csMat), value = TRUE)[1L]
  cells <- RSiena:::perturbationCells(
             setNames(list(1, 1), c(nm("transTrip"), nm("recip"))))

  sd_cells <- numeric(length(p0))
  for (cell in cells) {
    A  <- RSiena:::perturbationJacobian(cell$shifts, cs$density,
                                        cs$changeStatsMap)
    du <- RSiena:::utilityShift(A, th, length(p0))
    p  <- RSiena:::mlogit_update_r(p0, du, w$group_id, "alter")
    sd_cells <- sd_cells + cell$weight * p
  }

  ## marginalEffects drops rows with density == 0 (no-change alternatives).
  keep <- cs$density != 0L
  expect_equal(nrow(ref), sum(keep))
  expect_equal(ref$secondDiff, sd_cells[keep], tolerance = 1e-12)
})


# ── applyUtilityShift: one formula, checked against the paper ────────────────
#
# The paper's eq. (mlogitupdate), Ben-Akiva & Lerman 1985 p.114, is stated as
# valid for ANY linear-in-parameters preference shift:
#
#   p_ij(x; Df_i) = p_ij(x) exp(Df_ij) / sum_h p_ih(x) exp(Df_ih)
#
# The two code paths in the package are that formula specialised to the two
# shapes the shift takes -- a one-alternative perturbation (dyadic effect, n
# separate counterfactuals) and an ego-wide one (ego-specific effect, one
# counterfactual per ego). These tests check the implementation against the
# equation itself, not against the existing code, and then check that it
# reproduces the existing code where the existing code applies.

test_that("applyUtilityShift matches the paper's equation computed directly", {
  set.seed(21)
  n <- 15; g <- rep(1:3, each = 5)
  p <- as.vector(unlist(tapply(runif(n), g, function(z) z / sum(z))))
  shiftAll   <- rnorm(n, 0, 0.3)   # ego-wide part
  shiftFocal <- rnorm(n, 0, 0.3)   # focal-alternative part

  ## Direct transcription of eq. (mlogitupdate), one counterfactual per row:
  ## row j's counterfactual shifts every alternative by shiftAll, and j alone
  ## additionally by shiftFocal[j].
  direct <- numeric(n)
  for (j in seq_len(n)) {
    inGroup <- which(g == g[j])
    Df <- shiftAll[inGroup]
    Df[inGroup == j] <- Df[inGroup == j] + shiftFocal[j]
    direct[j] <- (p[j] * exp(shiftAll[j] + shiftFocal[j])) /
                 sum(p[inGroup] * exp(Df))
  }
  expect_equal(RSiena:::applyUtilityShift(p, shiftAll, shiftFocal, g),
               direct, tolerance = 1e-14)
})

test_that("applyUtilityShift reduces to both existing paths", {
  set.seed(5)
  n <- 12; g <- rep(1:3, each = 4)
  p <- as.vector(unlist(tapply(runif(n), g, function(z) z / sum(z))))
  d <- rnorm(n, 0, 0.5)

  ## Dyadic: the shift is on the focal alternative only, so the denominator
  ## collapses and no group sum is needed -- the "alter" path.
  expect_equal(RSiena:::applyUtilityShift(p, 0, d, g),
               RSiena:::mlogit_update_r(p, d, g, "alter"), tolerance = 1e-14)

  ## Ego-specific: the shift applies to every alternative of the ego -- the
  ## "ego" path, which needs the full denominator.
  expect_equal(RSiena:::applyUtilityShift(p, d, 0, g),
               RSiena:::mlogit_update_r(p, d, g, "ego"), tolerance = 1e-14)
})

test_that("an ego-wide shift that is uniform leaves probabilities unchanged", {
  ## A property of the multinomial the ego path must have and the collapsed
  ## form must not: adding the same constant to every alternative's utility
  ## cancels in the ratio. Catches a shift routed through the wrong part.
  set.seed(9)
  n <- 12; g <- rep(1:3, each = 4)
  p <- as.vector(unlist(tapply(runif(n), g, function(z) z / sum(z))))
  expect_equal(RSiena:::applyUtilityShift(p, rep(0.7, n), 0, g), p,
               tolerance = 1e-14)
  ## whereas the same constant as a focal shift genuinely moves probabilities
  expect_false(isTRUE(all.equal(
      RSiena:::applyUtilityShift(p, 0, rep(0.7, n), g), p)))
})

test_that("applyUtilityShift keeps each ego's probabilities a distribution", {
  ## Only meaningful for the ego-wide case: there is ONE counterfactual per
  ## ego, so its alternatives must still sum to 1. The dyadic case is n
  ## separate counterfactuals and carries no such constraint.
  set.seed(13)
  n <- 12; g <- rep(1:3, each = 4)
  p <- as.vector(unlist(tapply(runif(n), g, function(z) z / sum(z))))
  out <- RSiena:::applyUtilityShift(p, rnorm(n, 0, 0.4), 0, g)
  ## as.numeric(): tapply returns a 1-d array, which expect_equal would
  ## otherwise report as differing from a plain vector on structure.
  expect_equal(as.numeric(tapply(out, g, sum)), rep(1, 3), tolerance = 1e-14)
})


test_that("NA marks an unperturbed row and must survive the update", {
  ## A contrast leaves rows outside its range at NA, meaning "not perturbed",
  ## and those rows must DROP OUT of the aggregation. Treating NA as a zero
  ## shift returns p unchanged -- a valid-looking number that stays in the
  ## average and pulls it toward zero. The existing kernel propagates NA;
  ## so must this.
  p <- c(0.2, 0.3, 0.5); g <- c(1L, 1L, 1L); d <- c(0.4, NA, 0.1)
  new <- RSiena:::applyUtilityShift(p, 0, d, g)
  old <- RSiena:::mlogit_update_r(p, d, g, "alter")
  expect_true(is.na(new[2]))
  expect_equal(new, old, tolerance = 1e-14)

  ## and from the ego-wide part too
  expect_true(is.na(RSiena:::applyUtilityShift(p, c(0.1, NA, 0.2), 0, g)[2]))
})

test_that("the direction vector supplies the paper's creation/dissolution sign switch", {
  skip_on_cran()
  fx <- cs_fixture("ans", "mydata", "mymodel", "mynet")
  skip_if(is.null(fx), "base fixtures unavailable")

  ## Paper, ego-specific section: an ego-wide shift adds +beta*Ddelta to all
  ## CREATION alternatives, subtracts it from all DISSOLUTION alternatives,
  ## and leaves the no-change option at zero. In this implementation that
  ## sign structure comes from the `direction` multiplier, not from a
  ## separate rule -- so it has to be checked, or a sign error would be
  ## invisible.
  d     <- fx$cs$density
  csMap <- fx$cs$changeStatsMap
  eff   <- fx$resolve("transTrip")
  A <- RSiena:::perturbationJacobian(setNames(list(1), eff), d, csMap)
  expect_length(A, 1L)
  v <- A[[1]]$val

  expect_true(all(v[d ==  1L] > 0), info = "creation rows must shift up")
  expect_true(all(v[d == -1L] < 0), info = "dissolution rows must shift down")
  expect_true(all(v[d ==  0L] == 0), info = "no-change option must not shift")
  ## and symmetric: the same magnitude, opposite sign
  expect_equal(unique(v[d == 1L]), -unique(v[d == -1L]))
})


# ── Defects found by reading the logic, not by running it ────────────────────
#
# These four were found reviewing code that nothing calls yet, so no test could
# have caught them. Each is pinned here so the fix cannot regress.

test_that("an unreachable row stays NA through the shift construction", {
  ## .effectColumns used to zero NA shifts. NA means the perturbation does not
  ## reach that row (a contrast row outside its range); as a zero shift the row
  ## comes back unchanged and stays in the average, dragging it toward zero.
  csMap <- list(bases = c("x", "x"), types = c("eval", "creation"))
  A <- RSiena:::.effectColumns("x", c(1, NA, 1), c(1L, 1L, -1L), csMap, TRUE)
  expect_true(any(vapply(A, function(a) any(is.na(a$val)), logical(1L))))
  ## and it must survive all the way to the utility shift
  du <- RSiena:::utilityShift(A, c(0.5, 0.5), 3L)
  expect_true(is.na(du[2]))
})

test_that("a dependency shift that is entirely zero-or-NA does not crash", {
  ## `all(d == 0)` is NA as soon as one element is NA, and `if (NA)` errors
  ## rather than deciding. Reachable whenever a contrast leaves rows
  ## unperturbed and a dependency is declared.
  expect_no_error(
    RSiena:::perturbationCells(
      list(B = c(0, NA)),
      baseStats = list(B = c(1, 2), C = c(0, 4), A = c(0, 8)),
      dependencies = list(list(target = "A", terms = c("B", "C")))))
})

test_that("the unified update is the two existing kernels composed, in order", {
  ## No new C++ kernel: the closed form IS mlogit_update(ego) followed by
  ## mlogit_update(alter). Reusing them keeps both compiled paths and their
  ## exactness, and removes the hand-rolled group sum entirely (which had a
  ## silent failure mode on non-contiguous groups).
  set.seed(17); n <- 3000; g <- rep(seq_len(n / 30), each = 30)
  p <- as.vector(unlist(tapply(runif(n), g, function(z) z / sum(z))))
  sa <- rnorm(n, 0, 0.3); sf <- rnorm(n, 0, 0.3)

  composed <- RSiena:::mlogit_update_r(
                RSiena:::mlogit_update_r(p, sa, g, "ego"), sf, g, "alter")
  expect_equal(RSiena:::applyUtilityShift(p, sa, sf, g), composed,
               tolerance = 1e-14)

  ## The ORDER is forced by the algebra, not chosen: the ego-wide part is
  ## common to every alternative and so must sit inside every denominator.
  ## Reversed it is wrong by ~2e-2, i.e. not a rounding difference.
  reversed <- RSiena:::mlogit_update_r(
                RSiena:::mlogit_update_r(p, sf, g, "alter"), sa, g, "ego")
  expect_gt(max(abs(composed - reversed)), 1e-3)

  ## Single-shape shifts stay byte-exact against the kernel they reuse.
  expect_identical(RSiena:::applyUtilityShift(p, 0, sf, g),
                   RSiena:::mlogit_update_r(p, sf, g, "alter"))
  expect_identical(RSiena:::applyUtilityShift(p, sa, 0, g),
                   RSiena:::mlogit_update_r(p, sa, g, "ego"))
})

test_that("an all-NA direction does not turn a skip into an error", {
  ## any(mask) is NA when the direction is unknown, and `if (!NA)` errors.
  csMap <- list(bases = "x", types = "creation")
  expect_no_error(
    RSiena:::.effectColumns("x", 1, rep(NA_integer_, 3), csMap, TRUE))
})


# ── The analytic Jacobian must decline what it cannot do ─────────────────────
#
# firstDiffJacobian() is the derivative of the ONE-ALTERNATIVE update. An
# ego-wide perturbation renormalises over the ego's whole choice set and so has
# a different derivative, which is not implemented. Using the one-alternative
# form there gave an SE about twice the truth -- silently, because the point
# estimate DID use the ego update while the Jacobian did not.
#
# Pre-existing (HEAD's predictFirstDiffJac has no perturbType handling either).
# Finite differences are correct for this case, so the analytic path declines.

test_that("delta SE agrees with bootstrap for BOTH perturbation types", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("RSENA_FULL_TESTS"), "1"),
              "bootstrap comparison; RSENA_FULL_TESTS not set")
  ans <- load_fixture("ans"); mydata <- load_fixture("mydata")
  mymodel <- load_fixture("mymodel")
  skip_if(is.null(ans) || is.null(mydata) || is.null(mymodel),
          "base fixtures unavailable")

  for (pt in c("alter", "ego")) {
    mk <- function() {
      tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                                 type = "tieProb", level = "period",
                                 includeDefaults = FALSE)
      suppressMessages(set_target(tg, transTrip, diff = 1, perturbType = pt))
    }
    set.seed(101)
    bo <- suppressMessages(marginalEffects(ans, mydata, targets = mk(),
            control_uncertainty = set_postest_uncertainty_saom(
                mode = "bootstrap", nsim = 400),
            control_algo = set_postest_algo_saom(verbose = FALSE)))
    de <- suppressMessages(marginalEffects(ans, mydata, targets = mk(),
            control_uncertainty = set_postest_uncertainty_saom(
                mode = "delta", nsim = 1L),
            control_algo = set_postest_algo_saom(verbose = FALSE)))
    ratio <- bo$SE / de$SE
    ## Monte Carlo noise at nsim = 400 is a few percent; a wrong Jacobian form
    ## was off by a factor of two, so this band separates them cleanly.
    expect_true(all(ratio > 0.85 & ratio < 1.15),
      info = paste0("perturbType = ", pt, ": bootstrap/delta ratio ",
                    paste(round(ratio, 3), collapse = ", ")))
  }

  ## The same check for a second difference whose two sides are of DIFFERENT
  ## perturbation types.  A second difference composes two updates, and while
  ## the Jacobian applied one harness to every cell -- the ego one, whenever
  ## either side was ego -- the SE described a different counterfactual than
  ## the estimate.  On the Glasgow egoX x altX difference that landed at
  ## 0.31-0.70 of the bootstrap SE while finite differences sat at 0.78-1.03,
  ## and nothing in the suite noticed, because every second difference tested
  ## here had matching types on both sides.
  mk2 <- function() {
    tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                               type = "tieProb", level = "period",
                               includeDefaults = FALSE)
    suppressMessages(set_second_diff(tg, list(
        transTrip = list(diff = 1, perturbType = "ego"),
        recip     = list(contrast = c(0, 1), perturbType = "alter")),
        name = "mixed_sd"))
  }
  set.seed(202)
  bo2 <- suppressMessages(marginalEffects(ans, mydata, targets = mk2(),
          control_uncertainty = set_postest_uncertainty_saom(
              mode = "bootstrap", nsim = 400),
          control_algo = set_postest_algo_saom(verbose = FALSE)))
  de2 <- suppressMessages(marginalEffects(ans, mydata, targets = mk2(),
          control_uncertainty = set_postest_uncertainty_saom(
              mode = "delta", nsim = 1L),
          control_algo = set_postest_algo_saom(verbose = FALSE)))
  ratio2 <- bo2$SE / de2$SE
  expect_true(all(ratio2 > 0.85 & ratio2 < 1.15),
    info = paste0("mixed-type second difference: bootstrap/delta ratio ",
                  paste(round(ratio2, 3), collapse = ", ")))
})

test_that("an ego perturbation uses the analytic path, with the ego harness", {
  skip_on_cran()
  ans <- load_fixture("ans"); mydata <- load_fixture("mydata")
  mymodel <- load_fixture("mymodel")
  skip_if(is.null(ans) || is.null(mydata) || is.null(mymodel),
          "base fixtures unavailable")

  ## The split that makes this work: the CORE (A = dDu/dtheta -- which
  ## parameters a perturbation reaches, on which rows) is effect-dependent and
  ## identical for both perturbation types. Only the HARNESS differs, chaining
  ## A through the update's derivative; the ego-wide form carries a group sum
  ## because every alternative of the ego sits in the denominator.
  ##
  ## Before the ego harness existed, the one-alternative derivative was used
  ## for both and the ego SE came out about twice the truth.
  seen <- new.env(parent = emptyenv()); seen$ego <- 0L
  orig <- RSiena:::firstDiffJacobianEgo
  on.exit(assignInNamespace("firstDiffJacobianEgo", orig, ns = "RSiena"))
  assignInNamespace("firstDiffJacobianEgo",
    function(...) { seen$ego <- seen$ego + 1L; orig(...) }, ns = "RSiena")

  mk <- function(pt) {
    tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                               type = "tieProb", level = "period",
                               includeDefaults = FALSE)
    suppressMessages(set_target(tg, transTrip, diff = 1, perturbType = pt))
  }
  run <- function(pt) suppressMessages(marginalEffects(ans, mydata,
      targets = mk(pt),
      control_uncertainty = set_postest_uncertainty_saom(mode = "delta",
                                                         nsim = 1L),
      control_algo = set_postest_algo_saom(verbose = FALSE)))

  seen$ego <- 0L; a <- run("alter")
  expect_equal(seen$ego, 0L)      # dyadic: one-alternative harness
  seen$ego <- 0L; e <- run("ego")
  expect_gt(seen$ego, 0L)         # ego-wide: group-sum harness

  ## and the two must give different SEs -- if they agreed, the harness would
  ## not be doing anything, which is exactly the bug this replaces.
  expect_false(isTRUE(all.equal(a$SE, e$SE)))
})


# ── The invariant: point estimate and Jacobian share one Ddelta ──────────────
#
# The counterfactual is a modelling CHOICE (perturb in delta-space; propagate
# direct functional dependence; do not propagate indirect co-variation). So
# "correct" means CONSISTENT: the point estimate, the Jacobian and the
# delta-method SE must all describe the same counterfactual.
#
# If the Jacobian propagates a dependency the point estimate does not, the SE
# is the standard error of a DIFFERENT quantity -- and it looks entirely
# reasonable. That is precisely how the ego-perturbation SE came to be twice
# its true value.
#
# This checks the two derivations agree on Du itself, not merely that both
# produce numbers.

test_that("point estimate and Jacobian imply the SAME utility shift", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mydata2 <- load_fixture("mydata2")
  mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mydata2) || is.null(mymodel2),
          "ans2 fixtures unavailable")

  w  <- RSiena:::getStaticChangeContributions(ans = ans2, data = mydata2,
          effects = mymodel2, depvar = "mynet2", returnWide = TRUE)
  th <- coef(ans2)
  cs <- RSiena:::contribToChangeStats(w$contribMat, w$effectNames, theta = th)
  thEff <- RSiena:::buildThetaEff(th, cs$changeStatsMap)
  nm <- function(x)
    grep(paste0("(^|_)", x, "$"), colnames(cs$csMat), value = TRUE)[1L]
  B <- nm("transTrip"); C <- nm("recip"); A <- nm("transRecTrip")
  skip_if(any(is.na(c(A, B, C))), "model lacks the interaction trio")

  delta <- 1
  dep   <- list(list(target = A, terms = c(B, C)))

  ## Route 1 -- the shift set, then Du = Ddelta . theta.
  shifts <- RSiena:::.shiftSetFor(setNames(list(delta), B), dep, cs$csMat)
  Amat   <- RSiena:::perturbationJacobian(shifts, cs$density,
                                          cs$changeStatsMap)
  du_new <- RSiena:::utilityShift(Amat, th, nrow(cs$csMat))

  ## Route 2 -- shift the change statistics by hand and let the utility
  ## function say what happened.  Shares no code with route 1, which is what
  ## makes agreement evidence rather than tautology.
  ##
  ## A relation states how its target MOVES, not what its level is: A here is
  ## a real effect with its own change statistic, and the declaration says
  ## that statistic responds as a product would.  So A shifts by
  ## (B+delta)*C - B*C = delta*C, on top of whatever A already was.
  shifted <- cs$csMat
  shifted[, B] <- cs$csMat[, B] + delta
  shifted[, A] <- cs$csMat[, A] + delta * cs$csMat[, C]
  du_direct <-
    RSiena:::calculateUtility(shifted,   thEff, NULL, cs$densityIdx) -
    RSiena:::calculateUtility(cs$csMat, thEff, NULL, cs$densityIdx)

  expect_equal(du_new, du_direct, tolerance = 1e-12,
    info = paste("the Jacobian and the point estimate must describe the same",
                 "counterfactual; if they diverge the SE belongs to a",
                 "different quantity than the estimate"))
})


test_that("accumulated aggregation honours `condition` instead of dropping it", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("RSENA_FULL_TESTS"), "1"),
              "dynamic simulation; RSENA_FULL_TESTS not set")
  ans <- load_fixture("ans"); mydata <- load_fixture("mydata")
  mymodel <- load_fixture("mymodel"); mycontrols <- load_fixture("mycontrols")
  skip_if(is.null(ans) || is.null(mydata) || is.null(mymodel) ||
          is.null(mycontrols), "dynamic fixtures unavailable")

  ## Conditioning used to be discarded here with a warning, on the reasoning
  ## that accumulation sums over ministeps while a condition varies between
  ## them. But that IS what stratifying means: accumulate within each stratum.
  ## The columns only had to survive the ego normalisation and then join the
  ## grouping -- the same thing getGroupVars() does on the non-accumulated
  ## path.
  run <- function(acc) {
    tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
            type = "tieProb", level = "period", dynamic = TRUE,
            accumulated = acc, condition = "density",
            includeDefaults = FALSE)
    tg <- suppressMessages(set_target(tg, transTrip, diff = 1))
    set.seed(5)
    suppressMessages(marginalEffects(ans, mydata, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
      control_algo = set_postest_algo_saom(algorithm = mycontrols, n3 = 10,
                                           verbose = FALSE)))
  }
  plain <- expect_no_warning(run(FALSE))
  acc   <- expect_no_warning(run(TRUE))

  ## Same conditioning structure as the non-accumulated path.
  expect_true("mynet_density" %in% names(acc))
  expect_equal(nrow(acc), nrow(plain))
  ## and genuinely stratified, not one row repeated
  expect_gt(length(unique(acc$mynet_density)), 1L)
})
