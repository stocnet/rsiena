# test-postest-dependencies.R
#
# The effect-dependency layer: declaring how effects relate, so that
# perturbing one propagates to the change statistics that move with it.
#
# Iteration 1 supports multiplicative relationships only.  `A ~ B:C` means
# A's change statistic is the product of B's and C's, so perturbing B by delta
# also moves A by delta * s_C.  That is exactly what the existing compound
# effect path computes, so the declaration is resolved INTO
# intEffectNames/modEffectNames rather than given a separate numeric path.
#
# The decisive test is therefore an equivalence: a declared dependency must
# reproduce the hand-written interaction arguments bit for bit.  Anything else
# would mean the layer is computing something different from the machinery it
# is meant to be a front end for.

test_that("a dependency actually changes the result", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mydata2 <- load_fixture("mydata2")
  mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mydata2) || is.null(mymodel2),
          "ans2 fixtures unavailable")

  mk <- function(tg) marginalEffects(ans2, mydata2, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
      control_algo = set_postest_algo_saom(verbose = FALSE))

  base <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                               type = "tieProb", level = "period",
                               includeDefaults = FALSE)
  base <- suppressMessages(set_target(base, transTrip, diff = 1))
  withd <- suppressMessages(
    set_dependency(base, transRecTrip ~ transTrip:recip))

  ## Guards against a no-op layer: if declaring a dependency left the numbers
  ## untouched, the equivalence test above would still pass for the wrong
  ## reason (both paths doing nothing).
  ## A single target now comes back as the data.frame itself, not a
  ## one-element list, so no indexing.
  a <- mk(base); b <- mk(withd)
  expect_true(is.data.frame(a)); expect_true(is.data.frame(b))
  expect_false(isTRUE(all.equal(a$firstDiff, b$firstDiff)),
    info = "declaring a dependency must change the first differences")
})

test_that("a dependency reaches BOTH steps of a second difference", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mydata2 <- load_fixture("mydata2")
  mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mydata2) || is.null(mymodel2),
          "ans2 fixtures unavailable")

  ## A second difference COMPOSES two steps: effect2 moves from the base, then
  ## effect1 moves from there.  The relation has to be read on both, and the
  ## test for that needs no reference implementation -- P(AB) - P(A) - P(B) +
  ## P(base) is symmetric in the two effects, so swapping which one is
  ## "second" must not change the answer.  It did: while only one step read
  ## the relation, whichever effect sat in position 2 failed to move the
  ## dependent effect, and the two orders disagreed by the cross term.
  sd <- function(first, second) {
    tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                               type = "tieProb", level = "period",
                               includeDefaults = FALSE)
    tg <- suppressMessages(set_second_diff(tg,
            setNames(list(list(diff = 1), list(diff = 1)), c(first, second)),
            name = "sd"))
    tg <- suppressMessages(set_dependency(tg, transRecTrip ~ transTrip:recip))
    marginalEffects(ans2, mydata2, targets = tg,
        control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
        control_algo = set_postest_algo_saom(verbose = FALSE))$secondDiff
  }

  expect_equal(sd("transTrip", "recip"), sd("recip", "transTrip"),
               tolerance = 1e-10)
})

test_that("dependencies covering no target leave the result untouched", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mydata2 <- load_fixture("mydata2")
  mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mydata2) || is.null(mymodel2),
          "ans2 fixtures unavailable")

  mk <- function(tg) marginalEffects(ans2, mydata2, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
      control_algo = set_postest_algo_saom(verbose = FALSE))

  base <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                               type = "tieProb", level = "period",
                               includeDefaults = FALSE)
  base <- suppressMessages(set_target(base, density, contrast = c(-1, 1)))
  ## density appears nowhere in this relation, so it must be inert here.
  other <- suppressMessages(
    set_dependency(base, transRecTrip ~ transTrip:recip))

  expect_true(isTRUE(compare_me_output(mk(base), mk(other))))
})

test_that("unsupported dependency forms are rejected, not approximated", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mymodel2), "ans2 fixtures unavailable")

  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, transTrip, diff = 1))

  bad <- function(f) {
    t2 <- suppressMessages(set_dependency(tg, f))
    tryCatch({RSiena:::.targetsToEffectList(t2); "NO ERROR"},
             error = function(e) conditionMessage(e))
  }
  ## Only `a:b` is supported so far; everything else must error rather than
  ## be silently treated as a product.
  expect_match(bad(transRecTrip ~ transTrip == recip), "[Oo]nly pure-interaction")
  expect_match(bad(transRecTrip ~ gw(transTrip)),      "[Oo]nly pure-interaction")
  expect_match(bad(transRecTrip ~ transTrip + recip),  "[Oo]nly pure-interaction")
  expect_match(bad(transRecTrip ~ transTrip:transTrip), "both sides")

  ## A formula without a left-hand side names no dependent effect.
  expect_error(suppressMessages(set_dependency(tg, ~ transTrip * recip)),
               "left-hand side")
})

test_that("a derived effect is excluded by default but only warned about later", {
  skip_on_cran()
  ans3 <- load_fixture("ans3"); mymodel3 <- load_fixture("mymodel3")
  skip_if(is.null(ans3) || is.null(mymodel3), "ans3 fixtures unavailable")

  ## Declared at construction: the derived effect is not offered as a default
  ## target, because perturbing it holds its components fixed -- which the
  ## declaration says cannot happen.
  tg1 <- make_postest_targets(ans3, effects = mymodel3, depvar = "mynet3",
                              dependencies = list(unspInt ~ recip:inPop))
  expect_false("unspInt" %in% tg1$effectName1[tg1$include])
  expect_true(all(c("density", "recip", "inPop") %in%
                    tg1$effectName1[tg1$include]))

  ## Declared afterwards: warn, do NOT deselect.  Once the object exists a
  ## selection may be deliberate, and silently dropping it would override an
  ## explicit choice.
  tg2 <- make_postest_targets(ans3, effects = mymodel3, depvar = "mynet3")
  expect_true("unspInt" %in% tg2$effectName1[tg2$include])
  expect_warning(tg3 <- suppressMessages(
      set_dependency(tg2, unspInt ~ recip:inPop)), "selected target")
  expect_true("unspInt" %in% tg3$effectName1[tg3$include],
    info = "declaring a dependency must not silently drop a chosen target")

  ## Already deselected: nothing to warn about.
  tg4 <- suppressMessages(set_target(tg2, unspInt, include = FALSE))
  expect_no_warning(suppressMessages(
      set_dependency(tg4, unspInt ~ recip:inPop)))
})
