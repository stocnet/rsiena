# test-equivalence-harness.R
#
# Self-tests for the API-equivalence harness (helper-equivalence.R).
#
# Two jobs, both of which must hold BEFORE the harness is trusted to verify the
# config-object refactor:
#
#   (1) compare_me_output() must actually detect differences.  A comparison
#       function that returns TRUE for everything would report the refactor as
#       verified while checking nothing — the single worst outcome available
#       here.  Part A feeds it deliberately corrupted input and requires it to
#       fail in each case.
#
#   (2) Every corpus entry must be a valid call that produces output today.
#       Otherwise a broken entry is indistinguishable from a refactor
#       regression once the object form is added.  Part B runs them all.
#
# When the set_postest_*_saom() constructors land, Part C is added: run each
# corpus entry in both forms and require compare_me_output() to return TRUE.

# ── Part A: the comparator is not vacuous ────────────────────────────────────

.mk_df <- function() {
  data.frame(
    period    = c(1L, 1L, 2L, 2L),
    grp       = c("a", "b", "a", "b"),
    firstDiff = c(0.10, -0.20, 0.30, 0.40),
    SE        = c(0.01, 0.02, 0.03, 0.04),
    stringsAsFactors = FALSE
  )
}

test_that("compare_me_output accepts genuinely identical output", {
  a <- .mk_df()
  expect_true(isTRUE(compare_me_output(a, a)))
  expect_true(isTRUE(compare_me_output(a, .mk_df())))
})

test_that("compare_me_output detects a changed numeric value", {
  a <- .mk_df(); b <- .mk_df()
  b$firstDiff[3] <- b$firstDiff[3] + 1e-6      # small, but far above tolerance
  res <- compare_me_output(a, b)
  expect_false(isTRUE(res))
  expect_match(paste(res, collapse = " | "), "firstDiff")
})

test_that("compare_me_output detects a difference just above tolerance", {
  a <- .mk_df(); b <- .mk_df()
  b$SE[1] <- b$SE[1] + 1e-9
  expect_false(isTRUE(compare_me_output(a, b, tolerance = 1e-10)))
  ## ...and accepts it when the tolerance genuinely allows it
  expect_true(isTRUE(compare_me_output(a, b, tolerance = 1e-6)))
})

test_that("compare_me_output detects structural differences", {
  a <- .mk_df()

  b <- .mk_df(); b$SE <- NULL                       # dropped column
  expect_match(paste(compare_me_output(a, b), collapse = " | "), "only in flat")

  b <- .mk_df(); b$extra <- 1                       # added column
  expect_match(paste(compare_me_output(a, b), collapse = " | "), "only in object")

  b <- .mk_df()[, c("period", "grp", "SE", "firstDiff")]   # reordered columns
  expect_match(paste(compare_me_output(a, b), collapse = " | "), "different order")

  b <- .mk_df()[1:3, ]                              # dropped row
  expect_match(paste(compare_me_output(a, b), collapse = " | "), "row count")

  b <- .mk_df()[c(2, 1, 3, 4), ]                    # permuted rows
  expect_false(isTRUE(compare_me_output(a, b)))
})

test_that("compare_me_output detects character and NA differences", {
  a <- .mk_df()

  b <- .mk_df(); b$grp[2] <- "z"
  expect_match(paste(compare_me_output(a, b), collapse = " | "), "non-numeric")

  b <- .mk_df(); b$firstDiff[2] <- NA_real_
  expect_match(paste(compare_me_output(a, b), collapse = " | "), "NA pattern")

  ## NAs in the same places on both sides must NOT be reported as a difference
  a2 <- .mk_df(); a2$firstDiff[2] <- NA_real_
  b2 <- .mk_df(); b2$firstDiff[2] <- NA_real_
  expect_true(isTRUE(compare_me_output(a2, b2)))
})

test_that("compare_me_output rejects non-data.frame input", {
  expect_false(isTRUE(compare_me_output(.mk_df(), list(a = 1))))
  expect_false(isTRUE(compare_me_output(NULL, .mk_df())))
})

test_that("compare_me_output handles multi-effect list results", {
  a <- list(tt = .mk_df(), rp = .mk_df())

  expect_true(isTRUE(compare_me_output(a, a)))

  ## a difference buried inside ONE element must still be reported, named
  b <- list(tt = .mk_df(), rp = .mk_df())
  b$rp$firstDiff[1] <- b$rp$firstDiff[1] + 1e-6
  res <- compare_me_output(a, b)
  expect_false(isTRUE(res))
  expect_match(paste(res, collapse = " | "), "\\[rp\\]")
  expect_match(paste(res, collapse = " | "), "firstDiff")

  ## differing effect sets must be reported
  expect_false(isTRUE(compare_me_output(a, list(tt = .mk_df()))))
})

# ── Part B: every corpus entry is a valid call today ─────────────────────────
#
# Establishes the baseline half of the equivalence pair.  If an entry cannot
# run now, a failure after the refactor would be ambiguous between "the
# refactor broke it" and "the corpus entry was always wrong".

test_that("every equivalence-corpus entry runs and returns a data.frame", {
  fixtures <- me_corpus_fixtures()
  skip_if(is.null(fixtures$ans) || is.null(fixtures$mydata),
          "base fixtures unavailable")

  corpus <- me_corpus(fixtures)
  expect_gt(length(corpus), 0L)

  for (entry in corpus) {
    if (isTRUE(entry$slow) && !identical(Sys.getenv("RSENA_FULL_TESTS"), "1"))
      next
    if (any(vapply(fixtures[entry$needs], is.null, logical(1L))))
      next

    out <- tryCatch(do.call(marginalEffects, entry$args),
                    error = function(e) e)
    expect_false(inherits(out, "error"),
      info = paste0("corpus entry '", entry$name, "' failed: ",
                    if (inherits(out, "error")) conditionMessage(out) else ""))
    if (inherits(out, "error")) next

    ## Multi-effect calls return a named list of data.frames; single-effect
    ## calls return one frame.  Both are valid corpus results.
    frames <- if (is.data.frame(out)) list(out) else out
    expect_true(is.list(frames) && length(frames) > 0L,
      info = paste0("corpus entry '", entry$name,
                    "' must return a data.frame or a list of them"))
    expect_true(all(vapply(frames, is.data.frame, logical(1L))),
      info = paste0("corpus entry '", entry$name,
                    "' returned a non-data.frame element"))
    expect_true(all(vapply(frames, nrow, integer(1L)) > 0L),
      info = paste0("corpus entry '", entry$name, "' returned empty output"))

    ## Self-comparison must hold: a result is always equivalent to itself.
    ## Cheap, but it exercises the comparator against every real output shape
    ## the corpus produces — not just the synthetic frames in Part A.
    expect_true(isTRUE(compare_me_output(out, out)),
      info = paste0("self-comparison failed for '", entry$name, "'"))
  }
})

# ── Part C: flat form and object form must agree ─────────────────────────────
#
# The point of the whole harness.  For every corpus entry, run the call in its
# original flat form and in the config-object form derived mechanically by
# as_object_args(), and require the outputs to be identical.
#
# Scope at step 3: only the uncertainty and control domains move into objects;
# the effect specification is still passed flat.  Step 4 extends this to the
# model object.
#
# TRANSITIONAL.  This test compares old against new and therefore only works
# while both exist.  Step 5 removes the flat form, at which point Part C is
# deleted — the migrated test suite covers the object form and the step-1
# snapshots pin the numbers.  compare_me_output() is kept beyond that, since
# steps 6 and 6b need before/after comparison of internals.
#
# Skips automatically until the constructors exist.

test_that("flat and config-object forms produce identical output", {
  skip_if_not(exists("set_postest_uncertainty_saom", mode = "function") &&
              exists("set_postest_control_saom", mode = "function"),
              "config-object constructors not implemented yet (step 3)")

  fixtures <- me_corpus_fixtures()
  skip_if(is.null(fixtures$ans) || is.null(fixtures$mydata),
          "base fixtures unavailable")

  corpus  <- me_corpus(fixtures)
  checked <- 0L

  for (entry in corpus) {
    if (isTRUE(entry$slow) && !identical(Sys.getenv("RSENA_FULL_TESTS"), "1"))
      next
    if (any(vapply(fixtures[entry$needs], is.null, logical(1L))))
      next

    obj_args <- as_object_args(entry$args)
    skip_if(is.null(obj_args), "as_object_args() unavailable")

    ## Bootstrap draws are stochastic: seed both runs identically so any
    ## difference is attributable to the refactor rather than to RNG.
    set.seed(4242L); flat <- do.call(marginalEffects, entry$args)
    set.seed(4242L); objf <- do.call(marginalEffects, obj_args)

    res <- compare_me_output(flat, objf)
    expect_true(isTRUE(res),
      info = paste0("corpus entry '", entry$name, "': ",
                    paste(if (isTRUE(res)) "" else res, collapse = " | ")))
    checked <- checked + 1L
  }

  ## Guard against the harness silently checking nothing — e.g. if every entry
  ## skipped, or as_object_args() quietly returned inputs unchanged.
  expect_gt(checked, 0L)
})

test_that("as_object_args actually routes arguments into the objects", {
  skip_if_not(exists("set_postest_uncertainty_saom", mode = "function") &&
              exists("set_postest_control_saom", mode = "function"),
              "config-object constructors not implemented yet (step 3)")

  flat <- list(object = NULL, data = NULL, effectName1 = "transTrip",
               diff1 = 1, level = "period", type = "tieProb",
               uncertainty = TRUE, uncertaintyMode = "delta", nsim = 7L,
               dynamic = FALSE, n3 = 33, verbose = FALSE)
  obj <- as_object_args(flat)

  ## the moved arguments must NOT survive as flat arguments...
  for (nm in c("uncertainty", "uncertaintyMode", "nsim", "dynamic", "n3",
               "verbose"))
    if (nm %in% c("uncertaintyMode", "nsim", "dynamic", "n3", "verbose"))
      expect_false(nm %in% names(obj),
        info = paste0("'", nm, "' must be folded into a config object"))

  ## ...they must be inside the objects, with their values preserved
  expect_s3_class(obj$uncertainty, "sienaPostestUncertainty")
  expect_s3_class(obj$control, "sienaPostestControl")
  expect_equal(obj$uncertainty$mode, "delta")
  expect_equal(obj$uncertainty$nsim, 7L)
  expect_equal(obj$control$n3, 33L)
  expect_false(obj$control$dynamic)

  ## ...and the untouched ones must pass straight through
  expect_equal(obj$effectName1, "transTrip")
  expect_equal(obj$level, "period")
})
