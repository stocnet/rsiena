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
# Part C compared the flat and object forms; it was retired with the flat
# form in step 5c -- see the note where it stood.
# (historical) When the set_postest_*_saom() constructors landed, Part C ran each
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

# ── Part B: every corpus entry is a valid call ───────────────────────────────
#
# Originally the baseline half of the equivalence pair.  With Part C gone this
# stands on its own: a broad smoke check that every configuration the corpus
# describes still runs and returns something of the right shape.  Cheaper than
# the migrated suite's assertions and much wider in the combinations it
# covers, which is what makes it worth keeping.

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

    out <- tryCatch(do.call(marginalEffects, me_call(entry, fixtures)),
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

# ── Part C: REMOVED in step 5c ───────────────────────────────────────────────
#
# Part C ran every corpus entry through both the flat and the config-object
# form and required identical output.  Step 5c removed the flat form, so there
# is no longer a second implementation to compare against and the test cannot
# be written.
#
# What replaced it, and why the corpus is still here:
#
#   * tests/testthat/test-snapshots.R pins the NUMBERS -- migrated to the
#     object interface before the removal and verified to reproduce the cached
#     golden values, so the numeric reference survived the transition.
#   * tests/testthat/test-postest-behaviour-baseline.R pins the BEHAVIOURS
#     step 5d intends to change, so each change stays attributable.
#   * compare_me_output() is kept: steps 6 and 6b need before/after comparison
#     of internals when the analytic Jacobian lands.
#
# The corpus itself (me_corpus) is kept for Parts A and B, which do not need a
# second implementation.

# The test that stood here checked as_object_args() routed each argument into
# the right config object. It went with the translator in step 5c: there is no
# routing left to test, because there is only one way to pass these settings.
