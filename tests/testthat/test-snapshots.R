# test-snapshots.R
#
# Numerical regression tests for the postestimation pipeline.
# These tests pin down exact output values so that any refactor that silently
# breaks the math will fail here, even if structural tests still pass.
#
# Snapshot fixtures are built by helper-models.R under RSENA_REBUILD_MODELS=1
# and cached in tests/testthat/cache/.
#
# All deterministic snapshots use tolerance = 1e-10.
# Bootstrap sanity tests are stochastic; they use range / ordering checks only.

# ── Load base fixtures (needed for bootstrap sanity run) ─────────────────────
ans        <- load_fixture("ans")
mydata     <- load_fixture("mydata")
# ans2/mydata2 carry transRecTrip, needed by the interaction snapshot (F).
ans2       <- load_fixture("ans2")
mydata2    <- load_fixture("mydata2")
# The effects objects the targets are enumerated from.  Needed since these
# tests moved to the object interface: make_marginal_targets() takes the
# effects object rather than deriving it from the fit.
mymodel    <- load_fixture("mymodel")
mymodel2   <- load_fixture("mymodel2")

# One target, stated the same way in every snapshot.  These tests pin NUMBERS,
# so the call shape should be uniform and out of the way; the settings that
# vary (condition, and which effect is perturbed how) are the arguments.
#
# Every snapshot uses type = "tieProb" and level = "period", which were
# spelled out on each flat call; they are the defaults here for the same
# reason -- a snapshot that silently differed in either would be pinning a
# different quantity than its name says.
# The SE column of a snapshot, whatever vintage the cache is.  Fixtures are
# gitignored and only rebuilt on demand, so one built before the step-5d
# rename still carries "delta_se"; a fresh one carries "SE".  Same number.
snap_se <- function(d) if (!is.null(d[["SE"]])) d[["SE"]] else d[["delta_se"]]

snap_tg <- function(fit, eff, depvar, condition = NULL,
                    dependency = NULL, ...) {
  tg <- make_marginal_targets(fit, effects = eff, depvar = depvar,
                             type = "tieProb", level = "period",
                             condition = condition, includeDefaults = FALSE)
  perturb <- list(...)
  for (nm in names(perturb))
    tg <- suppressMessages(do.call(set_target,
            c(list(x = tg, shortNames = nm), perturb[[nm]])))
  if (!is.null(dependency))
    tg <- suppressMessages(set_dependency(tg, dependency))
  tg
}

# ── Load snapshot fixtures ───────────────────────────────────────────────────
snap_me_static_transTrip  <- load_fixture("snap_me_static_transTrip")
snap_me_static_secondDiff <- load_fixture("snap_me_static_secondDiff")
snap_me_delta_transTrip   <- load_fixture("snap_me_delta_transTrip")
snap_me_delta_interaction <- load_fixture("snap_me_delta_interaction")
snap_me_delta_density     <- load_fixture("snap_me_delta_density")
snap_predict_changeProb   <- load_fixture("snap_predict_changeProb")
snap_predict_tieProb      <- load_fixture("snap_predict_tieProb")

# ── Guard ─────────────────────────────────────────────────────────────────────
.all_snaps_present <- !any(
  vapply(list(snap_me_static_transTrip, snap_me_static_secondDiff,
              snap_me_delta_transTrip,
              snap_me_delta_interaction, snap_me_delta_density,
              snap_predict_changeProb, snap_predict_tieProb),
         is.null, logical(1L))
)

# ── (A) firstDiff static — golden values ─────────────────────────────────────

test_that("snapshot: marginalEffects static firstDiff golden values", {
  skip_if(!.all_snaps_present,
          "Snapshot cache missing — run with RSENA_REBUILD_MODELS=1")
  snap <- snap_me_static_transTrip
  expect_true(is.data.frame(snap))
  expect_true("firstDiff" %in% names(snap))
  # Re-run and compare exactly (same static contributions = deterministic)
  out <- marginalEffects(ans, mydata,
    targets = snap_tg(ans, mymodel, "mynet", condition = "density",
                      transTrip = list(diff = 1)),
    control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
    control_algo = set_postest_algo_saom(verbose = FALSE)
  )
  expect_equal(out$firstDiff, snap$firstDiff, tolerance = 1e-10,
               info = "firstDiff must match cached golden values exactly")
})

# ── (B) secondDiff static — golden values ────────────────────────────────────

test_that("snapshot: marginalEffects static secondDiff golden values", {
  skip_if(!.all_snaps_present,
          "Snapshot cache missing — run with RSENA_REBUILD_MODELS=1")
  snap <- snap_me_static_secondDiff
  tg <- snap_tg(ans, mymodel, "mynet", condition = "density")
  tg <- suppressMessages(set_second_diff(tg,
          list(transTrip = list(diff = 1),
               recip     = list(contrast = c(0, 1)))))
  out <- marginalEffects(ans, mydata, targets = tg,
    control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
    control_algo = set_postest_algo_saom(verbose = FALSE)
  )
  expect_equal(out$secondDiff, snap$secondDiff, tolerance = 1e-10,
               info = "secondDiff must match cached golden values exactly")
})

# ── (C) delta SE — golden values ─────────────────────────────────────────────
# Delta method is deterministic (no random draws) so the SE is exact.

test_that("snapshot: marginalEffects delta SE golden values", {
  skip_if(!.all_snaps_present,
          "Snapshot cache missing — run with RSENA_REBUILD_MODELS=1")
  snap <- snap_me_delta_transTrip
  expect_false(is.null(snap_se(snap)))
  out <- marginalEffects(ans, mydata,
    targets = snap_tg(ans, mymodel, "mynet", condition = "density",
                      transTrip = list(diff = 1)),
    control_uncertainty = set_postest_uncertainty_saom(
        mode = "delta", sd = TRUE, ci = FALSE, simMean = FALSE, nsim = 1L),
    control_algo = set_postest_algo_saom(verbose = FALSE)
  )
  expect_equal(out$firstDiff, snap$firstDiff, tolerance = 1e-10,
               info = "delta firstDiff must match cached golden values")
  expect_equal(out$SE, snap_se(snap), tolerance = 1e-10,
               info = "delta SE must match cached golden values exactly")
})

# ── (D) predict changeProb — golden values ───────────────────────────────────

test_that("snapshot: predict.sienaFit changeProb golden values", {
  skip_if(!.all_snaps_present,
          "Snapshot cache missing — run with RSENA_REBUILD_MODELS=1")
  snap <- snap_predict_changeProb
  out <- predict(
    object  = ans, newdata = mydata,
    type    = "changeProb",
    level   = "period", condition = "density",
    uncertainty = FALSE, verbose = FALSE
  )
  expect_equal(out$changeProb, snap$changeProb, tolerance = 1e-10,
               info = "changeProb must match cached golden values exactly")
})

# ── (E) predict tieProb — golden values ──────────────────────────────────────

test_that("snapshot: predict.sienaFit tieProb golden values", {
  skip_if(!.all_snaps_present,
          "Snapshot cache missing — run with RSENA_REBUILD_MODELS=1")
  snap <- snap_predict_tieProb
  out <- predict(
    object  = ans, newdata = mydata,
    type    = "tieProb",
    level   = "period", condition = "density",
    uncertainty = FALSE, verbose = FALSE
  )
  expect_equal(out$tieProb, snap$tieProb, tolerance = 1e-10,
               info = "tieProb must match cached golden values exactly")
})

# ── (F) Bootstrap sanity tests — range checks (stochastic, no golden values) ─
# These confirm SE > 0 and that bootstrap and delta SE are in the same ballpark.
# Uses set.seed for reproducibility of which rows are tested.

test_that("snapshot: bootstrap SE > 0 and finite for all period rows", {
  skip_slow()
  skip_if(is.null(ans), "ans fixture not available")
  set.seed(2026L)
  out_boot <- marginalEffects(ans, mydata,
    targets = snap_tg(ans, mymodel, "mynet", condition = "density",
                      transTrip = list(diff = 1)),
    control_uncertainty = set_postest_uncertainty_saom(
        mode = "bootstrap", nsim = 50),
    control_algo = set_postest_algo_saom(verbose = FALSE)
  )
  expect_true("SE" %in% names(out_boot))
  expect_true(all(is.finite(out_boot$SE)),
              info = "All bootstrap SE values must be finite")
  expect_true(all(out_boot$SE > 0),
              info = "All bootstrap SE values must be strictly positive")
})

test_that("snapshot: bootstrap SE and delta SE within 5x of each other", {
  skip_slow()
  skip_if(is.null(ans), "ans fixture not available")
  skip_if(!.all_snaps_present,
          "Snapshot cache missing — run with RSENA_REBUILD_MODELS=1")
  set.seed(2026L)
  out_boot <- marginalEffects(ans, mydata,
    targets = snap_tg(ans, mymodel, "mynet", condition = "density",
                      transTrip = list(diff = 1)),
    control_uncertainty = set_postest_uncertainty_saom(
        mode = "bootstrap", nsim = 100),
    control_algo = set_postest_algo_saom(verbose = FALSE)
  )
  snap_delta <- snap_me_delta_transTrip
  ## Compared directly rather than through merge() suffixes.  This used to
  ## read joined$SE_boot / joined$SE_delta -- names merge() only creates for
  ## columns present in BOTH frames under the same name.  The two frames
  ## called it "SE" and "delta_se", so both were NULL, the ratio was
  ## numeric(0), and all(numeric(0) > 0.2) is TRUE: the test asserted nothing
  ## for as long as it existed.
  expect_equal(nrow(out_boot), nrow(snap_delta),
               info = "bootstrap and delta runs must cover the same rows")
  boot_se  <- out_boot$SE
  delta_se <- snap_se(snap_delta)
  expect_false(is.null(boot_se) || is.null(delta_se))
  expect_gt(length(boot_se), 0L)
  ratio <- boot_se / delta_se
  expect_true(all(ratio > 0.2 & ratio < 5),
              info = paste("Bootstrap / delta SE ratio out of [0.2, 5] range:",
                           paste(round(ratio, 3), collapse = ", ")))
})

# ── (F) delta SE — interaction spec ──────────────────────────────────────────
#
# Tripwire for the analytic-Jacobian generalisation
# (postestimate_api_redesign.md Sec. 2.2).  calculateUtilityDiffJacobian()
# returns NULL for interaction = TRUE, so delta_se below is produced by the
# finite-difference fallback.  Sec. 2.2 replaces that with an analytic
# Jacobian; these numbers must NOT move when it does — the two paths compute
# the same derivative, one exactly and one by central differences.
#
# A mismatch here means either the analytic A-matrix is wrong, or the FD step
# was masking an error.  Tolerance is loosened to 1e-8 (not 1e-10) because
# central FD carries O(eps^2) truncation error, so exact bit-equality is not
# the right expectation once the analytic path lands.

test_that("snapshot: marginalEffects delta SE, interaction spec", {
  skip_if(!.all_snaps_present,
          "Snapshot cache missing — run with RSENA_REBUILD_MODELS=1")
  skip_if(is.null(ans2) || is.null(mydata2),
          "ans2/mydata2 fixtures unavailable")
  snap <- snap_me_delta_interaction
  expect_false(is.null(snap_se(snap)))
  out <- marginalEffects(ans2, mydata2,
    targets = snap_tg(ans2, mymodel2, "mynet2", condition = "recip",
                      dependency = transRecTrip ~ transTrip:recip,
                      transTrip = list(diff = 1)),
    control_uncertainty = set_postest_uncertainty_saom(
        mode = "delta", sd = TRUE, ci = FALSE, nsim = 1L),
    control_algo = set_postest_algo_saom(verbose = FALSE)
  )
  expect_equal(out$firstDiff, snap$firstDiff, tolerance = 1e-10,
               info = "interaction point estimate must not move")
  expect_equal(out$SE, snap_se(snap), tolerance = 1e-8,
               info = "interaction delta SE must not move (FD -> analytic)")
})

# ── (G) delta SE — density spec ──────────────────────────────────────────────
#
# Same tripwire.  Density returns NULL from calculateUtilityDiffJacobian too:
# its utility shift is -2 * changeUtil, so the A-matrix is dense over all K
# columns rather than sparse in one.  That makes it the case most likely to
# expose a mistake in the sparse-representation work (Sec. 2.5.3).

test_that("snapshot: marginalEffects delta SE, density spec", {
  skip_if(!.all_snaps_present,
          "Snapshot cache missing — run with RSENA_REBUILD_MODELS=1")
  snap <- snap_me_delta_density
  expect_false(is.null(snap_se(snap)))
  out <- marginalEffects(ans, mydata,
    targets = snap_tg(ans, mymodel, "mynet",
                      density = list(contrast = c(-1, 1))),
    control_uncertainty = set_postest_uncertainty_saom(
        mode = "delta", sd = TRUE, ci = FALSE, nsim = 1L),
    control_algo = set_postest_algo_saom(verbose = FALSE)
  )
  expect_equal(out$firstDiff, snap$firstDiff, tolerance = 1e-10,
               info = "density point estimate must not move")
  expect_equal(out$SE, snap_se(snap), tolerance = 1e-8,
               info = "density delta SE must not move (FD -> analytic)")
})
