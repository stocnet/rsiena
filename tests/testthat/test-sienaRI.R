# Tests for interpret_size (relative importance) and entropy
# Models: ans, mydata, mymodel, mycontrols — from helper-models.R

# ── interpret_size.sienaFit ──────────────────────────────────────────────────

test_that("interpret_size.sienaFit returns new-format sienaRI", {
  ri <- interpret_size(ans, data = mydata, depvar = "mynet",
                       uncertainty = FALSE)
  expect_s3_class(ri, "sienaRI")
  expect_true(!is.null(ri$data))
  expect_true(is.data.frame(ri$data))
  # RI columns for each non-rate effect
  riCols <- grep("^RI_", names(ri$data), value = TRUE)
  expect_true(length(riCols) >= 1)
  # sigma columns
  sigCols <- grep("^sigma_", names(ri$data), value = TRUE)
  expect_equal(length(sigCols), length(riCols))
  # effectNames populated
  expect_true(length(ri$effectNames) >= 1)
  expect_true(length(ri$shortEffectNames) >= 1)
  expect_equal(ri$dependentVariable, "mynet")
})

test_that("interpret_size.sienaFit RI values are in [0,1]", {
  ri <- interpret_size(ans, data = mydata, uncertainty = FALSE)
  riCols <- grep("^RI_", names(ri$data), value = TRUE)
  riDF <- as.data.frame(ri$data)
  riVals <- unlist(riDF[, riCols])
  riVals <- riVals[!is.na(riVals)]
  expect_true(all(riVals >= 0 & riVals <= 1))
})

test_that("interpret_size.sienaFit RI rows sum to ~1", {
  ri <- interpret_size(ans, data = mydata, uncertainty = FALSE)
  riCols <- grep("^RI_", names(ri$data), value = TRUE)
  riDF <- as.data.frame(ri$data)
  riMat <- as.matrix(riDF[, riCols])
  rowSums_ri <- rowSums(riMat, na.rm = TRUE)
  # Rows with non-NA values should sum close to 1
  valid <- !is.na(rowSums_ri) & rowSums_ri > 0
  expect_true(all(abs(rowSums_ri[valid] - 1) < 1e-10))
})

test_that("interpret_size.sienaFit with uncertainty adds uncertainty slot", {
  ri <- interpret_size(ans, data = mydata, uncertainty = TRUE,
                       nsim = 5, verbose = FALSE)
  expect_s3_class(ri, "sienaRI")
  expect_true(!is.null(ri$uncertainty))
  expect_true(is.data.frame(ri$uncertainty))
  expect_true("Mean" %in% names(ri$uncertainty))
})

test_that("interpret_size rejects endowment effects", {
  mymodel_endow <- getEffects(mydata)
  mymodel_endow <- includeEffects(mymodel_endow, transTrip, name = "mynet",
                                   type = "endow")
  expect_error(
    interpret_size(ans, data = mydata, effects = mymodel_endow,
                   uncertainty = FALSE),
    "endowment or creation"
  )
})

# ── relativeImportance.sienaFit ──────────────────────────────────────────────

test_that("relativeImportance.sienaFit returns sienaRI identical to interpret_size", {
  ri_is <- interpret_size(ans, data = mydata, depvar = "mynet",
                          uncertainty = FALSE)
  ri_ri <- relativeImportance(ans, data = mydata, depvar = "mynet",
                              uncertainty = FALSE)
  expect_s3_class(ri_ri, "sienaRI")
  expect_equal(ri_ri$data, ri_is$data)
  expect_equal(ri_ri$effectNames, ri_is$effectNames)
  expect_equal(ri_ri$shortEffectNames, ri_is$shortEffectNames)
  expect_equal(ri_ri$dependentVariable, ri_is$dependentVariable)
})

test_that("relativeImportance with getChangeStats attaches contributions", {
  ri <- relativeImportance(ans, data = mydata, depvar = "mynet",
                           uncertainty = FALSE, getChangeStats = TRUE)
  expect_s3_class(ri, "sienaRI")
  expect_true(!is.null(ri$changeStatistics))
  expect_true(!is.null(ri$changeStatistics$contribMat))
  expect_true(!is.null(ri$changeStatistics$effectNames))
})

test_that("interpret_size.sienaFit with getChangeStats=TRUE attaches contributions", {
  ri <- interpret_size(ans, data = mydata, depvar = "mynet",
                       uncertainty = FALSE, getChangeStats = TRUE)
  expect_s3_class(ri, "sienaRI")
  expect_true(!is.null(ri$changeStatistics))
})

test_that("relativeImportance rejects endowment effects", {
  mymodel_endow <- getEffects(mydata)
  mymodel_endow <- includeEffects(mymodel_endow, transTrip, name = "mynet",
                                   type = "endow")
  expect_error(
    relativeImportance(ans, data = mydata, effects = mymodel_endow,
                       uncertainty = FALSE),
    "endowment or creation"
  )
})

# ── interpret_size.sienaEffects (legacy path) ────────────────────────────────

test_that("interpret_size.sienaEffects returns new-format sienaRI", {
  ri_old <- interpret_size(mymodel, data = mydata, theta = ans$theta)
  expect_s3_class(ri_old, "sienaRI")
  # Now returns new format with $data data.frame (same as sienaFit path)
  expect_true(!is.null(ri_old$data))
  expect_true(is.data.frame(as.data.frame(ri_old$data)))
})

# ── print / summary / plot ───────────────────────────────────────────────────

test_that("print.sienaRI works for new format", {
  ri <- interpret_size(ans, data = mydata, uncertainty = FALSE)
  expect_output(print(ri), "Expected relative importance")
})

test_that("summary.sienaRI works for new format", {
  ri <- interpret_size(ans, data = mydata, uncertainty = FALSE)
  s <- summary(ri)
  expect_s3_class(s, "summary.sienaRI")
})

test_that("plot.sienaRI works for new format (no error)", {
  ri <- interpret_size(ans, data = mydata, uncertainty = FALSE)
  expect_no_error(plot(ri))
})

# ── entropy.sienaFit ─────────────────────────────────────────────────────────

test_that("entropy.sienaFit returns data.frame with R_entropy", {
  ent <- entropy(ans, data = mydata, depvar = "mynet",
                 uncertainty = FALSE, verbose = FALSE)
  expect_true(is.data.frame(ent))
  expect_true("R_entropy" %in% names(ent) || "Mean" %in% names(ent))
})

test_that("entropy.sienaFit with uncertainty returns MCSE/CI columns", {
  skip_slow()
  ent <- entropy(ans, data = mydata, nsim = 5,
                 uncertainty = TRUE, verbose = FALSE)
  expect_true(is.data.frame(ent))
  expect_true("Mean" %in% names(ent))
})

# ── entropy.default ──────────────────────────────────────────────────────────

test_that("entropy.default computes certainty correctly", {
  # Uniform distribution => certainty = 0
  p_uniform <- rep(1/4, 4)
  expect_equal(entropy(p_uniform), 0, tolerance = 1e-10)
  # Degenerate distribution => certainty = NA (single non-NA choice)
  p_one <- c(1, NA, NA)
  expect_true(is.na(entropy(p_one)))
  # Skewed distribution => certainty > 0
  p_skew <- c(0.9, 0.05, 0.05)
  expect_true(entropy(p_skew) > 0)
  expect_true(entropy(p_skew) < 1)
})

# ── getChangeStatistics ──────────────────────────────────────────────────────

test_that("getChangeStatistics static returns data.frame with effect columns", {
  cs <- getChangeStatistics(ans, data = mydata, depvar = "mynet")
  expect_true(is.data.frame(cs))
  expect_true("period" %in% names(cs))
  expect_true("ego" %in% names(cs))
  expect_true("choice" %in% names(cs))
  # Effect columns present
  effCols <- grep("_eval$", names(cs), value = TRUE)
  expect_true(length(effCols) >= 1)
})

test_that("getChangeStatistics static with flip=FALSE returns raw contributions", {
  cs_flip   <- getChangeStatistics(ans, data = mydata, depvar = "mynet",
                                    flip = TRUE)
  cs_noflip <- getChangeStatistics(ans, data = mydata, depvar = "mynet",
                                    flip = FALSE)
  # With flip, non-density columns on no-change rows get negated,
  # so the two should differ
  denCol <- grep("density", names(cs_flip), value = TRUE)[1]
  noChange <- cs_noflip[[denCol]] == -1
  if (any(noChange)) {
    nonDenCols <- setdiff(grep("_eval$", names(cs_flip), value = TRUE), denCol)
    if (length(nonDenCols) > 0) {
      expect_false(identical(cs_flip[[nonDenCols[1]]][noChange],
                              cs_noflip[[nonDenCols[1]]][noChange]))
    }
  }
})

test_that("getChangeStatistics dynamic returns data.frame", {
  skip_slow()
  cs_dyn <- getChangeStatistics(ans, data = mydata, depvar = "mynet",
                                 dynamic = TRUE, algorithm = mycontrols,
                                 n3 = 50)
  expect_true(is.data.frame(cs_dyn))
  expect_true("chain" %in% names(cs_dyn))
  expect_true("ministep" %in% names(cs_dyn))
})

# ── Rcpp helper correctness: compare against R reference implementation ───────

test_that("loo_change_probs matches R-loop reference", {
  contrib <- RSiena:::getStaticChangeContributions(
    ans = ans, data = mydata, effects = ans$requestedEffects,
    depvar = "mynet", returnWide = TRUE
  )
  theta_use <- ans$theta[seq_along(contrib$effectNames)]
  names(theta_use) <- contrib$effectNames
  cMat     <- contrib$contribMat
  gid      <- contrib$group_id
  nEff     <- ncol(cMat)

  # R reference: nEff separate matrix multiplies + grouped softmax
  util_full  <- as.numeric(cMat %*% theta_use)
  loo_ref    <- matrix(NA_real_, nrow = nrow(cMat), ncol = nEff)
  for (k in seq_len(nEff)) {
    th_k    <- theta_use; th_k[k] <- 0
    util_k  <- as.numeric(cMat %*% th_k)
    loo_ref[, k] <- as.numeric(
      RSiena:::softmax_arma_by_group(util_k, gid))
  }

  # Rcpp: one call
  loo_cpp <- RSiena:::loo_change_probs(cMat, theta_use, as.integer(gid))

  expect_equal(loo_cpp, loo_ref, tolerance = 1e-12)
})

test_that("l1d_grouped matches per-group R-loop reference", {
  contrib <- RSiena:::getStaticChangeContributions(
    ans = ans, data = mydata, effects = ans$requestedEffects,
    depvar = "mynet", returnWide = TRUE
  )
  theta_use <- ans$theta[seq_along(contrib$effectNames)]
  names(theta_use) <- contrib$effectNames
  cMat <- contrib$contribMat
  gid  <- contrib$group_id
  nEff <- ncol(cMat)

  util_full <- as.numeric(cMat %*% theta_use)
  fullProb  <- as.numeric(RSiena:::softmax_arma_by_group(util_full, gid))
  looProb   <- RSiena:::loo_change_probs(cMat, theta_use, as.integer(gid))

  uG  <- unique(gid)
  nG  <- length(uG)
  l1d_ref <- matrix(NA_real_, nrow = nG, ncol = nEff)
  for (i in seq_along(uG)) {
    rows <- which(gid == uG[i])
    if (sum(!is.na(fullProb[rows])) <= 1L) next
    for (k in seq_len(nEff)) {
      l1d_ref[i, k] <- sum(abs(fullProb[rows] - looProb[rows, k]),
                            na.rm = TRUE)
    }
  }

  l1d_cpp <- RSiena:::l1d_grouped(fullProb, looProb, as.integer(gid))
  expect_equal(l1d_cpp, l1d_ref, tolerance = 1e-12)
})

test_that("computeRelativeImportance L1D matches R-loop reference", {
  contrib <- RSiena:::getStaticChangeContributions(
    ans = ans, data = mydata, effects = ans$requestedEffects,
    depvar = "mynet", returnWide = TRUE
  )
  thetaHat <- RSiena:::nameThetaFromEffects(ans$theta, ans$requestedEffects)

  # New C++-backed path
  ri_new <- RSiena:::computeRelativeImportance(contrib, thetaHat,
                                               distFun = "L1D")
  ri_df  <- as.data.frame(ri_new)
  riCols <- grep("^RI_", names(ri_df), value = TRUE)
  ri_mat <- as.matrix(ri_df[, riCols])

  # R reference loop
  theta_use <- thetaHat[contrib$effectNames]
  cMat <- contrib$contribMat; gid <- contrib$group_id; nEff <- ncol(cMat)
  util_full <- as.numeric(cMat %*% theta_use)
  fullProb  <- as.numeric(RSiena:::softmax_arma_by_group(util_full, gid))
  loo_ref   <- matrix(NA_real_, nrow = nrow(cMat), ncol = nEff)
  for (k in seq_len(nEff)) {
    th_k <- theta_use; th_k[k] <- 0
    loo_ref[, k] <- as.numeric(
      RSiena:::softmax_arma_by_group(as.numeric(cMat %*% th_k), gid))
  }
  uG <- unique(gid); nG <- length(uG)
  ri_ref <- matrix(NA_real_, nrow = nG, ncol = nEff)
  for (i in seq_len(nG)) {
    rows <- which(gid == uG[i])
    ref  <- fullProb[rows]
    if (sum(!is.na(ref)) <= 1L) next
    l1d <- vapply(seq_len(nEff),
                  function(k) sum(abs(ref - loo_ref[rows, k]), na.rm = TRUE),
                  numeric(1))
    s <- sum(l1d, na.rm = TRUE)
    ri_ref[i, ] <- if (s > 0) l1d / s else NA_real_
  }

  # Values match within floating-point tolerance (strip dimnames before compare)
  valid <- !is.na(ri_ref[, 1L])
  expect_equal(unname(ri_mat[valid, , drop = FALSE]),
               unname(ri_ref[valid, , drop = FALSE]),
               tolerance = 1e-10)
})
