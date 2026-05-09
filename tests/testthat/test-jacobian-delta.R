## test-jacobian-delta.R
##
## Numeric regression for jacobianCondDelta with the new batched-perturbations
## API (Commit 2).
##
## Strategy: construct a mock estimator with a *known* analytical Jacobian,
## run jacobianCondDelta, and compare to both the analytical answer and to a
## naive one-at-a-time FD reference (what the old code did).
##
## The mock estimator satisfies the perturbations API:
##   estimator(theta, perturbations = NULL)  → named list of data.frames
##   estimator(theta, perturbations = list)  → list(hat=..., perturbations=list(...))
## It does NOT use any chains or siena07 — pure R arithmetic.

context("jacobianCondDelta — batched perturbations API")

make_mock_estimator <- function(f_list) {
  ## f_list: named list of functions theta → numeric vector (one per spec).
  ## Each function produces the 'outcome' column of that spec's data.frame.
  eff_names <- names(f_list)
  est <- function(theta, perturbations = NULL, useChangeContributions = FALSE) {
    eval_at <- function(th) {
      lapply(f_list, function(f) data.frame(outcome = f(th)))
    }
    hat <- eval_at(theta)
    if (is.null(perturbations)) return(hat)
    perts_out <- lapply(perturbations, eval_at)
    list(hat = hat, perturbations = perts_out)
  }
  attr(est, "eff_names") <- eff_names
  est
}

## ── Fixture ──────────────────────────────────────────────────────────────────
## f1(theta) = [theta1^2,  theta1*theta2,  theta2^2]    (3 rows in output)
## f2(theta) = [exp(theta1 + theta2)]                   (1 row)
##
## Analytical Jacobians at theta = (1, 2):
##   J1:  [[2*t1, 0], [t2, t1], [0, 2*t2]]  = [[2,0],[2,1],[0,4]]
##   J2:  [[exp(3), exp(3)]]                 = [[exp(3), exp(3)]]

thetaHat <- c(t1 = 1, t2 = 2)

f1 <- function(th) c(th["t1"]^2,  th["t1"] * th["t2"],  th["t2"]^2)
f2 <- function(th) exp(th["t1"] + th["t2"])

specs <- list(
  s1 = list(outcomeName = "outcome"),
  s2 = list(outcomeName = "outcome")
)

mock_est <- make_mock_estimator(list(s1 = f1, s2 = f2))

## ── 1. Analytical ground truth ───────────────────────────────────────────────
J1_analytic <- matrix(c(2, 0,
                         2, 1,
                         0, 4),
                       nrow = 3, ncol = 2, byrow = TRUE,
                       dimnames = list(NULL, c("t1", "t2")))

J2_analytic <- matrix(c(exp(3), exp(3)),
                       nrow = 1, ncol = 2, byrow = TRUE,
                       dimnames = list(NULL, c("t1", "t2")))

## ── 2. Naive one-at-a-time FD reference (old code behaviour) ─────────────────
naive_jac <- function(estimator, thetaHat, specs, eps = 1e-5) {
  nParams <- length(thetaHat)
  hat_val <- estimator(thetaHat)
  jac_list <- lapply(names(specs), function(specName) {
    oc    <- specs[[specName]]$outcomeName
    n_out <- length(hat_val[[specName]][[oc]])
    J     <- matrix(NA_real_, nrow = n_out, ncol = nParams,
                    dimnames = list(NULL, names(thetaHat)))
    for (k in seq_len(nParams)) {
      tp <- thetaHat; tp[k] <- tp[k] + eps
      tm <- thetaHat; tm[k] <- tm[k] - eps
      ep_oc <- estimator(tp)[[specName]][[oc]]
      em_oc <- estimator(tm)[[specName]][[oc]]
      J[, k] <- (ep_oc - em_oc) / (2 * eps)
    }
    J
  })
  names(jac_list) <- names(specs)
  jac_list
}

naive <- naive_jac(mock_est, thetaHat, specs, eps = 1e-5)

## ── 3. New batched jacobianCondDelta ──────────────────────────────────────────
batched <- jacobianCondDelta(mock_est, thetaHat, specs, eps = 1e-5)

## ── Tests ────────────────────────────────────────────────────────────────────

test_that("jacobianCondDelta returns hat for free", {
  expect_named(batched$hat, c("s1", "s2"))
  expect_equal(batched$hat$s1$outcome,
               unname(f1(thetaHat)), tolerance = 1e-12)
  expect_equal(batched$hat$s2$outcome,
               unname(f2(thetaHat)), tolerance = 1e-12)
})

test_that("jacobianCondDelta matches naive one-at-a-time FD (s1, 3 rows)", {
  expect_equal(batched$jac$s1, naive$s1, tolerance = 1e-10)
})

test_that("jacobianCondDelta matches naive one-at-a-time FD (s2, 1 row)", {
  expect_equal(batched$jac$s2, naive$s2, tolerance = 1e-10)
})

test_that("jacobianCondDelta matches analytical Jacobian for s1 (quadratic)", {
  eps_fd <- 1e-5
  ## FD error for quadratic is O(eps^2) ≈ 1e-10 — tight tolerance OK.
  expect_equal(batched$jac$s1, J1_analytic, tolerance = 1e-8)
})

test_that("jacobianCondDelta matches analytical Jacobian for s2 (exp)", {
  ## FD error for exp is O(eps^2 * exp(3)) ≈ 2e-9.
  expect_equal(batched$jac$s2, J2_analytic, tolerance = 1e-7)
})

test_that("seDeltaRows is consistent with seDelta on each row of s1 Jacobian", {
  covT <- diag(c(0.01, 0.04))
  J    <- batched$jac$s1
  se_rows <- seDeltaRows(J, covT)
  se_scalar <- vapply(seq_len(nrow(J)), function(i)
    seDelta(J[i, ], covT), numeric(1L))
  expect_equal(se_rows, se_scalar, tolerance = 1e-14)
})
