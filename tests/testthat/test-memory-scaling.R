# testthat::skip_on_cran()

# library(testthat)
# library(RSiena)

test_that("compute_memory_scale increases with actor/choice size", {
  d_small <- list(depvars = list(net = matrix(0, nrow = 20, ncol = 20)))
  d_large <- list(depvars = list(net = matrix(0, nrow = 100, ncol = 100)))

  s_small <- compute_memory_scale(d_small, depvar = "net", dynamic = FALSE)
  s_large <- compute_memory_scale(d_large, depvar = "net", dynamic = FALSE)

  expect_gte(s_small, 1L)
  expect_gt(s_large, s_small)
})

test_that("compute_memory_scale dynamic includes n3", {
  d <- list(depvars = list(net = matrix(0, nrow = 50, ncol = 50)))

  s_low_n3 <- compute_memory_scale(d, depvar = "net", dynamic = TRUE, n3 = 250)
  s_high_n3 <- compute_memory_scale(d, depvar = "net", dynamic = TRUE, n3 = 2000)

  expect_gte(s_low_n3, 1L)
  expect_gt(s_high_n3, s_low_n3)
})

test_that("resolve_batch_size shrinks with memory_scale and respects override", {
  b_low_scale <- resolve_batch_size(
    nsim = 120,
    nbrNodes = 2,
    batch_size = NULL,
    memory_scale = 1
  )

  b_high_scale <- resolve_batch_size(
    nsim = 120,
    nbrNodes = 2,
    batch_size = NULL,
    memory_scale = 6
  )

  b_override <- resolve_batch_size(
    nsim = 120,
    nbrNodes = 2,
    batch_size = 7,
    memory_scale = 999
  )

  expect_gt(b_low_scale, b_high_scale)
  expect_equal(b_override, 7L)
})
