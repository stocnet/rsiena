# Regression tests for chainStore_simulate batch seeding.
#
# chainStore_simulate varies the RNG seed per batch only when `seedBase` is
# supplied.  Its own comment says the per-batch seed exists so that "n3 above
# batchSize adds chains instead of repeating batch 1".  When the algorithm
# carries a fixed randomSeed -- which is ordinary practice, and what
# sienaAlgorithmCreate(seed = ...) produces -- an UNSEEDED store re-seeds
# identically on every batch and regenerates the SAME chains.  The effective
# number of chains then collapses to n3BatchSize, silently, while the result
# still reports n3.
#
# The controlling test is "batches differ": it is the property the streaming
# path depends on for its sample size to mean what it says.

# dynArgs as the streaming path builds them, with an algorithm that carries a
# FIXED seed -- that is the condition under which the defect appears, so the
# test pins it rather than inheriting whatever the fixture happens to have.
.seedTestDynArgs <- function() {
  ans    <- load_fixture("ans")
  mydata <- load_fixture("mydata")
  mymodel <- load_fixture("mymodel")
  algo   <- load_fixture("mycontrols")
  if (is.null(ans) || is.null(mydata) || is.null(mymodel) || is.null(algo))
    return(NULL)
  algo$randomSeed <- 12345L
  list(args  = list(ans = ans, data = mydata, algorithm = algo,
                    effects = mymodel, depvar = "mynet", returnWide = TRUE),
       theta = c(ans$rate, ans$theta))
}

# Compare two batches by their numeric payload rather than by identical() on
# the whole object, which can differ on attributes that carry no information.
.ccFingerprint <- function(b) {
  m <- if (!is.null(b$contribMat)) b$contribMat else b$changeStats$csMat
  c(nrow = nrow(m), sum = sum(as.numeric(m)),
    cs = sum(as.numeric(m) * seq_len(length(m))))
}

test_that("a seeded store regenerates a given batch identically", {
  skip_on_cran()
  dyn <- .seedTestDynArgs()
  skip_if(is.null(dyn), "fixture unavailable")

  s <- chainStore_simulate(dyn$args, batchSize = 2L, n3Total = 4L,
                           seedBase = 99L)
  a <- s$getBatch(1L, theta = dyn$theta)
  b <- s$getBatch(1L, theta = dyn$theta)
  expect_identical(.ccFingerprint(a), .ccFingerprint(b))
})

test_that("two stores built with the same seedBase agree batch for batch", {
  skip_on_cran()
  dyn <- .seedTestDynArgs()
  skip_if(is.null(dyn), "fixture unavailable")

  s1 <- chainStore_simulate(dyn$args, 2L, 4L, seedBase = 99L)
  s2 <- chainStore_simulate(dyn$args, 2L, 4L, seedBase = 99L)
  for (b in 1:2)
    expect_identical(.ccFingerprint(s1$getBatch(b, theta = dyn$theta)),
                     .ccFingerprint(s2$getBatch(b, theta = dyn$theta)))
})

test_that("different batches hold DIFFERENT chains", {
  # The defect: with a fixed algorithm$randomSeed and no per-batch seed,
  # batch 2 is a bit-for-bit copy of batch 1.
  skip_on_cran()
  dyn <- .seedTestDynArgs()
  skip_if(is.null(dyn), "fixture unavailable")

  s <- chainStore_simulate(dyn$args, 2L, 4L)        # as the streaming
                                                    # point-estimate path calls it
  expect_false(identical(.ccFingerprint(s$getBatch(1L, theta = dyn$theta)),
                         .ccFingerprint(s$getBatch(2L, theta = dyn$theta))))
})

test_that("raising n3 above the batch size adds chains", {
  # Equivalent statement of the same defect at the level a caller sees: the
  # streamed point estimate must not be invariant to n3.
  skip_on_cran()
  dyn <- .seedTestDynArgs()
  skip_if(is.null(dyn), "fixture unavailable")

  # Every batch must be a distinct draw, so four batches must yield four
  # distinct fingerprints. Comparing row counts alone would pass trivially,
  # because changing n3Total changes the NUMBER of batches either way.
  s <- chainStore_simulate(dyn$args, batchSize = 2L, n3Total = 8L)
  fps <- vapply(seq_len(s$nBatches),
                function(b) paste(.ccFingerprint(s$getBatch(b, theta = dyn$theta)),
                                  collapse = "|"),
                character(1L))
  expect_equal(length(unique(fps)), s$nBatches)
})
