\name{set_postest_control_saom}
\alias{set_postest_control_saom}
\alias{print.sienaPostestControl}
\alias{sienaPostestControl}
\title{Configure simulation and parallelisation control for SAOM postestimation}
\description{
  \code{set_postest_control_saom} constructs a \code{sienaPostestControl}
  object bundling the settings that control \emph{how} postestimation
  quantities (marginal effects, predictions) are computed by
  \code{\link{marginalEffects}} (and, prospectively,
  \code{\link{predict.sienaFit}}): static vs.\ dynamic simulation, chain
  storage and batching, parallelisation across a cluster, memory budgeting,
  and logging.
}
\usage{
set_postest_control_saom(
    dynamic                = FALSE,
    algorithm              = NULL,
    n3                     = 200,
    n3PointEst             = NULL,
    n3BatchSize            = 100L,
    chainStoreMode         = c("auto", "disk", "memory"),
    useChangeContributions = FALSE,
    chainStorePath         = NULL,
    useCluster             = FALSE,
    nbrNodes               = 1,
    clusterType            = c("PSOCK", "FORK"),
    cl                     = NULL,
    batchDir               = "temp",
    prefix                 = "simBatch_b",
    combineBatch           = TRUE,
    batchSize              = NULL,
    keepBatch              = FALSE,
    verbose                = TRUE,
    details                = FALSE,
    memoryScale            = NULL,
    batchUnitBudget        = 2.5e+08,
    dynamicMinistepFactor  = 10,
    saveDir                = NULL,
    gcEachBatch            = FALSE,
    gcEachSim              = FALSE
)

\method{print}{sienaPostestControl}(x, \dots)
}
\arguments{
  \item{dynamic}{Logical: if \code{TRUE}, use dynamic contributions from
      Markov chains simulated between observations (via
      \code{\link{siena07}}) rather than static change contributions
      extracted at the observation moments. Requires \code{algorithm}.}
  \item{algorithm}{A \code{sienaAlgorithm} object (see
      \code{\link{set_algorithm_saom}}), required when \code{dynamic =
      TRUE}. Supplying \code{dynamic = TRUE} without \code{algorithm} is
      an error raised immediately, before any simulation setup cost is
      incurred.}
  \item{n3}{Single positive numeric, coerced to integer: number of chains
      / theta draws simulated in phase 3 for dynamic predictions.}
  \item{n3PointEst}{Optional single positive numeric (coerced to integer):
      number of chains used for the point estimate specifically; defaults
      to \code{n3} when \code{NULL}.}
  \item{n3BatchSize}{Single positive numeric, coerced to integer: batch
      size for chain simulation under \code{dynamic = TRUE}.}
  \item{chainStoreMode}{Character: \code{"auto"}, \code{"disk"}, or
      \code{"memory"} — how simulated chains are cached for dynamic
      predictions.}
  \item{useChangeContributions}{Logical: reuse pre-stored change
      contributions (e.g.\ from a fitted \code{sienaFit} object) instead
      of recomputing them.}
  \item{chainStorePath}{Optional single character string: path for
      on-disk chain storage when \code{chainStoreMode = "disk"}.}
  \item{useCluster}{Logical: parallelise simulation across a cluster. If
      \code{FALSE} but \code{cl} is supplied, this is automatically set to
      \code{TRUE} (see Details) rather than silently running single-core.}
  \item{nbrNodes}{Single positive numeric, coerced to integer: number of
      cluster nodes to use. Overwritten by \code{length(cl)} when \code{cl}
      is supplied and \code{useCluster = FALSE} (see Details).}
  \item{clusterType}{Character: \code{"PSOCK"} or \code{"FORK"} — the
      \code{\link[parallel]{makeCluster}} cluster type.}
  \item{cl}{Optional existing cluster object (as created by
      \code{\link[parallel]{makeCluster}}) to reuse instead of creating a
      new one.}
  \item{batchDir}{Single non-NA character string: directory used for
      on-disk batch files.}
  \item{prefix}{Single non-NA character string: filename prefix for
      on-disk batch files.}
  \item{combineBatch}{Logical: combine batch results into a single object
      after simulation.}
  \item{batchSize}{Optional single positive numeric (coerced to integer):
      number of simulations per batch. When \code{NULL}, a batch size is
      planned automatically from \code{batchUnitBudget},
      \code{dynamicMinistepFactor}, and \code{memoryScale}.}
  \item{keepBatch}{Logical: keep on-disk batch files after combining
      instead of deleting them.}
  \item{verbose}{Print progress information. Unlike the other logical
      settings here, \code{verbose} is used as a \emph{level} elsewhere in
      the package (tested via \code{verbose >= 1}, \code{verbose >= 2}),
      so any single non-NA logical or numeric value is accepted and stored
      unchanged, rather than being restricted to \code{TRUE}/\code{FALSE}.}
  \item{details}{Logical: include detailed intermediate quantities in the
      output.}
  \item{memoryScale}{Optional single positive numeric: manual override for
      the per-unit memory scale used when auto-planning \code{batchSize}.
      \code{NULL} lets the planner estimate it.}
  \item{batchUnitBudget}{Single positive numeric: the memory budget (in
      abstract "units"; see \code{planBatch}) available for batching.
      Not coerced to integer.}
  \item{dynamicMinistepFactor}{Single positive numeric: multiplier used to
      estimate the number of ministeps per simulated chain when
      auto-planning \code{batchSize} under \code{dynamic = TRUE}. Not
      coerced to integer.}
  \item{saveDir}{Optional single non-NA character string: directory in
      which to save intermediate batch results for later inspection.}
  \item{gcEachBatch}{Logical: force garbage collection after each batch
      (trades speed for lower peak memory use).}
  \item{gcEachSim}{Logical: force garbage collection after each individual
      simulation (trades more speed for even lower peak memory use).}
  \item{x}{A \code{sienaPostestControl} object, as returned by
      \code{set_postest_control_saom}.}
  \item{\dots}{Further arguments, currently ignored (present for S3
      consistency with \code{\link{print}}).}
}
\value{
  \code{set_postest_control_saom} returns a list of class
  \code{"sienaPostestControl"} with elements \code{dynamic},
  \code{algorithm}, \code{n3}, \code{n3PointEst}, \code{n3BatchSize},
  \code{chainStoreMode}, \code{useChangeContributions},
  \code{chainStorePath}, \code{useCluster}, \code{nbrNodes},
  \code{clusterType}, \code{cl}, \code{batchDir}, \code{prefix},
  \code{combineBatch}, \code{batchSize}, \code{keepBatch}, \code{verbose},
  \code{details}, \code{memoryScale}, \code{batchUnitBudget},
  \code{dynamicMinistepFactor}, \code{saveDir}, \code{gcEachBatch}, and
  \code{gcEachSim}, and a \code{"version"} attribute giving the RSiena
  package version used to construct it.

  \code{print.sienaPostestControl} returns \code{x}, invisibly; it is
  called for its side effect of printing a human-readable summary.
}
\details{
  A \code{sienaPostestControl} object is passed to
  \code{\link{marginalEffects}} to control \emph{how} postestimation
  quantities are computed: static vs.\ dynamic simulation, chain storage
  and batching, parallelisation across a cluster, memory budgeting, and
  logging.

  \strong{Cluster normalisation}: passing \code{cl} without also setting
  \code{useCluster = TRUE} would otherwise run computations single-core
  silently. To make the intent explicit, if \code{useCluster = FALSE} and
  a non-empty \code{cl} is supplied, \code{useCluster} is set to
  \code{TRUE} and \code{nbrNodes} is set to \code{length(cl)}
  automatically.

  \strong{Dynamic predictions}: \code{dynamic = TRUE} simulates Markov
  chains between observations via \code{\link{siena07}} and requires an
  \code{algorithm} object (see \code{\link{set_algorithm_saom}});
  supplying \code{dynamic = TRUE} without \code{algorithm} is an error
  raised immediately, before any simulation setup cost is incurred.

  \strong{Cost-relevant settings}: \code{n3} (and, for dynamic
  predictions, \code{n3BatchSize}) control how many chains/draws are
  simulated and in what batch sizes; \code{useCluster}, \code{nbrNodes},
  and \code{clusterType} control parallelisation; \code{batchUnitBudget},
  \code{memoryScale}, and \code{dynamicMinistepFactor} control automatic
  batch-size planning when \code{batchSize = NULL}; and
  \code{gcEachBatch}/\code{gcEachSim} trade speed for lower peak memory
  use. These are the settings to check first when a postestimation call
  is slower or more memory-hungry than expected.
}
\seealso{
  \code{\link{set_postest_uncertainty_saom}} for the accompanying
  uncertainty-quantification object, \code{\link{marginalEffects}} for the
  function these settings configure, and \code{\link{set_algorithm_saom}}
  for the \code{algorithm} argument required when \code{dynamic = TRUE}.
}
\examples{
## Defaults: static predictions, single-core
c1 <- set_postest_control_saom()
print(c1)

## Dynamic predictions require an algorithm object
myalgo <- set_algorithm_saom(nsub = 1, n3 = 50, seed = 1293)
c2 <- set_postest_control_saom(dynamic = TRUE, algorithm = myalgo, n3 = 500)
print(c2)

## Parallel over 4 PSOCK workers, garbage-collect after each batch
c3 <- set_postest_control_saom(useCluster = TRUE, nbrNodes = 4,
                                gcEachBatch = TRUE)
print(c3)
}
\keyword{models}
