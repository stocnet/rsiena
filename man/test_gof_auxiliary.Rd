\name{test_gof_auxiliary}
\alias{test_gof_auxiliary}
\alias{sienaGOF-auxiliary}
\alias{sparseMatrixExtraction}
\alias{networkExtraction}
\alias{behaviorExtraction}
\alias{OutdegreeDistribution}
\alias{IndegreeDistribution}
\alias{BehaviorDistribution}
\alias{TriadCensus}
\alias{mixedTriadCensus}
\alias{dyadicCov}
\alias{egoAlterCombi}

\title{Auxiliary functions for goodness of fit assessment by \code{\link{test_gof}}}
\description{
 The functions given here are auxiliary to function \code{\link{test_gof}}
 which assesses goodness of fit for actor-oriented models.

 The auxiliary functions are, first, some functions of networks or behaviour
 (i.e., statistics) for which the simulated values for the fitted model
 are compared to the observed value; second, some extraction functions
 to extract the observed and simulated networks and/or behaviour
 from the \code{\link{sienaFit}} object produced by \code{\link{siena}}
 with \code{returnDeps=TRUE}.

 These functions are exported here mainly to enable users to write
 their own versions. At the end of this help page some non-exported functions
 are listed. These are not exported because they depend on packages that
 are not in the R base distribution; and to show templates for readers
 wishing to contruct their own functions.
 }

\usage{
OutdegreeDistribution(i, obsData, sims, period, groupName, varName,
         levls=0:8, cumulative=TRUE)

IndegreeDistribution(i, obsData, sims, period, groupName, varName,
         levls=0:8, cumulative=TRUE)

BehaviorDistribution(i, obsData, sims, period, groupName, varName,
         levls=NULL, cumulative=TRUE)

TriadCensus(i, obsData, sims, period, groupName, varName, levls=1:16)

mixedTriadCensus(i, obsData, sims, period, groupName, varName)

dyadicCov(i, obsData, sims, period, groupName, varName, dc)

egoAlterCombi(i, obsData, sims, period, groupName, varName, trafo=NULL)

sparseMatrixExtraction(i, obsData, sims, period, groupName, varName)

networkExtraction(i, obsData, sims, period, groupName, varName)

behaviorExtraction(i, obsData, sims, period, groupName, varName)
}
\arguments{
  \item{i}{ Index number of simulation to be extracted, ranging from 1 to
    \code{length(sims)}; if NULL, the data observation will be extracted.}
  \item{obsData}{ The observed data set to which the model was fitted; normally
    this is \code{x$f} where \code{x} is the \code{\link{sienaFit}}
    object for which the fit is being assessed.}
  \item{sims}{ The simulated data sets to be compared with the observed data;
    normally this is
    \code{x$sims} where \code{x} is the \code{\link{sienaFit}}
    object for which the fit is being assessed.}
  \item{period}{ Period for which data and simulations are used (may run from
    1 to number of waves - 1).}
  \item{groupName}{ Name of group; relevant for multi-group data sets;
    defaults in \code{\link{test_gof}} to "\code{Data1}".}
  \item{varName}{ Name of dependent variable. }
  \item{levls}{ Levels used as values of the auxiliary statistic.
    For \code{BehaviorDistribution}, this defaults to the observed
    range of values.}
  \item{cumulative}{ Are the distributions to be considered as raw or
    cumulative (<=) distributions?}
  \item{dc}{Dyadic covariate: either a matrix with dimensions
    \eqn{n \times n}{n * n}; or, as period-dependent values,
    an array with dimensions \eqn{n \times n \times (M-1)}{n * n * (M-1)};
    where \eqn{n} is the number of actors and \eqn{M} is the number of waves.
    There may be more time points, but those after \eqn{(M-1)} will not be used.}
  \item{trafo}{ For \code{egoAlterCombi}:
    transformation of the dependent behavior variable.
    Default: identity function. }
}
\details{
 The statistics should be chosen to represent features of the network
 that are not explicitly fit by the estimation procedure but can be
 considered important properties that the model at hand should represent well.
 The three given here are far from a complete set;
 they will be supplemented in due time by statistics depending on networks and
 behavior jointly.
 The examples below give a number of other statistics, using the packages
 \code{sna} and \code{igraph}.

 The \code{levls} parameter must be adapted to the range of values that is
 considered important. For indegrees and outdegrees, the whole range should
 usually be covered. If the range is large, which could be the case, e.g.,
 for indegrees of two-mode networks where the second mode has few nodes,
 think about the possibility of making a selection such as
 \code{levls=5*(0:20)} or \code{levls=c(0:4,5*(1:20))};
 which in most cases will make sense only if \code{cumulative=TRUE}.

 The method signature for the auxiliary statistics generally is\cr
 \code{function(i, obsData, sims, period, groupName, varName, \dots)}.
 For constructing new auxiliary statistics, it is helpful to
 study the code of \code{OutdegreeDistribution}, \cr
 \code{IndegreeDistribution},
 and \code{BehaviorDistribution} and of the example functions below.

 \code{TriadCensus} returns the distribution of the Holland-Leinhardt triad
 census according to the algorithm by Batagelj and Mrvar
 (implementation by Parimalarangan, Slota, and Madduri). An alternative is the
 \code{TriadCensus.sna} function mentioned below, from package \code{sna}, which
 gives the same results. Here the \code{levls} parameter can be used to exclude
 some triads, e.g., for non-directed networks.\cr
 The Batagelj-Mrvar algorithm is optimized for sparse, large graphs and may be
 much faster than the procedure implemented in \code{sna}. For dense graphs
 the \code{sna} procedure may be faster.

 \code{dyadicCov} assumes that \code{dc} is a categorical dyadic variable,
 and returns the frequencies of the non-zero values for realized ties.
 It is recommended to do this for a categorical covariate.
 Since zero values of dc are not counted, it may be advisable
 to code dc so that all non-diagonal values are non-zero,
 and all diagonal values are zero.

 \code{egoAlterCombi} returns the frequencies of realized ties for each
  ego-alter combination of values of the dependent behavior variable,
  transformed by trafo; if \code{NULL}, untransformed.\cr
  The range of this transformation should preferably have 5 or fewer values,
  and definitely not 10 or more, as this would lead to problems in
  the coding of ego-alter combinations.
    }

\value{
 \code{OutdegreeDistribution} returns a named vector, the distribution
  of the observed or simulated outdegrees for the values in \code{levls}.

 \code{IndegreeDistribution} returns a named vector, the distribution
  of the observed or simulated indegrees for the values in \code{levls}.

 \code{BehaviorDistribution} returns a named vector, the distribution of the
  observed or simulated behavioral variable for the values in \code{levls}.

 \code{TriadCensus} returns a named vector, the distribution of the
  Holland-Leinhardt triad census according to the algorithm by
  Batagelj and Mrvar.

 \code{mixedTriadCensus} returns a named vector, the distribution of the
  mixed triad census of Hollway, Lomi, Pallotti, and Stadtfeld (2017).
  See their Figure 1 for the meaning of the codes.
  In this figure, ties between the bottom nodes are for the first network,
  ties from the bottom to the top nodes are for the second network.
  The mixed triad census can be used for pairs of dependent networks of which
  the first must be one-mode and the second can be one-mode or two-mode.
  If the second is one-mode, the set of triads considered is only a subset
  of all mixed triads, and ties in the figure are directed upward;
  existence of other ties is not considered.

 \code{dyadicCov} returns a named vector, the frequencies of the non-missing
  non-zero values dc(ego,alter) of the observed or simulated (ego,alter) ties.
  
 \code{sparseMatrixExtraction} returns the simulated network as a
  \code{"TsparseMatrix"};
  this is the virtual class for sparse numeric matrices
  represented by triplets in the \code{Matrix} package.\cr
  Tie variables for ordered pairs with a missing value for wave=period or
  period+1 are zeroed; note that this also is done
  in \code{RSiena} for calculation of target statistics.
  Tie variables that are structurally determined at the beginning of a period
  are used to replace observed values at the end of the period;
  tie variables that are structurally determined at the end, but not
  the beginning, of a period are used to replace
  simulated values at the end of the period.
  \cr
  To treat the objects returned by this function as regular matrices,
  it is necessary to attach the \code{Matrix} package in your session.

\code{networkExtraction} returns the network as an edge list of class
  \code{network} according to the \code{network} package
  (used for package \code{sna}).
  Missing values and structural values are treated as in
  \code{sparseMatrixExtraction}, see above.

\code{behaviorExtraction} returns the dependent behavior
  variable as an integer vector.
  Values for actors with a missing value for wave=period or period+1 are
  transformed to \code{NA}.
}

\references{
   Batagelj, V., and Mrvar, A. (2001),
  A subquadratic triad census algorithm for large sparse networks with small
  maximum degree. \emph{Social Networks}, \bold{23}, 237--243.

   Holland, P.W., and Leinhardt, S. (1976),
  Local structure in social networks. \emph{Sociological Methodology}, \bold{6}, 1--45.

   Hollway, J., Lomi, A., Pallotti, F., and Stadtfeld, C. (2017),
   Multilevel social spaces: The network dynamics of organizational fields.
   \emph{Network Science}, \bold{5}, 187--212.

   Lospinoso, J.A. and Snijders, T.A.B. (2019),
  Goodness of fit for stochastic actor-oriented models.
  \emph{Methodological Innovations}, \bold{12}:2059799119884282.

   Parimalarangan S., Slota, G.M., and Madduri, K. (2017),
 Fast parallel graph triad census and triangle counting on shared-memory
  platforms, \emph{2017 IEEE International Parallel and Distributed Processing
  Symposium Workshops (IPDPSW)}, Lake Buena Vista, FL, pp. 1500-1509.
}
\author{Josh Lospinoso, Tom Snijders}
\seealso{\code{\link{siena}}, \code{\link{test_gof}} }
\examples{
### For use out of the box:

mynet1 <- as_dependent_rsiena(array(c(s501, s502), dim=c(50, 50, 2)))
mybeh <- as_dependent_rsiena(s50a[,1:2], type="behavior")
mydycov <- matrix(rep(1:5, 500), 50, 50) # artificial, just for trying
mydata <- make_data_rsiena(mynet1, mybeh)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(transTies, cycle3))
# Shorter phases 2 and 3, just for example:
myalgo <- set_algorithm_saom(nsub=1, n3=50, seed=122)
(ans <- siena(mydata, effects=myeff, control_algo=myalgo,
                                returnDeps=TRUE, silent=TRUE, batch=TRUE))
# NULL for the observations:
OutdegreeDistribution(NULL, ans$f, ans$sims, period=1, groupName="Data1",
  levls=0:7, varName="mynet1")
dyadicCov(NULL, ans$f, ans$sims, period=1, groupName="Data1",
  dc=mydycov, varName="mynet1")
# An arbitrary selection for simulation run i:
IndegreeDistribution(5, ans$f, ans$sims, period=1, groupName="Data1",
  varName="mynet1")
BehaviorDistribution(20, ans$f, ans$sims, period=1, groupName="Data1",
  varName="mybeh")
sparseMatrixExtraction(50, ans$f, ans$sims, period=1, groupName="Data1",
  varName="mynet1")
networkExtraction(40, ans$f, ans$sims, period=1, groupName="Data1",
  varName="mynet1")
behaviorExtraction(50, ans$f, ans$sims, period=1, groupName="Data1",
  varName="mybeh")

gofi <- test_gof(ans, IndegreeDistribution, verbose=TRUE, join=TRUE,
  varName="mynet1")
gofi
plot(gofi)

(gofo <- test_gof(ans, OutdegreeDistribution, verbose=TRUE, join=TRUE,
    varName="mynet1", cumulative=FALSE))
# cumulative is an example of "\dots".
plot(gofo)

traf <- function(z){pmax(pmin(z,4),2)}
(gofea <- test_gof(ans, egoAlterCombi, verbose=TRUE, join=TRUE,
    varName=c("mynet1","mybeh"), trafo=traf))
plot(gofea)
descriptives(gofea, showAll=TRUE)

(gofdc <- test_gof(ans, dyadicCov, verbose=TRUE, join=TRUE,
    dc=mydycov, varName="mynet1"))
plot(gofdc)

(gofb <- test_gof(ans, BehaviorDistribution, varName = "mybeh",
    verbose=TRUE, join=TRUE, cumulative=FALSE))
plot(gofb)

(goftc <- test_gof(ans, TriadCensus, verbose=TRUE, join=TRUE,
    varName="mynet1"))
plot(goftc, center=TRUE, scale=TRUE)
# For this type of auxiliary statistics
# it is advisable in the plot to center and scale.
# note the keys at the x-axis (widen the plot if they are not clear).
descriptives(goftc)

### The mixed triad census for co-evolution of one-mode and two-mode networks:
actors <- as_nodeset_rsiena(50, nodeSetName="actors")
activities <- as_nodeset_rsiena(20, nodeSetName="activities")
onemodenet <- as_dependent_rsiena(array(c(s501, s502), dim=c(50, 50, 2)),
                nodeSet="actors")
# Not meaningful, just for example:
twomodenet <- as_dependent_rsiena(
                array(c(s502[1:50, 1:20], s503[1:50, 1:20]), dim=c(50, 20, 2)),
                type= "bipartite", nodeSet=c("actors", "activities"))
twodata <- make_data_rsiena(onemodenet, twomodenet,
                nodeSets=list(actors, activities))
twoeff <- make_specification(twodata)
twoeff <- set_effect(twoeff, outActIntn, depvar="onemodenet",
                covar1="twomodenet")
twoeff <- set_effect(twoeff, outActIntn, depvar="twomodenet",
                covar1="onemodenet")
twoeff <- set_effect(twoeff, from, depvar="onemodenet", covar1="twomodenet")
twoeff <- set_effect(twoeff, to, depvar="twomodenet", covar1="onemodenet")
twoeff
# Shorter phases 2 and 3, just for example:

twoalgo <- set_algorithm_saom(nsub=1, n3=50, seed=1291)
(ans <- siena(twodata, effects=twoeff, control_algo=twoalgo,
                returnDeps=TRUE, silent=TRUE, batch=TRUE))
(gof.two <- test_gof(ans, mixedTriadCensus,
                varName=c("onemodenet", "twomodenet"), verbose=TRUE))
plot(gof.two, center=TRUE, scale=TRUE)

\dontrun{
### Here come some useful functions for building your own auxiliary statistics:
### First an extraction function.

# igraphNetworkExtraction extracts simulated and observed networks
# from the results of a siena run.
# It returns the network as an edge list of class "graph"
# according to the igraph package.
# Ties for ordered pairs with a missing value for wave=period or period+1
# are zeroed;
# note that this also is done in RSiena for calculation of target statistics.
# However, changing structurally fixed values are not taken into account.
igraphNetworkExtraction <- function(i, data, sims, period, groupName, varName) {
  require(igraph)
  dimsOfDepVar <- attr(data[[groupName]]$depvars[[varName]], "netdims")[1]
  missings <- is.na(data[[groupName]]$depvars[[varName]][,,period]) |
    is.na(data[[groupName]]$depvars[[varName]][,,period+1])
  if (is.null(i)) {
    # test_gof wants the observation:
    original <- data[[groupName]]$depvars[[varName]][,,period+1]
    original[missings] <- 0
    returnValue <- graph.adjacency(original)
  }
  else
  {
    missings <- graph_from_adjacency_matrix(missings)
    #test_gof wants the i-th simulation:
    returnValue <- difference(
      make_empty_graph(dimsOfDepVar) +
        edges(t(sims[[i]][[groupName]][[varName]][[period]][,1:2])),
      missings)
  }
  returnValue
}

### Then some auxiliary statistics.

# GeodesicDistribution calculates the distribution of non-directed
# geodesic distances; see ?sna::geodist
# The default for \code{levls} reflects that
# geodesic distances larger than 5 do not differ appreciably
# with respect to interpretation.
# Note that the levels of the result are named;
# these names are used in the \code{plot} method.
GeodesicDistribution <- function (i, data, sims, period, groupName,
  varName, levls=c(1:5,Inf), cumulative=TRUE, ...) {
  x <- networkExtraction(i, data, sims, period, groupName, varName)
  require(network)
  require(sna)
  a <- sna::geodist(symmetrize(x))$gdist
  if (cumulative)
  {
    gdi <- sapply(levls, function(i){ sum(a<=i) })
  }
  else
  {
    gdi <- sapply(levls, function(i){ sum(a==i) })
  }
  names(gdi) <- as.character(levls)
  gdi
}

# Holland and Leinhardt Triad Census from sna; see ?sna::triad.census.
# For undirected networks, call this with levls=1:4
TriadCensus.sna <- function(i, data, sims, period, groupName,
                            varName, levls=1:16){
  unloadNamespace("igraph") # to avoid package clashes
  require(network)
  require(sna)
  x <- networkExtraction(i, data, sims, period, groupName, varName)
  if (network.edgecount(x) <= 0){x <- symmetrize(x)}
  # because else triad.census(x) will lead to an error
  tc <- sna::triad.census(x)[levls]
  # names are transferred automatically
  tc
}

# Holland and Leinhardt Triad Census from igraph; see ?igraph::triad_census.
TriadCensus.i <- function(i, data, sims, period, groupName, varName){
  unloadNamespace("sna") # to avoid package clashes
  require(igraph)
  x <- igraphNetworkExtraction(i, data, sims, period, groupName, varName)
# suppressWarnings is used because else warnings will be generated
# when a generated network happens to be symmetric.
  setNames(suppressWarnings(triad_census(x)),
        c("003", "012", "102", "021D","021U", "021C", "111D", "111U",
          "030T", "030C", "201",  "120D", "120U", "120C", "210", "300"))
}

# CliqueCensus calculates the distribution of the clique census
# of the symmetrized network; see ?sna::clique.census.
CliqueCensus<-function (i, obsData, sims, period, groupName,
                        varName, levls = 1:5){
  require(sna)
  x <- networkExtraction(i, obsData, sims, period, groupName, varName)
  cc0 <- sna::clique.census(x, mode='graph', tabulate.by.vertex = FALSE,
    enumerate=FALSE)[[1]]
  cc <- 0*levls
  names(cc) <- as.character(levls)
  levels.used <- as.numeric(intersect(names(cc0), names(cc)))
  cc[levels.used] <- cc0[levels.used]
  cc
}

# Distribution of Bonacich eigenvalue centrality; see ?igraph::eigen_centrality.
EigenvalueDistribution <- function (i, data, sims, period, groupName,
                   varName, levls=c(seq(0,1,by=0.125)), cumulative=TRUE){
  require(igraph)
  x <- igraphNetworkExtraction(i, data, sims, period, groupName, varName)
  a <- igraph::eigen_centrality(x)$vector
  a[is.na(a)] <- Inf
  lel <- length(levls)
  if (cumulative)
  {
    cdi <- sapply(2:lel, function(i){sum(a<=levls[i])})
  }
  else
  {
    cdi <- sapply(2:lel, function(i){
      sum(a<=levls[i]) - sum(a <= levls[i-1])})
  }
  names(cdi) <- as.character(levls[2:lel])
  cdi
}

## Finally some examples
## of the three auxiliary statistics constructed above.
mynet1 <- as_dependent_rsiena(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- as_dependent_rsiena(s50a, type="behavior")
mydata <- make_data_rsiena(mynet1, mybeh)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(transTrip, cycle3))
myeff <- set_effect(myeff, outdeg, depvar="mybeh", covar1="mynet1")
myeff <- set_effect(myeff,  outdeg, depvar="mybeh", covar1="mynet1")
# Shorter phases 2 and 3, just for example:
myalgo <- set_algorithm_saom(nsub=1, n3=100, seed=1291)
(ans2 <- siena(mydata, effects=myeff, control_algo=myalgo,
           returnDeps=TRUE, batch=TRUE))
gofc <- test_gof(ans2, EigenvalueDistribution, varName="mynet1",
  verbose=TRUE, join=TRUE)
plot(gofc)
descriptives(gofc, showAll=TRUE)

goftc <- test_gof(ans2, TriadCensus, varName="mynet1",
           verbose=TRUE, join=TRUE)
plot(goftc, center=TRUE, scale=TRUE)
# For this type of auxiliary statistics
# it is advisable in the plot to center and scale.
# note the keys at the x-axis; these names are given by sna::triad.census
descriptives.test_gof(goftc)
round(descriptives.test_gof(goftc))

gofgd <- test_gof(ans2, GeodesicDistribution, varName="mynet1",
           verbose=TRUE, join=TRUE, cumulative=FALSE)
plot(gofgd)
# and without infinite distances:
gofgdd <- test_gof(ans2, GeodesicDistribution, varName="mynet1",
           verbose=TRUE, join=TRUE, levls=1:7, cumulative=FALSE)
plot(gofgdd)
}
}
\keyword{models, methods}
