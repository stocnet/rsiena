\name{sienaGOF-auxiliary}
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

\title{Auxiliary functions for goodness of fit assessment by \code{\link{sienaGOF}}}
\description{
 The functions given here are auxiliary to function \code{\link{sienaGOF}}  
 which assesses goodness of fit for actor-oriented models.

 The auxiliary functions are, first, some functions of networks or behaviour  
 (i.e., statistics) for which the simulated values for the fitted model  
 are compared to the observed value; second, some extraction functions  
 to extract the observed and simulated networks and/or behaviour  
 from the \code{\link{sienaFit}} object produced by \code{\link{siena07}}  
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
    defaults in \code{\link{sienaGOF}} to "\code{Data1}".}
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
 study the code of \code{OutdegreeDistribution}, \code{IndegreeDistribution},
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
 Since zero values of dc are not counted, it may be advisable
 to code dc so that all non-diagonal values are non-zero,
 and all diagonal values are zero.
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
\seealso{\code{\link{siena07}}, \code{\link{sienaGOF}} }
\examples{
### For use out of the box:

mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[,1:2], type="behavior")
mycov <- c(rep(1:3,16),1,2) # artificial, just for trying
mydycov <- matrix(rep(1:5, 500), 50, 50) # also artificial, just for trying
mydata <- sienaDataCreate(mynet1, mybeh)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTies, cycle3)
# Shorter phases 2 and 3, just for example:
myalgorithm <- sienaAlgorithmCreate(nsub=1, n3=50, seed=122, projname=NULL)
(ans <- siena07(myalgorithm, data=mydata, effects=myeff, returnDeps=TRUE,
   batch=TRUE))

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

gofi <- sienaGOF(ans, IndegreeDistribution, verbose=TRUE, join=TRUE,
  varName="mynet1")
gofi
plot(gofi)

(gofo <- sienaGOF(ans, OutdegreeDistribution, verbose=TRUE, join=TRUE,
    varName="mynet1", cumulative=FALSE))
# cumulative is an example of "\dots".
plot(gofo)

(gofdc <- sienaGOF(ans, dyadicCov, verbose=TRUE, join=TRUE,
    dc=mydycov, varName="mynet1"))
plot(gofdc)

# How to use dyadicCov for ego-alter combinations of a monadic variable:
mycov.egoalter <- outer(10*mycov, mycov ,'+')
diag(mycov.egoalter) <- 0
dim(mycov.egoalter) # 50 * 50 matrix
# This is a dyadic variable indicating ego-alter combinations of mycov.
# This construction works since mycov has integer values
# not outside the interval from 1 to 9 (actually, only 1 to 3).
# All cells of this matrix contain a two-digit number,
# left digit is row (ego) value, right digit is column (alter) value.
# See the top left part of the matrix:
mycov.egoalter[1:10,1:12]
# The number of values is the square of the number of values of mycov;
# therefore, unwise to do this for a monadic covariate with more than 5 values.
gof.mycov <- sienaGOF(ans, dyadicCov, verbose=TRUE, varName="mynet1",
    dc=mycov.egoalter)
plot(gof.mycov)
descriptives.sienaGOF(gof.mycov, showAll=TRUE)

(gofb <- sienaGOF(ans, BehaviorDistribution, varName = "mybeh",
    verbose=TRUE, join=TRUE, cumulative=FALSE))
plot(gofb)

(goftc <- sienaGOF(ans, TriadCensus, verbose=TRUE, join=TRUE,
    varName="mynet1"))
plot(goftc, center=TRUE, scale=TRUE)
# For this type of auxiliary statistics
# it is advisable in the plot to center and scale.
# note the keys at the x-axis (widen the plot if they are not clear).
descriptives.sienaGOF(goftc)

### The mixed triad census for co-evolution of one-mode and two-mode networks:
actors <- sienaNodeSet(50, nodeSetName="actors")
activities <- sienaNodeSet(20, nodeSetName="activities")
onemodenet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)),
                            nodeSet="actors")
# Not meaningful, just for example:
twomodenet <- sienaDependent(array(c(s502[1:50, 1:20], s503[1:50, 1:20]),
                                                        dim=c(50, 20, 2)),
                            type= "bipartite", nodeSet=c("actors", "activities"))
twodata <- sienaDataCreate(onemodenet, twomodenet,
                        nodeSets=list(actors, activities))
twoeff <- getEffects(twodata)
twoeff <- includeEffects(twoeff, outActIntn, name="onemodenet",
                            interaction1="twomodenet")
twoeff <- includeEffects(twoeff, outActIntn, name="twomodenet",
                            interaction1="onemodenet")
twoeff <- includeEffects(twoeff, from, name="onemodenet",
                            interaction1="twomodenet")
twoeff <- includeEffects(twoeff, to, name="twomodenet",
                            interaction1="onemodenet")
twoeff
# Shorter phases 2 and 3, just for example:
twoalgorithm <- sienaAlgorithmCreate(projname=NULL, nsub=1, n3=50,
                                     seed=5634)
(ans <- siena07(twoalgorithm, data=twodata, effects=twoeff, returnDeps=TRUE,
   batch=TRUE))
(gof.two <- sienaGOF(ans, mixedTriadCensus,
                        varName=c("onemodenet", "twomodenet"), verbose=TRUE))
plot(gof.two, center=TRUE, scale=TRUE)

\dontrun{
### Here come some useful functions for building your own auxiliary statistics:
### First an extraction function.

# igraphNetworkExtraction extracts simulated and observed networks
# from the results of a siena07 run.
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
    # sienaGOF wants the observation:
    original <- data[[groupName]]$depvars[[varName]][,,period+1]
    original[missings] <- 0
    returnValue <- graph.adjacency(original)
  }
  else
  {
    missings <- graph.adjacency(missings)
    #sienaGOF wants the i-th simulation:
    returnValue <- graph.difference(
      graph.empty(dimsOfDepVar) +
        edges(t(sims[[i]][[groupName]][[varName]][[period]][,1:2])),
      missings)
  }
  returnValue
}

### Then some auxiliary statistics.

# GeodesicDistribution calculates the distribution of non-directed
# geodesic distances; see ?sna::geodist
# The default for \code{levls} reflects that geodesic distances larger than 5
# do not differ appreciably with respect to interpretation.
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
TriadCensus.sna <- function(i, data, sims, period, groupName, varName, levls=1:16){
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
CliqueCensus<-function (i, obsData, sims, period, groupName, varName, levls = 1:5){
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

# Distribution of Bonacich eigenvalue centrality; see ?igraph::evcent.
EigenvalueDistribution <- function (i, data, sims, period, groupName, varName,
  levls=c(seq(0,1,by=0.125)), cumulative=TRUE){
  require(igraph)
  x <- igraphNetworkExtraction(i, data, sims, period, groupName, varName)
  a <- igraph::evcent(x)$vector
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

## Finally some examples of the three auxiliary statistics constructed above.
mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- sienaDependent(s50a, type="behavior")
mydata <- sienaDataCreate(mynet1, mybeh)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip, cycle3)
myeff <- includeEffects(myeff, outdeg, name="mybeh", interaction1="mynet1")
myeff <- includeEffects(myeff,  outdeg, name="mybeh", interaction1="mynet1")
# Shorter phases 2 and 3, just for example:
myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=200, seed=765, projname=NULL)
(ans2 <- siena07(myalgorithm, data=mydata, effects=myeff, returnDeps=TRUE,
   batch=TRUE))
gofc <- sienaGOF(ans2, EigenvalueDistribution, varName="mynet1",
  verbose=TRUE, join=TRUE)
plot(gofc)
descriptives.sienaGOF(gofc, showAll=TRUE)

goftc <- sienaGOF(ans2, TriadCensus, varName="mynet1", verbose=TRUE, join=TRUE)
plot(goftc, center=TRUE, scale=TRUE)
# For this type of auxiliary statistics
# it is advisable in the plot to center and scale.
# note the keys at the x-axis; these names are given by sna::triad.census
descriptives.sienaGOF(goftc)
round(descriptives.sienaGOF(goftc))

gofgd <- sienaGOF(ans2, GeodesicDistribution, varName="mynet1",
  verbose=TRUE, join=TRUE, cumulative=FALSE)
plot(gofgd)
# and without infinite distances:
gofgdd <- sienaGOF(ans2, GeodesicDistribution, varName="mynet1",
  verbose=TRUE, join=TRUE, levls=1:7, cumulative=FALSE)
plot(gofgdd)
}
}
\keyword{models}
