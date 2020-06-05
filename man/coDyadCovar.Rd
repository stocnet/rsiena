\name{coDyadCovar}
\alias{coDyadCovar}
\title{Function to create a constant dyadic covariate object.}
\description{
  This function creates a constant dyadic covariate object from a matrix.
}
\usage{
coDyadCovar(val, centered=TRUE, nodeSets=c("Actors", "Actors"),
    sparse=is(val,"dgTMatrix"), type=c("oneMode", "bipartite"))
}
\arguments{
  \item{val}{Matrix of covariate values. May be sparse, of type "dgTMatrix". }
  \item{centered}{Boolean: if TRUE, then the mean value is subtracted.}
  \item{nodeSets}{The name of the node sets with which this
    covariate is associated. If the entire data set contains more than
    one node set, then the node sets must be specified in all data objects.}
  \item{sparse}{Boolean: whether a sparse matrix or not.}
  \item{type}{oneMode or bipartite: whether the matrix refers to a one-mode
  or a bipartite (two-mode) network.}
}
\details{
  When part of a Siena data object, the covariate is assumed to be
  associated with the node sets named in nodeSets
  of the Siena data object. The name of the associated node sets will only
  be checked when the Siena data object is created.
}
\value{
   Returns the covariate as an object of class "coDyadCovar", in which
  form it can be used as an argument to \code{\link{sienaDataCreate}}.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{sienaDataCreate}}, \code{\link{varDyadCovar}},
 \code{\link{coCovar}}, \code{\link{varCovar}}}
\examples{

mydyadvar <- coDyadCovar(s503)

}
\keyword{classes}
