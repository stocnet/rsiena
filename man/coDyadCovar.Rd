\name{coDyadCovar}
\alias{coDyadCovar}
\title{Function to create a constant dyadic covariate object}
\description{
  This function creates a constant dyadic covariate object from a matrix.\cr
   This function now is superseded by function \code{\link{as_covariate_rsiena}}.
}
\usage{
coDyadCovar(val, centered=TRUE, nodeSets=c("Actors", "Actors"),
     warn=TRUE, sparse=inherits(val,"TsparseMatrix"), 
     type=c("oneMode", "bipartite"))
}
\arguments{
  \item{val}{Matrix of covariate values. May be sparse, of type
     \code{"TsparseMatrix"}. }
  \item{centered}{Boolean: if TRUE, then the mean value is subtracted.}
  \item{nodeSets}{The name of the node sets with which this
    covariate is associated. If the entire data set contains more than
    one node set, then the node sets must be specified in all data objects.}
  \item{warn}{Logical: is a warning given if all values are \code{NA},
  or all non-missing values are the same.}
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
  form it can be used as an argument to \code{\link{make_data_rsiena}}.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{as_covariate_rsiena}},
                    \code{\link{as_nodeset_rsiena}}}
\examples{
mydyadvar <- as_covariate_rsiena(s503, type="oneMode")
class(mydyadvar)
}
\keyword{classes}
