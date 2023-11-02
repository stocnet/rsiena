\name{sienaDependent}
\alias{sienaDependent}
\alias{sienaNet}
\title{Function to create a dependent variable for a Siena model}
\description{
  Creates a Siena dependent variable: either a network,  
  created from a matrix or array or list of sparse matrix of triples;   
  or a behavior variable, created from a matrix.\cr   
  \code{sienaDependent()} and \code{sienaNet()} are identical   
  functions; the second name was used from the start of the   
  \code{RSiena} package, but the first name    
  indicates more precisely the purpose of this function.
}
\usage{
sienaDependent(netarray, type=c("oneMode", "bipartite", "behavior", "continuous"),
nodeSet="Actors", sparse=is.list(netarray), allowOnly=TRUE, imputationValues=NULL)

sienaNet(netarray, type=c("oneMode", "bipartite", "behavior", "continuous"),
nodeSet="Actors", sparse=is.list(netarray), allowOnly=TRUE, imputationValues=NULL)
}
\arguments{
  \item{netarray}{\code{type="behavior"} or \code{"continuous"}: matrix
    (actors \eqn{\times}{*} waves).\cr
    \code{type="oneMode"} or \code{"bipartite"}: array of values or list of
    sparse matrices of type \code{"TsparseMatrix"},
	see the \code{\link{Matrix}} package;
    if an array is used, it should have dimensions:
    for a one-mode network, \eqn{n \times n \times M}{n * n * M},
    and for a two-mode network \eqn{n \times m \times M}{n * m * M},
    where \eqn{n} is the number of actors, \eqn{m} is the number
    of nodes in the second mode, and \eqn{M} is the number of waves.}
  \item{type}{type of dependent variable, default \code{oneMode}.}
  \item{nodeSet}{character string naming the appropriate node set.
    For a bipartite network, a vector containing 2 character strings:
    "rows" first, then "columns".}
  \item{sparse}{logical: TRUE indicates the data is in sparse matrix
    format, FALSE otherwise.}
  \item{allowOnly}{logical: If TRUE, it will be detected when between any
    two consecutive waves the changes are non-decreasing or non-increasing,
    and if this is the case, this will also be a constraint for the
    simulations between these two waves.
    This is done by means of the internal parameters \code{uponly} and
    \code{downonly}.
    If FALSE, the parameters \code{uponly} and \code{downonly} always are
    set to FALSE, and changes in dependent variables will not be constrained
    to be non-decreasing or non-increasing.
    This also will imply that some effects are excluded because they are
    superfluous in such constrained situations.
    This will be reported in the output of \code{\link{print01Report}}.

    For normal operation when this is the case for all periods, usually
    TRUE is the appropriate option. When it is only the case for some of the
    periods, and for data sets that will be part of a multi-group object
    created by \code{\link{sienaGroupCreate}}, FALSE usually is preferable.}
  \item{imputationValues}{for \code{behavior} or \code{continuous} dependent variables,
    a matrix with imputation values can be included that will be used instead
    of the default imputation values.}
}
\details{
  Adds attributes so that the array or list of matrices can be used in a
  Siena model fit.
}
\value{
  An object of class \code{sienaDependent}. An array or (networks only) a list of
  sparse matrices with attributes:
  \item{netdims}{Dimensions of the network or behavior
    variable: senders, receivers (1 for behavior), periods}
  \item{type}{oneMode, bipartite or behavior}
  \item{sparse}{Boolean: whether the network is given as a list of sparse
    matrices or not}
  \item{nodeSet}{Character string with name(s) of node set(s)}
  \item{allowOnly}{The value of the \code{allowOnly} parameter}
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/} .}
\author{Ruth Ripley and Tom A.B. Snijders}
\seealso{\code{\link{sienaDataCreate}}, \code{\link{sienaNodeSet}},
   \code{\link{sienaDataConstraint}}}
\examples{
mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- sienaDependent(s50a, type="behavior")
## note that the following example works although the node sets do not yet exist!
mynet3 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)),
       type="bipartite", nodeSet=c("senders", "receivers"))
## sparse matrix input
## To show this, we first go back from the adjacency matrices to edgelists.
## The manual shows one way to do this.
## Another way is to use the sparse matrix representation which internally
## indeed is an edge list:
library(Matrix)
sp501 <- as(Matrix(s501), "TsparseMatrix")
sp502 <- as(Matrix(s502), "TsparseMatrix")
sp503 <- as(Matrix(s503), "TsparseMatrix")
## If you are interested in the internal structure of these sparse matrices,
## you can request
str(sp501)
## Slot @i is the row, @j is the column, and @x the value;
## here the values all are 1.
## Slots @i and @j do not contain information about the number of nodes,
## so that is supplied additionally by @Dim.
mymatlist <- list(sp501, sp502, sp503)
mynet.sp <- sienaDependent(mymatlist)
}
\keyword{classes}
