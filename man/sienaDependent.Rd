\name{sienaDependent}
\alias{sienaDependent}
\alias{sienaNet}
\title{Function to create a dependent variable for a Siena model}
\description{
  Creates a Siena dependent variable: either a network, created from a matrix
  or array or list of sparse matrix of triples; or a behavior variable,
  created from a matrix.\cr
  \code{sienaDependent()} and \code{sienaNet()} are identical functions;
  the second name was used from the start of the \code{RSiena} package,
  but the first name indicates more precisely the purpose of this function.
}
\usage{
sienaDependent(netarray, type=c("oneMode", "bipartite", "behavior", "continuous"),
nodeSet="Actors", sparse=is.list(netarray), allowOnly=TRUE, imputationValues=NULL)

sienaNet(netarray, type=c("oneMode", "bipartite", "behavior", "continuous"),
nodeSet="Actors", sparse=is.list(netarray), allowOnly=TRUE, imputationValues=NULL)
}
\arguments{
  \item{netarray}{type="behavior" or "continuous": matrix 
    (actors \eqn{\times}{*} waves).\cr
	type="oneMode" or "bipartite": array of values or list of sparse matrices
    of type "dgTMatrix"; if an array is used, it should have
	dimensions: for a one-network, \eqn{n \times n \times M}{n * n * M},
	and for a two-mode network \eqn{n \times m \times M}{n * m * M}, 
	where \code{n} is the number of actors, \code{m} is the number
	of nodes in the second mode, and \code{M} is the number of waves.}
  \item{type}{type of dependent variable, default "oneMode".}
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
  \item{imputationValues}{for "behavior" or "continuous" dependent variables,
    a matrix with imputation values can be included that will be used instead 
    of the default imputation values.}
}
\details{
  Adds attributes so that the array or list of matrices can be used in a
  Siena model fit.
}
\value{
  An object of class "sienaDependent". An array or (networks only) a list of sparse
  matrices with attributes:
  \item{netdims}{Dimensions of the network or behavior
    variable: senders, receivers (1 for behavior), periods}
  \item{type}{oneMode, bipartite or behavior}
  \item{sparse}{Boolean: whether the network is given as a list of sparse matrices or not}
  \item{nodeSet}{Character string with name(s) of node set(s)}
  \item{allowOnly}{The value of the \code{allowOnly} parameter}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/} .}
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
tmps501 <- as(Matrix(s501), "dgTMatrix")
tmps502 <- as(Matrix(s502), "dgTMatrix")
tmps503 <- as(Matrix(s503), "dgTMatrix")
## If you are interested in the internal structure of these sparse matrices,
## you can request
str(tmps501)
## Slot @i is the row, @j is the column, and @x the value;
## here the values all are 1.
## The edge lists now are created as matrices in the following way,
## where the columns are sender, receiver, tie value, and wave:
mymat1 <- cbind(tmps501@i + 1, tmps501@j + 1, 1, 1)
mymat2 <- cbind(tmps502@i + 1, tmps502@j + 1, 1, 2)
mymat3 <- cbind(tmps503@i + 1, tmps503@j + 1, 1, 3)
mymat <- rbind(mymat1, mymat2, mymat3)
## mymat is an edge list and includes all 3 waves.
## The edge lists do not contain information about the number of nodes,
## so that has to be supplied additionally; here it is 50:
mymatlist <- by( mymat, mymat[, 4], function(x)
    spMatrix(50, 50, x[, 1], x[, 2], x[, 3]) )
mynet4 <- sienaDependent(mymatlist)
## or alternatively
mymat1 <- mymat[mymat[, 4] == 1, ]
mymat2 <- mymat[mymat[, 4] == 2, ]
mymat3 <- mymat[mymat[, 4] == 3, ]
mymat1s <- spMatrix(50, 50, mymat1[, 1], mymat1[, 2], mymat1[, 3])
mymat2s <- spMatrix(50, 50, mymat2[, 1], mymat2[, 2], mymat2[, 3])
mymat3s <- spMatrix(50, 50, mymat3[, 1], mymat3[, 2], mymat3[, 3])
mynet4 <- sienaDependent(list(mymat1s, mymat2s, mymat3s))
}
\keyword{classes}
