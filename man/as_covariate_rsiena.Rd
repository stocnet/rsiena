\name{as_covariate_rsiena}
\alias{as_covariate_rsiena}
\alias{covariate}
\title{Function to create a covariate object}
\description{
  This function creates a monadic or dyadic covariate object
}
\usage{
as_covariate_rsiena(val, type="monadic", centered=TRUE,
                    nodeSet="Actors", warn=TRUE,
                    imputationValues=NULL)
}
\arguments{
  \item{val}{Covariate values, given as a vector, or matrix 
  (may be sparse, of type \code{"TsparseMatrix"}), \cr or 3-dimensional array,
   or list of 3 sparse matrices of type \code{"TsparseMatrix"};\cr
   which of these structures is given
    determines the kind of covariate (see Details).}
  \item{type}{\code{"monadic"}, \code{"oneMode"} or \code{"bipartite"}:
    whether the  values refer to a monadic covariate, \cr
	and in case of a dyadic covariate, whether it is \cr
	 a one-mode or a bipartite (two-mode) dyadic covariate (see Details).}
  \item{centered}{Boolean: if TRUE, then the mean value is subtracted.}
  \item{nodeSet}{Name of node set: character string;\cr
    for \code{type="bipartite"}, referring to  bipartite (two-mode) dyadic covariates, \cr
	a vector of two node set names  (character strings).\cr
    If the entire data set contains more than one node set,\cr
    then the node sets must be specified in all data objects.}
  \item{warn}{Logical: determines whether a warning given if \cr
  all values are \code{NA}, or all non-missing values are the same.}
  \item{imputationValues}{For monadic covariates only:\cr
    Vector or matrix of covariate values of same length or dimensions as \code{val}, \cr
	to be used for imputation of \code{NA} values (if any) in \code{val}. \cr
	Must not contain  any \code{NA}.}
}
\details{This function determines the kind of covariate
  from the kind of arguments given:
  \itemize{
  \item if \code{val} is a vector,
         it will be a time-constant monadic (actor) covariate;
  \item if \code{val} is a matrix and \code{ type="monadic"},
     it will be a changing monadic (actor) covariate;
  \item if \code{val} is a matrix and
     \code{type="oneMode"} or \code{"bipartite"},
     it will be a time-constant dyadic covariate;
  \item if \code{val} is a  3-dimensional array or list of sparse matrices
    of type \code{"TsparseMatrix"}, it will be a changing dyadic covariate.
 }

  When part of a \code{siena} data object,
  the created covariate is associated with the node set
  \code{nodeSet} of this data object.\cr
  A dyadic covariate is a matrix where the rows correspond to one node set, e.g., \code{A}, \cr
  and the columns to a node set, e.g., \code{B};\cr
  if these are the same (\code{A=B}) the dyadic covariate is one-mode,\cr
  if they are different the dyadic covariate is bipartite (two-mode).
  
  For one-mode dyadic covariates one node set should be given for the \code{nodeSet} keyword, \cr
  and for two-mode dyadic covariates 
  two node sets should be given;\cr see the examples below.
  
  In practice, node sets needs to be specified only for use\cr
  in constructing data sets with at least one two-mode network.

  If there are any \code{NA} values in \code{val},\cr
  and \code{imputationValues} is given, and the covariate is monadic,\cr
  then the corresponding elements
  of \code{imputationValues} are used for imputation.\cr
  If \code{imputationValues} is \code{NULL}, imputation is by the mean value.  \cr
  In both cases, cases with imputed values\cr
  are not used for calculating target statistics (see the manual).

  For changing covariates (monadic or dyadic), the values for wave \code{m}\cr
  are supposed in the simulations to be valid in the whole period 
  from  wave \code{m} to wave \code{m+1}.\cr
  If the data set has \code{M} waves, this means that the values, if any, 
  for wave \code{M} will not be used.\cr
  These values may be given or omitted as part of \code{val}.

  There are \code{print} methods for classes \code{\link{coDyadCovar}} and
  \code{\link{varDyadCovar}}.

  This function replaces \code{\link{coCovar}}, \code{\link{varCovar}},
  \code{\link{coDyadCovar}}, \code{\link{varDyadCovar}},
  which now are obsolete.
}
\value{
  The covariate is returned as an object of class \code{\link{coCovar}},
  \code{\link{varCovar}},
  \code{\link{coDyadCovar}}, or  
  \code{\link{varDyadCovar}},\cr as the case may be.  This object
  can be used as an argument to \code{\link{make_data_rsiena}}.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Tom A.B. Snijders}
\seealso{\code{\link{make_data_rsiena}}, \code{\link{as_nodeset_rsiena}};
  this function replaces the functions \code{\link{coCovar}},
  \code{\link{varCovar}}, \code{\link{coDyadCovar}}, \code{\link{varDyadCovar}}.
}
\examples{
# constant monadic covariate:
myconstcovar <- as_covariate_rsiena(s50a[,1])
senders <- as_nodeset_rsiena(50, nodeSetName="senders")
receivers <- as_nodeset_rsiena(30, nodeSetName="receivers")
senders.attribute <- as_covariate_rsiena(rep(1:10, each=5),
                                         nodeSet="senders")
receivers.attribute <- as_covariate_rsiena(rep(1:5, each=6),
                                         nodeSet="receivers")
# varying monadic covariate:
myvarcovar <- as_covariate_rsiena(s50a, type="monadic")
senders.covariate <- as_covariate_rsiena(s50a, nodeSet="senders",
                                               type="monadic")
receivers.covariate <- as_covariate_rsiena(s50s[1:30,],
                            nodeSet="receivers", type="monadic")
# constant dyadic covariate:
mydyadvar <- as_covariate_rsiena(s503, type="oneMode")
# varying dyadic covariate:
mydyadvar <- as_covariate_rsiena(array(c(s501, s502), dim=c(50, 50, 2)),
                                type="oneMode")
mydyadvar
# using sparse matrices for a dyadic covariate:
require(Matrix)
sps501 <- as(s501,"TsparseMatrix")
sps502 <- as(s502,"TsparseMatrix")
sps501r <- as(s501[,1:30],"TsparseMatrix")
sps502r <- as(s502[,1:30],"TsparseMatrix")
spcDyadvar <- as_covariate_rsiena(sps501, type="oneMode")
spvDyadvar <- as_covariate_rsiena(list(sps501,sps502), type="oneMode")
spvDyadBipvar <- as_covariate_rsiena(list(sps501r,sps502r),
                   type="bipartite", nodeSet=c("senders","receivers"))
spvDyadvar
spvDyadBipvar
}
\keyword{classes}
