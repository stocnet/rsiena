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
    (may be sparse, of type \cr
    \code{"TsparseMatrix"}), or 3-dimensional array, 
    or list of 3 sparse matrices 
    of type \code{"TsparseMatrix"}; which of these structures
    determines the kind of covariate (see Details).}
  \item{type}{\code{"monadic"}, \code{"oneMode"} or \code{"bipartite"}: 
    whether the 
    values refer to a monadic covariate, and in case of a dyadic
    covariate, whether it is for a one-mode or a bipartite 
    (two-mode) network (see Details).}
  \item{centered}{Boolean: if TRUE, then the mean value is subtracted.}
  \item{nodeSet}{Name of node set: character string;
    for \code{type="bipartite"}, referring to
    bipartite (two-mode) dyadic covariates, a vector of two node set names
    (character strings).\cr  
    If the entire data set contains more than one node set,
    then the node sets must be specified in all data objects.}
  \item{warn}{Logical: determines whether a warning given if all values are 
   \code{NA}, or all non-missing values are the same.}
  \item{imputationValues}{For monadic covariates only:
    Vector or matrix of covariate values of same
    length or dimensions as \code{val}, to be used for imputation of
    \code{NA} values (if any) in \code{val}. Must not contain
    any \code{NA}.}
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
  \code{nodeSet} of this data object.
  In practice, the node set needs to be specified only in the case of
  the use of the covariate with a two-mode network.
  
  If there are any \code{NA} values in \code{val}, and
  \code{imputationValues} is given, and the covariate is monadic,
  then the corresponding elements
  of \code{imputationValues} are used for imputation.
  If \code{imputationValues} is \code{NULL}, imputation
  is by the mean value.  In both cases, cases with imputed values
  are not used for calculating target statistics (see the manual).
  
  For changing covariates (monadic or dyadic), the values for wave \code{m} 
  are supposed in the simulations to be valid in the whole period from 
  wave \code{m} to wave \code{m+1}.
  If the data set has \code{M} waves, this means that the values,
  if any, for wave \code{M} will not be used.  
  These values may be given or omitted as part of \code{val}.
  
  The function replaces \code{\link{coCovar}}, \code{\link{varCovar}},
  \code{\link{coDyadCovar}}, \code{\link{varDyadCovar}},
  which now are obsolete.
}
\value{
  The covariate is returned as an object of class \code{\link{coCovar}}, 
  \code{\link{varCovar}},
  \code{\link{coDyadCovar}}, or  \cr
  \code{\link{varDyadCovar}}, as the case may be.  This object
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
# using sparse matrices for a constant dyadic covariate:
require(Matrix)
sps501 <- as(s501,"TsparseMatrix")
sps502 <- as(s502,"TsparseMatrix")
sps502r <- as(s502[1:30,],"TsparseMatrix")
spcDyadvar <- as_covariate_rsiena(sps501, type="oneMode")
spvDyadvar <- as_covariate_rsiena(list(sps501,sps502), type="oneMode")
spvDyadBipvar <- as_covariate_rsiena(list(sps501,sps502), 
                   type="bipartite", nodeSet=c("senders","receivers"))
}
\keyword{classes}
