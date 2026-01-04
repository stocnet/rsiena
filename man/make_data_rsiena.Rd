\name{make_data_rsiena}
\alias{sienaDataCreate}
\alias{make_data_rsiena}
\alias{sienadata}
\title{Function to create a Siena data object}
\description{
    Creates a Siena data object from input dependent variables   
    (networks and possibly behavioural variables), covariates,
     and composition change objects.\cr
  The functions \code{make_data_rsiena} and \code{sienaDataCreate}
  are identical. The first name is preferred, 
  the other is kept for backward compatibility.
}
\usage{
make_data_rsiena(..., nodeSets=NULL)

sienaDataCreate(..., nodeSets=NULL, getDocumentation=FALSE)
}
\arguments{
  \item{\dots}{objects of class \code{\link{sienaDependent}},
    \code{\link{coCovar}}, \code{\link{varCovar}}, \code{\link{coDyadCovar}},
  \code{\link{varDyadCovar}}, and/or \code{\link{sienaCompositionChange}};
  or a list of such objects, of which the first element must not be
  a \code{\link{sienaCompositionChange}} object.
  There should be at least one \code{\link{sienaDependent}} object.\cr
  If there are one-mode as well as two-mode dependent networks,
  the one-mode networks should be mentioned first.}
  \item{nodeSets}{list of Siena node sets. Default is the single node
    set named \code{Actors}, length equal to the number of rows in the first
    object of class \code{sienaDependent}. If the entire data set contains more than
    one node set, then the node sets must have been specified in
    the creation of all data objects mentioned in \dots.}
  \item{getDocumentation}{Flag to allow documentation of internal
    functions, not for use by users.}
}
\details{
The function checks that the objects fit, that there is at least one dependent
variable, and adds various attributes to each variable describing the
data. If there is more than one nodeSet they must all be specified.\cr
Function \code{\link{print01Report}} will give a basic description of the data
object and is a check useful, e.g., for diagnosing problems.
}
\value{
An object of class \code{sienadata} which is designed to be used in a siena
model fit by \code{\link{siena}}. \cr
In \code{RSiena} versions prior to 1.6, the class name was \code{siena}. \cr
The components of the object are:
\item{nodeSets}{List of node sets involved}
\item{observations}{Integer indicating number of waves of data}
\item{depvars}{List of networks and behavior variables}
\item{cCovars}{List of constant covariates}
\item{vCovars}{List of changing covariates}
\item{dycCovars}{List of constant dyadic covariates}
\item{dyvCovars}{List of changing dyadic covariates}
\item{compositionChange}{List of composition change objects
  corresponding to the node sets}
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{
\code{\link{as_dependent_rsiena}}, \code{\link{as_covariate_rsiena}},
\code{\link{as_nodeset_rsiena}}, \cr
\code{\link{as_composition_rsiena}}, \code{\link{make_group_rsiena}},
\code{\link{make_constraint.sienadata}}, \cr
\code{\link{print01Report}}}
\examples{
mynet <- as_dependent_rsiena(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- as_dependent_rsiena(s50a, type="behavior")
mydata <- make_data_rsiena(mynet, mybeh)
# This gives the same result as
mydata <- make_data_rsiena(list(mynet, mybeh))
## And for a two-mode network
mynet1 <- as_dependent_rsiena(array(c(s501, s502), dim=c(50, 50, 2)), 
                              nodeSet="senders")
senders <- as_nodeset_rsiena(50, nodeSetName="senders")
receivers <- as_nodeset_rsiena(30, nodeSetName="receivers")
mynet2 <- as_dependent_rsiena(
            array(c(s501[,1:30], s502[,1:30]), dim=c(50, 30, 2)),
            nodeSet=c("senders", "receivers"))
(mydata <- make_data_rsiena(mynet1, mynet2, nodeSets=list(senders, receivers)))
\dontrun{
print01Report(mydata, modelname = "mydescription")
}
}
\keyword{classes}
