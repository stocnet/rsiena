\name{sienaGroupCreate}
\alias{sienaGroupCreate}
\alias{sienaGroup}
\title{Function to group together several Siena data objects }
\description{Creates an object of class "sienaGroup" from a list of  
  Siena data objects.
}
\usage{
sienaGroupCreate(objlist, singleOK = FALSE, getDocumentation=FALSE)
}
\arguments{
  \item{objlist}{ List of objects of class \code{siena}. }
  \item{singleOK}{ Boolean: is it OK to only have one object?}
  \item{getDocumentation}{Flag to allow documentation of internal
  functions, not for use by users.}
}
\details{

This function creates a Siena group object from several Siena data
objects ('groups'), all of which use networks, covariates and actor sets
with the same names. The variables must correspond exactly
between all data objects; the numbers of waves may differ.
It can be used as data input to \code{\link{siena07}} for the
multigroup option.
Also used internally for convenience with a single Siena data object.

Each covariate should either be centered in all groups, or non-centered
in all groups.
Centered actor covariates are re-centered at the overall mean.
This means that the original values are used, and the overall
mean of all non-missing observations is subtracted. Note that this implies
that group-dependent variables that are constant for all actors in each
group, can be used as centered actor covariates.

For combining two-wave with more-wave groups in one group object,
covariates that are changing covariates for the more-wave groups have to
be specified as changing covariates also for the two-wave groups.
This can be done by specifying them with values for the two waves;
for actor covariates this will be by using an \eqn{n \times 2}{n * 2} matrix,
for dyadic covariates an \eqn{n \times n \times 2}{n * n * 2} array
(or \eqn{n \times m \times 2}{n * m * 2} for the two-mode case).
The values for the second wave should be identical to those for the
first wave (they will be used only for centering operations).

For later use in \code{\link{siena07}},
it will often (but not always...) be helpful when creating the
Siena data objects in \code{objlist}
to use \code{allowOnly=FALSE} in the call of
\code{\link{sienaDependent}}; see the help page for this function.

If there are multiple dependent networks, it may be necessary to run
\code{\link{sienaDataConstraint}} before \code{sienaGroupCreate}
to ensure that these constraints are equal for all groups.
}

\value{
An object of class \code{sienaGroup}; this is a list
containing the input objects, with attributes:
  \item{netnames}{names of the dependent variables in each set}
  \item{symmetric}{vector of booleans, one for each dependent
    variable. TRUE if all occurrences of the network are symmetric.}
  \item{structural}{vector of booleans, indicating whether structurally
    fixed values occur in this network}
  \item{allUpOnly}{vector of booleans, indicating whether changes are
    all upwards in all the occurrences of this network}
  \item{allDownOnly}{similar to previous, but for downward changes}
  \item{anyUpOnly}{vector of booleans, indicating whether changes are
    all upwards in any of the occurrences of this network}
  \item{anyDownOnly}{similar to previous, but for downward changes}
  \item{types}{vector of network types of the dependent variables}
  \item{observations}{Total number of periods to process}
  \item{periodNos}{Sequence of numbers of periods which are not skipped
    in multigroup processing}
  \item{netnodeSets}{list of names of the node sets corresponding to the
    dependent variables}
  \item{cCovars}{names of the constant covariates, if any}
  \item{vCovars}{names of the changing  covariates, if any}
  \item{dycCovars}{names of the constant dyadic covariates, if any}
  \item{dyvCovars}{names of the changing dyadic covariates, if any}
  \item{ccnodeSets}{list of the names of the node sets corresponding to
    the constant covariates}
  \item{cvnodeSets}{list of the names of the node sets corresponding to
    the changing covariates}
  \item{dycnodeSets}{list of the names of the node sets corresponding to
    the constant dyadic covariates}
  \item{dyvcnodeSets}{list of the names of the node sets corresponding to
    the changing dyadic covariates}
  \item{compositionChange}{boolean: any composition change at all?}
  \item{exooptions}{named vector of composition change options for the
    node sets}
  \item{names}{Either from the input objects or "Data1", "Data2" etc}
  \item{class}{"sienaGroup" inheriting from "siena"}
  \item{balmean}{vector of means for balance calculations}
  \item{bRange}{vector of difference between maximum and minimum values
    for behavior variables, NA for other dependent variables}
  \item{behRange}{matrix of  maximum and minimum values
    for behavior variables, NA for other dependent variables}
  \item{bSim}{vector of similarity means for behavior variables, NA for
    other dependent variables}
  \item{bPoszvar}{vector of booleans indicating positive variance for
    behavior variables. NA for other dependent variables}
  \item{bMoreThan2}{vector of booleans indicating whether the behavior
    variables take more than 2 distinct values}
  \item{cCovarPoszvar}{vector of booleans indicating positive variance for
    constant covariates}
  \item{cCovarMoreThan2}{vector of booleans indicating whether the
    constant covariates take more than 2 distinct values}
  \item{cCovarRange}{vector of difference between maximum and minimum
    values for constant covariates}
  \item{cCovarRange2}{matrix of maximum and minimum values for constant
    covariates}
  \item{cCovarSim}{vector of similarity means for constant covariates}
  \item{cCovarMean}{vector of means for constant covariates}
  \item{vCovarRange}{vector of difference between maximum and minimum
    values for changing covariates}
  \item{vCovarSim}{vector of similarity means for changing covariates}
  \item{vCovarMoreThan2}{vector of booleans indicating whether the
    changing covariates take more than 2 distinct values}
  \item{vCovarPoszvar}{vector of booleans indicating positive variance for
    changing covariates}
  \item{vCovarMean}{vector of means for changing covariates}
  \item{dycCovarMean}{vector of means for constant dyadic covariates}
  \item{dycCovarRange}{vector of ranges for constant dyadic covariates}
  \item{dycCovarRange2}{matrix of maximum and minimum values
    for constant dyadic covariates}
  \item{dyvCovarRange}{vector of ranges for changing dyadic covariates}
  \item{dyvCovarMean}{vector of means for changing dyadic covariates}
  \item{anyMissing}{vector of booleans, one for each dependent variable,
    indicating the presence of any missing values}
  \item{netRanges}{matrix of maximum and minimum values for dependent
    networks, NA for behavior variables}
}
\references{See the Section on Multi-group Siena analysis in the manual
   available from \url{https://www.stats.ox.ac.uk/~snijders/siena/ }.}
\author{ Ruth Ripley, Modification by Tom Snijders }
\seealso{\code{\link{sienaDataCreate}}, \code{\link{sienaDataConstraint}} }
\examples{
Group1 <- sienaDependent(array(c(N3401, HN3401), dim=c(45, 45, 2)))
Group3 <- sienaDependent(array(c(N3403, HN3403), dim=c(37, 37, 2)))
Group4 <- sienaDependent(array(c(N3404, HN3404), dim=c(33, 33, 2)))
Group6 <- sienaDependent(array(c(N3406, HN3406), dim=c(36, 36, 2)))
# Illustration of the use of group-level variables:
# dum1 is a dummy variable for group 1,
# having constant value 1 in group 1, and constant value 0 in the other groups.
dum1.1  <- coCovar(c(rep(1,45)), warn = FALSE)
dum1.3  <- coCovar(c(rep(0,37)), warn = FALSE)
dum1.4  <- coCovar(c(rep(0,33)), warn = FALSE)
dum1.6  <- coCovar(c(rep(0,36)), warn = FALSE)
# In a similar way, dummies for the other groups can be defined.
dataset.1 <- sienaDataCreate(Friends = Group1, dum1 = dum1.1)
dataset.3 <- sienaDataCreate(Friends = Group3, dum1 = dum1.3)
dataset.4 <- sienaDataCreate(Friends = Group4, dum1 = dum1.4)
dataset.6 <- sienaDataCreate(Friends = Group6, dum1 = dum1.6)
(FourGroups <- sienaGroupCreate(list(dataset.1, dataset.3, dataset.4,
                                     dataset.6)))
class(FourGroups)
# The main effect of the group-level variable is the \code{egoX} effect:
myeff <- getEffects(FourGroups)
(myeff <- includeEffects(myeff, egoX, interaction1 = "dum1"))
}
\keyword{ classes }
