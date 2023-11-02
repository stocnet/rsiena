\name{includeEffects}
\alias{includeEffects}
\title{
  Function to include effects in a Siena model
}
\description{
  This function can be used for model specification  
  by modifying a Siena effects object.
}
\usage{
includeEffects(myeff, ..., include = TRUE, name = myeff$name[1], type = "eval",
 interaction1 = "", interaction2 = "", fix=FALSE, test=FALSE, character=FALSE,
 verbose = TRUE)
}
\arguments{
  \item{myeff}{a Siena effects object as created by \code{\link{getEffects}}
}
  \item{\dots}{
short names to identify the effects which should be included or excluded.}
  \item{include}{Boolean. default TRUE, but can be switched to FALSE to
    turn off an effect.}
  \item{name}{Name of dependent variable (network or behavior) for which effects
    are being included. Defaults to the first in the effects object.}
  \item{type}{
  Type of effects to be included: "eval", "endow", "creation", or "rate".}
  \item{interaction1}{
  Name of siena object where needed to completely identify the
  effects e.g. covariate name or behavior variable name.}
  \item{interaction2}{
  Name of siena object where needed to completely identify the
  effects e.g. covariate name or behavior variable name.}
  \item{fix}{
  Boolean. Are the effects to be fixed at the value stored in
  \code{myeff$initialValue} or not.}
  \item{test}{
  Boolean. Are the effects to be tested or not (requires \code{fix}).}
\item{character}{Boolean: are the effect names character strings or not.}
  \item{verbose}{Boolean: should the print of altered effects be produced.}
}
\details{Recall from the help page for \code{\link{getEffects}}
   that a Siena effects object (class \code{sienaEffects} or
   \code{sienaGroupEffects}) is a \code{\link{data.frame}};
   the rows in the data frame are the effects for this data set;
   some of the columns/variables of the data frame
   are used to identify the effect, other columns/variables define how this
   effect is used in the estimation. \cr   
   
  The function \code{includeEffects} operates as an interface
  setting the "include" column on selected rows of the effects object,
  to the value requested (\code{TRUE} or \code{FALSE}).
   The selected effects must be indicated by the arguments \ldots, \code{type},
  and (if necessary) \code{interaction1} and \code{interaction2}.
  The names \code{interaction1} and \code{interaction2} do not refer
  to interactions between effects, but to dependence of effects on
  other variables in the data set.
  The arguments should identify the effects completely.\cr
  The short names must not be set between quotes,
  unless you use \code{character=TRUE}.

  Note that the internal effect parameter has a default value which differs
  between effects. This can be set by function \code{\link{setEffect}}.
  Also the value of \code{myeff$initialValue} can be set by this function.
  The function \code{\link{setEffect}} operates on the effects object
  in a more detailed way, but applies to one effect at the time.

  Further information about Siena effects objects is given
  in the help page for \code{\link{getEffects}}.

  A list of all effects available in a given effects object (e.g.,
  \code{myeff}), including their names of dependent variables,
  effect names, short names,
  and values of \code{interaction1} and \code{interaction2} (if any),
  is obtained by executing \code{\link{effectsDocumentation}(myeff)}.

  The input names \code{interaction1} and \code{interaction2} do not themselves
  refer to created interactions, but to dependence of the base effects on
  other variables in the data set.
  They are used to completely identify the effects.
}
\value{
    An updated version of the input effects object, with the
    \code{include}, \code{test}, and \code{fix}
    columns for one or more rows updated.
    Details of the rows altered will be printed.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{
  \code{\link{getEffects}}, \code{\link{setEffect}},
  \code{\link{includeInteraction}},
  \code{\link{includeGMoMStatistics}},
  \code{\link{updateSpecification}},
  \code{\link{print.sienaEffects}}, \code{\link{effectsDocumentation}}
}

\examples{
mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh  <- sienaDependent(s50a, type="behavior")
mydata <- sienaDataCreate(mynet1, mybeh)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip, balance)
myeff <- includeEffects(myeff, avAlt, name="mybeh", interaction1="mynet1")
myeff
}
\keyword{ classes }
