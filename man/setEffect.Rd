\name{setEffect}
\alias{setEffect}
\title{
  Function to set various columns in an effects object in a Siena model.
}
\description{
  This function provides an interface to change various columns of
  a selected row of a Siena effects object.
}
\usage{
setEffect(myeff, shortName, parameter = NULL,
fix = FALSE, test = FALSE, random=FALSE, initialValue = 0, timeDummy = ",", include = TRUE,
name = myeff$name[1], type = "eval", interaction1 = "",
interaction2 = "", effect1=0, effect2=0, effect3=0,
period=1, group=1, character=FALSE, verbose = TRUE)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{myeff}{a Siena effects object as created by \code{\link{getEffects}}
}
\item{shortName}{A short name (all with or all without quotes)
  to identify the effect which should be changed.
}
  \item{parameter}{Value of internal effect parameter. If \code{NULL},
  no change is made.
}
  \item{fix}{For fixing effects. Boolean required. Default \code{FALSE}.
}
  \item{test}{For testing effects by score-type tests. Boolean required.
  Default \code{FALSE}.
}
  \item{random}{For specifying that effects will vary randomly;
  used only in function \code{sienaBayes} in package \code{multiSiena}.
  Not relevant for \code{RSiena} at this moment.
  Boolean required. Default \code{FALSE}.
}
  \item{initialValue}{Initial value for estimation. Default 0.
}
\item{timeDummy}{string: Comma delimited string of which periods to
dummy. Alternatively, use \code{\link{includeTimeDummy}}.
}
  \item{include}{Boolean. default \code{TRUE}, but can be switched to
    \code{FALSE} to turn off an effect.}
  \item{name}{Name of dependent variable (network or behavior)
    for which effects are being modified.
    Defaults to the first in the effects object.
}
  \item{type}{Character string indicating the
  type of the effect to be changed : currently "rate", "eval",
  "endow", or "creation". Default "eval".
}
  \item{interaction1}{
  Name of siena object where needed to completely identify the
  effect e.g. covariate name or behavior variable name.
}
  \item{interaction2}{
  Name of siena object where needed to completely identify the
  effect e.g. covariate name or behavior variable name.
}
\item{effect1}{
  Only for \code{shortName=unspInt, behUnspInt} or \code{contUnspInt},
  which means this is a user-defined interaction effect: \code{effect1} is
  the row number in \code{myeff} of the first component of
  the interaction effect.
}
\item{effect2}{
  See \code{effect1}: second component of interaction effect.
}
\item{effect3}{
  See \code{effect1}: third component of interaction effect.
}
\item{period}{
  Number of period if basic rate. Use numbering within groups.
}
\item{group}{
  Number of group if basic rate.
  Only relevant for \code{\link{sienaGroup}} data sets.
  }
\item{character}{Boolean: whether the short name is a character string.}
\item{verbose}{Boolean: should the print of altered effects be produced.}
}
\details{Recall from the help page for \code{\link{getEffects}}
   that a Siena effects object (class \code{sienaEffects} or
   \code{sienaGroupEffects}) is a \code{\link{data.frame}};
   the rows in the data frame are the effects for this data set;
   some of the columns/variables of the data frame
   are used to identify the effect, other columns/variables define how this
   effect is used in the estimation. \cr
   The function \code{\link{includeEffects}} can operate on several
   effects simultaneously, but in a less detailed way.
   The main use of \code{setEffect} is that it can change not only the
   value of the column \code{include}, but also
   those of \code{initialValue} and \code{parm}.
   The arguments \code{shortName}, \code{name}, \code{type},
   \code{interaction1}, \code{interaction2}, \code{effect1}, \code{effect2},
   \code{effect3}, \code{period}, and \code{group} should identify one
   effect completely. (Not all of them are needed;
   see \code{\link{getEffects}}.)\cr
   The call of \code{setEffect} will set, for this effect,
   the column elements of the resulting effects object for \code{parm},
   \code{fix}, \code{test}, \code{randomEffects}, \code{initialValue},
    \code{timeDummy}, and \code{include} to the values requested.\cr
    The \code{shortName} must not be set between quotes,
    unless you use \code{character=TRUE}.

   The input names \code{interaction1} and \code{interaction2} do not themselves
   refer to created interactions, but to dependence of the base effects on
   other variables in the data set.
   They are used to completely identify the effects.

   If a value for \code{parameter} is given,
   the occurrences of \code{#} in the original effect and function names
   are replaced by this value. If a value for \code{parameter} is not given,
   the default value of the internal effect parameter of this effect is used.
}
\value{
  An object of class \code{\link{sienaEffects}} or
  \code{\link{sienaGroupEffects}}. This will be an updated version of the
  input effects object, with one row updated.
  Details of the row altered will be printed, unless \code{verbose=FALSE}.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley, Tom Snijders}
\seealso{
  \code{\link{getEffects}}, \code{\link{includeEffects}},
  \code{\link{includeInteraction}}, \code{\link{includeGMoMStatistics}},
  \code{\link{updateSpecification}},
  \code{\link{print.sienaEffects}}, \code{\link{effectsDocumentation}}.
}

\examples{
mynet <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- sienaDependent(s50a, type="behavior")
mydata <- sienaDataCreate(mynet, mybeh)
myeff <- getEffects(mydata)
myeff <- setEffect(myeff, gwespFF)
# Specify an effect parameter:
myeff <- setEffect(myeff, outTrunc, parameter=1)
myeff
# Set the initial rate parameter for one period:
myeff <- setEffect(myeff, Rate, initialValue=1.5, name="mybeh",
                   type="rate", period=2)
myeff
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{ classes }
