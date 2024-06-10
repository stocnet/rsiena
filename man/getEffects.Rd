\name{getEffects}
\alias{getEffects}
\alias{sienaEffects}
\alias{sienaGroupEffects}
\title{Function to create a Siena effects object}
\description{
  Creates a basic list of effects for all dependent
  variables in the input siena object.
}
\usage{
getEffects(x, nintn = 10, behNintn=4, getDocumentation=FALSE, onePeriodSde=FALSE)
}
\arguments{
  \item{x}{an object of class `siena" or `sienaGroup"}
  \item{nintn}{Number of user-defined network interactions
       that can later be created.}
  \item{behNintn}{Number of user-defined behavior interactions
       that can later be created.}
  \item{getDocumentation}{Flag to allow documentation of internal
  functions, not for use by users.}
  \item{onePeriodSde}{Flag to indicate that the stochastic
       differential equation (SDE) model dZ(t) = [aZ(t) + b] dt + g dW(t)
       should be used, instead of the regular SDE with a scale parameter.
       This is only relevant in case the model includes a continuous
       dependent variable and one period is studied.}
}
\details{
Creates a data frame of effects for use in siena model estimation.
The regular way of changing this object is by the functions
\code{\link{includeEffects}}, \code{\link{setEffect}},
and \code{\link{includeInteraction}}.

Note that the class of the return object may be lost if the data.frame
is edited using \command{fix}.
See \code{\link{fix}} and \code{\link{edit.data.frame}}.
}
\value{An object of class \code{sienaEffects} or \code{sienaGroupEffects}:
this is a data frame of which the rows are the effects available
for data set \code{x}.\cr
  The effects object consists of consecutive parts, each of which relates
  to one dependent variable in the input object.
  The columns are:
  \item{name}{name of the dependent variable}
  \item{effectName}{name of the effect}
  \item{functionName}{name of the function}
  \item{shortName}{short name for the effect}
  \item{interaction1}{second variable to define the effect, if any}
  \item{interaction2}{third variable to define the effect, if any}
  \item{type}{"eval", "endow", "creation", "rate", or "gmm"}
  \item{basicRate}{boolean: whether a basic rate parameter}
  \item{include}{boolean: include in the model to be fitted or not}
  \item{randomEffects}{boolean: random or fixed effect. Currently not used.}
  \item{fix}{boolean: fix parameter value or not}
  \item{test}{boolean: test parameter value or not}
  \item{timeDummy}{comma separated list of periods, or "all",
    or "," for none --
    which time dummy interacted parameters should be included?}
  \item{initialValue}{starting value for estimation, also used for
    \code{fix} and \code{test}.}
  \item{parm}{internal effect parameter values}
  \item{functionType}{"objective" or "rate"}
  \item{period}{period for basic rate parameters}
  \item{rateType}{"Structural", "covariate", "diffusion"}
  \item{untrimmedValue}{Used to store initial values which could be
    trimmed}
  \item{effect1}{Used to indicate effect number in user-specified interactions}
  \item{effect2}{Used to indicate effect number in user-specified interactions}
  \item{effect3}{Used to indicate effect number in user-specified interactions}
  \item{interactionType}{Defines "dyadic" or "ego" or "OK" effects, used in
    \code{\link{includeInteraction}}}
  \item{local}{whether a local effect; used for the option \code{localML}
    in \code{\link{sienaAlgorithmCreate}}}
  \item{effectFn}{here NULL, but could be replaced by a function later}
  \item{statisticFn}{here NULL, but could be replaced by a function
    later}
  \item{netType}{Type of dependent variable:
    "oneMode", "behavior", or "bipartite"}
  \item{groupName}{name of relevant group data object}
  \item{group}{sequential number of relevant group data object in total}
  \item{effectNumber}{the sequence number of the row}
  The combination of \code{name}, \code{shortName}, \code{interaction1},
  \code{interaction2}, and \code{type} uniquely identifies any effect
  other than basic rate effects and user-specified interaction effects.
  For the latter, \code{effect1}, \code{effect2} and \code{effect3}
  are also required for the identification. The combination \code{name},
  \code{shortName}, \code{period} and \code{group} uniquely identifies a
  basic rate effect.

  The columns not used for identifying the effect define how the effect is used
  for the estimation.

  The columns \code{initialValue} and \code{parm} should not be confused:
  \code{initialValue} gives the initial value for the parameter to be estimated,
  indicated in the manual by \code{theta};
  \code{parm} gives the internal value of the parameter defining the effect,
  indicated in the manual (Chapter 12) by \code{p}, and is fixed
  during the estimation.

  Note that if an effects object is printed by \code{print(...)},
  by default only the included rows are printed.

  A list of all effects in a given effects object (e.g., \code{myeff}),
  including their names of dependent variables, effect names, short names,
  and values of interaction1 and interaction2 (if any),
  is obtained by executing \code{\link{effectsDocumentation}(myeff)}.

  As from version 1.3.24, effects object have a "version" attribute.
  Effects objects including interaction effects created by
  \code{\link{includeInteraction}} are not necessarily compatible
  between versions of \code{RSiena}. Therefore it is recommended,
  for effects objects including any interaction effects,
  to create them again when changing to a new version of \code{RSiena}.
  If an effects object including any interaction effects is used from an
  old version of \code{RSiena}, this will lead to a warning when running
  \code{\link{siena07}}. \cr
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley, Tom Snijders}
\seealso{\code{\link{sienaDataCreate}}, \code{\link{sienaGroupCreate}},
        \code{\link{includeEffects}}, \code{\link{setEffect}},
        \code{\link{includeGMoMStatistics}},
        \code{\link{updateSpecification}},
        \code{\link{print.sienaEffects}},
		\code{\link{effectsDocumentation}}}
\examples{
mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- sienaDependent(s50a, type="behavior")
mycovar <- coCovar(rnorm(50))
mydyadcovar <- coDyadCovar(matrix(as.numeric(rnorm(2500) > 2), nrow=50))
mydata <- sienaDataCreate(mynet1, mybeh, mycovar, mydyadcovar)
myeff <- getEffects(mydata)
myeff
}

% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{classes}
