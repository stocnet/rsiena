\name{includeInteraction}
\alias{includeInteraction}
\title{
Function to create user-specified interactions for a Siena model.
}
\description{This function allows the user to include or exclude
  an interaction effect in a Siena effects object.
}
\usage{
includeInteraction(myeff, ..., include = TRUE, name = myeff$name[1],
    type = "eval", interaction1 = rep("", 3), interaction2 = rep("", 3),
    fix = FALSE, test = FALSE, parameter = NULL, random = FALSE,
    character = FALSE, verbose = TRUE)
}
\arguments{
  \item{myeff}{a Siena effects object as created by \code{\link{getEffects}}}
  \item{\dots}{2 or 3 short names to identify the effects
    which should be interacted.}
  \item{include}{Boolean. default TRUE, but can be switched to FALSE to
    turn off an interaction.
  }
  \item{name}{Name of dependent variable (network or behavior)
    for which interactions are being defined.
    Defaults to the first in the effects object.
  }
  \item{type}{
    Type of effects to be interacted.
  }
  \item{interaction1}{
    Vector of Siena objects where needed to completely identify the
    effect e.g. covariate name or behavior variable name.
    Trailing blanks may be omitted.
  }
  \item{interaction2}{
    Vector of siena objects where needed to completely identify the
    effect e.g. covariate name or behavior variable name.
    Trailing blanks may be omitted.
  }
  \item{fix}{
    Boolean. Are the effects to be fixed at the value stored in
    \code{myeff$initialValue} or not.}
  \item{test}{
    Boolean. Are the effects to be tested or not (requires \code{fix}).}
  \item{parameter}{
    Value of internal effect parameter of this interaction effect.
	If NULL, no change is made.}
  \item{random}{For specifying that the interaction effect will vary randomly;
    not relevant for \code{RSiena} at this moment.
    Boolean required. Default FALSE.}
  \item{character}{Boolean: are the effect names character strings or not.}
  \item{verbose}{Boolean: should the print of altered effects be produced.}
}
\details{The details provided should uniquely identify up to three
  effects. If so, an interaction effect will be created and included or
  not in the model. \cr
  Whether interactions between two or three given effects
  can be created depends on their \code{interactionType} (which can be, for
  dependent network variables, empty, ego, or dyadic; and for dependent
  behavioral variables, empty or OK). Consult the section on Interaction
  Effects in the manual for this. The \code{interactionType} is shown
  in the list of effects obtained from the function
  \code{\link{effectsDocumentation}}.\cr
  The short names must not be set between quotes,
  unless you use \code{character=TRUE}.\cr
  From the point of view of model building it is usually advisable, when
  including an interaction effect in a model, also to include the
  corresponding main effects. This is however not enforced
  by \code{includeInteraction()}.

  Interaction effects are constructed from effects with shortName
  \code{unspInt} (for networks) and \code{behUnspInt} (for behavior)
  by specifying their elements \code{effect1} and \code{effect2},
  and possibly \code{effect3}.

  The number of possible user-specified interaction effects is limited
  by the parameters \code{nintn} (for dependent network variables)
  and \code{behNintn} (for dependent behavior variables)
  in the call of \code{\link{getEffects}}, which determine the numbers of
  effects with shortNames \code{unspInt} and \code{behUnspInt}.

  The input names \code{interaction1} and \code{interaction2} do not themselves
  refer to created interactions, but to dependence of the base effects on
  other variables in the data set.
  They are used to completely identify the effects.

  To set attributes of interaction effects in the effects object,
  function \code{\link{setEffect}} can also be used with short name
  \code{unspInt} or \code{behUnspInt}, and further using parameters
  \code{effect1} and \code{effect2}, and possibly \code{effect3}.

  Further information about Siena effects objects is given
  in the help page for \code{\link{getEffects}}.

  A list of all effects in a given effects object (e.g., \code{myeff}),
  including their names of dependent variables, effect names, short names,
  and values of interaction1 and interaction2 (if any),
  is obtained by executing \code{\link{effectsDocumentation}(myeff)}.
  }
\value{
  An updated version of the input effects object; if \code{include},
  containing the interaction effect between "effect1" and "effect2" and
  possibly "effect3"; if not, without this interaction effect.
  The shortName of the interaction effect is "unspInt" for network effects
  and "behUnspInt" for behavior effects. \cr
  If \code{verbose=TRUE}, details of the fields altered will be printed.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley, Tom Snijders}
\seealso{
  \code{\link{getEffects}}, \code{\link{setEffect}},
  \code{\link{includeEffects}}, \code{\link{effectsDocumentation}
}}
\examples{
mynet <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
alc <- varCovar(s50a)
mydata <- sienaDataCreate(mynet, alc)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
myeff <- includeInteraction(myeff, recip, outAct)
myeff <- includeEffects(myeff, egoX, altX, simX, interaction1="alc")
myeff <- includeInteraction(myeff, recip, simX, interaction1=c("", "alc"))
myeff
# How to set the effect parameter of an interaction:
myeff <- getEffects(mydata)
myeff <- setEffect(myeff, gwespFF, parameter=69)
myeff <- includeInteraction(myeff, recip, gwespFF, parameter=69)
myeff
}
\keyword{ classes }
