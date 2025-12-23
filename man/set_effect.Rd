\name{set_effect} 
\alias{set_effect}
\alias{set_effect.sienaEffects}
\alias{setEffect}
\title{
  Function to change the specification of a Siena model
  by setting various columns in an effects object.
}
\description{
  This function provides an interface to change the 
  specification of a Siena model by setting 
  various columns in an effects object.\cr
  The function \code{set_effect} combines the functionality
  of \code{setEffect} and \code{\link{includeEffects}}.
  \code{set_effect} now is preferred.
  The oher two are kept for backward compatibility.\cr
}
\usage{
\method{set_effect}{sienaEffects}(x, shortNames, type="eval",
            depvar=x$name[1], covar1="", covar2="",
            effect1=0, effect2=0, effect3=0, period=1, group=1, 
            parameter=NULL, initialValue = 0, random=FALSE,
            timeDummy = ",", 
            include=TRUE, fix=FALSE, test=FALSE, verbose=TRUE, ...)

setEffect(x, shortName, parameter = NULL,
    fix = FALSE, test = FALSE, random=FALSE, initialValue = 0,
    timeDummy = ",", include = TRUE,
    name = x$name[1], type = "eval", interaction1 = "",
    interaction2 = "", effect1=0, effect2=0, effect3=0,
    period=1, group=1, character=FALSE, verbose = TRUE)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{x}{a Siena effects object as created by \code{\link{make_specification}}
}
\item{shortNames}{A shortName or list of shortNames
  (all with or all without quotes)
  to identify the effect which should be changed.
}
\item{shortName}{A shortName to identify the effect 
    which should be changed.
}
  \item{depvar}{Name of dependent variable (network or behavior)
    for which interactions are being defined.
    Defaults to the first in the effects object.
  }
  \item{covar1}{
    Vector of Siena objects where needed to completely identify the
    effect e.g. covariate name or behavior variable name.
    See the examples below. Trailing blanks may be omitted.
  }
  \item{covar2}{
    Vector of Siena objects where needed to completely identify the
    effect e.g. covariate name or behavior variable name.
    Trailing blanks may be omitted.
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
    for which the effect is being modified.
    Defaults to the first in the effects object, which is the
   first dependent variable specified in \code{\link{make_data_rsiena}}.
}
  \item{type}{Character string indicating the
  type of the effect to be changed : "rate", "eval",
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
  the row number in \code{x} of the first component of
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
\item{...}{Additional arguments, not used now.}
}
\details{Recall from the help page for \code{\link{make_specification}}
   that a Siena effects object (class \cr
   \code{sienaEffects} or \code{sienaGroupEffects}) 
   is a \code{\link{data.frame}};
   the rows in the data frame are the effects for this data set;
   some of the columns/variables of the data frame
   are used to identify the effect, other columns/variables define how this
   effect is used in the estimation. \cr
   The function \code{set_effect} can operate on several
   effects simultaneously, but in a less detailed way.
   Using \code{set_effect} for only one effect can change not only the
   value of the column \code{include}, but also
   those of \code{initialValue} and \code{parm}.
   The arguments \code{shortName}, \code{name}, \code{type},
   \code{interaction1}, \code{interaction2}, \code{effect1}, \code{effect2},
   \code{effect3}, \code{period}, and \code{group} should identify one
   effect completely. \cr
   (Not all of them are needed; see \code{\link{make_specification}}.)\cr
   The call of \code{setEffect} will set, for this effect,
   the column elements of the resulting effects object for \code{parm},
   \code{fix}, \code{test}, \code{randomEffects}, \code{initialValue},
    \code{timeDummy}, and \code{include} to the values requested.\cr
    For \code{setEffect} the \code{shortName} must not 
    be set between quotes, unless \code{character=TRUE}.

   The input names \code{interaction1} and \code{interaction2}
   for \code{setEffect} do not themselves
   refer to created interactions, but to dependence of the base effects on
   other variables in the data set.
   They are used to completely identify the effects.

   A value for \code{parameter} should be given only if only one effect
   is specified in \code{shortNames}. If it is given,
   the occurrences of \code{#} in the original effect and function names
   are replaced by this value. If a value for \code{parameter} is not given,
   the current value of the internal effect parameter of this effect is used.
}
\value{
  An object of class \code{\link{sienaEffects}} or
  \code{\link{sienaGroupEffects}}. This will be an updated version of the
  input effects object, with one row updated.
  Details of the rows altered will be printed, unless \code{verbose=FALSE}.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley, Tom Snijders}
\seealso{
  \code{\link{make_specification}}, 
  \code{\link{set_interaction}}, \code{\link{includeGMoMStatistics}}, \cr
  \code{\link{update_specification}},
  \code{\link{print.sienaEffects}}, \code{\link{effectsDocumentation}}.
}

\examples{
mynet <- as_dependent_rsiena(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- as_dependent_rsiena(s50a, type="behavior")
mydata <- make_data_rsiena(mynet, mybeh)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, gwespFF)
myeff <- set_effect(myeff, list(inPop, outAct))
myeff <- set_effect(myeff, "gwespFB")
myeff
# Specify an effect parameter:
myeff <- set_effect(myeff, outTrunc, parameter=1)
myeff
# Set the initial rate parameter for one period:
myeff <- set_effect(myeff, Rate, initialValue=1.5, depvar="mybeh",
                   type="rate", period=2)
myeff
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{classes, methods}
