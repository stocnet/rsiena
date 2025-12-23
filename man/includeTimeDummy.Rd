\name{includeTimeDummy}
\alias{includeTimeDummy}
\title{Function to include time dummy effects in a Siena model}
\description{
  This function specifies time heterogeneity for selected effects   
  in a Siena model, by interacting them with time dummies,   
  without explicitly using time-dependent covariates.
 }

\usage{
includeTimeDummy(myeff, ..., timeDummy="all", depvar=myeff$name[1], 
       type="eval", covar1="", covar2="", include=TRUE, character=FALSE)
}
\arguments{
  \item{myeff}{A Siena effects object as created by 
              \code{\link{make_specification}}.}
  \item{\dots}{
   Short names to identify the effects for which interactions with time dummies
   should be included or excluded.
   This function cannot be used for regular interaction effects.}
  \item{timeDummy}{Character string. Either "all" or the periods for which
   to create dummies (from 1 to (number of waves - 1)), space delimited.}
  \item{include}{Boolean. default TRUE, but can be switched to FALSE to
    turn off an effect.}
  \item{depvar}{Name of dependent network or behavioral variable for which
    effects are being included. Defaults to the first in the effects object.}
  \item{type}{Type of dummy effects to be interacted.}
  \item{covar1}{
   Name of variable where needed to completely identify the
   effects e.g. covariate name or behavior variable name.}
  \item{covar2}{
   Name of variable where needed to completely identify the
   effects e.g. covariate name or behavior variable name.}
  \item{character}{Boolean: are the effect names character strings or not}
}
\details{
  The arguments (\ldots, \code{name, covar1, covar2})
  should identify the effects completely.\cr
  See \code{\link{make_specification}} and \code{\link{effectsDocumentation}}
  for more information about this.

  This function operates by setting the \code{timeDummy} column
  on selected rows of a Siena effects object, thereby specifying
  interactions of the specified effect or effects with dummy variables
  for the specified periods.
  The \code{timeDummy} column of \code{myeff} will be set to include the
  values requested if \code{include=TRUE}, and to exclude them for
  \code{include=FALSE}.

  For an effects object in which the \code{timeDummy} column
  of some of the included effects includes some or all period numbers,
  interactions of those effects with ego effects of time dummies for the
  indicated periods will also be estimated by \code{\link{siena}}.
  For the outdegree effect this is just the ego effect of the time dummies.
  If \dots does not include the outdegree effect, then still this
  ego effect will be created, but its parameter will be fixed to 0.

  An alternative to the use of \code{\link{includeTimeDummy}}
  is to define time-dependent actor covariates (dummy variables or other
  functions of wave number that are the same for all actors),
  include these in the data set through \code{\link{make_data_rsiena}},
  and include interactions of other effects with ego effects of these
  time-dependent actor covariates by \code{\link{set_interaction}}.
  This is illustrated in an example in the help file for
  \code{\link{test_time}}.
  Using \code{\link{includeTimeDummy}} is easier; on the other hand,
  using self-defined interactions with time-dependent variables
  gives more control (e.g., it will allow to specify linear time dependence
  and test time heterogeneity for interaction effects).
}
\value{
    An updated version of \code{myeff}, with the \code{timeDummy}
    column for one or more rows updated.
    Details of the rows altered will be printed.
}

\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}
  for general information on RSiena.

}
\author{Josh Lospinoso}
\seealso{\code{\link{sienaTimeTest}}, \code{\link{make_specification}},
  \code{\link{includeEffects}}, \code{\link{effectsDocumentation}}.}
\examples{
\dontrun{
## Estimate a restricted model
mynet <- as_dependent_rsiena(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mydata <- make_data_rsiena(mynet)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(outAct, transTrip))
(ans <- siena(mydata, effects=myeff))

## Conduct the score type test to assess whether heterogeneity is present.
tt <- sienaTimeTest(ans)
summary(tt)

## Suppose that we wish to include a time dummy.
## Since there are three waves, the number of periods is two.
## This means that only one time dummy can be included for
## the interactions. The default is for period 2;
## an equivalent model, but with different parameters
## (that can be transformed into each other) is obtained
## when the dummies are defined for period 1.
myeff <- includeTimeDummy(myeff, density, recip, timeDummy="2")
myeff       # Note the \code{timeDummy} column.
(ans2 <- siena(mydata, effects=myeff))

## Re-assess the time heterogeneity
tt2 <- sienaTimeTest(ans2)
summary(tt2)

## And so on..
}
}
\keyword{models}
