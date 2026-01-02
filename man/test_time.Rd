\name{test_time}
\alias{test_time}
\alias{test_time.sienaFit}
\alias{sienaTimeTest}
\alias{print.sienaTimeTest}
\title{Function to assess and account for time heterogeneity of parameters}
\description{
 Takes a \code{sienaFit} object estimated by Method of Moments,
 and tests for time heterogeneity by the addition of interactions
 with time dummy variables at waves \code{m=2...(M-1)}.
 The test used is the score-type test of Schweinberger (2012).

 Tests for joint significance, parameter-wise significance,
 period-wise significance, individual significance,
 and one-step estimates of the unrestricted model parameters
 are returned in a list.\cr
  \code{test_time()} and \code{sienaTimeTest()}
  are identical functions. The first name is preferred,
  the other is kept for backward compatibility.
 }
\usage{
\method{test_time}{sienaFit}(x, effects=NULL,
              excludedEffects=NULL, condition=FALSE, ...)

sienaTimeTest(x, effects=NULL,
              excludedEffects=NULL, condition=FALSE)
}
\arguments{
  \item{x}{ A \code{sienaFit} object
      returned by \code{\link{siena}}. }
  \item{effects}{Optional vector of effect numbers to test. Use the
    numbering on the print of the sienaFit object.}
  \item{excludedEffects}{Optional vector of effect numbers
    for which time heterogeneity is not to be tested. Use the
    numbering on the print of the sienaFit object.}
  \item{condition}{ Whether to orthogonalize effect-wise score-type tests
    and individual significance tests against estimated effects
    and un-estimated dummy terms, or just against estimated effects.}
  \item{\dots}{For other arguments. Not used currently.  }
}
\details{
  This test follows the score type test of Schweinberger (2012) as
 elaborated by Lospinoso et al. (2011) by using statistics already
 calculated at each wave to obtain vectors of partitioned moment
 functions corresponding to a restricted model (the model in the
 \code{sienaFit} object; used as null hypothesis) and an
 unrestricted model (which contains dummies for waves
 \code{m=2...(M-1)}; used as alternative hypothesis).

 \code{condition=TRUE} leads to a rough-and-easy approximation to
 controlling the mentioned tests also for the unestimated effects.

 After assessing time heterogeneity, effects objects can be modified by
 adding numbers of all or some periods to the \code{timeDummy} column.
 This is facilitated by the \code{\link{includeTimeDummy}} function.
 For an effects object in which the \code{timeDummy} column
 of some of the included effects includes some or all period numbers,
 interactions of those effects  with time dummies for the indicated
 periods will also be estimated.

 An alternative to the use of \code{\link{includeTimeDummy}}
 is to define time-dependent actor covariates (dummy variables or other
 functions of wave number that are the same for all actors),
 include these in the data set through \code{\link{make_data_rsiena}},
 and include interactions of other effects with ego effects of these
 time-dependent actor covariates by \code{\link{set_interaction}}.
 This is illustrated in an example below.
 Using \code{\link{includeTimeDummy}} is easier;
 using self-defined interactions with time-dependent variables
 gives more control.

 If you wish to use this function with \code{sienaFit} objects that use the
 finite differences method of derivative estimation, or which use maximum
 likelihood estimation, you must request the derivatives to be returned
 by wave using \code{byWave=TRUE} option in the call of \code{\link{siena07}}.

 Effects leading to dummy interactions that are collinear with the
 model originally fitted, after excluding the effects mentioned,
 will be automatically excluded from the time heterogeneity testing.

 If \code{test_time} gives errors that there are too
 many collinear effects, run it with a smaller set of effects as specified
 by the \code{effects} parameter.
 For example, if the model has 40 effects of which the
 first 8 are rate parameters and therefore uninteresting,
 and there is such an error message, try \code{effects=9:30};
 if that still does not work, decrease the upper limit of 30, if it does work
 increase it, to find the largest possible set of effects for which
 heterogeneity assessment still is possible; then as a next step
 try the remaining effects in a similar way.

 Also if the execution is time-consuming, e.g., for a multi-group
 \code{sienaFit} object with many groups and many effects,
 it can be helpful to carry out the function in smaller subsets of effects.
}
\value{
  \code{test_time} returns a list containing many items,
    including the following:
  \item{JointTest }{A chi-squared test for joint significance of the dummies.}
  \item{EffectTest }{A chi-squared test for joint significance
      across dummies for each separate effect.}
  \item{GroupTest }{A chi-squared test for joint significance
      across dummies; if \code{sienaFit} is a fit for a multi-group object
      then these refer to each group; else they refer to each period.}
  \item{IndividualTest }{A matrix displaying initial estimates, one-step
    estimates, and \eqn{p}-values for the individual interactions.}
}

\references{
J.A. Lospinoso, M. Schweinberger, T.A.B. Snijders, and R.M. Ripley
(2011), Assessing and Accounting for Time Heterogeneity in Stochastic
Actor Oriented Models. \emph{Advances in Data Analysis and Computation},
\bold{5}, 147--176.

M. Schweinberger (2012), Statistical modeling of network panel data:
Goodness-of-fit. \emph{British Journal of Statistical
and Mathematical Psychology} \bold{65}, 263--281.
}

\author{Josh Lospinoso, Tom Snijders}
\seealso{\code{\link{siena}}, \code{\link{plot.sienaTimeTest}},
  \code{\link{includeTimeDummy}}}
\examples{
## Estimate a restricted model
myalgorithm <- set_algorithm_saom(nsub=1, n3=50)
# Short estimation not for practice, just for having a quick demonstration
mynet1 <- as_dependent_rsiena(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mydata <- make_data_rsiena(mynet1)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, transTrip)
ans <- siena(mydata, effects=myeff, batch=TRUE, control_algo=myalgorithm)

## Conduct the score-type test to assess whether heterogeneity is present.
tt <- test_time(ans)
summary(tt)

## Suppose that we wish to include time dummies.
## Add them in the following way:
myeff <- includeTimeDummy(myeff, recip, transTrip, timeDummy="2")
ans2 <- siena(mydata, effects=myeff, batch=TRUE, control_algo=myalgorithm)

## Re-assess the time heterogeneity
(tt2 <- test_time(ans2))

## And so on..

## A demonstration of the plotting facilities, on a larger dataset:
## (Of course pasting these identical sets of three waves after each other
##  in a sequence of six is not really meaningful. It's just a demonstration.)

myalgo <- set_algorithm_saom(nsub=1, n3=50, seed=654)
mynet1 <- as_dependent_rsiena(array(c(s501, s502, s503, s501, s503, s502),
                               dim=c(50, 50, 6)))
mydata <- make_data_rsiena(mynet1)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(transTrip, outAct))
myeff <- includeTimeDummy(myeff, density, timeDummy="all")
myeff <- includeTimeDummy(myeff, recip, timeDummy="2,3,5")
myeff <- includeTimeDummy(myeff, outAct, timeDummy="4")
\dontrun{
(ansp <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE))
ttp <- test_time(ansp, effects=2:3)
# Using effects=1:4 leads to an informative error message.

## Pairwise plots show
plot(ttp, pairwise=TRUE)
}

## Instead of working with includeTimeDummy,
## you can also define time dummies explicitly;
## this may give more control and more clarity:
dum2 <- matrix(c(0,1,0,0,0), nrow=50, ncol=5, byrow=TRUE)
dum3 <- matrix(c(0,0,1,0,0), nrow=50, ncol=5, byrow=TRUE)
dum4 <- matrix(c(0,0,0,1,0), nrow=50, ncol=5, byrow=TRUE)
dum5 <- matrix(c(0,0,0,0,1), nrow=50, ncol=5, byrow=TRUE)
time2 <- as_covariate_rsiena(dum2, type="monadic")
time3 <- as_covariate_rsiena(dum3, type="monadic")
time4 <- as_covariate_rsiena(dum4, type="monadic")
time5 <- as_covariate_rsiena(dum5, type="monadic")
mydata <- make_data_rsiena(mynet1, time2, time3, time4, time5)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(transTrip, outAct))
## corresponding to includeTimeDummy(myeff, density, timeDummy="all"):
myeff <- set_effect(myeff, egoX, covar1='time2')
myeff <- set_effect(myeff, egoX, covar1='time3')
myeff <- set_effect(myeff, egoX, covar1='time4')
myeff <- set_effect(myeff, egoX, covar1='time5')
## corresponding to myeff <- includeTimeDummy(myeff, recip, timeDummy="2,3,5"):
myeff <- set_interaction(myeff, list(egoX, recip), covar1=c('time2', ''))
myeff <- set_interaction(myeff, list(egoX, recip), covar1=c('time3', ''))
myeff <- set_interaction(myeff, list(egoX, recip), covar1=c('time5', ''))
## corresponding to myeff <- includeTimeDummy(myeff, outAct, timeDummy="4"):
myeff <- set_interaction(myeff, list(egoX, outAct), covar1=c('time4', ''))
\dontrun{
(anspp <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE))
## anspp contains identical results as ansp above.
}

## A demonstration of RateX heterogeneity.
\dontrun{
mynet1 <- as_dependent_rsiena(array(c(s501, s502, s503), dim=c(50, 50, 3)))
myccov <- as_covariate_rsiena(s50a[,1], type="monadic")
mydata <- make_data_rsiena(mynet1, myccov)
myeff <- make_specification(mydata)
myeff <-  set_effect(myeff, list(transTrip, outAct))
myeff <- includeTimeDummy(myeff, RateX, type="rate", covar1="myccov")
(ans <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE))
}
}
\keyword{tests,methods}
