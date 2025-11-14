\name{sienaRIDynamics}
\alias{sienaRIDynamics}
\alias{plot.sienaRIDynamics}
\title{Functions to assess the dynamics of relative importance of effects
       between observation moments}
\description{
   The function \code{sienaRIDynamics} returns the relative importance of
   effects in a SAOM according to the measure of relative importance described
   in Indlekofer and Brandes (2013).
   More precisely, it returns aggregates of relative importances over
   simulated micro-steps between observation moments as described in Section
   3.2 of this paper.
   The measure is based on the influence of effects on simulated
   micro-steps and takes the complete model specification into account.
   Therefore, required arguments are the analysed data given as a
   \code{siena} data object (necessary for the micro-steps simulations)
   as well as the complete model specification represented either by an
   estimated \code{sienaFit} object or by the triple consisting of a
   suitable parameter vector \code{theta} and the corresponding
   \code{sienaAlgorithm} and \code{sienaEffects} objects.
 }
\usage{
  sienaRIDynamics(data, ans=NULL, theta=NULL, algorithm=NULL, effects=NULL,
     depvar=NULL, intervalsPerPeriod=NULL, n3 = NULL, useChangeContributions=FALSE)
  \method{plot}{sienaRIDynamics}(x, staticRI=NULL, col=NULL,
     ylim=NULL, width=NULL, height=NULL, legend=TRUE, legendColumns=NULL,
     legendHeight=NULL, cex.legend=NULL, \dots)
}
\arguments{
  \item{data}{\code{siena} data object representing the analyzed data and
  resulting from a call to \code{\link{sienaDataCreate}}.}
  \item{ans}{\code{sienaFit} object resulting from a call to
    \code{\link{siena07}}.
    The \code{sienaFit} object contains all necessary information on the model
    specification, in particular, the vector of parameter values
    \code{ans$theta}, the used algorithm for estimation \code{ans$x},
    and information on included model effects \code{ans$effects}.
    If \code{ans} is a valid \code{sienaFit} object, the calculations of
    relative importances are based on \code{ans$theta}, \code{ans$x},
    and \code{ans$effects}.
    Alternatively, the necessary information can be given directly
    as a suitable parameter vector \code{theta}, a \code{sienaAlgorithm} object,
    and a \code{sienaEffects} object.
    In this case, \code{ans} has to be unspecified (i.e., \code{ans=NULL}).}
  \item{theta}{Vector of parameter values of effects included in the model.
    Length of \code{theta} has to be equal to the number of included effects
    (rate effects as well as evaluation effects).}
  \item{algorithm}{\code{sienaAlgorithm} object as resulting from a call to
    \code{\link{sienaAlgorithmCreate}}.
    Works only for estimation by Method of Moments (i.e., if
    \code{maxlike = FALSE}).}
  \item{effects}{\code{sienaEffects} object specifying which effects are
    included the model.}
  \item{depvar}{If the model contains more than one dependent variable,
    it has to be specified for which dependent variable the relative importances
    of associated effects should be determined.
    Only those simulated micro-steps that refer to the selected dependent
    variable are considered in the aggregated values of relative importance
    of the corresponding effects.
    If \code{depvar=NULL}, relative importances will be determined for
    the first element in the list of dependent variables of \code{data}
    (\code{attr(data$depvars,"names")[1]}).}
  \item{intervalsPerPeriod}{For analyzing the dynamics of relative
    importances between observation moments, the time interval between
    observation moments is divided into \code{intervalsPerPeriod} subintervals.
    The values of relative importance are aggregated over all micro-steps
    of one subinterval. Default value is 10 subintervals per period.}
  \item{n3}{Integer specifying the number of chains to be simulated. 
    Defaults to the value in \code{algorithm} or what was set when estimating \code{ans}.
    A very large value will result in a long computation time and large memory usage.}
  \item{useChangeContributions}{Logical specifying whether previously extracted change 
    contributions should be used instead of rerunning phase 3 simulations.}
  \item{x}{\code{sienaRIDynamics} object resulting from a call to
    \code{sienaRIDynamics}.}
  \item{staticRI}{\code{\link{sienaRI}} object corresponding to the same model
    and data as \code{x}.
    If \code{staticRI} is specified, the expected relative importances of
    effects at observation moment are indicated in the plot by circlets
    at the beginning of each period.}
  \item{col}{Colors used in the plot. If \code{col=NULL} a default color
    scheme is used.}
  \item{ylim}{A vector of two numbers specifying the maximum and minimum
    value of the y-range visible in the plot. If \code{ylim=NULL} default
    values are used.}
  \item{width}{Width of the plot. If \code{width=NULL} a default value is used.}
  \item{height}{Height of the plot. If \code{height=NULL} a
    default value is used.}
  \item{legend}{Boolean: if \code{TRUE} a legend is added to the plot.
    \code{x$effectNames} are used as labels.}
  \item{legendColumns}{Number of columns in legend. If
    \code{legendColumns=NULL} a default value is used.
    Only effective if \code{legend=TRUE}.}
  \item{legendHeight}{Height of legend. If \code{legendHeight=NULL}
    a default value is used. Only effective if \code{legend=TRUE}.}
  \item{cex.legend}{Specifies the relative font size of legend labels.}
  \item{\dots}{Other arguments.}
}
\details{
\code{sienaRIDynamics} takes the data as well as the complete model
specification into account.
Therefore, necessary arguments are the analyzed data given as a
\code{siena} data object as well as the complete model specification
represented either by an estimated \code{sienaFit} object or by the
triple consisting of a suitable parameter vector \code{theta}
and the corresponding \code{sienaAlgorithm} and \code{sienaEffects} objects.

Note that \code{sienaRIDynamics} works only with Method of Moments
(i.e., for \code{sienaAlgorithm} objects with \code{maxlike = FALSE}).

If \code{ans} is a valid \code{sienaFit} object that resulted from a
call to \code{\link{siena07}} with unconditional estimation (i.e.,
for \code{sienaAlgorithm} objects with \code{cond = FALSE}),
the calculations of relative importances are based on \code{ans$theta},
\code{ans$x}, and \code{ans$effects}.
(Note that if the estimation of \code{ans} has been conditional on a
dependent variable, it is necessary to give the associated
\code{sienaEffects} object explicitly to \code{sienaRIDynamics}).

Alternatively, the necessary information can be given directly as a
suitable parameter vector \code{theta} (containing necessary rate
and evaluation parameters), a \code{sienaAlgorithm} object, and a
\code{sienaEffects} object. In this case, \code{ans} has to be
unspecified, i.e., \code{ans=NULL}.

If the sienaFit object \code{ans} was already estimated with
\code{returnChangeContributions=TRUE}, the change contributions
are stored in \code{ans$changeContributions} and can be used to
estimate the relative importances of effects in \code{sienaRIDynamics}
without running any additional simulations.
In this case, \code{useChangeContributions} has to be set to \code{TRUE}
and \code{n3} is ignored, since the number of chains is determined by the
previous \code{siena07} run.
If \code{useChangeContributions=FALSE}, new chains are simulated internally
with \code{siena07} and the change contributions are calculated from these
chains. The number of simulations is determined by the value of \code{n3} or 
the value of \code{algorithm$n3} if \code{algorithm} is specified. Otherwise,
it defaults to the value of \code{ans$n3} if \code{ans} is specified.
}
\value{
\code{sienaRIDynamics} returns an object of class \code{sienaRIDynamics}.

 A returned \code{sienaRIDynamics} object stores the relative importances
   of effects in simulated micro-steps aggregated to \code{intervalsPerPeriod}
   averages per period. For details, see Section 3.2 of
   Indlekofer and Brandes (2013).

A \code{sienaRIDynamics} object consists of following components:
    \describe{
    \item{\code{intervalValues}}{the relative importances of included
          effects for simulated micro-steps of the considered dependent
          variable aggregated over subintervals.}
    \item{\code{dependentVariable}}{the name of the dependent variable
          the results refer to.}
    \item{\code{effectNames}}{the names of considered effects.}
    }
}
\seealso{\code{\link{sienaRI} and code{\link{siena07}}}}
\references{
  Indlekofer, N. and Brandes, U. (2013), 
  Relative Importance of Effects in Stochastic Actor-oriented Models.
  \emph{Network Science}, \bold{1 (3)}, 275--297.
}
\author{Natalie Indlekofer}
\examples{
mycontrols <- sienaAlgorithmCreate(nsub=2, n3=50, cond=FALSE)
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, density, recip, transTies, nbrDist2)
ans <- siena07(mycontrols, data=mydata, effects=mymodel, batch=TRUE,
               returnChangeContributions=TRUE)

RIDynamics1 <- sienaRIDynamics(mydata, ans=ans, useChangeContributions=TRUE)
RIDynamics1
plot(RIDynamics1)


mycontrols2 <- sienaAlgorithmCreate(nsub=2, n3=50, cond=TRUE)
ans2 <- siena07(mycontrols2, data=mydata, effects=mymodel, batch=TRUE)

RIDynamics2 <- sienaRIDynamics(mydata, ans=ans2, effects=mymodel)
RIDynamics2

RIDynamics3 <- sienaRIDynamics(data=mydata, theta=c(ans2$rate,ans2$theta),
             algorithm=mycontrols2, effects=mymodel, intervalsPerPeriod=4)
RIDynamics3

mycontrols3 <- sienaAlgorithmCreate(nsub=2, n3=50)
mynet2 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- sienaDependent(s50a, type="behavior")
mydata2 <- sienaDataCreate(mynet2, mybeh)
mymodel2 <- getEffects(mydata2)
mymodel2 <- includeEffects(mymodel2, density, recip, transTies, transTrip, nbrDist2)
ans3 <- siena07(mycontrols3, data=mydata2, effects=mymodel2, batch=TRUE)

RIDynamics4 <- sienaRIDynamics(mydata2, ans=ans3, depvar="mybeh")
RIDynamics4

RIDynamics5 <- sienaRIDynamics(mydata2, ans=ans3, depvar="mynet2", n3 = 100)
RIDynamics5


RI5 <- sienaRI(mydata2, ans=ans3)
plot(RIDynamics5, staticRI=RI5[[1]])

}
\keyword{models}