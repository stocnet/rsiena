\name{siena08}
\alias{siena08}
\alias{sienaMeta}
\title{
Function to perform a meta analysis of a collection of Siena fits.
}
\description{
Estimates a meta analysis based on a collection of Siena fits.
}
\usage{
siena08(..., projname = "sienaMeta", bound = 5, alpha = 0.05, maxit=20)
}
\arguments{
  \item{\dots}{names of \code{\link{sienaFit}} objects, returned from
    \code{\link{siena07}}.
    They will be renamed if entered in format \code{newname=oldname}.
    It is also allowed to give for \dots a list of \code{\link{sienaFit}}
    objects.
}
\item{projname}{
  Base name of report file if required
}
\item{bound}{
  Upper limit of standard error for inclusion in the meta analysis.
}
\item{alpha}{
  1 minus confidence level of confidence intervals.
}
\item{maxit}{
  Number of iterations of iterated least squares procedure.
}
}
\details{
  A meta analysis is performed as described
  in the Siena manual, section "Meta-analysis of Siena results". This
  consists of three parts: an iterated weighted least squares (IWLS)
  modification of the method described in the reference below;
  maximum likelihood estimates and confidence intervals based on
  profile likelihoods under normality assumptions;
  and Fisher combinations of left-sided and right-sided \eqn{p}-values.
  These are produced for all effects separately.

  Note that the corresponding effects must have the same effect name in each
  model fit. This implies that at least covariates and behavior variables must
  have the same name in each model fit.
}
\value{An object of class \code{\link{sienaMeta}}.  There are \code{print},
  \code{summary} and \code{plot} methods for this class.
  This object contains at least the following.
  \item{thetadf}{Data frame containing the coefficients, standard errors
    and score test results}
  \item{projname}{Root name for any output file to be produced by the print method}
  \item{bound}{Estimates with standard error above this value were
    excluded from the calculations}
  \item{scores}{Object of class \code{by} indicating, for each effect in
    the models, whether score test  information was present.}
  \item{requestedEffects}{The \code{requestedEffects} component of the first
    \code{sienaFit} object in ....}
  \item{muhat}{The vector of IWLS estimates.}
  \item{se.muhat}{The vector of standard errors of the IWLS estimates.}
  \item{theta}{The vector of ML estimates \code{mu.ml} (see below).}
  \item{se}{The vector of standard errors of the ML estimates \code{mu.ml.se}
    (see below).}

  Then for each effect, there is a list with at least the following.

  \item{cor.est}{Spearman rank correlation coefficient between estimates
                 and their standard errors.}
  \item{cor.pval}{p-value for above}
  \item{regfit}{Part of the result of the fit of \code{\link{iwlsm}}.}
  \item{regsummary}{The summary of the fit, which includes the
                    coefficient table.}
  \item{Tsq}{test statistic for effect zero in every model}
  \item{pTsq}{p-value for above}
  \item{tratio}{test statistics that mean effect is 0}
  \item{ptratio}{p-value for above}
  \item{Qstat}{Test statistic for variance of effects is zero}
  \item{pttilde}{p-value for above}
  \item{cjplus}{Test statistic for at least one theta strictly greater than 0}
  \item{cjminus}{Test statistic for at least one theta strictly less than 0}
  \item{cjplusp}{p-value for \code{cjplus}}
  \item{cjminusp}{p-value for \code{cjminus}}
  \item{mu.ml}{ML estimate of population mean}
  \item{mu.ml.se}{standard error of ML estimate of population mean}
  \item{sigma.ml}{ML estimate of population standard deviation}
  \item{mu.confint}{confidence interval for population mean
                    based on profile likelihood}
  \item{sigma.confint}{confidence interval for population standard deviation
                    based on profile likelihood}
  \item{n1}{Number of fits on which the meta analysis is based}
  \item{cjplus}{Test statistic for combination of right one-sided
  Fisher combination tests}
  \item{cjminus}{Test statistic for combination of left one-sided
  Fisher combination tests}
  \item{cjplusp}{p-value for \code{cjplus}}
  \item{cjminusp}{p-value for \code{cjminus}}
  \item{scoreplus}{Test statistic for combination of right one-sided
  \eqn{p}-values from score tests}
  \item{scoreminus}{Test statistic for combination of left one-sided
  \eqn{p}-values from score tests}
  \item{scoreplusp}{p-value for \code{scoreplus}}
  \item{scoreminusp}{p-value for \code{scoreminus}}
  \item{ns}{Number of fits on which the score test analysis is based}
}
\references{
  Snijders, T.A.B, and Baerveldt, C. (2003),
  A Multilevel Network Study of the Effects of Delinquent Behavior
  on Friendship Evolution.
  \emph{Journal of Mathematical Sociology} \bold{27}, 123--151.

  See also the manual (Section 11.2)
  and \url{https://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Ruth Ripley, Tom Snijders}

\seealso{\code{\link{print.sienaMeta}}, \code{\link{funnelPlot}},
  \code{\link{meta.table}}, \code{\link{iwlsm}}, \code{\link{siena07}}
}
\examples{
\dontrun{
# A meta-analysis for three groups does not make much sense
# for generalizing to a population of networks,
# but the Fisher combinations of p-values are meaningful.
# However, using three groups does show the idea.

Group1 <- sienaDependent(array(c(N3401, HN3401), dim=c(45, 45, 2)))
Group3 <- sienaDependent(array(c(N3403, HN3403), dim=c(37, 37, 2)))
Group4 <- sienaDependent(array(c(N3404, HN3404), dim=c(33, 33, 2)))
dataset.1 <- sienaDataCreate(Friends = Group1)
dataset.3 <- sienaDataCreate(Friends = Group3)
dataset.4 <- sienaDataCreate(Friends = Group4)
OneAlgorithm <- sienaAlgorithmCreate(projname = "SingleGroups", seed=128)
effects.1 <- getEffects(dataset.1)
effects.3 <- getEffects(dataset.3)
effects.4 <- getEffects(dataset.4)
effects.1 <- includeEffects(effects.1, transTrip)
effects.1 <- setEffect(effects.1, transRecTrip, fix=TRUE, test=TRUE)
effects.3 <- includeEffects(effects.3, transTrip)
effects.3 <- setEffect(effects.3, transRecTrip, fix=TRUE, test=TRUE)
effects.4 <- includeEffects(effects.4, transTrip)
effects.4 <- setEffect(effects.4, transRecTrip, fix=TRUE, test=TRUE)
ans.1 <- siena07(OneAlgorithm, data=dataset.1, effects=effects.1, batch=TRUE)
ans.3 <- siena07(OneAlgorithm, data=dataset.3, effects=effects.3, batch=TRUE)
ans.4 <- siena07(OneAlgorithm, data=dataset.4, effects=effects.4, batch=TRUE)
ans.1
ans.3
ans.4
(meta <- siena08(ans.1, ans.3, ans.4))
plot(meta, which=2:3, layout = c(2,1))
# For specifically presenting the Fisher combinations:
# First determine the components of meta with estimated effects:
which.est <- sapply(meta, function(x){ifelse(is.list(x),!is.null(x$cjplus),FALSE)})
Fishers <- t(sapply(1:sum(which.est),
        function(i){c(meta[[i]]$cjplus, meta[[i]]$cjminus,
                        meta[[i]]$cjplusp, meta[[i]]$cjminusp, 2*meta[[i]]$n1 )}))
Fishers <- as.data.frame(Fishers, row.names=names(meta)[which.est])
names(Fishers) <- c('Fplus', 'Fminus', 'pplus', 'pminus', 'df')
Fishers
round(Fishers,4)
}
}
\keyword{models}
