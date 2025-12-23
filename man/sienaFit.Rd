\name{sienaFit.methods}
\alias{sienaFit.methods}
\alias{print.summary.sienaFit}
\alias{sienaFit}
\alias{xtable.sienaFit}
\alias{print.xtable.sienaFit}
\alias{siena_table}
\alias{print.sienaFit}
\alias{summary.sienaFit}
\alias{write_result.sienaFit}
\alias{write_result}
\alias{coef}
\alias{coef.sienaFit}
\alias{vcov}
\alias{vcov.sienaFit}
\title{Methods for processing sienaFit objects, 
         produced by \code{\link{siena}}
}
\usage{
\method{coef}{sienaFit}(object, dropRates=TRUE, shortenNames=TRUE, ...)

\method{vcov}{sienaFit}(x, dropRates=TRUE, shortenNames=TRUE, ...)

\method{print}{sienaFit}(x, tstat=TRUE, ...)

\method{summary}{sienaFit}(object, ...)

\method{print}{summary.sienaFit}(x, matrices=TRUE, ...)

\method{write_result}{sienaFit}(x, type="tex",
    fileName=paste(deparse(substitute(x)),".",type,sep=""),
    vertLine=TRUE, tstatPrint=FALSE, sig=FALSE, d=3, ...)

\method{xtable}{sienaFit}(x, caption = NULL, label = NULL, align = NULL,
    digits = NULL, display = NULL, ...)

siena_table(x, type="tex", 
    fileName=paste(deparse(substitute(x)), ".", type,sep=""),
    vertLine=TRUE, tstatPrint=FALSE, sig=FALSE, d=3, nfirst=NULL)
}
\arguments{
  \item{object}{An object of class \code{\link{sienaFit}}, produced by
    \code{\link{siena}}. For \code{siena_table}, objects of class
    \code{sienaBayes} are also permitted.}
  \item{x}{An object of class \code{sienaFit}, or
    \code{summary.sienaFit} as appropriate. For \code{siena_table},
    objects of class \code{sienaBayes} are also permitted.}
  \item{dropRates}{Logical. Should the results for the basic rate parameters
    also be given?}
  \item{shortenNames}{Logical. Should the effect names
    be shortened?}
  \item{matrices}{Boolean: whether also to print in the summary the
   covariance matrix of the estimates, the derivative matrix of expected
   statistics X by parameters, and the covariance matrix of the statistics.}
  \item{tstat}{Boolean: if this is NULL, the t-statistics for convergence
     will not be added to the report.}
  \item{type}{Type of output to produce; must be either \code{"tex"} or
    \code{"html"}.}
  \item{fileName}{Name of the file; defaults to the name of the \code{sienaFit}
    object. "" indicates output to the console.}
  \item{vertLine}{Boolean: add vertical lines separating the columns in
    \code{siena_table}.}
  \item{tstatPrint}{Boolean: add a column of significance t values (parameter
    estimate/standard error estimate) to \code{siena_table}.}
  \item{sig}{Boolean: adds symbols (daggers and asterisks) indicating
    significance levels for the parameter estimates to \code{siena_table}.}
  \item{d}{The number of decimals places used in \code{siena_table}.}
  \item{caption}{
    See documentation for \code{\link[xtable]{xtable}}.}
  \item{label}{
    See documentation for \code{\link[xtable]{xtable}}.}
  \item{align}{
    See documentation for \code{\link[xtable]{xtable}}.}
  \item{digits}{
    See documentation for \code{\link[xtable]{xtable}}.}
  \item{display}{
    See documentation for \code{\link[xtable]{xtable}}}
  \item{nfirst}{Only relevant for the \code{multiSiena} package.}
  \item{\dots}{
    Add extra parameters for  \code{\link[xtable]{print.xtable}}
    here. e.g.
    \code{type}, \code{file}.}
}
\description{
  \code{print}, \code{summary}, and \code{xtable} methods for
  \code{sienaFit} objects.
}
\value{
  Function \code{coef.sienaFit} gives the estimated parameters,
  and \code{vcov.sienaFit} their estimated variance-covariance matrix.
  
  The function \code{print.sienaFit} prints a table containing estimated
  parameter values, standard errors and (optionally)
  t-statistics for convergence.

  The function \code{summary.sienaFit} prints a table containing
  estimated parameter values, standard errors and t-statistics for
  convergence together with the covariance matrix of the estimates, the
  derivative matrix of expected statistics \code{X} by parameters, and the
  covariance matrix of the expected statistics \code{X}.

  The function \code{xtable.sienaFit} creates an object of class
  \code{xtable.sienaFit} which inherits from class \code{xtable} and
  passes an extra arguments to the \code{print.xtable}.

  The function \code{siena_table} outputs a latex or html table of
  the estimates and standards errors of a \code{\link{sienaFit}} object.
  The table will be written to a file in the current directory and
  has a footnote reporting the maximum of the convergence t-ratios.
  Endowment or creation effects will be denoted, respectively,
  by 'maintenance' or 'creation'.

  See the manual for how to import the html tables easily into MS-Word.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Ruth Ripley, Charlotte Greenan, Tom Snijders}

\seealso{\code{\link[xtable]{xtable}},
 \code{\link[xtable]{print.xtable}}, \code{\link{siena}}}
\examples{
myalgorithm <- set_algorithm_saom(nsub=2, n3=100)
mynet <-  as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- make_data_rsiena(mynet)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, outAct)
ans <- siena(mydata, effects=myeff,
                     control_algo=myalgorithm, batch=TRUE)
ans
coef(ans)
round(vcov(ans), 4)
round(sqrt(diag(vcov(ans))), 4)
summary(ans)
\dontrun{
xtable(ans, type="html", file="ans.html")
write_result(ans, type="html", tstat=TRUE, d=2)
}
}
\keyword{methods}
