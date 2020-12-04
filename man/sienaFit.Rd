\name{sienaFit.methods}
\alias{sienaFit.methods}
\alias{print.summary.sienaFit}
\alias{sienaFit}
\alias{xtable.sienaFit}
\alias{print.xtable.sienaFit}
\alias{siena.table}
\alias{print.sienaFit}
\alias{summary.sienaFit}
\title{Methods for processing sienaFit objects, produced by \code{\link{siena07}}.
}
\usage{
\method{print}{sienaFit}(x, tstat=TRUE, ...)

\method{summary}{sienaFit}(object, ...)

\method{print}{summary.sienaFit}(x, matrices=TRUE, ...)

\method{xtable}{sienaFit}(x, caption = NULL, label = NULL, align = NULL,
                digits = NULL, display = NULL, ...)

siena.table(x, type="tex", file=paste(deparse(substitute(x)), ".", type,sep=""),
            vertLine=TRUE, tstatPrint=FALSE, sig=FALSE, d=3, nfirst=NULL)
}
\arguments{
  \item{object}{An object of class \code{\link{sienaFit}}, produced by
    \code{\link{siena07}}.}
  \item{x}{An object of class \code{sienaFit}, or
    \code{summary.sienaFit} as appropriate.}
  \item{matrices}{Boolean: whether also to print in the summary the 
   covariance matrix of the estimates, the derivative matrix of expected 
   statistics X by parameters, and the covariance matrix of the statistics.}
  \item{tstat}{Boolean: if this is NULL, the t-statistics for convergence
     will not be added to the report.}
  \item{type}{Type of output to produce; must be either \code{"tex"} or
    \code{"html"}.}
  \item{file}{Name of the file; defaults to the name of the \code{sienaFit}
    object. "" indicates output to the console.}
  \item{vertLine}{Boolean: add vertical lines separating the columns in
    \code{siena.table}.}
  \item{tstatPrint}{Boolean: add a column of significance t values (parameter
    estimate/standard error estimate) to \code{siena.table}.}
  \item{sig}{Boolean: adds symbols (daggers and asterisks) indicating
    significance levels for the parameter estimates to \code{siena.table}.}
  \item{d}{The number of decimals places used in \code{siena.table}.}
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
  \item{nfirst}{Only relevant for the \code{RSienaTest} package.}
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

  The function \code{siena.table} outputs a latex or html table of
  the estimates and standards errors of a \code{\link{sienaFit}} object.
  The table will be written to a file in the current directory and
  has a footnote reporting the maximum of the convergence t-ratios.

  See the manual for how to import these tables easily into MS-Word.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Ruth Ripley, Charlotte Greenan, Tom Snijders}

\seealso{\code{\link[xtable]{xtable}}, \code{\link[xtable]{print.xtable}},
  \code{\link{siena07}}}
\examples{
myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=100, projname=NULL)
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)
ans
summary(ans)
\dontrun{
xtable(ans, type="html", file="ans.html")
siena.table(ans, type="html", tstat=TRUE, d=2)
}
}
\keyword{methods}
