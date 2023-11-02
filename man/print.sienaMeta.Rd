\name{print.sienaMeta}
\alias{print.sienaMeta}
\alias{summary.sienaMeta}
\alias{plot.sienaMeta}
\alias{print.summary.sienaMeta}
\alias{meta.table}
\title{Methods for processing sienaMeta objects}
\usage{
\method{print}{sienaMeta}(x, file=FALSE, reportEstimates=FALSE, ...)

\method{summary}{sienaMeta}(object, file=FALSE, extra=TRUE, ...)

\method{print}{summary.sienaMeta}(x, file=FALSE, extra=TRUE, ...)

\method{plot}{sienaMeta}(x, ..., which = 1:length(x$theta),
                         useBound=TRUE, layout = c(2,2))
meta.table(x, d=3, option=2,
    filename=paste(deparse(substitute(x)),'_global.tex',sep=""), align=TRUE)
}
\description{
  \code{print}, \code{summary}, and \code{plot} methods for   
   \code{sienaMeta} objects; and a function to write a LaTeX table.
}
\arguments{
  \item{object}{An object of class \code{sienaMeta}.}
  \item{x}{An object of class \code{sienaMeta}, or
    \code{summary.sienaMeta} as appropriate.}
  \item{file}{Boolean: if TRUE, sends output to file named
    \code{x$projname.txt}. If FALSE, output is to the terminal.}
  \item{reportEstimates}{Boolean: whether to report all estimates
    and standard errors.}
  \item{extra}{Boolean: if TRUE, prints more information.}
  \item{which}{Set of effects contained in the plot
    (given by sequence numbers).}
  \item{useBound}{Boolean: whether to restict plotted symbols to the bound
    used in the call of \code{sienaMeta}.}
  \item{layout}{Vector giving number of rows and columns
    in the arrangement of the several panels in a
    rectangular array, possibly spanning multiple pages.}
  \item{d}{Number of decimals to be used in table.}
  \item{option}{1: results without normality assumptions;
    2: results with normality assumptions, with confidence intervals;
    3: results with normality assumptions, with standard errors.}
  \item{filename}{filename for output; if "", printed to the console.}
  \item{align}{Whether to align numbers at the decimal point.}
  \item{\dots}{For extra arguments (none used at present).}
}
\value{
 The function \code{print.sienaMeta} prints details of the merged
 estimates of the meta-analysis carried out by \code{\link{siena08}},
 with test statistics. See the help page for \code{\link{siena08}}
 for what is produced by this function.

 The function \code{summary.sienaMeta} prints details as for the
 \code{print} method, but also details of the \code{sienaFit} objects
 included.

 Output from either can be directed to a file by using the argument
 \code{file}. It will be appended to any existing file of the same
 name: \code{projname.txt}  where \code{projname} is the value of the
 argument to \code{siena08}.

 The function \code{meta.table} writes a combined table of results
 for all parameters to a LaTeX file in (as default) the current working
 directory. This table is a more compact version of
 the results presented by \code{print.sienaMeta}.

 The function \code{plot.sienaMeta} plots estimates against standard
 errors for each effect, with reference lines added at the two-sided
 significance threshold 0.05. It returns an object of class \code{trellis},
 of the \pkg{lattice} package. Effects for which a score test
 was requested are not plotted.\cr
 Another funnel plot, not using \code{\link{siena08}},
 is available as \code{\link{funnelPlot}}.
}
\references{ 
Snijders, T.A.B, and Baerveldt, C. (2003),
A Multilevel Network Study of the Effects of Delinquent Behavior
on Friendship Evolution.
\emph{Journal of Mathematical Sociology} \bold{27}, 123--151.

  See also the Siena manual and
  \url{https://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Ruth Ripley, Tom Snijders
}
\seealso{\code{\link{siena08}}}

\examples{
\dontrun{
# A meta-analysis for three groups does not make much sense
# for generalizing to a population of networks,
# but it the Fisher combinations of p-values are meaningful.
# But using three groups shows the idea.

Group1 <- sienaDependent(array(c(N3401, HN3401), dim=c(45, 45, 2)))
Group3 <- sienaDependent(array(c(N3403, HN3403), dim=c(37, 37, 2)))
Group4 <- sienaDependent(array(c(N3404, HN3404), dim=c(33, 33, 2)))
dataset.1 <- sienaDataCreate(Friends = Group1)
dataset.3 <- sienaDataCreate(Friends = Group3)
dataset.4 <- sienaDataCreate(Friends = Group4)
OneAlgorithm <- sienaAlgorithmCreate(projname = "SingleGroups")
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
meta <- siena08(ans.1, ans.3, ans.4)
print(meta, reportEstimates=FALSE)
print(meta)
summary(meta)
# For specifically presenting the Fisher combinations:
# First determine the number of estimated effects:
(neff <- sum(sapply(meta, function(x){ifelse(is.list(x),
        !is.null(x$cjplus),FALSE)})))
Fishers <- t(sapply(1:neff,
        function(i){c(meta[[i]]$cjplus, meta[[i]]$cjminus,
                        meta[[i]]$cjplusp, meta[[i]]$cjminusp, 2*meta[[i]]$n1 )}))
Fishers <- as.data.frame(Fishers, row.names=names(meta)[1:neff])
names(Fishers) <- c('Fplus', 'Fminus', 'pplus', 'pminus', 'df')
Fishers
# For plotting:
plo <- plot(meta, layout = c(3,1))
plo
plo[3]
# Show effects of bound (bounding at 0.4 is not reasonable, just for example)
meta <- siena08(ans.1, ans.3, ans.4, bound=0.4)
plot(meta, which=c(2,3), layout=c(2,1))
plot(meta, which=c(2,3), layout=c(2,1), useBound=FALSE)
meta.table(meta, option=3, file='')
}}
\keyword{print}
