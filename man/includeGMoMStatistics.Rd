\name{includeGMoMStatistics}
\alias{includeGMoMStatistics}
\title{
  Function to include GMoM statistics in a Siena model
}
\description{
  This function can be used for including one or more GMoM statistics   
  by modifying a Siena effects object.
}
\usage{
  includeGMoMStatistics(myeff, ..., include=TRUE, 
                        depvar=myeff$name[1], covar1="", covar2="",
                        character=FALSE, verbose=TRUE)
}
\arguments{
  \item{myeff}{a Siena effects object as created by 
      \code{\link{make_specification}}.
  }
  \item{\dots}{
    short names to identify the GMoM statistics which should be included or
    excluded.}
  \item{include}{Boolean; default TRUE, but can be switched to FALSE to
    turn off an effect.}
  \item{depvar}{Name of dependent variable (network or behavior) for which
    statistics are being included. Defaults to the first in the effects object.}
  \item{covar1}{
    Name of siena object where needed to completely identify the
    effects e.g. covariate name or behavior variable name.}
  \item{covar2}{
    Name of siena object where needed to completely identify the
    effects e.g. covariate name or behavior variable name.}
  \item{character}{Boolean: are the statistic names character strings or not.}
  \item{verbose}{Boolean: should the print of altered statistic be produced.}
}
\details{The names \code{covar1} and \code{covar2} refer to
   the dependence of the GMoM statistics on other variables in the data set.
  The arguments should identify the GMoM statistic completely.
  The type does not have to be specified, as it is \code{gmm} for all
  GMoM statistiscs in the effects object. \cr
  The short names must not be set between quotes,
  unless you use \code{character=TRUE}.

  The function \code{includeGMoMStatistics} operates as an interface
  setting the "include" column on selected rows of the effects object,
  to the value requested (\code{TRUE} or \code{FALSE}).
}
\value{
  An updated version of the input effects object, with the \code{include}
  column for one or more rows updated.
  Details of the rows altered will be printed.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Viviana Amati.
}
\seealso{
  \code{\link{make_specification}},
  \code{\link{set_effect}},
  \code{\link{set_interaction}},
  \code{\link{print.sienaEffects}}
}

\examples{
  mynet1 <- as_dependent_rsiena(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  mybeh  <- as_dependent_rsiena(s50a, type="behavior")
  mydata <- make_data_rsiena(mynet1, mybeh)
  myeff <- make_specification(mydata)
  myeff <- includeGMoMStatistics(myeff, egoX_gmm, covar1="mybeh")
  myeff
}
\keyword{classes}
