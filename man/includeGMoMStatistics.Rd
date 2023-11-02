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
  includeGMoMStatistics(myeff, ..., include=TRUE, name=myeff$name[1],
                        interaction1="", interaction2="",
                        character=FALSE, verbose=TRUE)
}
\arguments{
  \item{myeff}{a Siena effects object as created by \code{\link{getEffects}}
  }
  \item{\dots}{
    short names to identify the GMoM statistics which should be included or
    excluded.}
  \item{include}{Boolean; default TRUE, but can be switched to FALSE to
    turn off an effect.}
  \item{name}{Name of dependent variable (network or behavior) for which
    statistics are being included. Defaults to the first in the effects object.}
  \item{interaction1}{
    Name of siena object where needed to completely identify the
    effects e.g. covariate name or behavior variable name.}
  \item{interaction2}{
    Name of siena object where needed to completely identify the
    effects e.g. covariate name or behavior variable name.}
  \item{character}{Boolean: are the statistic names character strings or not.}
  \item{verbose}{Boolean: should the print of altered statistic be produced.}
}
\details{The names \code{interaction1} and \code{interaction2} refer to
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
  \code{\link{getEffects}},
  \code{\link{includeEffects}},
  \code{\link{setEffect}},
  \code{\link{includeInteraction}},
  \code{\link{print.sienaEffects}}
}

\examples{
  mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  mybeh  <- sienaDependent(s50a, type="behavior")
  mydata <- sienaDataCreate(mynet1, mybeh)
  myeff <- getEffects(mydata)
  myeff <- includeGMoMStatistics(myeff, egoX_gmm, interaction1="mybeh")
  myeff
}
\keyword{ classes }
