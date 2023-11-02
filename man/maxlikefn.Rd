\name{maxlikefn}
\alias{maxlikefn}
\title{A ML version of FRAN}
\description{
 A function to be called as "FRAN". All in R. very slow.   
 work in progress.
}
\usage{
maxlikefn(z, x, INIT = FALSE, TERM = FALSE, data,
effects = NULL, nstart = 1000, pinsdel = 0.6,
pperm = 0.3, prelins = 0.1, multfactor=2, promul = 0.1,
promul0 = 0.5, pdiaginsdel = 0.1, fromFiniteDiff = FALSE,
noSamples = 1, sampInterval = 50, int = 1)
}
\arguments{
  \item{z}{control object, passed in automatically in siena07}
  \item{x}{model object, passed in automatically in siena07}
  \item{INIT}{if TRUE, do initial processing. May be required to set up z}
  \item{TERM}{if TRUE, do end processing. }
  \item{data}{A siena object }
  \item{effects}{list of data frames as returned by getEffects}
  \item{nstart}{Number of MH steps at the start, after making the chain}
  \item{pinsdel}{Probability of insert/delete step
}
  \item{pperm}{
Probability of permutation step. (set to zero in startup phase.)
}
  \item{prelins}{Insertion probability in InsDelPermute
  }
  \item{multfactor}{Factor controlling number of MH steps. Will be read
	from the model in preference, and that is easier to alter! But I
	don"t want to alter that program yet..
  }

  \item{promul}{Probability of choosing a random single multiple in
	InsDelPermute in start up phase.
}
  \item{promul0}{Probability of choosing a random single multiple in
	InsDelPermute not in startup phase
}
  \item{pdiaginsdel}{Probability of insertion or deletion of a diagonal
	link.
}
  \item{fromFiniteDiff}{Should always be FALSE}
  \item{noSamples}{
	Number of chains to be returned
}
\item{sampInterval}{
  If multiple chains are returned, the number of steps between each
}
  \item{int}{
	Number of parallel MCMC chains to pursue.
  }
}
\details{
  This can be used for the element FRAN of the model object.
  The arguments with no defaults must be passed in on the call to
  siena07.
  Also you must set the option \code{maxlike=TRUE} in the call to
  \code{sienaAlgorithmCreate()}
}
\value{Depends on the call. If \code{INIT} or \code{initC} or
  \code{TERM}
  are true, returns \code{z}, the control object. Otherwise,
  returns a list containing:
  \item{fra}{Simulated scores}
  \item{dff}{2nd deriv, not phase 2}
  \item{OK}{could be set to FALSE if serious error has occurred}
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{\code{\link{siena07}} }
\examples{
\dontrun{
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff<- getEffects(mydata)
myalgor<- sienaAlgorithmCreate(nsub=2, n3=100, maxlike=TRUE)
ans<- siena07(myalgor, data=mydata, effects=myeff,  batch=TRUE)
}}
