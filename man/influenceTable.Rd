\name{influenceTable}
\alias{influenceTable}
\alias{influencetable}
\alias{influence.Table}
\alias{influence.table}
\title{Function to construct influence tables for SAOMs}
\description{
 The function \code{influenceTable} constructs influence tables
 which may be helpful for the interpretation of results
 for network and behavior dynamics,
 for an estimation result represented by a
 \code{sienaFit} object created by \code{\link{siena07}}.
}

\usage{
influenceTable(x, xd, netname, behname,
    as.matrix=FALSE, levls=NULL, levls.alt=levls,
    out.ego=1, silent=FALSE, nfirst=x$nwarm+1,
    include.endow=FALSE, include.creation=FALSE)
}
\arguments{
  \item{x}{An object of class \code{\link{sienaFit}}, produced by a
    call to \code{\link{siena07}} for a model including
    a behavioral dependent variable;
    or an object of class \code{\link{sienaMeta}}, produced by a
    call to \code{\link{siena08}};
    or an object of class \code{sienaBayesFit}, produced by a
    call to \code{sienaBayes}.  }
  \item{xd}{If \code{x} is of class \code{\link{sienaFit}},
    a \code{\link{siena}} data set used to produce \code{x}.\cr
    If \code{x} is of class \code{\link{sienaMeta}} or \code{sienaBayesFit},
    one of the data sets used to produce \code{x}; preferable a representative
    one.}
  \item{netname}{character string: name of network dependent variable.}
  \item{behname}{character string: name of behavior dependent variable.}
  \item{as.matrix}{boolean: will the table be returned as a matrix.}
  \item{levls}{levels for ego.}
  \item{levls.alt}{levels for alter.}
  \item{out.ego}{presumed outdegree of ego for effects \code{totSim} and
       \code{totAlt}.}
  \item{silent}{boolean; if \code{FALSE}, the parameters taken from \code{x}
     are reproduced at the console.}
  \item{nfirst}{If \code{x} is of class \code{sienaBayesFit}: first run in
       posterior sample used for constructing the table.}
  \item{include.endow}{boolean, used only for models utilizing influence and
        endowment effects: by specifying \code{include.endow=TRUE},                                         #
        the sum of evaluation and endowment effects is used.}
  \item{include.creation}{boolean, used only for models utilizing influence and
         creation effects: by specifying \code{include.creation=TRUE},                                         #
        the sum of evaluation and creation effects is used.}
}
\details{
 This functions is used for constructing influence tables
 for the interpretation of results for network and behavior dynamics
 obtained with the \code{RSiena} or \code{multiSiena} packages.

 In matrix form, each row corresponds to a given average behavior of the
 alters, to whom the focal actor (ego) is connected by an outgoing tie;
 the columns are the different potential values of ego’s own behavior.
 The rows are for the values in \code{levls.alt}, the columns
 for the values in \code{levls}.

 If \code{levls} is \code{NULL} (the default), the levels of
 ego's behavior are taken as the integer range of
 the dependent actor variable.
 In most applications \code{levls.alt} will be the same as \code{levls},
 in which case it does not have to be specified.

 The data set \code{xd} is only used to get means and ranges
 which are used somewhere in the effects.
 
 This function has its own print method.
}

\value{
    An object of class \code{influenceTable}.\cr
  This represents the joint contribution to the evaluation function
  of effects \code{"linear"}, \code{"quad"}, \code{"avAlt"}, \code{"avSim"},
  \code{"totAlt"}, \code{"totSim"},  \code{"avAttHigher"}, \code{"avAttLower"},
  \code{"threshold"}, \code{"threshold2"}, \code{"threshold3"},
  and \code{"threshold4"}.\cr
  If \code{as.matrix=TRUE}, the object is a matrix;
  if \code{as.matrix=FALSE} (the default), it is a data frame, where the rows
  are all combinations of \code{levls} and \code{levls.alt}, and the columns:
  \item{alter}{alter's values as a factor;}
  \item{zalter}{the values \code{levls.alt}}
  \item{zego}{the values \code{levls}}
  \item{select}{the value of the influence table.}
}

\references{
  \url{https://www.stats.ox.ac.uk/~snijders/siena/}\cr
  See the manual, Sections 13.2 and 13.4.
}
\author{Tom Snijders}
\seealso{\code{\link{siena07}}, \code{\link{siena08}} }
\examples{
mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mybeh  <- sienaDependent(s50a[,1:2], type="behavior")
mydata <- sienaDataCreate(mynet, mybeh)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, avAlt, name="mybeh", interaction1="mynet")
myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=100, seed=1291)
# nsub=2, n3=100 is used here for having a brief computation, not for practice.
ans <- siena07(myalgorithm, data=mydata, effects=myeff, silent=TRUE, batch=TRUE)
influenceTable(ans, mydata, "mynet", "mybeh")
}

\keyword{models}
