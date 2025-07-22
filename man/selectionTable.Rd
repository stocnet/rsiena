\name{selectionTable}
\alias{selectionTable}
\alias{selectiontable}
\alias{selection.Table}
\alias{selection.table}
\title{Function to construct selection tables for SAOMs}
\description{
 The function \code{selectionTable} constructs a selection table
 which may be helpful for the interpretation of results
 for network dynamics,
 for an estimation result represented by a
 \code{sienaFit} object created by \code{\link{siena07}}.
}

\usage{
selectionTable(x, xd, name, vname,
    as.matrix=FALSE,levls=NULL, levls.alt=levls,
    nfirst=x$nwarm+1, multiplier=1,
    include.endow=FALSE, include.creation=FALSE,
    silent=FALSE)
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
  \item{name}{character string: name of network dependent variable.}
  \item{vname}{character string: name of actor covariate (should be centered!).}
  \item{as.matrix}{boolean: will the table be returned as a matrix.}
  \item{levls}{levels for ego.}
  \item{levls.alt}{levels for alter.}
  \item{nfirst}{If of class \code{sienaBayesFit}: first run in
       posterior sample used for constructing the table.}
  \item{multiplier}{multiplier for the range of the actor covariate.}
  \item{include.endow}{boolean, used only for models utilizing selection and
        endowment effects: by specifying \code{include.endow=TRUE},
        the sum of evaluation and endowment effects is used.}
  \item{include.creation}{boolean, used only for models utilizing selection and
         creation effects: by specifying \code{include.creation=TRUE},
        the sum of evaluation and creation effects is used.}
  \item{silent}{boolean; if \code{FALSE}, some information
     is reproduced at the console.}
}
\details{
 This functions is used for constructing selection tables
 for the interpretation of results for network dynamics
 obtained with the \code{RSiena} or \code{multiSiena} packages.

 In matrix form, each row corresponds to ego's value of the covariate,
 and each column to alter's value.
 The table entries are the joint contributions of covariate effects
 to the objective function, for this covariate and for
 the combinations of ego and alter values.
  Effects currently implemented are \code{"altX"}, \code{"altSqX"},
  \code{"egoX"}, \code{"egoSqX"},
  \code{"egoXaltX"}, \code{"simX"}, \code{"diffX"}, \code{"diffSqX"},
  \code{"higher"}, \code{"sameX"}, \code{"egoDiffX"}, and
  \code{"egoPlusAltX"}.\cr

 If \code{levls} is \code{NULL} (the default), the levels of
 ego's behavior are taken as the integer range of
 the dependent actor variable.
 In most applications \code{levls.alt} will be the same as \code{levls},
 in which case it does not have to be specified.

 The \code{multiplier} is used in case the variable \code{vname} has a different
 natural scale, and the values for ego and alter as reported
 should be multiplied.
 The values given for \code{levls} and \code{levls.alt} are before this
 multiplication, i.e., they are the values as in the data set.

 The data set \code{xd} is only used to get the means and similarity means
 which are subtracted somewhere in the effects.

 This function has its own print method.
}

\value{
  An object of class  \code{selectionTable}.\cr
  If \code{as.matrix=TRUE}, this object is a matrix;
  if \code{as.matrix=FALSE} (the default), it is a data frame, where the rows
  are all combinations of \code{levls} and \code{levls.alt}, and the columns:
  \item{ego}{ego's values as a factor;}
  \item{vego}{the values \code{levls}}
  \item{valter}{the values \code{levls.alt}}
  \item{select}{the value of the selection table.}
}

\references{
  \url{https://www.stats.ox.ac.uk/~snijders/siena/}\cr
  See the manual, Sections 13.2 and 13.4.
}
\author{Tom Snijders}
\seealso{\code{\link{siena07}}, \code{\link{siena08}} }
\examples{
mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mycov  <- coCovar(s50a[,1])
mydata <- sienaDataCreate(mynet, mycov)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, simX, interaction1="mycov")
myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=100, seed=1291)
# nsub=2, n3=100 is used here for having a brief computation, not for practice.
ans <- siena07(myalgorithm, data=mydata, effects=myeff, silent=TRUE, batch=TRUE)
selectionTable(ans, mydata, "mynet", "mycov")
}

\keyword{models}
