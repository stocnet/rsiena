\name{n3401}
\alias{n3401}
\alias{N3401}
\alias{N3403}
\alias{N3404}
\alias{N3406}
\docType{data}
\title{Network data: excerpt from "Dutch Social Behavior Data Set"
of Chris Baerveldt.
}
\description{
Matrices N3401, N3403, N3404, N3406, and HN3401, HN3403, HN3404,  
HN3406 are two waves of networks for four schools 
(numbered 1, 3, 4, 6).  
}
\format{
 Adjacency matrices for the network at two time points.
 The matrices with name N... are the first wave,
 those with name HN... are the second wave.

 There is a tie from pupil i to pupil j if i says that he/she   
 receives and/or gives emotional support from/to pupil j.  
 The data are part of a larger data set (see source below)  
 and were collected under the direction of Chris Baerveldt.
}
\source{
\url{https://www.stats.ox.ac.uk/~snijders/siena/CB_data.zip}
}
\references{
Houtzager, B. and Baerveldt, C. (1999), Just like Normal.
A Social Network Study of the Relation between Petty Crime
and the Intimacy of Adolescent Friendships.
\emph{Social Behavior and Personality} \bold{27}, 177--192.

Snijders, Tom A.B, and Baerveldt, C. (2003),
A Multilevel Network Study of the Effects of Delinquent Behavior
on Friendship Evolution.
\emph{Journal of Mathematical Sociology} \bold{27}, 123--151.

See \url{https://www.stats.ox.ac.uk/~snijders/siena/BaerveldtData.html}
}
\examples{
mynet <- sienaDependent(array(c(N3401, HN3401), dim=c(45, 45, 2)))
mydata <- sienaDataCreate(mynet)
}
\keyword{datasets}
