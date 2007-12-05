\name{getKnots}
\alias{getKnots}

\title{Spline-knots for plw and lmw}
\usage{
getKnots(x,nKnots=10,nOut=2000,nIn=4000)
}
\description{
Computes a set of nKnots interior knots( if possible) plus 2 boundary knots so that:

1) the nOut smallest and highest data points (in x) lies below and above the lower and upper boundary knots respectively.

2) there is at least nIn data points between all knots.

}
\arguments{
  \item{x}{Data vector}
  \item{nKnots}{Number of interior knots}
  \item{nOut}{Number of data points below and above the lower and upper boundary knots respectively.}
  \item{nIn}{Number of data points between knots.}
}
\details{
See the definition (R-code) for details.
}
\value{
A vector of knots.
}
\author{Magnus \eqn{\mbox{\AA}}{A}strand}
\seealso{plw, lmw, estimateSigmaMVbeta}
\keyword{univar}
\keyword{manip}

