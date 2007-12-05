\name{logitT}
\alias{logitTTransform}
\alias{logitTStat}
\alias{studenttTTest}

\title{logit-t and t-test by row}
\usage{

logitTTransform(pm)
logitTStat(affy.batch,group)
studenttTTest(x, group)

}
\description{
Functions for the logit-t test (Lemon et al. 2003) and the ordinary t-test computed for each row of an matrix.
}
\arguments{
  \item{pm}{A matrix of Pm intensities}
  \item{affy.batch}{An AffyBatch object}
  \item{group}{A group indicator vector, should have values 1 and 2 only.}
  \item{x}{A matrix}
}
\details{
See the definition (R-code) of each function for details.
}
\value{

logitTTransform returns a matrix

logitTStat returns a vector with the logit-t statistic for each probe set.

studenttTTest returns a vector with t-statistic for each row of x.

}
\references{
Lemon et al. (2003). A high performance test of differential gene expression for oligonucleotide arrays. Genome Biol. 2003; 4(10):R67
}
\author{Magnus \eqn{\mbox{\AA}}{A}strand}

\examples{
# ------------------------------------------
# Example analyzing the 6 arrays in the 
# AffySpikeU95Subset data set

# Loading the data
data(AffySpikeU95Subset)

# Vector with groups assignment 
group<-factor(rep(1:2,each=3))

# logit-T statistic
logitT<-logitTStat(AffySpikeU95Subset,
                   as.numeric(group))

# Computing RMA expression index
data.rma<-exprs(rma(AffySpikeU95Subset))


# Ordinary t-test by row/gene
studentT<-studenttTTest(data.rma, as.numeric(group))

# Comparing genes ranked top-20
logitTTop20  <- rank(-abs(logitT)) < 21
studentTTop20<- rank(-abs(studentT)) < 21
table(logitTTop20,studentTTop20)

}
\keyword{univar}
\keyword{manip}
