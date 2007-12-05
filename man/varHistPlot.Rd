\name{varHistPlot}
\alias{varHistPlot}

\title{Variance histogram and density}
\usage{
varHistPlot(model,main="Histogram variance estimators",
           histCol=8,densityCol=2,drawLegend=TRUE)
}
\description{
Will produce a histogram of observed variance estimators (logged) together with the fitted density.
}
\arguments{
  \item{model}{On object obtained from the function plw or lmw.}
  \item{main}{Main title of plot.}
  \item{histCol}{Color for histogram bars.}
  \item{densityCol}{Color for density function.}
  \item{drawLegend}{To draw a legend or not.}
}
\author{Magnus \eqn{\mbox{\AA}}{A}strand}
\seealso{plw, lmw}
\examples{

# ------------------------------------------
# Example using the result of the analysis of
# the 6 arrays in the AffySpikeU95Subset data set

# Loading the data
data(AffySpikeU95Subset)

# Defining design and contrast matrix
group<-factor(rep(1:2,each=3))
design<-model.matrix(~group-1)
contrast<-matrix(c(1,-1),1,2)

# Analyzing using plw
model1<-plw(AffySpikeU95Subset,design=design,contrast=contrast,
            epsilon=0.01)

## Look at fitted vs observed density for log(s2)
varHistPlot(model1)

}
\keyword{univar}
\keyword{manip}

