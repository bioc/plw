\name{scaleParameterPlot}
\alias{scaleParameterPlot}

\title{Scale parameter plotted against mean intensity}
\usage{
scaleParameterPlot(model,main="Scale parameter curve",
           col=1,pch='.',lty=1,curveCol=2,knotsPch=19,knotsCol=3)
}
\description{
Will produce a scatter plot of variance estimators (logged) for each probe (probe set) against the corresponding mean intensity
together with the fitted scale-parameter curve and points showing the knots of the used spline. 
}
\arguments{
  \item{model}{On object obtained from the function plw or lmw.}
  \item{main}{Main title of plot.}
  \item{col}{Color for individual points (mean,logs2).}
  \item{pch}{Plot symbol for individual points (mean,logs2).}
  \item{lty}{Line type for fitted scale parameter curve.}
  \item{curveCol}{Line color for fitted scale parameter curve.}
  \item{knotsPch}{Plot symbol for spline knots.}
  \item{knotsCol}{Plot color for spline knots.}
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

## Look at fitted curve for scale parameter
scaleParameterPlot(model1)

}
\keyword{univar}
\keyword{manip}

