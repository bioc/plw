\name{topRankSummary}
\alias{topRankSummary}
\alias{plotSummaryLog2FC}
\alias{plotSummaryT}


\title{Return or plots analysis result for top ranking or selected probe sets}
\usage{

plotSummaryLog2FC(model,nGenes=50,genesOfRank=1:nGenes,genes=NULL)
plotSummaryT(model,nGenes=50,genesOfRank=1:nGenes,genes=NULL)
topRankSummary(model,nGenes=50,genesOfRank=1:nGenes,genes=NULL)
}
\description{
Returns (or plots) t-statistic and/or log2FC for each probe and median for each probe set.
P
}
\arguments{
  \item{model}{On object obtained from the function plw.}
  \item{nGenes}{Gives summary for the nGenes top ranking genes}
  \item{genesOfRank}{Gives summary for genes ranked genesOfRank}
  \item{genes}{Gives summary for specific genes.}
}
\author{Magnus \eqn{\mbox{\AA}}{A}strand}
\seealso{plw}
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

## Selecting top genes
topRankSummary(model1,nGenes=10)

## Plotting t-statistics and log2FC for top genes
par(mfrow=c(1,2))
plotSummaryT(model1,nGenes=20)
plotSummaryLog2FC(model1,nGenes=20)
}
\keyword{univar}
\keyword{manip}
