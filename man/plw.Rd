\name{plw}
\alias{plw}
\title{Probe level Locally moderated Weighted median-t.}
\description{
Computes locally moderated weighted median t-test for microarray data.
}
\usage{
plw(x,design=rep(1,ncol(x)),contrast=matrix(1),
    probenames = unlist(ifelse(class(x) == "AffyBatch", 
                                list(p = probeNames(x)),
                                list(p = NULL))),
    maxIter = 200, epsilon = 1e-06, verbose = TRUE,
    nknots = 10, nOut = 2000, nIn = 4000, knots = NULL,
    checkRegulation = TRUE)
}
\arguments{
  \item{x}{Data, log2(PM) intensities or an AffyBatch object, see details}
  \item{design}{design matrix}
  \item{contrast}{contrast matrix}
  \item{probenames}{If not null, it is used to group PM probes into probe sets, see details.}
  \item{maxIter}{maximum number of iterations}
  \item{epsilon}{convergence criteria}
  \item{verbose}{print computation info or not}
  \item{nknots}{Number of knots of spline for \eqn{\nu}{v}}
  \item{nOut}{Parameter for calculating knots, see getKnots}
  \item{nIn}{Parameter for calculating knots, see getKnots}
  \item{knots}{Knots, if not NULL it overrides nknots, nOut and nIn}
  \item{checkRegulation}{If TRUE, data is checked for a correct specified contrast (see details)}
}
\details{
This function computes the Probe level Locally moderated Weighted median-t statistic (PLW) described in 
\eqn{\mbox{\AA}}{A}strand (2007b),
specially design for Affymetrix type data, or other microarray data with multiple probes.

The data object x should be either a matrix of perfect match (PM) intensities, or an object of class AffyBatch.
When x is a matrix of PM intensities, the intensities should be background corrected, normalized, and logged (with base 2).
If x is an AffyBatch object, the default background correction and normalization of RMA is applied to x.

When probenames is not null, it should be a vector of length equal to the number rows in the matrix x, giving the probe-set identity for each PM probe.
Use the function probeNames in the affy package to get probenames when x is a matrix of log2(PM) intensities.

Inference is done for each PM probe, thus moderated t-statistic, p-value and log2(FC) is calculated for each probe.
The median t-statistics for each probe-set is also computed.

Each PM probe g (row of x) is modeled as:
\deqn{y_g|c_g \sim N(\mu_g,c_g\Sigma)}{y_g|c_g ~ N(mu_g,c_g*Sigma)}
\deqn{c_g \sim \mbox{InvGamma}(m/2,m\nu/2)}{c ~ InvGamma(m/2,m*v/2)}
where \eqn{\nu}{v} is function of the mean intensity: \eqn{\nu(\bar\mu_g)}{v(mean(mu_g))},
\eqn{N}{N} denotes a multivariate normal distribution, 
\eqn{\Sigma}{Sigma} is a covariance matrix and 
\eqn{\mbox{InvGamma}(a,b)}{InvGamma(a,b)} is the inverse-gamma distribution with density function
\deqn{f(x)=(b)^{a} \exp\{-b/x\} x^{-a-1}/\Gamma(a)}{
      f(x)=b^a exp\{-b/x\} x^{-a-1} /Gamma(a)}

Given the design matrix D, \eqn{\mu_g}{mu_g} equals \eqn{D\gamma_g}{D*gamma_g}, 
and given the contrast matrix C the hypothesis \eqn{C\gamma_g=0}{C*gamma_g=0} is tested.
C should be a one row matrix of same length as the column vector \eqn{\gamma_g}{gamma_g}.

See examples on how to specify the design and contrast matrices.

A cubic spline is used to parameterize the smooth function \eqn{\nu(x)}{v(x)} 
\deqn{\nu(x) = \exp\{ H(x)^T \beta \}}{v(x)=exp{H(x)^T beta}}
where 
\eqn{H:R \to R^{2p-1}}{H:R->R^(2p-1)} 
is a set B-spline basis functions for a given set of p interior spline-knots, see chapter 5 of Hastie et al. (2001).

The parameter estimation procedure is based on the assumption that the specified contrast is close to zero
for most genes, or at least that the median contrast over all genes is close to zero. A check is run on data
to validate this assumptions. If the checking fails, with the error message "warning: most genes appears to be regulated..."
and if YOU ARE SURE that the design and contrast is correct, use checkRegulation=FALSE.


}
\value{
  \item{Sigma}{Estimated covariance matrix for \eqn{y=P^Tx}{y=P' x}}
  \item{m}{Estimated shape parameter for inverse-gamma prior for probe variances} 
  \item{v}{Estimated scale parameter curve for inverse-gamma prior for probe variances}
  \item{converged}{T if the EM algorithms converged}
  \item{iter}{Number of iterations}
  \item{modS2}{Moderated estimator of probe-specific variances}
  \item{histLogS2}{Histogram of log(s2) where s2 is the ordinary variance estimator}
  \item{fittedDensityLogS2}{The fitted density for log(s2)}
  \item{logs2}{Variance estimators, logged with base 2.}
  \item{medianT}{Median moderated t-statistic for each probe-set}
  \item{t}{Moderated t-statistic for each PM probe}
  \item{coefficients}{Estimated contrast for each PM probe}
  \item{p.value}{P-value from the moderated t-statistic for each PM probe}
  \item{dfT}{Degrees of freedom of the moderated t-statistic}
  \item{weights}{Weights for estimating the contrast}
  \item{P}{Transformation matrix}
  \item{beta}{Estimated parameter vector \eqn{\beta}{beta} of spline for \eqn{\nu(x)}{v(x)} }
  \item{knots}{The knots used in spline for \eqn{\nu(x)}{v(x)}}
  \item{x}{The input vector covariate vector x}
}
\references{

Hastie, T., Tibshirani, R., and Friedman, J. (2001). The Elements of Statistical Learning, volume 1. Springer, first edition.

Kristiansson, E., Sj\eqn{\mbox{\"o}}{o}gren, A., Rudemo, M., Nerman, O. (2005). Weighted Analysis of Paired Microarray Experiments. Statistical Applications in Genetics and Molecular Biology 4(1)

\eqn{\mbox{\AA}}{A}strand, M. et al. (2007a). Improved covariance matrix estimators for weighted analysis of microarray data. Journal of Computational Biology, Accepted.

\eqn{\mbox{\AA}}{A}strand, M. et al. (2007b). Empirical Bayes models for multiple-probe type arrays at the probe level. Bioinformatics, Submitted 1 October 2007.
}
\author{Magnus \eqn{\mbox{\AA}}{A}strand}
\seealso{estimateSigmaMVbeta, lmw}
\examples{

# ------------------------------------------
# Example analyzing the 6 arrays in the 
# AffySpikeU95Subset data set

# Loading the data
data(AffySpikeU95Subset)

# Defining design and contrast matrix
group<-factor(rep(1:2,each=3))
design<-model.matrix(~group-1)
contrast<-matrix(c(1,-1),1,2)

# Analyzing with an AffyBatch object as input
model1<-plw(AffySpikeU95Subset,design=design,contrast=contrast,
            epsilon=0.01)

## Look at fitted vs observed density for log(s2)
varHistPlot(model1)

## Look at fitted curve for scale parameter
scaleParameterPlot(model1)

## Selecting top genes
topRankSummary(model1,nGenes=10)

## Plotting t-statistics and log2FC for top genes
par(mfrow=c(1,2))
plotSummaryT(model1,nGenes=20)
plotSummaryLog2FC(model1,nGenes=20)


###---------------------------------------
# Analyzing with BG-adjusted and normalized PM data
pm1<-pm(bg.correct.rma(AffySpikeU95Subset, bgtype = 2))
pm2<-matrix(.C("qnorm_c", as.double(as.vector(pm1)), 
               as.integer(nrow(pm1)), 
               as.integer(ncol(pm1)))[[1]],
               nrow(pm1),ncol(pm1))
data<-log2(pm2)

probenames<-probeNames(AffySpikeU95Subset)
model2<-plw(data,design=design,contrast=contrast,
            probenames=probenames,epsilon=0.01)

###---------------------------------------
# Model1 and model2 should give identical result
# For example identical top ranking:
range(topRankSummary(model1)$t-
      topRankSummary(model2)$t,na.rm=TRUE)

}
\keyword{htest}
\keyword{models}
\keyword{design}
