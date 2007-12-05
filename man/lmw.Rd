\name{lmw}
\alias{lmw}
\title{Locally Moderated Weighted-t.}
\description{
Computes Locally Moderated Weighted t-test for microarray data.
}
\usage{
lmw(x,design=rep(1,ncol(x)),contrast=matrix(1), meanX=NULL,
    maxIter = 200, epsilon = 1e-06, verbose = TRUE,
    nknots = 10, nOut = 2000, nIn = 4000, knots = NULL,
    checkRegulation = TRUE)
}
\arguments{
  \item{x}{Data, log2 expression indexes.}
  \item{design}{design matrix}
  \item{contrast}{contrast matrix}
  \item{meanX}{Covariate used to model scale parameter, default=NULL (see details)}
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
This function computes the Locally Moderated Weighted-t statistic (LMW) described in \eqn{\mbox{\AA}}{A}strand (2007b), thus calculating
locally moderated weighted t-statistic, p-value and log2(FC) for each row of the data matrix x.

Each gene g (row of x) is modeled as:
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

For details about the model see Kristiansson et al. (2005), \eqn{\mbox{\AA}}{A}strand et al. (2007a,2007b).

As specified above, \eqn{\nu}{v} is modeled as a function of mean intensity: \eqn{\nu(\bar\mu_g)}{v(mean(mu_g))}.
If the parameter meanX is not NULL, meanX is used instead of the mean intensity when modeling \eqn{\nu}{v}.
Thus, if meanX is not NULL, meanX must be a vector of length equal to the number of rows of the data matrix x.

The parameter estimation procedure is based on the assumption that the specified contrast is close to zero
for most genes, or at least that the median contrast over all genes is close to zero. A check is run on data
to validate this assumptions. If the checking fails, with the error message "warning: most genes appears to be regulated..."
and if YOU ARE SURE that the design and contrast is correct, use checkRegulation=FALSE.

}
\value{
  \item{Sigma}{Estimated covariance matrix for \eqn{y=P^Tx}{y=P' x}}
  \item{m}{Estimated shape parameter for inverse-gamma prior for gene variances} 
  \item{v}{Estimated scale parameter curve for inverse-gamma prior for gene variances}
  \item{converged}{T if the EM algorithms converged}
  \item{iter}{Number of iterations}
  \item{modS2}{Moderated estimator of gene-specific variances}
  \item{histLogS2}{Histogram of log(s2) where s2 is the ordinary variance estimator}
  \item{fittedDensityLogS2}{The fitted density for log(s2)}
  \item{logs2}{Variance estimators, logged with base 2.}
  \item{t}{Moderated t-statistic}
  \item{coefficients}{Estimated contrast}
  \item{p.value}{P-value from the moderated t-statistic}
  \item{dfT}{Degrees of freedom of the moderated t-statistic}
  \item{weights}{Weights for estimating the contrast}
  \item{P}{Transformation matrix}
  \item{beta}{Estimated parameter vector \eqn{\beta}{beta} of spline for \eqn{\nu(x)}{v(x)} }
  \item{knots}{The knots used in spline for \eqn{\nu(x)}{v(x)}}
}
\references{

Hastie, T., Tibshirani, R., and Friedman, J. (2001). The Elements of Statistical Learning, volume 1. Springer, first edition.

Kristiansson, E., Sj\eqn{\mbox{\"o}}{o}gren, A., Rudemo, M., Nerman, O. (2005). Weighted Analysis of Paired Microarray Experiments. Statistical Applications in Genetics and Molecular Biology 4(1)

\eqn{\mbox{\AA}}{A}strand, M. et al. (2007a). Improved covariance matrix estimators for weighted analysis of microarray data. Journal of Computational Biology, Accepted.strand

\eqn{\mbox{\AA}}{A}strand, M. et al. (2007b). Empirical Bayes models for multiple-probe type arrays at the probe level. Bioinformatics, Submitted 1 October 2007.
}
\author{Magnus \eqn{\mbox{\AA}}{A}strand}
\seealso{estimateSigmaMVbeta, plw}
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

# Computing RMA expression index
data.rma<-exprs(rma(AffySpikeU95Subset))

# Analyzing
model1<-lmw(data.rma,design=design,contrast=contrast,epsilon=0.01)

## Look at fitted vs observed density for log(s2)
varHistPlot(model1)

## Look at fitted curve for scale parameter
scaleParameterPlot(model1)

}
\keyword{htest}
\keyword{models}
\keyword{design}
