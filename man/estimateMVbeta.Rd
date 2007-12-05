\name{estimateMVbeta}
\alias{estimateMVbeta}
\title{Zero mean multivariate t-dist. with covariate dependent scale.}
\description{
Estimate the parameters m and \eqn{\nu}{v} of the multivariate t-distribution with zero expectation,
where \eqn{\nu}{v} is modeled as smooth function of a covariate. 
The covariance matrix \eqn{\Sigma}{Sigma} is assumed to be known.
}
\usage{
estimateMVbeta(y, x, Sigma, maxIter = 200, epsilon = 1e-06,
    verbose = FALSE, nknots = 10, nOut = 2000, nIn = 4000,
    iterInit = 3, br = NULL)
}
\arguments{
  \item{y}{Data matrix}
  \item{x}{Covariate vector}
  \item{Sigma}{Covariance matrix}
  \item{maxIter}{Maximum number of iterations}
  \item{epsilon}{Convergence criterion}
  \item{verbose}{Print computation info or not}
  \item{nknots}{Number of knots of spline for \eqn{\nu}{v}}
  \item{nOut}{Parameter for calculating knots, see getKnots}
  \item{nIn}{Parameter for calculating knots, see getKnots}
  \item{iterInit}{Number of iteration in when initiating \eqn{\Sigma}{Sigma}}
  \item{br}{Knots, overrides nknots, n.out and n.in}
}
\details{
The multivariate t-distribution is parametrized as:

\deqn{y|c \sim N(\mu,c\Sigma)}{y|c ~ N(mu,c*Sigma)}
\deqn{c \sim \mbox{InvGamma}(m/2,m\nu/2)}{c ~ InvGamma(m/2,m*v/2)}
where \eqn{\nu}{v} is function of the covariate x: \eqn{\nu(x)}{v(x)} and
\eqn{N}{N} denotes a multivariate normal distribution, 
\eqn{\Sigma}{Sigma} is a covariance matrix and 
\eqn{\mbox{InvGamma}(\alpha,\beta)}{InvGamma(a,b)} is the inverse-gamma distribution with density function
\deqn{f(x)=(\beta)^{\alpha} \exp\{-\beta/x\} x^{-\alpha-1}/\Gamma(\alpha)}{
      f(x)=b^a exp\{-b/x\} x^{-a-1} /Gamma(a)}


A cubic spline is used to parameterize the smooth function \eqn{\nu(x)}{v(x)} 
\deqn{\nu(x) = \exp\{ H(x)^T \beta \}}{v(x)=exp{H(x)^T beta}}
where 
\eqn{H:R \to R^{2p-1}}{H:R->R^(2p-1)} 
is a set B-spline basis functions for a given set of p interior spline-knots, see chapter 5 of Hastie (2001).
In this application \eqn{\mu}{mu} equals zero, and m is the degrees of freedom.
}
\value{
  \item{Sigma}{The input covariance matrix for y}
  \item{m}{Estimated shape parameter for inverse-gamma prior for gene variances} 
  \item{v}{Estimated scale parameter curve for inverse-gamma prior for gene variances}
  \item{converged}{TRUE if the EM algorithms converged}
  \item{iter}{Number of iterations}
  \item{modS2}{Moderated estimator of gene-specific variances}
  \item{histLogS2}{Histogram of log(s2) where s2 is the ordinary variance estimator}
  \item{fittedDensityLogS2}{The fitted density for log(s2)}
  \item{logs2}{Variance estimators, logged with base 2.}
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
\seealso{plw, lmw, estimateSigmaMVbeta}
\keyword{htest}
\keyword{models}


