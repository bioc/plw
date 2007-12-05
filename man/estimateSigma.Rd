\name{estimateSigma}
\alias{estimateSigma}
\title{Fit zero mean multivariate t-distribution, known df}
\usage{
estimateSigma(y, m, v, maxIter = 100, epsilon = 1e-06, verbose = FALSE)
}
\description{
Estimate the covariance matrix \eqn{\Sigma}{Sigma} of the multivariate t-distribution with zero expectation assuming the degrees of freedom is known.
}
\arguments{
  \item{y}{data matrix}
  \item{m}{degrees of freedom}
  \item{v}{scale parameter}
  \item{maxIter}{maximum number of iterations}
  \item{epsilon}{convergence criteria}
  \item{verbose}{print computation info or not}
}
\details{
The multivariate t-distribution is parametrized as:
\deqn{y|c \sim N(\mu,c\Sigma)}{y|c ~ N(mu,c*Sigma)}
\deqn{c \sim \mbox{InvGamma}(m/2,m\nu/2)}{c ~ InvGamma(m/2,m*v/2)}
Here \eqn{N}{N} denotes a multivariate normal distribution, 
\eqn{\Sigma}{Sigma} is a covariance matrix and 
\eqn{\mbox{InvGamma}(\alpha,\beta)}{InvGamma(a,b)} is the inverse-gamma distribution with density function
\deqn{f(x)=(\beta)^{\alpha} \exp\{-\beta/x\} x^{-\alpha-1}/\Gamma(\alpha)}{
      f(x)=b^a exp\{-b/x\} x^{-a-1} /Gamma(a)}
In this application \eqn{\mu}{mu} equals zero, and m is the degrees of freedom.
}
\value{
  \item{Sigma}{Estimated covariance matrix for y}
  \item{iter}{Number of iterations}
}
\references{

Hastie, T., Tibshirani, R., and Friedman, J. (2001). The Elements of Statistical Learning, volume 1. Springer, first edition.

Kristiansson, E., Sj\eqn{\mbox{\"o}}{o}gren, A., Rudemo, M., Nerman, O. (2005). Weighted Analysis of Paired Microarray Experiments. Statistical Applications in Genetics and Molecular Biology 4(1)

\eqn{\mbox{\AA}}{A}strand, M. et al. (2007a). Improved covariance matrix estimators for weighted analysis of microarray data. Journal of Computational Biology, Accepted.

\eqn{\mbox{\AA}}{A}strand, M. et al. (2007b). Empirical Bayes models for multiple-probe type arrays at the probe level. Bioinformatics, Submitted 1 October 2007.
}
\author{Magnus \eqn{\mbox{\AA}}{A}strand}
\seealso{estimateSigmaMV}
\keyword{htest}
\keyword{models}
\keyword{design}

