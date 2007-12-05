\name{statByIndex}
\alias{madByIndex}
\alias{meanByIndex}
\alias{medianByIndex}
\alias{orderStatByIndex}
\alias{sdByIndex}
\alias{meanSdByRow}

\title{Computes statistics by index or by row}
\usage{
madByIndex(x,index)
meanByIndex(x,index)
medianByIndex(x,index)
orderStatByIndex(x,index,orderStat)
sdByIndex(x,index)
meanSdByRow(mat)
}
\description{
These function give the same result as 
   by(x,index,mad)
   by(x,index,mean)
   by(x,index,median)
but are much faster.
NOTE: The index vector is assumed to be SORTED and should contain INTEGER values only.

The function meanSdByRow computes mean and standard deviation for each row of the matrix mat.
A list with mean and sd is returned and gives the the same result as:

    list(mean=apply(mat,1,mean),sd=apply(mat,1,sd))
}
\arguments{
  \item{x}{Data vector}
  \item{index}{Index vector}
  \item{orderStat}{Which order statistic to compute}
  \item{mat}{Matrix}
}
\details{
See the definition (R-code) of each function for details.
}
\value{
 All but the last function:
 A vector with the statistic for each level if index.
 meanSdByRow:
 A list with items mean and sd. 
}
\author{Magnus \eqn{\mbox{\AA}}{A}strand}
\seealso{by, apply}
\examples{

## Example 1
## Computing, mad, mean and median by index.
## Compares with the result obtained using by(...) 

n<-10000
x<-rnorm(n)
index<-sort(round(runif(n,0.5,10.5)))

mad1<-madByIndex(x,index)
mad2<-by(x,index,mad)

mean1<-meanByIndex(x,index)
mean2<-by(x,index,mean)

median1<-medianByIndex(x,index)
median2<-by(x,index,median)

par(mfrow=c(2,2),mar=c(4,4,1.5,.5),mgp=c(1.5,.25, 0))
plot(mad1,mad2,main="Comparing mad",pch=19)
abline(a=0,b=1,col=2)
plot(mean1,mean2,main="Comparing mean",pch=19)
abline(a=0,b=1,col=2)
plot(median1,median2,main="Comparing median",pch=19)
abline(a=0,b=1,col=2)


## Example 2
## Computing, median by index
## Compares with the running time when using by(...)
n<-200000
x<-rnorm(n)
index<-sort(round(runif(n,0.5,10.5)))

system.time(median1<-medianByIndex(x,index))

system.time(median2<-by(x,index,median))


## Example 3
## Computing, mean and sd by row 
## Compares with using apply
nrow<-5000
ncol<-20
mat<-matrix(rnorm(ncol*nrow),nrow,ncol)

system.time(res1<-meanSdByRow(mat))
system.time(res2<-list(mean=apply(mat,1,mean),sd=apply(mat,1,sd)))

par(mfrow=c(1,2),mar=c(4,4,1.5,.5),mgp=c(1.5,.25, 0))
plot(res1$mean,res2$mean,pch='.')
plot(res1$sd,res2$sd,pch='.')

}
\keyword{univar}
\keyword{manip}

