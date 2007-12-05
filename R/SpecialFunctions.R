## --------------------------------------------------------------------------
SSdensity<- function(x, m, v, df)
{
    t <- (1/m)^(-m/2) * v^(m/2) * exp(log(x)* ((df-2)/2) -log(m*v+x) * 
            ((m+df)/2) )  * gamma((m+df)/2) / ( gamma(m/2) * gamma(df/2) )
    return(t)
}

## --------------------------------------------------------------------------
logSSdensity<- function(x, m, v, df)
{
    return(exp(x)*SSdensity(exp(x), m, v, df))
}
## --------------------------------------------------------------------------
meanSdByRow<-function(mat){
    n<-nrow(mat)
    m<-ncol(mat)
    aux <- .C("MeanAndSd", x=as.double(t(mat)), n=as.integer(n),
             m=as.integer(m), mean=double(length=n), sd2=double(length=n),
             PACKAGE="plw")
    ii<-aux$sd2 < 0
    aux$sd2[ii]<-0
    return(list(mean=aux$mean, sd=sqrt(aux$sd2)))
}

## --------------------------------------------------------------------------
ReduzeToFullRank<-function(M){
    if (length(M) == 0)
        return(0)
    n<-nrow(M)
    i<-findLinearDepependentRow(M)
    while(n>1&i>0){
        M<-M[-i, ]
        n<-nrow(M)
        i<-findLinearDepependentRow(M)
    }
    return(M)
}

## --------------------------------------------------------------------------
findLinearDepependentRow<-function(M){
    i<-0
    n<-nrow(M)
    if(n==1)
        return(0)
    rankM<-computeRank(M)
    found<-FALSE
    while(!found & i<n){
        i<-i+1
        ri<-computeRank(M[-i, ])
        if(ri==rankM) found<-TRUE
   }
   return(i*found)
}

## --------------------------------------------------------------------------
computeRank<-function (M, eps = 10^(-5))
{
    if (length(M) == 0)
        return(0)
    return(sum(abs(svd(M)$d) > eps))
}

## --------------------------------------------------------------------------
## kovarians skattning givet mean=mu==0
## OBS y skall vara nxp, och sc 1xn
covZeroMeanScaledData<-function(y, sc){
    n<-nrow(y)
    p<-ncol(y)
    aux<-.C("cov_zero_mean_scaled_data", y=as.double(y), sc=as.double(sc), 
            n=as.integer(n), m=as.integer(p), covmat=double(length=(p*p)), 
            PACKAGE="plw")$covmat
    return(matrix(aux, p, p))
}

## --------------------------------------------------------------------------
## kovarians skattning för givet mean=mu
covFixedMean<-function(y, mu=rep(0,ncol(y))){
    y<-scale(y, center=mu, scale=FALSE)
    p<-ncol(y)
    n<-nrow(y)
    est1<-cov(y)
    est2<-rep(1/n, n)%*%y
    est3<-est1*(n-1)/n+matrix(est2, ncol=1)%*%est2
    return(est3)
}

## --------------------------------------------------------------------------
## kovarians skattning för givet mean=mu
weightedCovFixedMean<-function(y, mu=rep(0, ncol(y)), weights=rep(1, nrow(y))  ){
    weights<-weights/mean(weights)
    y<-scale(y, center=mu,scale=FALSE)
    p<-ncol(y)
    res<-matrix(0, p, p)
    for(i in 1:p) { for(j in i:p) {   res[i, j]<-res[j, i]<-mean(y[, i]*y[, j]*weights) }}
    return(res)
}


## --------------------------------------------------------------------------
## maximering mha golden ratio metoden, bra när mål funktionen är tung att beräkna
## Antar att max för funktionen f ligger mellan a och b
goldenRatio<-function(a, b, f, eps){
    r<-(sqrt(5)-1)/2
    r<-c(1-r, r)
    x<-a+r*(b-a)
    f1<-f(x[1])
    f2<-f(x[2])
    while(abs(a-b)>eps){
        if(f2<f1){
            b<-x[2]
            x<-a+r*(b-a)
            f2<-f1
            f1<-f(x[1])
        }
        else{
            a<-x[1]
            x<-a+r*(b-a)
            f1<-f2
            f2<-f(x[2])
        }
    }       
    return(mean(x))
}


## --------------------------------------------------------------------------
matrixSquareRot<-function(x){
    temp<-eigen(x)
    if(!all(temp$values>0)){
        cat("Matrix is not positive definite\n")
        return(NULL)
    }
    else{
        sigma.5<-t(temp$vectors%*%diag(sqrt(temp$values)))
        return(sigma.5)
    }
}



## --------------------------------------------------------------------------
estimateMV<-function(SS, dfss, max.iter=100, epsilon=1e-06, verbose=FALSE){
## Givet SS, som är kvadrat summa med df=dfss, anpassa en prior

    ## Tar fram skattningar MHA moment 
    ii<- SS>0 & is.finite(SS)
    temp<-var(log(SS[ii]/dfss))-trigamma(dfss/2)
    if(temp>0){
        m0<-2*triGammaInverse(temp)
    }
    else{
        m0<-dfss*5
    }
    v0<-exp(mean(log(SS[ii]/dfss))-log(m0)+log(dfss)+
                        digamma(m0/2)-digamma(dfss/2))

    ## Tar fram ML skattningar
    if(length(SS[ii])>5000){
        h1<-hist(log(SS[ii]), breaks=1000, plot=FALSE)
    }
    else{
        n<-length(SS[ii])
        h1<-list(mids=log(SS[ii]), density=rep(1/n, n))
    }

    f1<-function(x, m, v){ return(   df(x/(dfss*v), dfss, m)/(dfss*v)  ) }
    f<-function(mv){ 
        m<-exp(mv[1])
        v<-exp(mv[2])
        return(-(sum(((h1$mids)+log(f1(exp(h1$mids), m, v)))*h1$density)))
    }
    opt.logmv<-optim(log(c(m0, v0)), f)
    m0<-exp(opt.logmv$par[1])
    v0<-exp(opt.logmv$par[2])
    return(list(m=m0, v=v0, histSS=h1, fitted.density.SS=list(x=h1$mids,
                y=exp(h1$mids)*f1(exp(h1$mids), m0, v0))))
}

## --------------------------------------------------------------------------
triGammaInverse<-function(x)
{
    if (!is.numeric(x))
        stop("Non-numeric argument to mathematical function")
    if (length(x) == 0)
        return(numeric(0))
    omit <- is.na(x)
    if (any(omit)) {
        y <- x
        if (any(!omit))
            y[!omit] <- Recall(x[!omit])
        return(y)
    }
    omit <- (x < 0)
    if (any(omit)) {
        y <- x
        y[omit] <- NaN
        warning("NaNs produced")
        if (any(!omit))
            y[!omit] <- Recall(x[!omit])
        return(y)
    }
    omit <- (x > 1e+07)
    if (any(omit)) {
        y <- x
        y[omit] <- 1/sqrt(x[omit])
        if (any(!omit))
            y[!omit] <- Recall(x[!omit])
        return(y)
    }
    omit <- (x < 1e-06)
    if (any(omit)) {
        y <- x
        y[omit] <- 1/x[omit]
        if (any(!omit))
            y[!omit] <- Recall(x[!omit])
        return(y)
    }
    y <- 0.5 + 1/x
    iter <- 0
    repeat {
        iter <- iter + 1
        tri <- trigamma(y)
        dif <- tri * (1 - tri/x)/psigamma(y, deriv = 2)
        y <- y + dif
        if (max(-dif/y) < 1e-08)
            break
        if (iter > 50) {
            warning("Iteration limit exceeded")
            break
        }
    }
    y
}
