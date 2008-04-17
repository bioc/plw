## ---------------------------------------------------------------------------
madByIndex<-function(x, index)
{
    if(length(x)!=length(index)){
        cat("error in madByIndex:  arguments x and index must", 
            "be of equal length", "\n")
        return (x)
    }
    aux1<-.C("mad_by_index", x=as.double(x), index=as.integer(index), 
             n=as.integer(length(x)), res=integer(length=1), PACKAGE="plw")

    return( aux1$x[1:aux1$res]/qnorm(3/4) )
}

## ---------------------------------------------------------------------------
medianByIndex<-function(x, index)
{
    if(length(x)!=length(index)){
        cat("error in medianByIndex:  arguments x and index must", 
            "be of equal length", "\n")
        return (x)
    }
    aux1<-.C("median_by_index", x=as.double(x), index=as.integer(index), 
             n=as.integer(length(x)), res=integer(length=1), PACKAGE="plw")

    return( aux1$x[1:aux1$res] )
}

## ---------------------------------------------------------------------------
meanByIndex<-function(x, index)
{
    if(length(x)!=length(index)){
        cat("error in meanByIndex:  arguments x and index must", 
            "be of equal length", "\n")
        return (x)
    }
    aux1<-.C("mean_by_index", x=as.double(x), index=as.integer(index), 
             n=as.integer(length(x)), res=integer(length=1), PACKAGE="plw")

    return( aux1$x[1:aux1$res] )
}

## ---------------------------------------------------------------------------
sdByIndex<-function(x, index)
{
    if(length(x)!=length(index)){
        cat("error in sdByIndex:  arguments x and index must", 
            "be of equal length", "\n")
        return (x)
    }
    aux1<-.C("sd_by_index", x=as.double(x), index=as.integer(index), 
             n=as.integer(length(x)), res=integer(length=1), PACKAGE="plw")

    return( aux1$x[1:aux1$res] )
}

## ---------------------------------------------------------------------------
orderStatByIndex<-function(x, index, orderStat)
{
    if(length(x)!=length(index)){
        cat("error in orderStatByIndex:  arguments x and index must", 
            "be of equal length", "\n")
        return (x)
    }
    aux1<-.C("order_stat_by_index", x=as.double(x), index=as.integer(index), 
             q=as.integer(orderStat-1), n=as.integer(length(x)), 
             res=integer(length=1), PACKAGE="plw")
    return( aux1$x[1:aux1$res] )
}


## ---------------------------------------------------------------------------
getKnots<-function(x, nKnots=10, nOut=2000, nIn=4000)
{
    if( length(x)< (nOut*2+nIn) ){ ## Not possible to set knots as specified
        br<-range(x)
        br<-seq(br[1], br[2], length=5)[2:4]
        cat(" Not possible to set knots as specified.", "\n", 
            "3 equally spaced knots will be used", "\n")
        return(br)
    }

    xSorted <- sort(x)
    n <- length(xSorted)
    bLow <- xSorted[c(nOut, n - nOut)]
    bUp <- bLow[2]
    bLow <- bLow[1]
    ii <- xSorted > bLow[1] & xSorted < bUp[1]
    xii <- xSorted[ii]
    br <- seq(bLow[1], bUp[1], length = nKnots + 2)
    h1 <- hist(xii, breaks = br, plot = FALSE)
    while (min(h1$counts) < nIn & nKnots > 2 & sum(h1$counts) > nIn) {
        i <- which.min(h1$counts)
        if (i == 1) {
            bLow <- c(xii[nIn], bLow)
            nKnots <- nKnots - 1
        }
        if (i == length(h1$counts)) {
            bUp <- c(xii[sum(h1$counts) - nIn], bUp)
            nKnots <- nKnots - 1
        }
        if (i > 1 & i < length(h1$counts)) {
            nKnots <- nKnots - 1
        }
        ii <- xSorted > bLow[1] & xSorted < bUp[1]
        xii <- xSorted[ii]
        br <- seq(bLow[1], bUp[1], length = nKnots + 2)
        h1 <- hist(xii, breaks = br, plot = FALSE)
    }
    if (sum(h1$counts) <= nIn) {
        br <- range(br)
    }
    br <- sort(c(bLow, br[-c(1, length(br))], bUp))
    return(br)
}

## ---------------------------------------------------------------------------
getSubset<-function(x, brr, nn, nMin)
{
    n<-length(x) 
    subSet<-NULL
    for(i in 1:length(nn)){
        ii<-(1:n)[x>=brr[i] & x<=brr[i+1]]
        if(nn[i]>nMin){
            ii<-sort(sample(ii, nMin))
        }
        subSet<-c(subSet, ii)
    }
    return(subSet)
}
    


## ---------------------------------------------------------------------------
estimateSigmaMVbeta<-function(y, x, maxIter=200, epsilon=0.000001, 
                              verbose=FALSE, nknots=10, nOut=2000, nIn=4000,
                              iterInit=3, br=NULL)
{
    p<-ncol(y)
    n<-nrow(y)
    sigma0<-matrix(0, p, p)
    diag(sigma0)<-1

    orderX<-order(x)

    ## Computing knots, inner, outer + boundary
    ## creates X matrix for the spline
    if(length(br)==0){
        br<-getKnots(x, nKnots=nknots, nOut=nOut, nIn=nIn)
        br<-sort(c(mean(br[1:2]), br))
    }
    mat1<-ns(x, knots=br[-c(1, length(br))], 
             Boundary.knots=br[c(1, length(br))])
    mat2<-cbind(rep(1, n), mat1)
    meanMat2<-c(1, rep(1, n)%*%mat1/n)

    est0<-NULL
    if(iterInit>0){
        if(verbose) cat("Init sigma using estimateSigmaMV\n")
        est0<-estimateSigmaMV(y, maxIter=iterInit, verbose=verbose)
        sigma0<-est0$Sigma   
    }
    ty<-t(y)

    if(verbose) cat("Init m and spline for v\n")
    ## Setting start values for m0 and v0, using moments for log(s2hat)
    s2hat <-.C("getSS", A=as.double(ginv(sigma0)), y=as.double(ty), 
               n=as.integer(n), p=as.integer(p), res=double(length=n),
               PACKAGE="plw")$res/p
    ii<- (s2hat > 0 & is.finite(s2hat))
    s2hatMin<-min(s2hat[ii])
    s2hat[!ii]<-s2hatMin 
    m1<-lm(log(s2hat)~mat1)
    v0<-exp(m1$fitted.values)
    beta0<-m1$coefficients
    est1<-estimateMV(s2hat*p/(v0), p)
    v0<-v0*est1$v
    beta0[1]<-beta0[1]+log(est1$v)
    m0<-est1$m

    iter<-0
    converged<-FALSE

    if(verbose) cat("Iterating EM-steps\n iter      sigma      m, v\n")

    while(iter<maxIter & !converged){
        iter<-iter+1
        temp2 <-.C("getSS", A=as.double(ginv(sigma0)), y=as.double(ty),
                   n=as.integer(n), p=as.integer(p), res=double(length=n),
                   PACKAGE="plw")$res

        ii<- (temp2 > 0 & is.finite(temp2))
        s2hatMin<-min(temp2[ii])
        temp2[!ii]<-s2hatMin

        temp2 <- temp2+m0*v0

        ## Optimizing beta 
        delta<-(m0+p)/temp2
        aux<-specialOptim(mat2, delta, beta0, meanMat2, 
                          reltol=sqrt(.Machine$double.eps)/
                                 ifelse(iter<10, 1, 100))
        betaHat<-aux$par

        vHat<-exp(betaHat[1]+mat1%*%betaHat[-1])
    
        ## Optimizing m 
        temp1<-beta0[1]+mat1%*%beta0[-1]
        temp1<-mean(temp1)-mean(log(temp2/2))+digamma((m0+p)/2)-
               mean(exp(temp1)*(m0+p)/temp2)

        f<-function(m){
            return( (m/2*(log(m/2)+temp1) - log(gamma(m/2))) )
        }
        mHat<-0.2+(0:1)
        logLik<-f(mHat)
        while(diff(logLik)>0){
            mHat<-mHat+1;
            logLik<-f(mHat)
        }
        if(mHat[1]>0.2){
            mHat[1]<-mHat[1]-1
        }

        ## Optimal m is between  mHat[1] and mHat[2]
        ## use golden ratio to optimize m
        mHat<-goldenRatio(mHat[1], mHat[2], f, epsilon/100)

        ## Optimize sigma
        sigmaHat<-covZeroMeanScaledData(y, (m0+p)/temp2)  

        ## Scaling v and sigma
        betaHat[1]<-betaHat[1]+log(mean(diag(sigmaHat)))
        sigmaHat<-sigmaHat/mean(diag(sigmaHat))

        ## Difference from previous iteration
        diffSigma<-as.vector(sigmaHat-sigma0)
        diffMV<-c(m0, log(v0)/10)-c(mHat, log(vHat)/10)

        if(verbose){
            cat(sprintf("%5d      %1.2e   %1.2e", iter, max(abs(diffSigma)),
                max(abs(diffMV))), "\n")
        }
        absDiff<-max(abs(c(diffSigma, diffMV)))
   
        converged<- (absDiff<epsilon)
        ##Uppdating estimate
        m0<-mHat
        beta0<-betaHat
        sigma0<-sigmaHat
        v0<-vHat
   }

   if(verbose) {
        cat("Done!\nSummarizing\n")
   }
   temp<-ginv(sigma0)%*%ty 
   logs2 <- log(.C("DiagAtimesBv2", a=as.double(ty), b=as.double(temp), 
                   n=as.integer(n), m=as.integer(p), diag=double(length=n),
                   PACKAGE="plw")$diag/p  )

   meanLogs2<-log(v0*m0)-log(p)+digamma(p/2)-digamma(m0/2)
   modS2<-temp2/(m0+p)

   dHat<-density(logs2-log(v0))
   d0<-list(x=dHat$x, y=logSSdensity(dHat$x+log(p), m0, 1, p))
   dHat<-hist(logs2-log(v0), breaks=40, plot=FALSE)

   return(list(Sigma=sigma0, m=m0, v=v0, converged=converged, 
               iter=iter, modS2=modS2, histLogS2=dHat,
               fittedDensityLogS2=d0, logs2=logs2, beta=beta0,
               knots=br, x=x))
}



## ---------------------------------------------------------------------------
estimateMVbeta<-function(y, x, Sigma, maxIter=200, epsilon=0.000001,
                         verbose=FALSE, nknots=10, nOut=2000, nIn=4000, 
                         iterInit=3, br=NULL)
{
    p<-ncol(y)
    n<-nrow(y)

    orderX<-order(x)

    ## Computing knots, inner and outer + boundary
    ## creates X matrix for the spline
    if(length(br)==0){
        br<-getKnots(x, nKnots=nknots, nOut=nOut, nIn=nIn)
        br<-sort(c(mean(br[1:2]), br))
    }
    mat1<-ns(x, knots=br[-c(1, length(br))],
             Boundary.knots=br[c(1, length(br))])
    mat2<-cbind(rep(1, n), mat1)
    meanMat2<-c(1, rep(1, n)%*%mat1/n)

    ty<-t(y)

    if(verbose) cat("Init m and spline for v\n")
    ## Setting start values for m0 and v0, using moments for log(s2hat)
    s2hat <-.C("getSS", A=as.double(ginv(Sigma)), y=as.double(ty), 
               n=as.integer(n), p=as.integer(p), res=double(length=n),
               PACKAGE="plw")$res/p

    ii<- (s2hat > 0 & is.finite(s2hat))
    s2hatMin<-min(s2hat[ii])
    s2hat[!ii]<-s2hatMin 
    m1<-lm(log(s2hat)~mat1)
    v0<-exp(m1$fitted.values)
    beta0<-m1$coefficients

    est1<-estimateMV(s2hat*p/(v0), p)
    v0<-v0*est1$v
    beta0[1]<-beta0[1]+log(est1$v)
    m0<-est1$m

    iter<-0
    converged<-FALSE

    if(verbose) cat("Iterating EM-steps\n iter      sigma      m, v\n")

    while(iter<maxIter & !converged){
        iter<-iter+1
        temp2 <-.C("getSS", A=as.double(ginv(Sigma)), y=as.double(ty),
                   n=as.integer(n), p=as.integer(p), res=double(length=n),
                   PACKAGE="plw")$res

        ii<- (temp2 > 0 & is.finite(temp2))
        s2hatMin<-min(temp2[ii])
        temp2[!ii]<-s2hatMin

        temp2 <- temp2+m0*v0

        ## Optimizing beta 
        delta<-(m0+p)/temp2
        aux<-specialOptim(mat2, delta, beta0, meanMat2, 
                          reltol=sqrt(.Machine$double.eps)/
                                 ifelse(iter<10, 1, 100))
        betaHat<-aux$par

        vHat<-exp(betaHat[1]+mat1%*%betaHat[-1])
    
        ## Optimizing m 
        temp1<-beta0[1]+mat1%*%beta0[-1]
        temp1<-mean(temp1)-mean(log(temp2/2))+digamma((m0+p)/2)-
               mean(exp(temp1)*(m0+p)/temp2)

        f<-function(m){
            return( (m/2*(log(m/2)+temp1) - log(gamma(m/2))) )
        }
        mHat<-0.2+(0:1)
        logLik<-f(mHat)
        while(diff(logLik)>0){
            mHat<-mHat+1
            logLik<-f(mHat)
        }
        if(mHat[1]>0.2){
            mHat[1]<-mHat[1]-1
        }

        ## Optimal m is between  mHat[1] and mHat[2]
        ## use golden ratio to optimize m
        mHat<-goldenRatio(mHat[1], mHat[2], f, epsilon/100)

        ## Difference from previous iteration
        diffMV<-c(m0, log(v0)/10)-c(mHat, log(vHat)/10)

        if(verbose){
            cat(sprintf("%5d      %1.2e   %1.2e", iter, 0, 
                max(abs(diffMV))), "\n")
        }
        absDiff<-max(abs(c(diffMV)))
   
        converged<- (absDiff<epsilon)
        ##Uppdating estimates
        m0<-mHat
        beta0<-betaHat
        v0<-vHat
    }

    if(verbose) {
        cat("Done!\nSummarizing\n")
    }
    temp<-ginv(Sigma)%*%ty 
    logs2 <- log(.C("DiagAtimesBv2", a=as.double(ty), b=as.double(temp), 
                    n=as.integer(n), m=as.integer(p), diag=double(length=n),
                    PACKAGE="plw")$diag/p  )

    meanLogs2<-log(v0*m0)-log(p)+digamma(p/2)-digamma(m0/2)

    modS2<-temp2/(m0+p)

    dHat<-density(logs2-log(v0))
    d0<-list(x=dHat$x, y=logSSdensity(dHat$x+log(p), m0, 1, p))
    dHat<-hist(logs2-log(v0), breaks=40, plot=FALSE)

    return(list(Sigma=Sigma, m=m0, v=v0, converged=converged, iter=iter,
                modS2=modS2, histLogS2=dHat, fittedDensityLogS2=d0, 
                logs2=logs2, beta=beta0, knots=br, x=x))
}

## ---------------------------------------------------------------------------
lmw <- function(x, design=rep(1, ncol(x)), contrast=matrix(1), meanX=NULL,
                maxIter=200, epsilon=0.000001, verbose=TRUE, nknots=10,
                nOut=2000, nIn=4000, knots=NULL, checkRegulation=TRUE)
{
    cl <- match.call()

    if(!plwCheckInput(x, design, contrast, checkRegulation=checkRegulation)){
        cat("lmw aborted", "\n"); return(NULL);
    }
    if(nrow(contrast)>1) {   contrast<-contrast[1, ]  }

    P<-plwGetTransformationMatrix(design, contrast)

    p<-ncol(x)
    if(is.null(meanX)){
       meanX<-x%*%rep(1/p, p)
    }

    geneNames<-rownames(x)
    arrayNames<-colnames(x)

    y<-x%*%P
    p<-ncol(y)

    fit1<-estimateSigmaMVbeta(y[, -p], meanX, maxIter=maxIter, 
                              epsilon=epsilon, verbose=verbose,
                              nknots=nknots, nOut=nOut, nIn=nIn, br=knots)

    fit2<-estimateSigma(y, fit1$m, fit1$v, maxIter=maxIter, epsilon=epsilon,
                        verbose=verbose)

    fit1$Sigma<-fit2$Sigma
    fit1$iter<-c(fit1$iter, fit2$iter) 
    fit1$converged<-all(fit1$iter<maxIter)

    D<-matrix((1:p)==p)*1
    gammaHat<- ginv(t(D)%*%ginv(fit2$Sigma)%*%D )%*%t(D)%*%ginv(fit2$Sigma)

    weights<-as.vector(gammaHat%*%t(P))

    gammaHat<- as.vector(y%*%t(gammaHat))
    varGammaHat<-ginv(t(D)%*%ginv(fit2$Sigma)%*%D )

    modT<-as.vector(gammaHat)/sqrt(as.vector(fit1$modS2)*varGammaHat)
    dfT<-fit1$m+ncol(P)-1
    modP<- 2-pt(abs(modT), dfT)*2

    names(fit1$modS2)<-names(gammaHat)<-geneNames
    names(modP)<-names(modT)<-geneNames
    rownames(P)<-names(weights)<-arrayNames
    colnames(P)<-c(paste("A", 1:(ncol(P)-1), sep=""), "B")

    fit1<-append(fit1, list(call=cl, t=modT, p.value=modP,
                            coefficients=gammaHat, P=P, dfT=dfT,
                            weights=weights))
  
    class(fit1)<-"plw"
    fit1
}

## ---------------------------------------------------------------------------
plwCheckInput<-function(x, design, contrast, checkRegulation=TRUE)
{
    check <- is.matrix(x)
    okInput <- check
    okReturn <- ifelse(check, "", "error: x must be a matrix\n")
    check <- is.matrix(design)
    okInput <- check & okInput
    okReturn <- cat(okReturn, ifelse(check, "",
                                       "error: design must be a matrix\n"), 
                     sep = "")
    check <- is.matrix(contrast)
    okInput <- check & okInput
    okReturn <- cat(okReturn, ifelse(check, "", 
                                       "error: contrast must be a matrix\n"), 
                     sep = "")
    if (okInput) {
        check <- nrow(design) == ncol(x) & ncol(design) == ncol(contrast)
        okInput <- check & okInput
        okReturn <- cat(okReturn, ifelse(check, "", 
            paste("error: wrong dimensions in input,",
                  " the must fulfill\n       nrow(design)",
                  "==ncol(x) and ncol(design)==ncol(contrast)\n", sep="")
                                          ), sep = "")
        check <- nrow(contrast) == 1
        okReturn <- cat(okReturn, ifelse(check | !okInput,"",
           "warning: contrast has multiple rows, will use first row only\n"), 
            sep = "")
        contrast <- matrix(contrast[1, ], nrow = 1)
    }
    if (okInput & checkRegulation) {
        temp <- t(contrast %*% ginv(t(design) %*% design) %*% t(design))
        temp <- x %*% temp
        check <- abs(median(temp)/mad(temp)) < 1
        okInput <- check & okInput 
        okReturn <- cat(okReturn, ifelse(check, "",
           paste("warning: most genes appears to be regulated.",
                                "Check your contrast matrix or use",
                                "checkRegulation=FALSE.\n")), sep = "")
    }
    if (length(okReturn) > 0)
        cat(okReturn)
    return(okInput)
}


## ---------------------------------------------------------------------------
plw <- function(x, design=rep(1, ncol(x)), contrast=matrix(1), 
                probenames=unlist(ifelse(class(x)=="AffyBatch",
                                          list(p=probeNames(x)),
                                          list(p=NULL))), 
                maxIter=200, epsilon=0.000001, verbose=TRUE, nknots=10,
                nOut=2000, nIn=4000, knots=NULL, checkRegulation=TRUE)
{
    cl <- match.call()
    if(class(x)=="AffyBatch"){
        pm1<-pm(x)
        if(!plwCheckInput(pm1, design, contrast, 
                          checkRegulation=FALSE)){
            cat("plw aborted.", "\n"); return(NULL);
        }
        if(length(probenames)!=nrow(pm1)){
            cat("The vector probenames is not of correct length.", "\n", 
                "plw aborted.", "\n");return(NULL);
        }
        if(verbose){
            cat("Background correcting and normalizing.", "\n")
        }
        pm1<-pm(bg.correct.rma(x, bgtype = 2))
        rows <- nrow(pm1)
        cols <- ncol(pm1)
        pm1<-log2(matrix(.C("qnorm_c", as.double(as.vector(pm1)),
                            as.integer(rows), as.integer(cols),
                            PACKAGE="affy")[[1]], rows, cols))
        rownames(pm1)<-probenames

        if(!plwCheckInput(pm1, design, contrast, 
                          checkRegulation=checkRegulation)){
            cat("plw aborted.", "\n"); return(NULL);
        }
    }
    else{
        pm1<-x
        if(!plwCheckInput(pm1, design, contrast,
                          checkRegulation=checkRegulation)){
            cat("plw aborted.", "\n")
        }
        if(length(probenames)!=nrow(pm1)){
            cat("The vector probenames is not of correct length.", "\n", 
                "plw aborted.", "\n")
        }
    }

    ii<-order(as.numeric(factor(probenames)))
    probenames<-probenames[ii]
    pm1<-pm1[ii, ]

    probeId<-as.numeric(factor(probenames))
    geneNames<-levels(factor(probenames))
    arrayNames<-colnames(pm1)

    if(nrow(contrast)>1) {   contrast<-contrast[1, ]  }

    P<-plwGetTransformationMatrix(design, contrast)

    p<-ncol(pm1)
    meanX<-pm1%*%rep(1/p, p)


    y<-pm1%*%P
    p<-ncol(y)
 
    fit1<-estimateSigmaMVbeta(y[, -p], meanX, maxIter=maxIter,
                              epsilon=epsilon, verbose=verbose, nknots=nknots, 
                              nOut=nOut, nIn=nIn, br=knots)

    fit2<-estimateSigma(y, fit1$m, fit1$v, maxIter=maxIter, epsilon=epsilon,
                        verbose=verbose)

    fit1$Sigma<-fit2$Sigma
    fit1$iter<-c(fit1$iter, fit2$iter) 
    fit1$converged<-all(fit1$iter<maxIter)

    D<-matrix((1:p)==p)*1
    gammaHat<- ginv(t(D)%*%ginv(fit2$Sigma)%*%D )%*%t(D)%*%ginv(fit2$Sigma)

    weights<-as.vector(gammaHat%*%t(P))

    gammaHat<- as.vector(y%*%t(gammaHat))
    varGammaHat<-ginv(t(D)%*%ginv(fit2$Sigma)%*%D )

    modT<-as.vector(gammaHat)/sqrt(as.vector(fit1$modS2)*varGammaHat)
    dfT<-fit1$m+ncol(P)-1
    modP<- 2-pt(abs(modT), dfT)*2

    medianT<-medianByIndex(modT, probeId)

    names(fit1$modS2)<-names(gammaHat)<-probenames
    names(modP)<-names(modT)<-probenames
    names(medianT)<-geneNames
 
    rownames(P)<-names(weights)<-arrayNames
    colnames(P)<-c(paste("A", 1:(ncol(P)-1), sep=""), "B")

    fit1<-append(fit1, list(call=cl, medianT=medianT, t=modT, 
                            p.value=modP, coefficients=gammaHat, 
                            P=P, dfT=dfT, weights=weights))
    class(fit1)<-"plw"
    fit1
}


## ---------------------------------------------------------------------------
## Print function for plw objects
print.plw<-function (x, ...)
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("Number of arrays     :", length(x$weights), "\n")
    if(!is.null(x$medianT)){
        cat("Number of probe-sets :", length(x$medianT), "\n")
        cat("Number of PM probes  :", length(x$t), "\n")
    }
    else{
        cat("Number of probe-sets :", length(x$t), "\n")
    }
    cat("Number of knots for v:", length(x$knots), "\n")
    cat("m parameter          :", round(x$m, 3), "\n")
    cat("Df for probe t-stat. :", round(x$dfT, 1), "\n")
    cat("Convergence status   :", x$converged, "\n")
    cat("Number of iterations :", x$iter, "\n")
}


## ---------------------------------------------------------------------------
## Plot showing fitted vs observed density for log(s2)
varHistPlot<-function(model, main="Histogram variance estimators", histCol=8,
                      densityCol=2, drawLegend=TRUE)
{
    temp<-c(length(model$fittedDensityLogS2$x),
            length(model$fittedDensityLogS2$y))
    if(class(model$histLogS2)=="histogram" &  
       class(model$fittedDensityLogS2)=="list" &
       sd(temp)==0 & min(temp)>1)
    {
        plot(1, xlim=range(model$histLogS2$breaks), 
             ylim=range(c(model$histLogS2$density, 
                          model$fittedDensityLogS2$y)),
             xlab="log(s2)", ylab="", main=main)
        lines(model$histLogS2, freq=FALSE, col=histCol)
        lines(model$fittedDensityLogS2, col=densityCol)
        if(drawLegend){
            legend("topleft", c("Observed density", "Fitted density"), 
                   inset=0.03, col=c(histCol, densityCol), lty=c(0, 1),
                   pch=c(15, -1), pt.cex=c(2.5, 1))
            legend("topleft", c("Observed density", "Fitted density"), 
                   inset=0.03, col=c(1, densityCol), lty=c(0, 1), 
                   pch=c(22, -1), pt.cex=c(2.5, 1))
        }
    }
    else{
        cat("Error: The model argument must be an",
            "output object from plw or lwp.", "\n")
    }
}

## ---------------------------------------------------------------------------
## Plot showing fitted curve for scale parameter  
scaleParameterPlot<-function(model, main="Scale parameter curve", col=1,
                             pch='.', lty=1, curveCol=2, knotsPch=19,
                             knotsCol=3)
{
    temp<-c(length(model$x), length(model$logs2), length(model$v), 
           length(model$knots))
    if(sd(temp[-4])==0 & min(temp[-4])>10 & temp[4]>1){
        plot(model$x, model$logs2, pch=pch, col=col, xlab="Mean intensity", 
             ylab="log(s2)", ylim=quantile(model$logs2, 
                                  c(5/length(model$logs2),
                                  1-5/length(model$logs2))),
             main=main)
        lines(model$x[order(model$x)], log(model$v[order(model$x)]), 
              col=curveCol, lty=lty)
        for(i in 1:length(model$knots)){
            ii<-which.min(abs(model$knots[i]-model$x))
            points(model$x[ii], log(model$v[ii]), pch=knotsPch, col=knotsCol)
        }
        legend("topright", inset=0.03, lty=c(lty, 0),
               col=c(curveCol, knotsCol), c("log(v)", "knots"), bg="white",
               pch=c(NA, knotsPch))
    }
    else{
        cat("Error: The model argument must be an",
            "output object from plw or lwp.", "\n")
    }
}


#    temp<-c(length(model$t), length(model$coefficients),
#            length(model$medianT))
#    if(sd(temp[-3])!=0 | min(temp)==0 ){
#        cat("Error: The model argument must be an",
#            "output object from plw.", "\n")
#        return(NULL)
#    }

## ---------------------------------------------------------------------------
## Compute summary tables for top ranking or seletced genes
topRankSummary<-function(model, nGenes=50, genesOfRank=1:nGenes,
                         genes=NULL)
{

    temp <- c(length(model$t), length(model$coefficients), length(model$medianT))
    if(sd(temp[-3]) != 0){
        cat("Error: The model argument must be an", "output object from plw or lmw.",
            "\n")
        return(NULL)
    }
    if (sd(temp[-3]) == 0 & min(temp) == 0) {
        return( topRankSummaryLMW(model, nGenes = nGenes, genesOfRank = genesOfRank, genes = genes) )
    }

    nGenes<-length(genesOfRank)
    ii<-order(-abs(model$medianT))[genesOfRank]
    medianT<-model$medianT[ii]

    if(!is.null(genes)){
        ii<-match(genes, names(model$medianT))
        if(sum(is.na(ii))==0){
         medianT<-model$medianT[ii]
         genesOfRank<-rank(-abs(model$medianT), ties.method="random")[ii]
         nGenes<-length(medianT)
        }
        else{
         cat("Argument genes has at least one element with no",
             "matching genes in data. Argument genes ignored.", "\n")
        } 
    }  

    ii<-names(model$t) %in% names(medianT)
    probeT<-model$t[ii]
    probeEst<-model$coefficients[ii]
    probeId<-match(names(probeT), names(medianT))
    ii<-order(probeId)
    probeT<-probeT[ii]
    probeEst<-probeEst[ii]
    probeId<-probeId[ii]
    medianEst<-medianByIndex(probeEst, probeId)
    maxNProbe<-max(table(probeId))

    summaryT<-matrix(NA, nGenes, maxNProbe+2)
    colnames(summaryT)<-c("Rank", "Median-t", 
                           paste("t-stat p", 1:maxNProbe, sep=""))
    rownames(summaryT)<-names(medianT)
    summaryLog2FC<-summaryT
    colnames(summaryLog2FC)<-c("Rank", "Median log2FC",
                                paste("log2FC p", 1:maxNProbe, sep=""))

    summaryT[, 2]<-medianT
    summaryLog2FC[, 2]<-medianEst

    summaryLog2FC[, 1]<-summaryT[, 1]<-genesOfRank

    noProbes<-table(probeId)
    for(i in 1:length(unique(noProbes))){
        k<-unique(noProbes)[i]
        rows<-noProbes==k
        ids<-probeId %in% (1:nGenes)[rows] 
        summaryT[rows, (1:k)+2]<-matrix(probeT[ids], ncol=k, byrow=TRUE)
        summaryLog2FC[rows, (1:k)+2]<-matrix(probeEst[ids], ncol=k,
                                              byrow= TRUE)
    }
    res<-list(t=summaryT, log2FC=summaryLog2FC) 
    class(res)<-"topRankSummary"
#    return(list(t=summaryT, log2FC=summaryLog2FC))
    return(res)
}


## ---------------------------------------------------------------------------
topRankSummaryLMW<-function (model, nGenes, genesOfRank, genes)
{
    temp <- c(length(model$t), length(model$coefficients))
    if(sd(temp) != 0){
        cat("Error: The model argument must be an", "output object from lmw.",
            "\n")
        return(NULL)
    }
    nGenes <- length(genesOfRank)
    ii <- order(-abs(model$t))[genesOfRank]
    tStat <- model$t[ii]
    geneEst<-model$coefficients[ii]
    if (!is.null(genes)) {
        ii <- match(genes, names(model$t))
        if (sum(is.na(ii)) == 0) {
            tStat <- model$t[ii]
            geneEst<-model$coefficients[ii]
            genesOfRank <- rank(-abs(model$t), ties.method = "random")[ii]
            nGenes <- length(tStat)
        }
        else {
            cat("Argument genes has at least one element with no",
                "matching genes in data. Argument genes ignored.",
                "\n")
        }
    }
    summaryT <- cbind(genesOfRank,tStat,geneEst)
    colnames(summaryT) <- c("Rank", "t-statistic","Estimate")
    rownames(summaryT) <- names(tStat)
    return(summaryT)
}



## ---------------------------------------------------------------------------
print.topRankSummary<-function (x, digits=3, ...){
    qT <-t(apply(x$t[, -(1:2)], 1, quantile, c(0.25, 0.75),
                 na.rm = TRUE, names = FALSE))
    qFC<-t(apply(x$log2FC[, -(1:2)], 1, quantile, c(0.25, 0.75),
                 na.rm = TRUE, names = FALSE))
    temp<-cbind(x$t[, 1:2], qT, x$log2FC[, 2])
    colnames(temp)<-c("Rank", "Median t", "  Q1-t", "  Q3-t", "Med. log2FC")
    print(temp, digits=digits)
}

## ---------------------------------------------------------------------------
## Plot summary of t-statistics for top ranking or seletced genes
plotSummaryT<-function(model, nGenes=50, genesOfRank=1:nGenes,
                       genes=NULL)
{
    data<-topRankSummary(model, nGenes, genesOfRank, genes)$t
    if(length(data)==0){
        return(NULL)
    }
    par0<-par()
    par(mar=c(4, max(nchar(rownames(data)))/2+1, 1, 1), las=1)
    medianT<-data[, 2]
    nGenes<-length(medianT)
    probeId<-rep(1:nrow(data), ncol(data)-2)
    probeT<-as.vector(data[, -(1:2)])
    ii<-is.na(probeT)
    probeId<-probeId[!ii]
    probeT<-probeT[!ii]
    plot(1, xlim=range(probeT), ylim=range(probeId), col=0, xlab="t-statistic",
         ylab="", axes=FALSE)
    abline(h=length(medianT):1, col=8, lty=2)
    abline(v=0, col=8, lty=2)
    points(probeT, nGenes-probeId+1, pch=19, cex=0.5, col=2)
    points(medianT, length(medianT):1, pch=19)
    axis(1)
    axis(2, length(medianT):1, names(medianT))
    box()
    par(mar=par0$mar, las=par0$las)
}

## ---------------------------------------------------------------------------
## Plot summary of log2FC for top ranking or seletced genes
plotSummaryLog2FC<-function(model, nGenes=50, genesOfRank=1:nGenes,
                            genes=NULL)
{
    data<-topRankSummary(model, nGenes, genesOfRank, genes)$log2FC
    if(length(data)==0){
        return(NULL)
    }
    par0<-par()
    par(mar=c(4, max(nchar(rownames(data)))/2+1, 1, 1), las=1)
    medianEst<-data[, 2]
    nGenes<-length(medianEst)
    probeId<-rep(1:nrow(data), ncol(data)-2)
    probeEst<-as.vector(data[, -(1:2)])
    ii<-is.na(probeEst)
    probeId<-probeId[!ii]
    probeEst<-probeEst[!ii]
    plot(1, xlim=range(probeEst), ylim=range(probeId), col=0,
         xlab="log2 FC", ylab="", axes=FALSE)
    abline(h=length(medianEst):1, col=8, lty=2)
    abline(v=0, col=8, lty=2)
    points(probeEst, nGenes-probeId+1, pch=19, cex=0.5, col=2)
    points(medianEst, length(medianEst):1, pch=19)
    axis(1)
    axis(2, length(medianEst):1, names(medianEst))
    box()
    par(mar=par0$mar, las=par0$las)
}


## ---------------------------------------------------------------------------
estimateSigma<-function(y, m, v, maxIter=100, epsilon=0.000001, 
                        verbose=FALSE)
{
    p<-ncol(y)
    n<-nrow(y)
    sigma0<-diag(rep(1, p))
    ty<-t(y)

    iter<-0
    converged<-FALSE

    if(verbose) cat("Iterating EM-steps\n iter      sigma      m, v\n")
    #######  Em algo
    while(iter<maxIter & !converged){
       iter<-iter+1
       temp2 <-.C("getSS", A=as.double(ginv(sigma0)), y=as.double(ty), 
                  n=as.integer(n), p=as.integer(p), res=double(length=n), 
                  PACKAGE="plw")$res +m*v

       #######  Sedan på sigma0
       sigmaHat<-covZeroMeanScaledData(y, (m+p)/(temp2))  

       #### Kollar på diff mot förra iterationen
       diffSigma<-as.vector(sigmaHat-sigma0)

       if(verbose){
           cat(sprintf("%5d      %1.2e", iter, max(abs(diffSigma))), "\n")
       }
       absDiff<-max(abs(c(diffSigma)))
  
       converged<- (absDiff<epsilon)

       ######Er sätter nu m0, v0, sigma0 med respektive hattar
       sigma0<-sigmaHat
   }

    return(list(Sigma=sigma0, iter=iter))
}

## ---------------------------------------------------------------------------
estimateSigmaMV<-function(y, maxIter=100, epsilon=0.000001, verbose=FALSE)
{
    p<-ncol(y)
    n<-nrow(y)
    sigma0<-diag(rep(1, p))
    ty<-t(y)
 
    ### Sätter start värden på m0 och v0, med hjälp av moment for log(s2hat)
    s2hat<-meanSdByRow(y)$sd^2
    ii<-s2hat>0
    temp<-var(log(s2hat[ii]))-trigamma(p/2)
    m0<-101
    if(temp>0.01){
        m0<-2*triGammaInverse(temp)
    }
    v0<-exp(mean(log(s2hat[ii]))-log(m0)+log(p)+digamma(m0/2)-digamma(p/2))

    iter<-0
    converged<-FALSE
 
    if(verbose) cat("Iterating EM-steps\n iter      sigma      m, v\n")
    #######  Em algo
    while(iter<maxIter & !converged){
        iter<-iter+1
        temp2 <-.C("getSS", A=as.double(ginv(sigma0)), y=as.double(ty), 
                   n=as.integer(n), p=as.integer(p), res=double(length=n), 
                   PACKAGE="plw")$res +m0*v0

        s1<-mean( log( temp2/2 )) - digamma((m0+p)/2) 
        s2<-mean( ( (m0+p) / ( temp2 )  ) )

        vHat<-1/s2
 
        ## Maximera map m, steg 1: öka m tills loglik minskar
        mHat<-0.2+(0:1)
        logLik<-mHat/2*(log(mHat)-log(2*s2)-s1-1)-log(gamma(mHat/2))
        while(diff(logLik)>0){
            mHat<-mHat+1
            logLik<-mHat/2*(log(mHat)-log(2*s2)-s1-1)-log(gamma(mHat/2))
        }
        if(mHat[1]>0.2){ mHat[1]<-mHat[1]-1 }

        f<-function(m) { m/2*(log(m)-log(2*s2)-s1-1)-log(gamma(m/2)) }
        mHat<-goldenRatio(mHat[1], mHat[2], f, epsilon/100)

        #######  Sedan på sigma0
        sigmaHat<-covZeroMeanScaledData(y, (m0+p)/(temp2))  

        ##### Skalar v och sigma
        vHat<-vHat*mean(diag(sigmaHat))
        sigmaHat<-sigmaHat/mean(diag(sigmaHat))

        #### Kollar på diff mot förra iterationen
        diffSigma<-as.vector(sigmaHat-sigma0)
        diffMV<-c(m0, v0)-c(mHat, vHat)

        if(verbose){
            cat(sprintf("%5d      %1.2e   %1.2e", iter, max(abs(diffSigma)),
                        max(abs(diffMV))), "\n")
        }
        absDiff<-max(abs(c(diffSigma, diffMV)))
   
        converged<- (absDiff<epsilon)

        ######Er sätter nu m0, v0, sigma0 med respektive hattar
        m0<-mHat
        v0<-vHat
        sigma0<-sigmaHat
    }
    if(verbose) cat("Done!\nSummarizing\n")

    modS2<-(temp2)/(m0+p)

    temp2<-(temp2-m0*v0)/p
    temp2<-temp2[temp2>0]
    dHat<-density(log(temp2))
    d0<-list(x=dHat$x, y=logSSdensity(dHat$x+log(p), m0, v0, p))
    dHat<-hist(log(temp2), breaks=60, plot=FALSE)

    return(list(Sigma=sigma0, m=m0, v=v0, converged=converged, iter=iter, 
                modS2=modS2, histLogS2=dHat, fittedDensityLogS2=d0))
}

## ---------------------------------------------------------------------------
plwGetTransformationMatrix<-function(design, contrast){
    p<-nrow(design)
    temp<-design%*%ginv(t(design)%*%design)
    A<-diag(rep(1, p))-temp%*%t(design)
    A<-t(ReduzeToFullRank(A))
    B<-temp%*%t(contrast)
    P<-cbind(A, B)
    return(P)
}


## ---------------------------------------------------------------------------
specialOptim<-function(mat2, delta, beta, meanMat2, maxit=1000, trace=0, 
                       reltol=sqrt(.Machine$double.eps))
{
    n<-nrow(mat2)
    n0<-ncol(mat2)
    aux<-.C("SpecialOptim", n0=as.integer(n0), par=as.double(beta), 
            value=double(1), mat2=as.double(mat2), delta=as.double(delta), 
            meanMat2=as.double(meanMat2), nn=as.integer(n), 
            maxit=as.integer(maxit), trace=as.integer(trace), 
            abstol=as.double(-100000), reltol=as.double(reltol), 
            fncount=as.integer(1), grcount=as.integer(1), 
            convergence=as.integer(1), PACKAGE="plw")
    return(list(par=aux$par, value=aux$value, 
                counts=c(aux$fncount, aux$grcount), 
                convergence=aux$convergence, message=NULL)) 
}
