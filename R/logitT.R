## The logit transformation as defined in Lemon et al (2003)
logitTTransform<-function(pm){
    for(i in 1:ncol(pm)){
        range.pm<-range(pm[, i])
        A<-range.pm[2]+0.001*diff(range.pm)
        N<-range.pm[1]-0.001*diff(range.pm)
        temp<-log(pm[, i]-N)-log(A-pm[, i])
        m<-mean(temp)
        std<-sd(temp)
        pm[, i]<-(temp-m)/std
    }
    return(pm)
}

logitTStat<-function(affy.batch, group){
    pm.transformed<-logitTTransform(pm(affy.batch))
    probe.t<-studenttTTest(pm.transformed, group)
    probe.names<-factor(probeNames(affy.batch))
    res<- medianByIndex(probe.t, as.numeric(probe.names))
    names(res)<-levels(probe.names)
    return(res)
}

studenttTTest<-function (x, group){
    group1<-group==1
    n1<-sum(group1)
    n2<-sum(!group1)
    if(n1<2 | n2<2){
        stop("increase sample size!")
    }
    m1<-meanSdByRow(x[, group1])
    m2<-meanSdByRow(x[, !group1])
    v.pooled = ((n1-1)*m1$sd^2 + (n2-1)*m2$sd^2)/(n1+n2 - 2)
    diff = m1$mean - m2$mean
    std = sqrt((1/n1 + 1/n2) * v.pooled)
    ii<-std>0
    t.stat<-diff
    t.stat[ii]<-diff[ii]/std[ii]
    max.t<-max(abs(t.stat[ii]))
    t.stat[!ii] = sign(diff[!ii])*(max.t+1)
    return(t.stat)
}
