#' Methods for partially matched samples
#'
#' This function allows you to use weighted z-test or modified test statistics approach to deal with partially matched samples.
#' @param x a (non-empty) numeric vector of data values in tumor samples
#' @param y a (non-empty) numeric vector of data values in normal samples
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param method method which is implemented for partially matched samples, "modified.t" is modified t-statistic of Kim, "corrected.z" is corrected Z-test of Looney and Jones, "MLE.homo" is MLE based test of Ekbohm under homoscedasticity, "MLE.hetero" is MLE based test of Lin and Stivers under heteroscedasticity, "weighted.z" is weighted Z-test combination (default).
#' @param data an optional matrix or data frame containing the variables in the formula formula. By default the variables are taken from environment(formula).
#' @details 
#' Partially matched samples can be viewd as data genereated from two experiemntal designs: (i) n1 mathced 
#' pairs, and (ii) independent groups with n2 and n3 per group, where both designs intend to estimate the 
#' same parameter. 
#' 
#' Here, we implement modified test statistics approach or weighted Z-test for partially mateched samples.
#' According to the False Discovery Rate (FDR), False Non-discovery Rate (FNR), and receiver operating characteristic curve (ROC),
#' the results suggest that weighting by the square root of sample size under the heteroscedascticity (weigthed z-test) is robust and performs well among these five methods. 
#' 
#' See \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3717400/} for more details.
#' 
#' To see the 'description file', please click the 'Index' at the bottom of this page.
#' @return \code{ts}: the value of the test statistic
#' @return \code{pvalue}: the p-value for the test
#' @export
#' @examples
#' ## Example 1: generate random normal data
#' x = rnorm(100,0,1)
#' y = rnorm(100,0.5,1)
#' x.with.na = c(x[1:20],NA,x[21:50],NA,NA,NA,x[51:100])
#' y.with.na = c(NA,NA,y[1:80],NA,y[81:100],NA)
#' PMS.test(x=x.with.na,y=y.with.na,method="weighted.z",alternative="two.sided")
#' 
#' ## Example 2: a dataset with missing values
#' library(datasets)
#' data<-matrix(data=presidents,30,4,byrow = FALSE)
#' data1<-data.frame(data)
#' x<-data1$X1
#' y<-data1$X4
#' PMS.test(x=x,y=y,method="modified.t",alternative="greater")

PMS.test=function(x,y,alternative="two.sided",method="weighted.z",data=NULL){
  if(!is.null(data)){
    attach(data)
  }
  na_x=which(!is.na(x)) #position for non-null tumor samples
  na_y=which(!is.na(y)) #position for non-null normal samples
  n_1=intersect(na_x,na_y) #position for paired samples
  n_2=setdiff(na_x,na_y) #position for tumor unmatched samples
  n_3=setdiff(na_y,na_x) #position for normal unmatched sample
  n1=length(n_1)
  n2=length(n_2)
  n3=length(n_3)
  # modified t-statistic of Kim
  if(method=="modified.t"){
    D=x[n_1]-y[n_1]
    nH=2/(1/n2+1/n3)
    ts=(n1*mean(D)+nH*(mean(x[n_2])-mean(y[n_3])))/sqrt(n1*var(D)+nH^2*(var(y[n_3])/n3+var(x[n_2])/n2))
    if(alternative=="two.sided"){
      pvalue=2*pnorm(abs(ts),lower.tail=F)
    }
    else if(alternative=="greater"){
      pvalue=pnorm(ts,lower.tail=F)
    }
    else if(alternative=="less"){
      pvalue=pnorm(ts,lower.tail=T)
    }
    cat('\n')
    cat(sprintf('     %s',"Modified t-statistic of Kim"))
    cat('\n\n')
  }
  # corrected Z-test of Looney and Jones
  else if(method=="corrected.z"){
    n12=n1+n2
    n13=n1+n3
    ts=(mean(x[c(n_1,n_2)])-mean(y[c(n_1,n_3)]))/sqrt(var(x[c(n_1,n_2)])/n12+var(y[c(n_1,n_3)])/n13-2*n1*cov(x[n_1],y[n_1])/(n12*n13))
    if(alternative=="two.sided"){
      pvalue=2*pnorm(abs(ts),lower.tail=F)
    }
    else if(alternative=="greater"){
      pvalue=pnorm(ts,lower.tail=F)
    }
    else if(alternative=="less"){
      pvalue=pnorm(ts,lower.tail=T)
    }
    cat('\n')
    cat(sprintf('      %s',"Corrected Z-test of Looney and Jonesm"))
    cat('\n\n')
  }
  # MLE based test of Ekbohm under homoscedasticity
  else if(method=="MLE.homo"){
    ST1square=var(x[n_1])
    SN1square=var(y[n_1])
    r=cov(x[n_1],y[n_1])/sqrt(ST1square*SN1square)
    rsquare=r^2
    r2p1=1+rsquare
    n12=n1+n2
    n13=n1+n3
    n23=n2+n3
    n2m3=n2*n3
    n1m1=n1-1
    n2m1=n2-1
    n3m1=n3-1
    fgd=n12*n13-n2m3*rsquare
    Tbar=mean(x[n_2])
    Nbar=mean(y[n_3])
    ts=(n1*(n13+n2*r)/fgd*(mean(x[n_1])-Tbar)-n1*(n12+n3*r)/fgd*(mean(y[n_1])-Nbar)+Tbar-Nbar)/sqrt((var(x[n_1])*n1m1+var(y[n_1])*n1m1+r2p1*(var(x[n_2])*n2m1+var(y[n_3])*n3m1))/(2*n1m1+r2p1*(n2m1+n3m1))*(2*n1*(1-r)+n23*(1-rsquare))/(n12*n13-n2m3*rsquare))
    if(alternative=="two.sided"){
      pvalue=2*pt(abs(ts),df=n1,lower.tail=F)
    }
    else if(alternative=="greater"){
      pvalue=pt(ts,df=n1,lower.tail=F)
    }
    else if(alternative=="less"){
      pvalue=pt(ts,df=n1,lower.tail=T)
    }
    cat('\n')
    cat(sprintf('      %s',"MLE based test of Ekbohm under homoscedasticity"))
    cat('\n\n')
  } 
  # MLE based test of Lin and Stivers under heteroscedasticity
  else if(method=="MLE.hetero"){
    ST1square=var(x[n_1])
    SN1square=var(y[n_1])
    STN1=cov(x[n_1],y[n_1])
    r=STN1/sqrt(ST1square*SN1square)
    rsquare=r^2
    n12=n1+n2
    n13=n1+n3
    n1m1=n1-1
    fgd=n12*n13-n2*n3*rsquare
    Tbar=mean(x[n_2])
    Nbar=mean(y[n_3])
    f=n1*(n13+n2*STN1/ST1square)/fgd
    g=n1*(n12+n3*STN1/SN1square)/fgd
    ts=(f*(mean(x[n_1])-Tbar)-g*(mean(y[n_1])-Nbar)+Tbar-Nbar)/sqrt(((f^2/n1+(1-f)^2/n2)*ST1square*n1m1+(g^2/n1+(1-g)^2/n3)*SN1square*n1m1-2*f*g*STN1*n1m1/n1)/n1m1)
    if(alternative=="two.sided"){
      pvalue=2*pt(abs(ts),df=n1,lower.tail=F)
    }
    else if(alternative=="greater"){
      pvalue=pt(ts,df=n1,lower.tail=F)
    }
    else if(alternative=="less"){
      pvalue=pt(ts,df=n1,lower.tail=T)
    }
    cat('\n')
    cat(sprintf('      %s',"MLE based test of Lin and Stivers under heteroscedasticity"))
    cat('\n\n')
  }
  # weighted Z-test combination
  else if(method=="weighted.z"){
    D=x[n_1]-y[n_1]
    n23=n2+n3
    if(alternative=="two.sided"){
      p1=pnorm(mean(D)/sd(D)*sqrt(n1),lower.tail=F)
      z1=qnorm(1-p1)
      p2=pnorm((mean(x[n_2])-mean(y[n_3]))/sqrt(var(x[n_2])/n2+var(y[n_3])/n3),lower.tail=F)
      z2=qnorm(1-p2)
      pc=1-pnorm((sqrt(n1)*z1+sqrt(n23)*z2)/sqrt(n1+n23))
      ts=qnorm(1-pc)
      if(pc<0.5){
        pvalue=2*pc
      }
      else{
        pvalue=2*(1-pc)
      }
    }
    else if(alternative=="greater"){
      p1=pnorm(mean(D)/sd(D)*sqrt(n1),lower.tail=F)
      z1=qnorm(1-p1)
      p2=pnorm((mean(x[n_2])-mean(y[n_3]))/sqrt(var(x[n_2])/n2+var(y[n_3])/n3),lower.tail=F)
      z2=qnorm(1-p2)
      pvalue=1-pnorm((sqrt(n1)*z1+sqrt(n23)*z2)/sqrt(n1+n23))
      ts=qnorm(1-pvalue)
    }
    else if(alternative=="less"){
      p1=pnorm(mean(D)/sd(D)*sqrt(n1),lower.tail=T)
      z1=qnorm(1-p1)
      p2=pnorm((mean(x[n_2])-mean(y[n_3]))/sqrt(var(x[n_2])/n2+var(y[n_3])/n3),lower.tail=T)
      z2=qnorm(1-p2)
      pvalue=1-pnorm((sqrt(n1)*z1+sqrt(n23)*z2)/sqrt(n1+n23))
      ts=qnorm(pvalue)
    }
    cat('\n')
    cat(sprintf('     %s',"Weighted Z-test Combination"))
    cat('\n\n')
  }
  cat(sprintf('%s: %f','Test-Statistic',ts))
  cat('\n')
  cat(sprintf('%s: %f','P-value',pvalue))    
  cat('\n\n')
  cat(sprintf('%s',paste('## We recommend weighted z-test combination method, since results suggest that ',
                        'weighting by the square root of sample size under the heteroscedascticity ', 
                         'assumption is robust and performs well in various scenarios.',sep='')))
  if(!is.null(data)){
    detach(data)
  }
}


