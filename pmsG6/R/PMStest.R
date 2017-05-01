#' Methods for partially matched samples
#'
#' This function allows you to use weighted Z-test or modified test statistics approach to deal with partially matched samples.
#' @param x a (non-empty) numeric vector of data values with missing values.
#' @param y a (non-empty) numeric vector of data values with different missing values of \code{x}.
#' @param alternative a character string specifying the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.
#' @param method method which is implemented for partially matched samples, \code{"weighted.z"} is weighted Z-test combination (default), \code{"modified.t"} is modified t-statistic of Kim et al, \code{"corrected.z"} is corrected Z-test of Looney and Jones, \code{"MLE.homo"} is MLE based test of Ekbohm under homoscedasticity, \code{"MLE.hetero"} is MLE based test of Lin and Stivers under heteroscedasticity.
#' @param data an optional matrix or data frame containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula)}.
#' @details 
#' \code{x} and \code{y} both must have at least two numeric values and two missing values. The missing values of \code{y} must be in different position with the missing values of \code{x}. If not, there is no need to use this function, you can simply clean your data and try regular \code{"t.test"}.
#' 
#' \code{alternative = "greater"} is the alternative that x has a larger mean than y.
#' 
#' Partially matched samples can be viewd as data genereated from two experiemntal designs: (i) n1 mathced 
#' pairs, and (ii) independent groups with n2 and n3 per group, where both designs intend to estimate the 
#' same parameter. 
#' 
#' Here, we implement modified test statistics approach or weighted Z-test for partially mateched samples.
#' According to the False Discovery Rate (FDR), False Non-discovery Rate (FNR), and receiver operating characteristic curve (ROC),
#' the results suggest that weighting by the square root of sample size (weighted Z-test) performs well in various scenarios.
#' 
#' This function does not calculate the confidence interval, because the confidence interval has no direct relationship with mean difference in all of the five methods.
#' 
#' See \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3717400/} for more details.
#' @return A list with class \code{"htest"} containing the following components:
#' @return \item{statistic}{the value of the test statistic. Either Z-statistic or t-statistic.}
#' @return \item{parameter}{the degrees of freedom for the t-statistic.}
#' @return \item{p.value}{the p-value for the test.}
#' @return \item{estimate}{the estimated difference in means.}
#' @return \item{null.value}{the specified hypothesized value of the mean difference.}
#' @return \item{alternative}{a character string describing the alternative hypothesis.}
#' @return \item{method}{a character string indicating what type of test was performed.}
#' @return \item{data.name}{a character string giving the names of the data.}
#' @author Weifan Han <\email{woailm1993@@yahoo.com}>
#' @author Fan Zhang <\email{frank08081993@@gmail.com}>
#' @author Peishan Han <\email{peishan.han@@stonybrook.edu}>
#' @author Yudong Li <\email{yudong.li@@stonybrook.edu}>
#' @author Elizabeth Pemberton <\email{pembertonei@gmail.com}>
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

PMS.test=function(x,y,alternative=c("two.sided","less","greater"),method=c("weighted.z","modified.t","corrected.z","MLE.homo","MLE.hetero"),data=NULL){
  alternative=alternative[1]
  method=method[1]
  if(!is.null(data)){
    attach(data)
  }
  x_name=deparse(substitute(x))
  y_name=deparse(substitute(y))
  data_name=paste(paste(x_name,collapse="\n"),"and",paste(y_name,collapse="\n"))
  null_value=c(0)
  names(null_value)=c("difference in means")
  na_x=which(!is.na(x))
  na_y=which(!is.na(y))
  esti=c(mean(x[na_x]),mean(y[na_y]))
  names(esti)=c(paste("mean of ",x_name),paste("mean of ",y_name))
  n_1=intersect(na_x,na_y)
  n_2=setdiff(na_x,na_y)
  n_3=setdiff(na_y,na_x)
  n1=length(n_1)
  n2=length(n_2)
  n3=length(n_3)
  if(n2<=1|n3<=1){
    stop("Your data does not have enough missing values, or your x and y have missing values in identical rows, please try to remove missing values and then use regular t.test.\n")
  }
  if(method=="modified.t"){
    method_output="Modified t-statistic of Kim et al"
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
    names(ts)="Z"
    output_list=list(statistic=ts,p.value=pvalue,estimate=esti,null.value=null_value,alternative=alternative,method=method_output,data.name=data_name)
  }
  else if(method=="corrected.z"){
    method_output="Corrected Z-test of Looney and Jones"
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
    names(ts)="Z"
    output_list=list(statistic=ts,p.value=pvalue,estimate=esti,null.value=null_value,alternative=alternative,method=method_output,data.name=data_name)
  }
  else if(method=="MLE.homo"){
    method_output="MLE Based Test of Ekbohm under Homoscedasticity"
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
    param=n1
    names(param)="df"
    names(ts)="t"
    output_list=list(statistic=ts,parameter=param,p.value=pvalue,estimate=esti,null.value=null_value,alternative=alternative,method=method_output,data.name=data_name)
  } 
  else if(method=="MLE.hetero"){
    method_output="MLE Based Test of Lin and Stivers under Heteroscedasticity"
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
    param=n1
    names(param)="df"
    names(ts)="t"
    output_list=list(statistic=ts,parameter=param,p.value=pvalue,estimate=esti,null.value=null_value,alternative=alternative,method=method_output,data.name=data_name)
  }
  else if(method=="weighted.z"){
    method_output="Weighted Z-test Combination"
    n23=n2+n3
    if(var.test(x[n_2],y[n_3])$p.value<0.05){
      partial_var=0
    }
    else{
      partial_var=1
    }
    if(alternative=="two.sided"){
      p1=t.test(x[n_1],y[n_1],alternative="greater",paired=T)$p.value
      z1=qnorm(1-p1)
      p2=t.test(x[n_2],y[n_3],alternative="greater",var.equal=partial_var)$p.value
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
      p1=t.test(x[n_1],y[n_1],alternative="greater",paired=T)$p.value
      z1=qnorm(1-p1)
      p2=t.test(x[n_2],y[n_3],alternative="greater",var.equal=partial_var)$p.value
      z2=qnorm(1-p2)
      pvalue=1-pnorm((sqrt(n1)*z1+sqrt(n23)*z2)/sqrt(n1+n23))
      ts=qnorm(1-pvalue)
    }
    else if(alternative=="less"){
      p1=t.test(x[n_1],y[n_1],alternative="less",paired=T)$p.value
      z1=qnorm(1-p1)
      p2=t.test(x[n_2],y[n_3],alternative="less",var.equal=partial_var)$p.value
      z2=qnorm(1-p2)
      pvalue=1-pnorm((sqrt(n1)*z1+sqrt(n23)*z2)/sqrt(n1+n23))
      ts=qnorm(pvalue)
    }
    names(ts)="Z"
    output_list=list(statistic=ts,p.value=pvalue,estimate=esti,null.value=null_value,alternative=alternative,method=method_output,data.name=data_name)
  }
  names(pvalue)="p-value"
  if(method!="weighted.z"){
    warning("We recommend weighted z-test combination method, since our experience show that weighting by the square root of sample size performs well in various scenarios.")
  }
  whole_var=var.test(x[na_x],y[na_y])$p.value
  if(whole_var<0.05&method!="MLE.hetero"){
    warning("You can try MLE based test of Lin and Stivers under heteroscedasticity, since your data shows heteroscedasticity.")
  }
  if(whole_var>=0.05&method!="MLE.homo"){
    warning("You can try MLE based test of Ekbohm under homoscedasticity, since your data shows homoscedasticity.")
  }
  if(!is.null(data)){
    detach(data)
  }
  structure(output_list,class="htest")
}