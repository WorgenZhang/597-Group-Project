
PMS.test=function(x,y,alternative="two.sided",method="5",conf.level=0.95,data=NULL){
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
  if(method=="1"){
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
  }
  # corrected Z-test of Looney and Jones
  else if(method=="2"){
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
  }
  # MLE based test of Ekbohm under homoscedasticity
  else if(method=="3"){
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
  }
  # MLE based test of Lin and Stivers under heteroscedasticity
  else if(method=="4"){
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
  }
  # weighted Z-test combination
  else if(method=="5"){
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
  }
  print(c(ts,pvalue))
  if(!is.null(data)){
    detach(data)
  }
}
x=c(0,1,2,3,NA,NA,NA,9)
y=c(3,9,NA,NA,4,5,NA,NA)
PMS.test(x=x,y=y,method="5",alternative="two.sided")

#for method 4
set.seed(123)
rho <- 0.5
s1 = s2 <- sqrt(2)
mu1 <- 2
mu2 <- -2
U <- rnorm(50)
V <- rnorm(50)
X <- mu1 + s1*U
Y <- mu2 + s2*(rho*U + sqrt(1 - rho^2)*V)
data1=data.frame(X,Y)
data1[1,1]=NA
data1[2,2]=NA
data1[3,1]=NA
data1[4,2]=NA
