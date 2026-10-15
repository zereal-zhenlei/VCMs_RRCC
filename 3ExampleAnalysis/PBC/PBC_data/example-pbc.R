   library(SparseM)
   library(quantreg)

   # analysis the pbc data with censored quantile regression 
   data1<-read.table("D:\\11.txt",header=F)  # read the data 
  
   
   data=data1[(1:312),]
   
  # data=as.numeric(data)
   data=as.matrix(data,ncol=312,20)
   y=as.numeric(data[,2])
   y=log(y)
   age=as.numeric(data[,5])
   alb=as.numeric(data[,13])
   alk=as.numeric(data[,15])
   asc=as.numeric(data[,7])
   bil=as.numeric(data[,11])
   ede=as.numeric(data[,10])
   hep=as.numeric(data[,8])
   pro=as.numeric(data[,19])
   sex=as.numeric(data[,6])
   spi=as.numeric(data[,9])
   # Take the natural logarithmic scale four variables
   alb1=as.numeric(alb[1:312])
   alk1=as.numeric(alk[1:312])
   bil1=as.numeric(bil[1:312])
   pro1=as.numeric(pro[1:312])

   alb1=log(alb1)
   alk1=log(alk1)
   bil1=log(bil1)
   pro1=log(pro1)
   x=cbind(age,alb1,alk1,asc,bil1,ede,hep,pro,sex,spi)
   x=scale(x)
   y=scale(y)
   p=dim(x)[2]
   n=dim(x)[1]
   lam=diag(dim(x)[2])
   y0=rep(0,p)
   betastart=rq(y~x-1,0.5)$coef
   print(betastart)
   x0=lam*log(n)/abs(betastart)
   yp=c(y,y0)
   xp=rbind(x,x0)
   betaf=rq(yp~xp-1,0.5)$coef
   betaf[abs(betaf)<1e-06]=0
   print(betaf)

