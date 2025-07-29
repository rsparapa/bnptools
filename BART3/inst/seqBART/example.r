############################################################
############################################################
##File to illustrate use of R functions
############################################################
############################################################
library(MASS)
library(Rlab)
library(msm)

# CHANGE DIRECTORY --------------------------------
source('DIRECTORY/load-pbart.R')
#--------------------------------------------------
#simulated data
set.seed(21)
  n=150
  pp=4
  xx<-matrix(NA,n,pp)
  varm1<-matrix(0.8,pp,pp)
  diag(varm1)<-1
  xx<-mvrnorm(n,rep(0,pp),varm1)
  xx[xx[,3]>0,3]<-1
  xx[xx[,3]<=0,3]<-0
  yy<-rowMeans(xx)+rnorm(n)
  for(i in 2:4){
    ran<-runif(n)
    xx[ran<0.1,i]=NA
  }


# run the sequential bart program
datatype=c(0,0,1,0)
type=1

impute<-BartMI(bartmpicstem,xx=xx,yy=yy,datatype=datatype,
               type=type,numskip=199,burn=1000)

# five imputed datasets
imputed1<-impute$imputed1
imputed2<-impute$imputed2
imputed3<-impute$imputed3
imputed4<-impute$imputed4
imputed5<-impute$imputed5

