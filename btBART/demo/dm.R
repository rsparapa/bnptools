## original file is/was ...
## /home/btuyishimire/allfiles/RESEARCH/NEWBART/TEST/A1CDATA/a1cdata_full.R
library(btBART)
library(BART3)
library(BayesTree)
library(lmeVarComp)
source("lpml.R")
##source("/home/btuyishimire/allfiles/RESEARCH/Rfunctions/lpml.R")

#A1C DATA FROM THE PHI DIABETES STUDY
##path <- "/home/btuyishimire/Documents/PHI/TARIMA/manuscript3/illustrative_example/"
##data <- read.csv(paste(path,"monthly_data12.csv",sep=""),header=TRUE)
data <- read.csv("monthly_data12.csv",header=TRUE)
data <- data[!is.na(data$A1C6),]

#-------  1. Check additivity of health and demographic factors ------------#

data_full <- data
n_full <- nrow(data_full)
# data_notrt <- data[data$newcurr_insulin=="No or Unknown" & data$newcurr_met == "No or Unknown",]
# n_notrt <- nrow(data_notrt)

demograph <- c( "age","Female","race_cat","new_insurance","marital_status_cat", "religion_cat")
health <- c("A1C0","DEPRESS","HYPERTSN", "bmicat","hdlcat","ldlcat","systolicbpcat","diastolicbpcat")
trt <- c("newcurr_insulin","neworder_insulin","newdisp_insulin","newcurr_met", "neworder_met","newdisp_met" )

ndpost <- 1000; burn <- 100

## A1C
#----- First do analysis using data without treatment 

xtrain1 <- makeind(data_full[,demograph])
xtrain2 <- makeind(data_full[,health])
xtrain_all <- cbind(xtrain1,xtrain2)
ytrain <- data_full$A1C6

post_1b <- wbart( x.train  = xtrain_all,
                  y.train  = ytrain,
                  ndpost   = ndpost,
                  nskip    =  burn)
 
mse1 <- mean((ytrain-post_1b$yhat.train.mean)^2)

sigs1 <- post_1b$sigma[(burn+1):(burn+ndpost)]
mus1 <- post_1b$yhat.train
out1 <- lpml(ytrain,mus1,sigs1)

#TWO BARTS

post_2b <- twbarts( x.train1 = xtrain1,
                    x.train2 = xtrain2,
                    y.train  = ytrain,   
                    ndpost   = ndpost, 
                    nskip    =  burn)
 
sigs2 <- post_2b$sigma[(burn+1):(burn+ndpost)]  
mus2  <-  post_2b$yhat.train
out2  <-  lpml(ytrain,mus2,sigs2) 
mse2 <- mean((ytrain-post_2b$yhat.train.mean)^2)

pbf_orig   <- exp(out1$LPML-out2$LPML)
pbf_orig

## log A1C
#----- First do analysis using data without treatment 

xtrain1 <- makeind(data_full[,demograph])
xtrain2 <- makeind(data_full[,health])
xtrain_all <- cbind(xtrain1,xtrain2)
ytrain <- data_full$logA1C6

post_1b <- wbart( x.train  = xtrain_all,
                  y.train  = ytrain,
                  ndpost   = ndpost,
                  nskip    =  burn)
 
mse1 <- mean((ytrain-post_1b$yhat.train.mean)^2)

sigs1 <- post_1b$sigma[(burn+1):(burn+ndpost)]
mus1 <- post_1b$yhat.train
out1 <- lpml(ytrain,mus1,sigs1)

#TWO BARTS

post_2b <- twbarts( x.train1 = xtrain1,
                    x.train2 = xtrain2,
                    y.train  = ytrain,   
                    ndpost   = ndpost, 
                    nskip    =  burn)
 
sigs2 <- post_2b$sigma[(burn+1):(burn+ndpost)]  
mus2  <-  post_2b$yhat.train
out2  <-  lpml(ytrain,mus2,sigs2) 
mse2 <- mean((ytrain-post_2b$yhat.train.mean)^2)

pbf_log   <- exp(out1$LPML-out2$LPML)
pbf_log

## A1C
ytrain <- data_full$A1C6

post_1b <- wbart( x.train  =  data_full[,c("age","A1C0")] ,
                  y.train  = ytrain,
                  ndpost   = ndpost,
                  nskip    =  burn)
 
sigs1 <- post_1b$sigma[(burn+1):(burn+ndpost)]
mus1 <- post_1b$yhat.train
out1 <- lpml(ytrain,mus1,sigs1)

#TWO BARTS

post_2b <- twbarts( x.train1 = data_full[,c("age")],
                    x.train2 = data_full[,c("A1C0")],
                    y.train  = ytrain,   
                    ndpost   = ndpost, 
                    nskip    =  burn)  

sigs2 <- post_2b$sigma[(burn+1):(burn+ndpost)]   
mus2  <-  post_2b$yhat.train 
out2  <-  lpml(ytrain,mus2,sigs2) 
mse2 <- mean((ytrain-post_2b$yhat.train.mean )^2)

pbf_ageA1C0_ori   <- exp(out1$LPML-out2$LPML)

#Normalize data
agen<- (data_full$age-min(data_full$age))/(max(data_full$age)-min(data_full$age))
A1C0n<- (data_full$A1C0-min(data_full$A1C0))/(max(data_full$A1C0)-min(data_full$A1C0))

xx <- cbind(agen,A1C0n)

addtest_ageA1C0_ori <- test.additivity(xx, ytrain)
addtest_ageA1C0_ori 
pbf_ageA1C0_ori 

## log A1C

ytrain <- data_full$logA1C6

post_1b <- wbart( x.train  =  data_full[,c("age","A1C0")] ,
                  y.train  = ytrain,
                  ndpost   = ndpost,
                  nskip    =  burn)
 
sigs1 <- post_1b$sigma[(burn+1):(burn+ndpost)]
mus1 <- post_1b$yhat.train
out1 <- lpml(ytrain,mus1,sigs1)

#TWO BARTS

post_2b <- twbarts( x.train1  = data_full[,c("age")],
                    x.train2  =  data_full[,c("A1C0")],
                    y.train  = ytrain,   
                    ndpost   = ndpost, 
                    nskip    =  burn)  

sigs2 <- post_2b$sigma[(burn+1):(burn+ndpost)]   
mus2  <-  post_2b$yhat.train 
out2  <-  lpml(ytrain,mus2,sigs2) 
mse2 <- mean((ytrain-post_2b$yhat.train.mean )^2)

pbf_ageA1C0_log   <- exp(out1$LPML-out2$LPML)


#Normalize data
agen<- (data_full$age-min(data_full$age))/(max(data_full$age)-min(data_full$age))
A1C0n<- (data_full$A1C0-min(data_full$A1C0))/(max(data_full$A1C0)-min(data_full$A1C0))

xx <- cbind(agen,A1C0n)

addtest_ageA1C0_log <- test.additivity(xx, ytrain)
addtest_ageA1C0_log 
pbf_ageA1C0_log 

## plotting

par(mfrow=c(2,2))   
  
  plot(data_full$age,data_full$A1C6,
       ylab = "A1C6",
       xlab = "age",cex=0.5,pch=20,
       main="Age vs A1C6")
 
  lines(lowess(data_full$age,data_full$A1C6), col="red",lwd=3)
  
  plot(data_full$age,data_full$logA1C6,
       ylab = "logA1C6",
       xlab = "age",cex=0.5,pch=20,
       main="Age vs logA1C6")
 
  lines(lowess(data_full$age,data_full$logA1C6), col="red",lwd=3)
   
   plot(data_full$A1C0,data_full$A1C6,
       ylab = "A1C6",
       xlab = "A1C0",cex=0.5,pch=20,
       main="A1C0 vs A1C6")
 
 lines(lowess(data_full$A1C0,data_full$A1C6), col="red",lwd=3)
  
  plot(data_full$A1C0,data_full$logA1C6,
       ylab = "logA1C6",
       xlab = "A1C0",cex=0.5,pch=20,
       main="A1C0 vs logA1C6")
 
  lines(lowess(data_full$A1C0,data_full$logA1C6), col="red",lwd=3)

par(mfrow=c(1,1)) 

