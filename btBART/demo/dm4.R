## original file is/was ...
## ~btuyishimire/allfiles/RESEARCH/NEWBART/TEST/A1CDATA/a1cdata_full.R
## compare results with 
## ~btuyishimire/allfiles/RESEARCH/NEWBART/TEST/A1CDATA/a1cdata_full.pdf
## library(btBART)
library(BART4)
library(BART3)
library(BayesTree)
library(lmeVarComp)
library(mgcv)
##source("~btuyishimire/allfiles/RESEARCH/Rfunctions/lpml.R")
source("../R/lpml.R")

##A1C DATA FROM THE PHI DIABETES STUDY
##path="~btuyishimire/Documents/PHI/TARIMA/manuscript3/illustrative_example/"
##path <- system.file('dm', package='btBART')
path <- '../inst/dm'
data <- read.csv(paste(path,"monthly_data12.csv",sep="/"),header=TRUE)
demograph <- c("age","Female","race_cat","new_insurance",
               "marital_status_cat","religion_cat")
health <- c("A1C0","DEPRESS","HYPERTSN","bmicat","hdlcat","ldlcat",
            "systolicbpcat","diastolicbpcat")
trt <- c("newcurr_insulin","neworder_insulin","newdisp_insulin",
         "newcurr_met","neworder_met","newdisp_met")
## data <- data[!is.na(data$A1C6), c('A1C6', 'logA1C6',
##   demograph, health, trt)]
## data$age <- pmax(21, pmin(floor(data$age), 85))
## write.csv(data, file='../inst/dm/monthly_data12.csv', row.names=FALSE)

##-------  1. Check additivity of health and demographic factors ------------

data_full <- data
n_full <- nrow(data_full)
## data_notrt <- data[data$newcurr_insulin=="No or Unknown" & data$newcurr_met == "No or Unknown",]
## n_notrt <- nrow(data_notrt)

ndpost <- 1000; burn <- 100

## A1C
##----- First do analysis using data without treatment 

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

##TWO BARTS

post_2b <- twbarts( x.train1 = xtrain1,
                    x.train2 = xtrain2,
                    y.train  = ytrain,   
                    ndpost   = ndpost, 
                    nskip    = burn)
 
sigs2 <- post_2b$sigma[(burn+1):(burn+ndpost)]  
mus2  <- post_2b$yhat.train
out2  <- lpml(ytrain,mus2,sigs2) 
mse2  <- mean((ytrain-post_2b$yhat.train.mean)^2)

(pbf_orig <- exp(out1$LPML-out2$LPML))

## log A1C
##----- First do analysis using data without treatment 

xtrain1 <- makeind(data_full[,demograph])
xtrain2 <- makeind(data_full[,health])
xtrain_all <- cbind(xtrain1,xtrain2)
ytrain <- data_full$logA1C6

post_1b <- wbart( x.train  = xtrain_all,
                  y.train  = ytrain,
                  ndpost   = ndpost,
                  nskip    = burn)
 
mse1 <- mean((ytrain-post_1b$yhat.train.mean)^2)

sigs1 <- post_1b$sigma[(burn+1):(burn+ndpost)]
mus1 <- post_1b$yhat.train
out1 <- lpml(ytrain,mus1,sigs1)

##TWO BARTS

post_2b <- twbarts( x.train1 = xtrain1,
                    x.train2 = xtrain2,
                    y.train  = ytrain,   
                    ndpost   = ndpost, 
                    nskip    = burn)
 
sigs2 <- post_2b$sigma[(burn+1):(burn+ndpost)]  
mus2  <- post_2b$yhat.train
out2  <- lpml(ytrain,mus2,sigs2) 
mse2  <- mean((ytrain-post_2b$yhat.train.mean)^2)

(pbf_log <- exp(out1$LPML-out2$LPML))

## A1C
ytrain <- data_full$A1C6

post_1b <- wbart( x.train  = data_full[,c("age","A1C0")] ,
                  y.train  = ytrain,
                  ndpost   = ndpost,
                  nskip    = burn)
 
sigs1 <- post_1b$sigma[(burn+1):(burn+ndpost)]
mus1 <- post_1b$yhat.train
out1 <- lpml(ytrain,mus1,sigs1)

##TWO BARTS

post_2b <- twbarts( x.train1 = data_full[,c("age")],
                    x.train2 = data_full[,c("A1C0")],
                    y.train  = ytrain,   
                    ndpost   = ndpost, 
                    nskip    = burn)  

sigs2 <- post_2b$sigma[(burn+1):(burn+ndpost)]   
mus2  <- post_2b$yhat.train 
out2  <- lpml(ytrain,mus2,sigs2) 
mse2  <- mean((ytrain-post_2b$yhat.train.mean )^2)

##Normalize data
agen <- (data_full$age-min(data_full$age))/(max(data_full$age)-min(data_full$age))
A1C0n <- (data_full$A1C0-min(data_full$A1C0))/(max(data_full$A1C0)-min(data_full$A1C0))

xx <- cbind(agen,A1C0n)

(addtest_ageA1C0_ori <- test.additivity(xx, ytrain))
(pbf_ageA1C0_ori <- exp(out1$LPML-out2$LPML))

ga1_s <- gam(ytrain~s(data_full$age,data_full$A1C0,bs="tp"), method="REML")
ga2_s <- gam(ytrain~s(data_full$age,bs="tp")+ s(data_full$A1C0,bs="tp"),method="REML")
ga1_te <- gam(ytrain~te(data_full$age,data_full$A1C0,bs="tp"), method="REML")
ga2_te <- gam(ytrain~te(data_full$age,bs="tp")+ te(data_full$A1C0,bs="tp"),method="REML")

AIC(ga1_s,ga2_s)
AIC(ga1_te,ga2_te)
anova.gam(ga1_s,ga2_s,test="F")
anova(ga1_te,ga2_te,test="F")

## log A1C

ytrain <- data_full$logA1C6

post_1b <- wbart( x.train  = data_full[,c("age","A1C0")] ,
                  y.train  = ytrain,
                  ndpost   = ndpost,
                  nskip    = burn)
 
sigs1 <- post_1b$sigma[(burn+1):(burn+ndpost)]
mus1 <- post_1b$yhat.train
out1 <- lpml(ytrain,mus1,sigs1)

##TWO BARTS

post_2b <- twbarts( x.train1 = data_full[,c("age")],
                    x.train2 = data_full[,c("A1C0")],
                    y.train  = ytrain,   
                    ndpost   = ndpost, 
                    nskip    = burn)  

sigs2 <- post_2b$sigma[(burn+1):(burn+ndpost)]   
mus2  <- post_2b$yhat.train 
out2  <- lpml(ytrain,mus2,sigs2) 
mse2  <- mean((ytrain-post_2b$yhat.train.mean )^2)

##Normalize data
agen <- (data_full$age-min(data_full$age))/(max(data_full$age)-min(data_full$age))
A1C0n <- (data_full$A1C0-min(data_full$A1C0))/(max(data_full$A1C0)-min(data_full$A1C0))

xx <- cbind(agen,A1C0n)

(addtest_ageA1C0_log <- test.additivity(xx, ytrain))
(pbf_ageA1C0_log <- exp(out1$LPML-out2$LPML))

ga1_s <- gam(ytrain~s(data_full$age,data_full$A1C0,bs="tp"), method="REML")
ga2_s <- gam(ytrain~s(data_full$age,bs="tp")+ s(data_full$A1C0,bs="tp"),method="REML")
ga1_te <- gam(ytrain~te(data_full$age,data_full$A1C0,bs="tp"), method="REML")
ga2_te <- gam(ytrain~te(data_full$age,bs="tp")+ te(data_full$A1C0,bs="tp"),method="REML")

AIC(ga1_s,ga2_s)
AIC(ga1_te,ga2_te)
anova.gam(ga1_s,ga2_s,test="F")
anova.gam(ga1_te,ga2_te,test="F")

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

##dev.copy2pdf(file='dm.pdf')
