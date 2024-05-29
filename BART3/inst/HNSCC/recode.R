
## hnscc = read.csv("hnscc.csv", stringsAsFactors=TRUE)
## hnscc <- data.frame(hnscc[,c("Age", "Sex", "Race", "Smoking", 
## "Alcohol", "Prior_Chemo", "Primary_Disease_Site", "P16", "ECS",
## "Positive_M", "Close_M", "PNI", "LVI", "Risk_Stratification", 
## "DFS", "Lymph_Nodes_Involved_2", "DFS_months", "study",
## "Initial_Diagnosis_Date", "Surgery_date")])
## save(hnscc, file = '~/git/bnptools/BART3/data/hnscc.rda')

library(BART3)
data(hnscc)

hnscc$Study =  ifelse(hnscc$study=="new",1,0)
hnscc$White = ifelse(hnscc$Race=="White",1,0)
hnscc$Black = ifelse(hnscc$Race=="Black",1,0)
hnscc$Asian = ifelse(hnscc$Race=="Asian",1,0)
hnscc$Male = ifelse(hnscc$Sex=="Male",1,0)
hnscc$Smoker = ifelse(hnscc$Smoking=="Yes",1,0)
hnscc$Alcohol = ifelse(hnscc$Alcohol=="Yes",1,0)
##hnscc$Alcohol = hnscc$Alcoholic
hnscc$PriorChemo = ifelse(hnscc$Prior_Chemo=="Yes",1,0)
hnscc$Larynx = ifelse(hnscc$Primary_Disease_Site=="Larynx",1,0)
hnscc$OralCavity = ifelse(hnscc$Primary_Disease_Site=="Oral Cavity",1,0)
hnscc$Oropharynx = ifelse(hnscc$Primary_Disease_Site=="Oropharynx",1,0)
hnscc$P16M = ifelse(hnscc$P16=="",1,0)
hnscc$P16P = ifelse(hnscc$P16=="Positive",1,0)
hnscc$P16N = ifelse(hnscc$P16=="Negative",1,0)
hnscc$ECSM = ifelse(hnscc$ECS=="",1,0)
hnscc$ECSY = ifelse(hnscc$ECS=="Yes",1,0)
hnscc$ECSN = ifelse(hnscc$ECS=="No",1,0)
hnscc$Positive_MM = ifelse(hnscc$Positive_M=="",1,0)
hnscc$Positive_MY = ifelse(hnscc$Positive_M=="Yes",1,0)
hnscc$Positive_MN = ifelse(hnscc$Positive_M=="No",1,0)
hnscc$Close_MM = ifelse(hnscc$Close_M=="",1,0)
hnscc$Close_MY = ifelse(hnscc$Close_M=="Yes",1,0)
hnscc$Close_MN = ifelse(hnscc$Close_M=="No",1,0)
hnscc$PNIM = ifelse(hnscc$PNI=="",1,0)
hnscc$PNIY = ifelse(hnscc$PNI=="Yes",1,0)
hnscc$PNIN = ifelse(hnscc$PNI=="No",1,0)
hnscc$LVIM = ifelse(hnscc$LVI=="",1,0)
hnscc$LVIY = ifelse(hnscc$LVI=="Yes",1,0)
hnscc$LVIN = ifelse(hnscc$LVI=="No",1,0)
hnscc$LN2M = ifelse(hnscc$Lymph_Nodes_Involved_2=="",1,0)
hnscc$LN2Y = ifelse(hnscc$Lymph_Nodes_Involved_2=="Yes",1,0)
hnscc$LN2N = ifelse(hnscc$Lymph_Nodes_Involved_2=="No",1,0)
hnscc$RSM = ifelse(hnscc$Risk_Stratification=="",1,0)
hnscc$RSY = ifelse(hnscc$Risk_Stratification=="High Risk",1,0)
hnscc$RSN = ifelse(hnscc$Risk_Stratification=="Not High Risk",1,0)

hnscc$month_diag_surg = (as.Date(hnscc$Surgery_date) - as.Date(hnscc$Initial_Diagnosis_Date))/30.25

ddate = hnscc$Initial_Diagnosis_Date
hnscc$dyear = as.numeric(substr(ddate, 1, 4))

sdate = hnscc$Surgery_date
hnscc$syear = as.numeric(substr(sdate, 1, 4))

hnscc$diag_year = ifelse(hnscc$dyear>2001 & hnscc$dyear<2011,1,2)
hnscc$diag_year = ifelse(is.na(hnscc$diag_year)==TRUE, 2, hnscc$diag_year)

times <- hnscc$DFS_months
delta <- hnscc$DFS

hnscc <- data.frame(hnscc[,c("DFS_months", "DFS", 
"Study", "diag_year", "Age", "Male", 
"White", "Black", "Asian", "Smoker", 
"Alcohol", "PriorChemo", "Oropharynx", "Larynx", "OralCavity",
"P16M", "P16P", "P16N", "ECSM", "ECSY", "ECSN",
"Positive_MM", "Positive_MY", "Positive_MN", 
"Close_MM", "Close_MY", "Close_MN", 
"PNIM", "PNIY", "PNIN", "LVIM", "LVIY", "LVIN", 
"RSM", "RSY", "RSN",  
"LN2M", "LN2Y", "LN2N")])

