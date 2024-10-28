#Libraries
library(dplyr)
library(lmtest)
library(oglmx)
library(multiverse)
library(lmtest)
library(sandwich)
library(tidyr)
library(readr)

#Set seed
set.seed(123)

#Set WD
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

#Load data
mydata=readRDS("Data/treatedData.rds")

#Set dependent variable
depVar="Q4b" #"Q2" #Q1 #Q5 #Q6a #Q6a

if(depVar=="Q1"){
  mydata$dependentVarLoop=mydata$Q1_factorVar
  mydata$timeSubmitDependentLoop=mydata$Time_Q1_Page.Submit
}

if(depVar=="Q2"){
  mydata$dependentVarLoop=mydata$Q2_factorVar
  mydata$timeSubmitDependentLoop=mydata$Time_Q2_Page.Submit
}
if(depVar=="Q4a"){
  mydata$dependentVarLoop=mydata$Q4a_factorVar
  mydata$timeSubmitDependentLoop=mydata$Time_Q4a_Page.Submit
}
if(depVar=="Q4b"){
  mydata$dependentVarLoop=mydata$Q4b_factorVar
  mydata$timeSubmitDependentLoop=mydata$Time_Q4b_Page.Submit
}
if(depVar=="Q5"){
  mydata$dependentVarLoop=mydata$Q5_factorVar
  mydata$timeSubmitDependentLoop=mydata$Time_Q5_Page.Submit
}
if(depVar=="Q6a"){
  mydata$dependentVarLoop=mydata$Q6a_factorVar
  mydata$timeSubmitDependentLoop=mydata$Time_Q6a_Page.Submit
}
if(depVar=="Q6b"){
  mydata$dependentVarLoop=mydata$Q6b_factorVar
  mydata$timeSubmitDependentLoop=mydata$Time_Q6b_Page.Submit
}

#Create multiverse
Ms=multiverse()

inside(Ms, {
  df <- mydata %>%
    dplyr::filter(branch(excl_totalTime,
                         "none" ~ TRUE,
                         "lessThanThreeMin" ~ TotalTimeAW > 180,
                         "lessThanFiveMin" ~ TotalTimeAW > 300,
                         "lessThanFifthQuant" ~ TotalTimeAW > quantile(TotalTimeAW,c(0.05)),
                         "lessThanTenthQuant" ~ TotalTimeAW > quantile(TotalTimeAW,c(0.10)))) %>%
    dplyr::filter(branch(excl_screenTime,
                         "none" ~ TRUE,
                         "lessThanFiveSeconds" ~ timeSubmitDependentLoop > 5,
                         "lessThanTenSeconds" ~ timeSubmitDependentLoop > 10,
                         "lessThanFifthQuantScreen" ~ timeSubmitDependentLoop > quantile(timeSubmitDependentLoop,c(0.05)),
                         "lessThanTenthQuantScreen" ~ timeSubmitDependentLoop > quantile(timeSubmitDependentLoop,c(0.10)))) %>%
    dplyr::filter(branch(excl_dietInconsistent,
                         "none" ~ TRUE,
                         "Strong" ~ inconsistentDietIdentityStrong == 0 ,
                         "Weak" ~ inconsistentDietIdentityWeak == 0))
})

inside(Ms, {
  formula_Ms= as.numeric(dependentVarLoop) ~  branch(ageVar,"none" ~  NULL,"linear"~age,"polyn"~age+ageSqrd) +
    branch(femaleVar,"none"~NULL,"dummy"~female) +
    branch(political_orientation,
           "none"~NULL,
           "continuous" ~ politicalNADummy+politicalNARecoded,
           "smallCenter" ~ relevel(politicalCategories1, ref="Center"),
           "largeCenter" ~ relevel(politicalCategories2, ref="Center")) +
    branch(incomeVar,"none"~NULL,"income"~relevel(Income_factor, ref="1000 à 2000 euros")) +
    branch(foodVar,"none"~NULL,"2ndDimOnly"~ABC,"2Dims"~ABC+pcaFoodDim1)+
    branch(religionVar,"none"~NULL,"linear"~religionNARecoded, "polyn"~religionNARecoded+religionNARecodedSqrd)
})

inside(Ms, {
  model_type <- branch(model_type,
                       lm ="lm",
                       ologit ="ologit",
                       oprobit ="oprobit")
})

inside(Ms, {
  if(model_type=="lm") fit_Ms<-lm(formula_Ms, data = df) %>% coeftest(vcov=vcovHC)
  if(model_type=="ologit") fit_Ms<-summary(oglmx(formula_Ms, data = df, link="logit", constantMEAN=FALSE, constantSD=FALSE, delta=0, threshparam=NULL, robust=TRUE))$estimate
  if(model_type=="oprobit") fit_Ms<-summary(oglmx(formula_Ms, data = df, link="probit", constantMEAN=FALSE, constantSD=FALSE, delta=0, threshparam=NULL, robust=TRUE))$estimate
})

multiverse::expand(Ms)
execute_multiverse(Ms, parallel=TRUE, progress = TRUE)


MsExp=multiverse::expand(Ms)
length(unique(MsExp$.universe))
MsExp$tmpVar=sapply(1:dim(MsExp)[1], function(i){ MsExp[i,]$.results[[1]]$fit_Ms})
MsExp$Count=sapply(1:dim(MsExp)[1], function(i){ dim(MsExp[i,]$tmpVar[[1]])[1]})
MsExp=MsExp[!sapply(1:dim(MsExp)[1], function(i){ is.null(MsExp[i,]$Count[[1]]) }),]
MsExp$VarForUnnest=sapply(1:dim(MsExp)[1], function(row) seq(1, unlist(MsExp[row,"Count"])))
MsExp=MsExp %>% unnest(VarForUnnest)
MsExp$ResultsVar=sapply(1:dim(MsExp)[1], function(i){MsExp[i,]$tmpVar[[1]][1:dim(MsExp[i,]$tmpVar[[1]])[1],]})

MsExp[,c("Term","Estimate","StdErr","tvalue","pvalue")]=t(sapply(1:dim(MsExp)[1], function(i){
  if(MsExp[i,]$Count==1) var=MsExp[i,]$ResultsVar[[1]]
  if(MsExp[i,]$Count>1) var=MsExp[i,]$ResultsVar[[1]][MsExp[i,]$VarForUnnest,]
  term_tmp=ifelse(MsExp[i,]$Count==1,"Intercept",rownames(MsExp[i,]$ResultsVar[[1]])[MsExp[i,]$VarForUnnest])
  return(c(term_tmp,var[1],var[2],var[3],var[4]))
}))
MsExp$Estimate=as.numeric(MsExp$Estimate)
MsExp$StdErr=as.numeric(MsExp$StdErr)
MsExp$tvalue=as.numeric(MsExp$tvalue)
MsExp$pvalue=as.numeric(MsExp$pvalue)

#View(MsExp)

#Prepare results Table
MsExp$TermRecoded=MsExp$Term
unique(MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="Intercept","",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="Threshold (1->2)","",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="Threshold (2->3)","",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="Threshold (3->4)","",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="Threshold (4->5)","",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="Threshold (5->6)","",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="(Intercept)","",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="politicalNADummy","",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="religionNARecoded" & MsExp$religionVar=="linear","Religion (linear)",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="religionNARecoded" & MsExp$religionVar=="polyn","Religion (poly)",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="religionNARecodedSqrd","Religion Sqrd (poly)",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="ABC","Animal-based consumption",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="pcaFoodDim1","",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=='relevel(Income_factor, ref = \"1000 à 2000 euros\")Moins de 1000 euros',"Below 1000 EUR",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=='relevel(Income_factor, ref = \"1000 à 2000 euros\")2001 à 3000 euros',"Between 2000 and 3000 EUR",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=='relevel(Income_factor, ref = \"1000 à 2000 euros\")Ne souhaite pas répondre',"",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=='relevel(Income_factor, ref = \"1000 à 2000 euros\")Plus de 3001 euros',"More than 3000 EUR",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="politicalNARecoded","Political Rightism",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="relevel(politicalCategories1, ref = \"Center\")Extreme-left","Politics - Extreme Left",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="relevel(politicalCategories1, ref = \"Center\")Left","Politics - Left",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="relevel(politicalCategories1, ref = \"Center\")Right","Politics - Right",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="relevel(politicalCategories1, ref = \"Center\")Extreme-right","Politics - Extreme Right",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="relevel(politicalCategories1, ref = \"Center\")No answer","",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="relevel(politicalCategories2, ref = \"Center\")Extreme-left","Politics - Extreme Left",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="relevel(politicalCategories2, ref = \"Center\")Left","Politics - Left",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="relevel(politicalCategories2, ref = \"Center\")Right","Politics - Right",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="relevel(politicalCategories2, ref = \"Center\")Extreme-right","Politics - Extreme Right",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="relevel(politicalCategories2, ref = \"Center\")No answer","",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="female","Female",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="age" & MsExp$ageVar=="linear","Age (linear)",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="age" & MsExp$ageVar=="polyn","Age (polyn)",MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="ageSqrd","Age Sqrd (polyn)",MsExp$TermRecoded)
setdiff(MsExp$TermRecoded,"") 
vecTermSelected=setdiff(MsExp$TermRecoded,"") 

MsExp$pvalueDecision=ifelse(MsExp$pvalue<0.01,"p≤0.01",
                            ifelse(MsExp$pvalue<0.05,"p≤0.05",
                                   ifelse(MsExp$pvalue<0.1,"p≤0.10","p>0.10")))

MsExp$positiveSign=ifelse(MsExp$Estimate>0,1,0)

saveRDS(MsExp,file=paste0("ResultsMultiverse/MsExp_",depVar))

vecOrderVar=c("Female","Age (linear)","Age (polyn)","Age Sqrd (polyn)","Animal-based consumption",
              "Below 1000 EUR","Between 1000 and 2000 EUR","Between 2000 and 3000 EUR","More than 3000 EUR",
              "Political Rightism","Politics - Extreme Left","Politics - Left","Politics - Center","Politics - Right","Politics - Extreme Right","Religion (linear)","Religion (poly)","Religion Sqrd (poly)")

matResults=matrix(data=NA,nrow=length(vecOrderVar),ncol=9)
rownames(matResults)=vecOrderVar
colnames(matResults)=c("Negative p≤0.01","Negative p≤0.05","Negative p≤0.10","p>0.10","Positive p≤0.10","Positive p≤0.05","Positive p≤0.01","Share sig. neg.","Share sig. pos.")
counterVar=1
for(i in vecOrderVar){
  print(i)
  matResults[counterVar,1]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==0 & MsExp$pvalueDecision=="p≤0.01",])
  matResults[counterVar,2]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==0 & MsExp$pvalueDecision=="p≤0.05",])
  matResults[counterVar,3]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==0 & MsExp$pvalueDecision=="p≤0.10",])
  matResults[counterVar,4]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$pvalueDecision=="p>0.10",])
  matResults[counterVar,5]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==1 & MsExp$pvalueDecision=="p≤0.10",])
  matResults[counterVar,6]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==1 & MsExp$pvalueDecision=="p≤0.05",])
  matResults[counterVar,7]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==1 & MsExp$pvalueDecision=="p≤0.01",])
  matResults[counterVar,1:7]=matResults[counterVar,1:7]/nrow(MsExp[MsExp$TermRecoded==i,])
  counterVar=counterVar+1
}
matResults[,8]=rowSums(matResults[,1:3])
matResults[,9]=rowSums(matResults[,5:7])
matResults


write.csv(as.data.frame(matResults), file=paste0("Export/",depVar,".csv"))