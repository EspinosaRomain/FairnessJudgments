#Libraries
library(dplyr)
library(lmtest)
library(oglmx)
library(multiverse)
library(lmtest)
library(sandwich)
library(tidyr)
library(readr)
library(nnet)

#Set seed
set.seed(123)

#Set WD
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

#Load data
mydata=readRDS("Data/treatedData.rds")

#Set dependent variable
depVar="Q6b"

if(depVar=="Q3"){
  mydata$dependentVarLoop=relevel(mydata$Q3_factorVar, ref = "Equal weight")
  mydata$timeSubmitDependentLoop=mydata$Time_Q3_Page.Submit
}
if(depVar=="Q5"){
  mydata$dependentVarLoop=relevel(mydata$Q5_factorVar, ref = "Equal weight")
  mydata$timeSubmitDependentLoop=mydata$Time_Q5_Page.Submit
}
if(depVar=="Q6a"){
  mydata$dependentVarLoop=relevel(mydata$Q6a_factorVar, ref = "Neither")
  mydata$timeSubmitDependentLoop=mydata$Time_Q6a_Page.Submit
}
if(depVar=="Q6b"){
  mydata$dependentVarLoop=relevel(mydata$Q6b_factorVar, ref = "Neither")
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
  fit_Ms<-summary(multinom(formula_Ms, data = df))
})

multiverse::expand(Ms)
execute_multiverse(Ms, parallel=TRUE, progress = TRUE)


MsExp=multiverse::expand(Ms)
length(unique(MsExp$.universe))
MsExp$tmpVar=sapply(1:dim(MsExp)[1], function(i){ MsExp[i,]$.results[[1]]$fit_Ms})
MsExp$Count=sapply(1:dim(MsExp)[1], function(i){ length(MsExp[i,]$tmpVar[[1]]$coefficients)/2})
MsExp=MsExp[!sapply(1:dim(MsExp)[1], function(i){ is.null(MsExp[i,]$Count[[1]]) }),]
MsExp$VarForUnnest=sapply(1:dim(MsExp)[1], function(row) seq(1, unlist(MsExp[row,"Count"])))
MsExp=MsExp %>% unnest(VarForUnnest)
MsExp$Answer=lapply(1:dim(MsExp)[1], function(row) return(1:2))
MsExp=MsExp %>% unnest(Answer)
MsExp$EstimatesVar=lapply(1:dim(MsExp)[1], function(i){MsExp[i,]$tmpVar[[1]]$coefficients})
MsExp$StdErrorVar=lapply(1:dim(MsExp)[1], function(i){MsExp[i,]$tmpVar[[1]]$standard.errors})

MsExp[,c("Term","Estimate","StdErr","tvalue","pvalue")]=t(sapply(1:dim(MsExp)[1], function(i){
  var1=MsExp[i,]$EstimatesVar[[1]][MsExp[i,]$Answer,MsExp[i,]$VarForUnnest]
  var2=MsExp[i,]$StdErrorVar[[1]][MsExp[i,]$Answer,MsExp[i,]$VarForUnnest]
  term_tmp=ifelse(MsExp[i,]$Count==1,"Intercept",colnames(MsExp[i,]$EstimatesVar[[1]])[MsExp[i,]$VarForUnnest])
  var3=var1/var2
  var4=(1 - pnorm(abs(var3), 0, 1)) * 2
  return(c(term_tmp,var1,var2,var3,var4))
}))
MsExp$Estimate=as.numeric(MsExp$Estimate)
MsExp$StdErr=as.numeric(MsExp$StdErr)
MsExp$tvalue=as.numeric(MsExp$tvalue)
MsExp$pvalue=as.numeric(MsExp$pvalue)


#Prepare results Table
MsExp$TermRecoded=MsExp$Term
unique(MsExp$TermRecoded)
MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="Intercept","",MsExp$TermRecoded)
# MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="Threshold (1->2)","",MsExp$TermRecoded)
# MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="Threshold (2->3)","",MsExp$TermRecoded)
# MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="Threshold (3->4)","",MsExp$TermRecoded)
# MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="Threshold (4->5)","",MsExp$TermRecoded)
# MsExp$TermRecoded=ifelse(MsExp$TermRecoded=="Threshold (5->6)","",MsExp$TermRecoded)
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

saveRDS(MsExp,file=paste0("ResultsMultiverse/MsExpMULTINOMIAL_",depVar))

vecOrderVar=c("Female","Age (linear)","Age (polyn)","Age Sqrd (polyn)","Animal-based consumption",
              "Below 1000 EUR","Between 1000 and 2000 EUR","Between 2000 and 3000 EUR","More than 3000 EUR",
              "Political Rightism","Politics - Extreme Left","Politics - Left","Politics - Center","Politics - Right","Politics - Extreme Right","Religion (linear)","Religion (poly)","Religion Sqrd (poly)")


#FIRST ANSWER
matResultsAnswer1=matrix(data=NA,nrow=length(vecOrderVar),ncol=9)
rownames(matResultsAnswer1)=vecOrderVar
colnames(matResultsAnswer1)=c("Negative p≤0.01","Negative p≤0.05","Negative p≤0.10","p>0.10","Positive p≤0.10","Positive p≤0.05","Positive p≤0.01","Share sig. neg.","Share sig. pos.")
counterVar=1
for(i in vecOrderVar){
  print(i)
  matResultsAnswer1[counterVar,1]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==0 & MsExp$pvalueDecision=="p≤0.01" & MsExp$Answer==1,])
  matResultsAnswer1[counterVar,2]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==0 & MsExp$pvalueDecision=="p≤0.05" & MsExp$Answer==1,])
  matResultsAnswer1[counterVar,3]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==0 & MsExp$pvalueDecision=="p≤0.10" & MsExp$Answer==1,])
  matResultsAnswer1[counterVar,4]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$pvalueDecision=="p>0.10" & MsExp$Answer==1,])
  matResultsAnswer1[counterVar,5]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==1 & MsExp$pvalueDecision=="p≤0.10" & MsExp$Answer==1,])
  matResultsAnswer1[counterVar,6]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==1 & MsExp$pvalueDecision=="p≤0.05" & MsExp$Answer==1,])
  matResultsAnswer1[counterVar,7]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==1 & MsExp$pvalueDecision=="p≤0.01" & MsExp$Answer==1,])
  matResultsAnswer1[counterVar,1:7]=matResultsAnswer1[counterVar,1:7]/nrow(MsExp[MsExp$TermRecoded==i & MsExp$Answer==1,])
  counterVar=counterVar+1
}
matResultsAnswer1[,8]=rowSums(matResultsAnswer1[,1:3])
matResultsAnswer1[,9]=rowSums(matResultsAnswer1[,5:7])
matResultsAnswer1

write.csv(as.data.frame(matResultsAnswer1), file=paste0("Export/Multinom_Answer1_",depVar,".csv"))

#SECOND ANSWER
matResultsAnswer2=matrix(data=NA,nrow=length(vecOrderVar),ncol=9)
rownames(matResultsAnswer2)=vecOrderVar
colnames(matResultsAnswer2)=c("Negative p≤0.01","Negative p≤0.05","Negative p≤0.10","p>0.10","Positive p≤0.10","Positive p≤0.05","Positive p≤0.01","Share sig. neg.","Share sig. pos.")
counterVar=1
for(i in vecOrderVar){
  print(i)
  matResultsAnswer2[counterVar,1]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==0 & MsExp$pvalueDecision=="p≤0.01" & MsExp$Answer==2,])
  matResultsAnswer2[counterVar,2]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==0 & MsExp$pvalueDecision=="p≤0.05" & MsExp$Answer==2,])
  matResultsAnswer2[counterVar,3]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==0 & MsExp$pvalueDecision=="p≤0.10" & MsExp$Answer==2,])
  matResultsAnswer2[counterVar,4]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$pvalueDecision=="p>0.10" & MsExp$Answer==2,])
  matResultsAnswer2[counterVar,5]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==1 & MsExp$pvalueDecision=="p≤0.10" & MsExp$Answer==2,])
  matResultsAnswer2[counterVar,6]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==1 & MsExp$pvalueDecision=="p≤0.05" & MsExp$Answer==2,])
  matResultsAnswer2[counterVar,7]=nrow(MsExp[MsExp$TermRecoded==i & MsExp$positiveSign==1 & MsExp$pvalueDecision=="p≤0.01" & MsExp$Answer==2,])
  matResultsAnswer2[counterVar,1:7]=matResultsAnswer2[counterVar,1:7]/nrow(MsExp[MsExp$TermRecoded==i & MsExp$Answer==2,])
  counterVar=counterVar+1
}
matResultsAnswer2[,8]=rowSums(matResultsAnswer2[,1:3])
matResultsAnswer2[,9]=rowSums(matResultsAnswer2[,5:7])
matResultsAnswer2


write.csv(as.data.frame(matResultsAnswer2), file=paste0("Export/Multinom_Answer2_",depVar,".csv"))
