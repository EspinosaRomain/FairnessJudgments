#Load libraries
library(ggplot2)
library(pastecs)
library(tidyverse)
library(reshape2)
library("corrplot")
library("factoextra")
library("FactoMineR")
library(dplyr)
library(ggridges) 
library(hrbrthemes) 
library(memisc)

#Set WD
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

#Load data
mydata=readRDS("Data/anonymizedData.csv")

#Age 
table(mydata$Age)
mydata$age=mydata$Age
mydata$age=ifelse(mydata$age=="75 or more","75",mydata$age)
table(mydata$age)
mydata$age=as.numeric(mydata$age)
mydata$ageSqrd=mydata$age^2

#Gender
#Issue of recoding with Qualtrics
mydata$Gender
mydata$Gender=ifelse(mydata$Gender=="Célibataire","Homme",mydata$Gender)
mydata$Gender=ifelse(mydata$Gender=="En couple, mais pas marié","Femme",mydata$Gender)
mydata$Gender=ifelse(mydata$Gender=="Marié","Autre",mydata$Gender)
table(mydata$Gender)

#Political self-position
mydata$political=as.numeric(mydata$Political.ideology_1)
table(mydata$political)

#Income levels
table(mydata$Income)

#Marital status
table(mydata$Marital.status)

#Demographics

#1. Gender
mydata$female=ifelse(mydata$Gender=="Femme",1,0)
round(table(mydata$female)/length(mydata$female),4)*100

#2. Age
mydata$ageCategory=NA
mydata$ageCategory=ifelse(mydata$age<=29,"18 to 29 y.o.",mydata$ageCategory)
mydata$ageCategory=ifelse(mydata$age>=30 & mydata$age<=44,"30 to 44 y.o.",mydata$ageCategory)
mydata$ageCategory=ifelse(mydata$age>=45 & mydata$age<=59,"45 to 59 y.o.",mydata$ageCategory)
mydata$ageCategory=ifelse(mydata$age>=60,"60+ y.o.",mydata$ageCategory)
round(table(mydata$ageCategory)/sum(table(mydata$ageCategory)),4)*100

#3. Income
table(mydata$Income)
vecRevenus=round(table(mydata$Income)/length(mydata$Income),4)*100
vecRevenus

#Religion
mydata$religionNADummy=ifelse(mydata$Religion=="",1,0)
mydata$religionNARecoded=cases(
  1<-mydata$Religion=="1. Non, pas du tout",
  2<-mydata$Religion=="2",
  3<-mydata$Religion=="3",
  4<-mydata$Religion=="4",
  5<-mydata$Religion=="5",
  6<-mydata$Religion=="6",
  7<-mydata$Religion=="7. Oui, beaucoup",
  .default=0
)
mydata$religionNARecoded=ifelse(mydata$religionNARecoded==0,mean(mydata[mydata$religionNARecoded>0,]$religionNARecoded),mydata$religionNARecoded)
table(mydata$religionNARecoded)
mydata$religionNARecodedSqrd=mydata$religionNARecoded^2

#Quality check
#Short answers: time<10seconds per screen on SW questions
shortPage1=mydata$Time_Q1_Page.Submit<10
shortPage2=mydata$Time_Q2_Page.Submit<10
shortPage3=mydata$Time_Q3_Page.Submit<10
shortPage4a=mydata$Time_Q4a_Page.Submit<10
shortPage4b=mydata$Time_Q4b_Page.Submit<10
shortPage5=mydata$Time_Q5_Page.Submit<10
shortPage6a=mydata$Time_Q6a_Page.Submit<10
shortPage6b=mydata$Time_Q6b_Page.Submit<10
mydata$totalShort=shortPage1+shortPage2+shortPage3+shortPage4a+shortPage4b+shortPage5+shortPage6a+shortPage6b
table(mydata$totalShort)

#Generate Control Variables
#Income
uniqueIncome=unique(mydata$Income)
mydata$Income_factor=factor(ifelse(mydata$Income=="","Ne souhaite pas répondre",mydata$Income), levels=c(uniqueIncome[5],uniqueIncome[1],uniqueIncome[3],uniqueIncome[4],uniqueIncome[2]))

#Food:
# Food.consumption_1_num: Red meat
# Food.consumption_2_num: White meat
# Food.consumption_3_num: Fish
# Food.consumption_4_num: Eggs
# Food.consumption_5_num: Dairy
# Food.consumption_6_num: Vegetables
# Food.consumption_7_num: Legumes
# Food.consumption_8_num: Fruits
# Food.consumption_9_num: Starchy products
vecFoods=c("Food.consumption_1","Food.consumption_2","Food.consumption_3","Food.consumption_4","Food.consumption_5","Food.consumption_6","Food.consumption_7","Food.consumption_8","Food.consumption_9")
for(i in vecFoods){
  mydata[[paste0(i,"_num")]]=ifelse(mydata[[i]]=="Jamais",0,NA)
  mydata[[paste0(i,"_num")]]=ifelse(mydata[[i]]=="Quelques fois par an",1,mydata[[paste0(i,"_num")]])
  mydata[[paste0(i,"_num")]]=ifelse(mydata[[i]]=="Quelques fois par mois",2,mydata[[paste0(i,"_num")]])
  mydata[[paste0(i,"_num")]]=ifelse(mydata[[i]]=="Quelques fois par semaine",3,mydata[[paste0(i,"_num")]])
  mydata[[paste0(i,"_num")]]=ifelse(mydata[[i]]=="Presque à chaque repas",4,mydata[[paste0(i,"_num")]])
}
pca_food=prcomp(mydata[,paste0(vecFoods,"_num")], center = TRUE,scale. = TRUE)
summary(pca_food)
var_pca_food<- get_pca_var(pca_food)
var_pca_food$cor[,1:3]
mydata$ABC=get_pca_ind(pca_food)$coord[,2]
mydata$pcaFoodDim1=get_pca_ind(pca_food)$coord[,1]

#Political orientation:
#Three ways to classify political orientation:
#1. Continuous variable with a dummy for missing values. politicalNADummy+politicalNARecoded
#2. Categorical variable: Extreme-left: 0,1,2; Left:3,4; Center: 5; Right: 6,7; Extreme-Right: 8,9,10: politicalCategories1
#3. Categorical variable: Extreme-left: 0,1; Left:2,3; Center: 4,5,6; Right: 7,8; Extreme-Right: 9,10: politicalCategories2
mydata$politicalNADummy=ifelse(is.na(mydata$political),1,0)
mydata$politicalNARecoded=ifelse(is.na(mydata$political),0,mydata$political)
mydata$politicalCategories1=cases(
      "Extreme-left"=mydata$political<3,
      "Left"=mydata$political<5,
      "Center"=mydata$political<6,
      "Right"=mydata$political<8,
      "Extreme-right"=mydata$political<11,
      "No answer"=is.na(mydata$political)
      )
table(mydata$political,mydata$politicalCategories, useNA = "always")
mydata$politicalCategories2=cases(
  "Extreme-left"=mydata$political<2,
  "Left"=mydata$political<4,
  "Center"=mydata$political<7,
  "Right"=mydata$political<9,
  "Extreme-right"=mydata$political<11,
  "No answer"=is.na(mydata$political)
)
table(mydata$political,mydata$politicalCategories2, useNA = "always")

#Time spent
mydata$TotalTimeAW=mydata$Time_Q1_Page.Submit+mydata$Time_Q2_Page.Submit+mydata$Time_Q3_Page.Submit+
  mydata$Time_Q4a_Page.Submit+mydata$Time_Q4b_Page.Submit+mydata$Time_Q5_Page.Submit+
  mydata$Time_Q6a_Page.Submit+mydata$Time_Q6b_Page.Submit
stat.desc(mydata$TotalTimeAW)

#Exclusion rules
#Dietary identity
table(mydata$Diet.identity, useNA="always")
uniqueDiet=unique(mydata$Diet.identity)
mydata$inconsistentDietIdentityStrong=cases(
  1<-mydata$Diet.identity=="Vegan" & (mydata$Food.consumption_1!="Jamais" | 
                                       mydata$Food.consumption_2!="Jamais"| 
                                       mydata$Food.consumption_3!="Jamais" |
                                       mydata$Food.consumption_4!="Jamais" |
                                     mydata$Food.consumption_5!="Jamais"),
  1<-mydata$Diet.identity=="Végétarien" & (mydata$Food.consumption_1!="Jamais" | 
                                       mydata$Food.consumption_2!="Jamais"| 
                                       mydata$Food.consumption_3!="Jamais"),
    .default=0
)
stat.desc(mydata$inconsistentDietIdentityStrong)
mydata$inconsistentDietIdentityWeak=cases(
  1<-mydata$Diet.identity=="Vegan" & (mydata$Food.consumption_1_num>1 | 
                                        mydata$Food.consumption_2_num>1| 
                                        mydata$Food.consumption_3_num>1|
                                        mydata$Food.consumption_4_num>1|
                                        mydata$Food.consumption_5_num>1),
  1<-mydata$Diet.identity=="Végétarien" & (mydata$Food.consumption_1_num>1| 
                                             mydata$Food.consumption_2_num>1| 
                                             mydata$Food.consumption_3_num>1),
  .default=0
)
stat.desc(mydata$inconsistentDietIdentityWeak)


#Recode dependent variables
#Q1
uniqueQ1=unique(mydata$Q1_SW)
table(mydata$Q1_SW)
mydata$Q1_factorVar=ifelse(mydata$Q1_SW==uniqueQ1[1],"More weight","")
mydata$Q1_factorVar=ifelse(mydata$Q1_SW==uniqueQ1[2],"Equal weight",mydata$Q1_factorVar)
mydata$Q1_factorVar=ifelse(mydata$Q1_SW==uniqueQ1[3],"Much less weight",mydata$Q1_factorVar)
mydata$Q1_factorVar=ifelse(mydata$Q1_SW==uniqueQ1[4],"A bit less weight",mydata$Q1_factorVar)
mydata$Q1_factorVar=ifelse(mydata$Q1_SW==uniqueQ1[5],"No weight",mydata$Q1_factorVar)
mydata$Q1_factorVar=ifelse(mydata$Q1_SW==uniqueQ1[6],"Instrumental weight",mydata$Q1_factorVar)
mydata$Q1_factorVar=factor(mydata$Q1_factorVar, levels = c("No weight","Instrumental weight","Much less weight", "A bit less weight", "Equal weight", "More weight"))
table(mydata$Q1_factorVar)

#Q2
uniqueQ2=unique(mydata$Q2_SW)
table(mydata$Q2_SW)
mydata$Q2_factorVar=ifelse(mydata$Q2_SW==uniqueQ2[1],"Direct and indirect","")
mydata$Q2_factorVar=ifelse(mydata$Q2_SW==uniqueQ2[2],"Direct only",mydata$Q2_factorVar)
mydata$Q2_factorVar=factor(mydata$Q2_factorVar, levels = c("Direct only", "Direct and indirect"))

#Q3
uniqueQ3=unique(mydata$Q3_SW)
table(mydata$Q3_SW)
mydata$Q3_factorVar=ifelse(mydata$Q3_SW==uniqueQ3[1],"Equal weight","")
mydata$Q3_factorVar=ifelse(mydata$Q3_SW==uniqueQ3[2],"More weight to cow",mydata$Q3_factorVar)
mydata$Q3_factorVar=ifelse(mydata$Q3_SW==uniqueQ3[3],"Less weight to cow",mydata$Q3_factorVar)
mydata$Q3_factorVar=factor(mydata$Q3_factorVar, levels = c("Less weight to cow", "Equal weight", "More weight to cow"))

#Q4a
uniqueQ4a=unique(mydata$Q4a_SW)
table(mydata$Q4a_SW)
mydata$Q4a_factorVar=ifelse(mydata$Q4a_SW==uniqueQ4a[1],"Avoid this life","")
mydata$Q4a_factorVar=ifelse(mydata$Q4a_SW==uniqueQ4a[2],"Against avoiding",mydata$Q4a_factorVar)
mydata$Q4a_factorVar=factor(mydata$Q4a_factorVar, levels = c("Against avoiding", "Avoid this life"))

#Q4b
uniqueQ4b=unique(mydata$Q4b_SW)
table(mydata$Q4b_SW)
mydata$Q4b_factorVar=ifelse(mydata$Q4b_SW==uniqueQ4b[1],"Create this life","")
mydata$Q4b_factorVar=ifelse(mydata$Q4b_SW==uniqueQ4b[2],"Against creating",mydata$Q4b_factorVar)
mydata$Q4b_factorVar=factor(mydata$Q4b_factorVar, levels = c("Against creating","Create this life"))

#Q5
uniqueQ5=unique(mydata$Q5_SW)
table(mydata$Q5_SW)
uniqueQ5
mydata$Q5_factorVar=ifelse(mydata$Q5_SW==uniqueQ5[1],"Equal weight","")
mydata$Q5_factorVar=ifelse(mydata$Q5_SW==uniqueQ5[2],"More weight for captive animal",mydata$Q5_factorVar)
mydata$Q5_factorVar=ifelse(mydata$Q5_SW==uniqueQ5[3],"More weight for wild animal",mydata$Q5_factorVar)
mydata$Q5_factorVar=factor(mydata$Q5_factorVar, levels = c("More weight for captive animal","Equal weight","More weight for wild animal"))

#Q6a
uniqueQ6a=unique(mydata$Q6a_SW)
table(mydata$Q6a_SW)
mydata$Q6a_factorVar=ifelse(mydata$Q6a_SW==uniqueQ6a[1],"Desirable birth","")
mydata$Q6a_factorVar=ifelse(mydata$Q6a_SW==uniqueQ6a[2],"Neither",mydata$Q6a_factorVar)
mydata$Q6a_factorVar=ifelse(mydata$Q6a_SW==uniqueQ6a[3],"Non-desirable birth",mydata$Q6a_factorVar)
mydata$Q6a_factorVar=factor(mydata$Q6a_factorVar, levels = c("Non-desirable birth", "Neither", "Desirable birth"))

#Q6b
uniqueQ6b=unique(mydata$Q6b_SW)
table(mydata$Q6b_SW)
mydata$Q6b_factorVar=ifelse(mydata$Q6b_SW==uniqueQ6b[1],"Neither","")
mydata$Q6b_factorVar=ifelse(mydata$Q6b_SW==uniqueQ6b[2],"Desirable inequality reduction",mydata$Q6b_factorVar)
mydata$Q6b_factorVar=ifelse(mydata$Q6b_SW==uniqueQ6b[3],"Non-desirable inequality reduction",mydata$Q6b_factorVar)
mydata$Q6b_factorVar=factor(mydata$Q6b_factorVar, levels = c("Non-desirable inequality reduction",
                                                             "Neither",
                                                             "Desirable inequality reduction"))


#Treated data
saveRDS(mydata,"Data/treatedData.rds")