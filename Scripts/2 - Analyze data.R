#Load libraries
library(dplyr)
library(cluster)
library(NbClust)
library(asbio)

#Set WD
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

#Load data
mydata=readRDS("Data/treatedData.rds")

#Descriptive statistics
t1=table(mydata$Q1_SW)
t1
write.table(t1,"Export/tableQ1_SW.csv", sep=";")
as.matrix(table(mydata$Q1_SW)/sum(table(mydata$Q1_SW)))

t2=table(mydata$Q2_SW)
t2
write.table(t2,"Export/tableQ2_SW.csv", sep=";")
as.matrix(table(mydata$Q2_SW)/sum(table(mydata$Q2_SW)))

t3=table(mydata$Q3_SW)
t3
write.table(t3,"Export/tableQ3_SW.csv", sep=";")
as.matrix(table(mydata$Q3_SW)/sum(table(mydata$Q3_SW)))

t4a=table(mydata$Q4a_SW)
t4a
write.table(t4a,"Export/tableQ4a_SW.csv", sep=";")
as.matrix(table(mydata$Q4a_SW)/sum(table(mydata$Q4a_SW)))

t4b=table(mydata$Q4b_SW)
t4b
write.table(t4b,"Export/tableQ4b_SW.csv", sep=";")
as.matrix(table(mydata$Q4b_SW)/sum(table(mydata$Q4b_SW)))

t5=table(mydata$Q5_SW)
t5
write.table(t5,"Export/tableQ5_SW.csv", sep=";")
as.matrix(table(mydata$Q5_SW)/sum(table(mydata$Q5_SW)))

t6a=table(mydata$Q6a_SW)
t6a
write.table(t6a,"Export/tableQ6a_SW.csv", sep=";")
as.matrix(table(mydata$Q6a_SW)/sum(table(mydata$Q6a_SW)))

t6b=table(mydata$Q6b_SW)
t6b
write.table(t6b,"Export/tableQ6b_SW.csv", sep=";")
as.matrix(table(mydata$Q6b_SW)/sum(table(mydata$Q6b_SW)))

#Repugnant Conclusion
mydata$RC=ifelse(mydata$Q6a_SW=="La naissance de cet animal serait souhaitable" & mydata$Q6b_SW=="Cette politique de réduction des inégalités est souhaitable",1,0)
round(table(mydata$RC)/sum(table(mydata$RC))*100,3)

#Clustering
#Number of possible groups:
length(levels(mydata$Q1_factorVar))*length(levels(mydata$Q2_factorVar))*length(levels(mydata$Q3_factorVar))*length(levels(mydata$Q4a_factorVar))*length(levels(mydata$Q4b_factorVar))*length(levels(mydata$Q5_factorVar))*length(levels(mydata$Q6a_factorVar))*length(levels(mydata$Q6b_factorVar))

diss_matrix <- daisy(mydata[,c("Q1_factorVar","Q2_factorVar","Q3_factorVar",
                               "Q4a_factorVar","Q4b_factorVar","Q5_factorVar",
                               "Q6a_factorVar","Q6b_factorVar"),], metric = "gower", type = list(ordered=1:8))

hclust_result <- hclust(diss_matrix, method = "ward.D2")

resClust=NbClust(diss=diss_matrix,distance=NULL,method = "ward.D2", index = "silhouette",  min.nc=2, max.nc=10)
print(resClust$All.index)
print(resClust$Best.nc)

# Cut dendrogram into k clusters
clusters=cutree(hclust_result, k = 2)

# Add cluster membership to the original data
mydata$Cluster <- clusters
table(clusters)
table(clusters)/sum(table(clusters))

#Look at the differences with k=2 clusters
listVar=c("Q1_factorVar","Q2_factorVar","Q3_factorVar",
          "Q4a_factorVar","Q4b_factorVar","Q5_factorVar",
          "Q6a_factorVar","Q6b_factorVar")
loop_counter=1
for(i in listVar){
  tmp=table(mydata[[i]], mydata$Cluster)
  tmp=sapply(1:2, function(j){round(tmp[,j]/sum(tmp[,j])*100,1)})
  tmp=cbind(tmp,NA)
  for(j in 1:dim(tmp)[1]){
    tmp[j,3]=round(prop.test(n = c(table(clusters)[1], table(clusters)[2]), x = c(table(mydata[[i]], mydata$Cluster)[j,1],  table(mydata[[i]], mydata$Cluster)[j,2]))$p.value,4)
  }
  if(loop_counter==1) savedMat=tmp
  if(loop_counter>1) savedMat=rbind(savedMat,tmp)
  loop_counter=loop_counter+1
}
savedMat

#Linear regression
resLM=lm(Cluster ~ relevel(Q1_factorVar, ref="Equal weight")+ relevel(Q2_factorVar, ref="Direct and indirect")+relevel(Q3_factorVar, ref="Equal weight")+relevel(Q4a_factorVar, ref="Avoid this life")+relevel(Q4b_factorVar, ref="Create this life")+relevel(Q5_factorVar, ref="Equal weight")+relevel(Q6a_factorVar, ref="Desirable birth")+relevel(Q6b_factorVar, ref="Neither"), data = mydata)
sumLM=summary(resLM)

savedMat=cbind(savedMat,NA)
savedMat[1:4,4]=round(sumLM$coefficients[2:5,4],4)
savedMat[6:7,4]=round(sumLM$coefficients[6:7,4],4)
savedMat[9,4]=round(sumLM$coefficients[8,4],4)
savedMat[11:12,4]=round(sumLM$coefficients[9:10,4],4)
savedMat[14,4]=round(sumLM$coefficients[11,4],4)
savedMat[16,4]=round(sumLM$coefficients[12,4],4)
savedMat[18:20,4]=round(sumLM$coefficients[13:15,4],4)
savedMat[22,4]=round(sumLM$coefficients[16,4],4)
savedMat[24,4]=round(sumLM$coefficients[17,4],4)
savedMat

savedMat=cbind(savedMat,NA)
savedMat[1,5]=round(chisq.test(table(mydata$Q1_factorVar,mydata$Cluster))$p.value,4)
savedMat[7,5]=round(chisq.test(table(mydata$Q2_factorVar,mydata$Cluster))$p.value,4)
savedMat[9,5]=round(chisq.test(table(mydata$Q3_factorVar,mydata$Cluster))$p.value,4)
savedMat[12,5]=round(chisq.test(table(mydata$Q4a_factorVar,mydata$Cluster))$p.value,4)
savedMat[14,5]=round(chisq.test(table(mydata$Q4b_factorVar,mydata$Cluster))$p.value,4)
savedMat[16,5]=round(chisq.test(table(mydata$Q5_factorVar,mydata$Cluster))$p.value,4)
savedMat[19,5]=round(chisq.test(table(mydata$Q6a_factorVar,mydata$Cluster))$p.value,4)
savedMat[22,5]=round(chisq.test(table(mydata$Q6b_factorVar,mydata$Cluster))$p.value,4)

lmWithoutQ1=lm(Cluster ~ relevel(Q2_factorVar, ref="Direct and indirect")+relevel(Q3_factorVar, ref="Equal weight")+relevel(Q4a_factorVar, ref="Avoid this life")+relevel(Q4b_factorVar, ref="Create this life")+relevel(Q5_factorVar, ref="Equal weight")+relevel(Q6a_factorVar, ref="Desirable birth")+relevel(Q6b_factorVar, ref="Neither"), data = mydata)
lmWithoutQ2=lm(Cluster ~ relevel(Q1_factorVar, ref="Equal weight")+relevel(Q3_factorVar, ref="Equal weight")+relevel(Q4a_factorVar, ref="Avoid this life")+relevel(Q4b_factorVar, ref="Create this life")+relevel(Q5_factorVar, ref="Equal weight")+relevel(Q6a_factorVar, ref="Desirable birth")+relevel(Q6b_factorVar, ref="Neither"), data = mydata)
lmWithoutQ3=lm(Cluster ~ relevel(Q1_factorVar, ref="Equal weight")+relevel(Q2_factorVar, ref="Direct and indirect")+relevel(Q4a_factorVar, ref="Avoid this life")+relevel(Q4b_factorVar, ref="Create this life")+relevel(Q5_factorVar, ref="Equal weight")+relevel(Q6a_factorVar, ref="Desirable birth")+relevel(Q6b_factorVar, ref="Neither"), data = mydata)
lmWithoutQ4a=lm(Cluster ~ relevel(Q1_factorVar, ref="Equal weight")+relevel(Q2_factorVar, ref="Direct and indirect")+relevel(Q3_factorVar, ref="Equal weight")+relevel(Q4b_factorVar, ref="Create this life")+relevel(Q5_factorVar, ref="Equal weight")+relevel(Q6a_factorVar, ref="Desirable birth")+relevel(Q6b_factorVar, ref="Neither"), data = mydata)
lmWithoutQ4b=lm(Cluster ~ relevel(Q1_factorVar, ref="Equal weight")+relevel(Q2_factorVar, ref="Direct and indirect")+relevel(Q3_factorVar, ref="Equal weight")+relevel(Q4a_factorVar, ref="Avoid this life")+relevel(Q5_factorVar, ref="Equal weight")+relevel(Q6a_factorVar, ref="Desirable birth")+relevel(Q6b_factorVar, ref="Neither"), data = mydata)
lmWithoutQ5=lm(Cluster ~ relevel(Q1_factorVar, ref="Equal weight")+relevel(Q2_factorVar, ref="Direct and indirect")+relevel(Q3_factorVar, ref="Equal weight")+relevel(Q4a_factorVar, ref="Avoid this life")+relevel(Q4b_factorVar, ref="Create this life")+relevel(Q6a_factorVar, ref="Desirable birth")+relevel(Q6b_factorVar, ref="Neither"), data = mydata)
lmWithoutQ6a=lm(Cluster ~ relevel(Q1_factorVar, ref="Equal weight")+relevel(Q2_factorVar, ref="Direct and indirect")+relevel(Q3_factorVar, ref="Equal weight")+relevel(Q4a_factorVar, ref="Avoid this life")+relevel(Q4b_factorVar, ref="Create this life")+relevel(Q5_factorVar, ref="Equal weight")+relevel(Q6b_factorVar, ref="Neither"), data = mydata)
lmWithoutQ6b=lm(Cluster ~ relevel(Q1_factorVar, ref="Equal weight")+relevel(Q2_factorVar, ref="Direct and indirect")+relevel(Q3_factorVar, ref="Equal weight")+relevel(Q4a_factorVar, ref="Avoid this life")+relevel(Q4b_factorVar, ref="Create this life")+relevel(Q5_factorVar, ref="Equal weight")+relevel(Q6a_factorVar, ref="Desirable birth"), data = mydata)

#partial_r2(resLM, covariates = 'relevel(Q1_factorVar, ref = "Equal weight")No weight')

savedMat=cbind(savedMat,NA)
savedMat[1,6]=round(asbio::partial.R2(lmWithoutQ1,resLM),4)
savedMat[7,6]=round(asbio::partial.R2(lmWithoutQ2,resLM),4)
savedMat[9,6]=round(asbio::partial.R2(lmWithoutQ3,resLM),4)
savedMat[12,6]=round(asbio::partial.R2(lmWithoutQ4a,resLM),4)
savedMat[14,6]=round(asbio::partial.R2(lmWithoutQ4b,resLM),4)
savedMat[16,6]=round(asbio::partial.R2(lmWithoutQ5,resLM),4)
savedMat[19,6]=round(asbio::partial.R2(lmWithoutQ6a,resLM),4)
savedMat[22,6]=round(asbio::partial.R2(lmWithoutQ6b,resLM),4)

colnames(savedMat)=c("Share Group 1","Share Group 2", "Pvalue - propTest", "Pvalue - OLS","Chi2-pvalue","PartialR2")
savedMat
write.table(savedMat,"Export/tableClusterAnalysis.csv", sep=";")

#Linear regression with demographics
resLMClust=lm(Cluster~female+age+ageSqrd+relevel(politicalCategories2, ref="Center")+relevel(Income_factor, ref="1000 à 2000 euros")+ABC+pcaFoodDim1+religionNARecoded, data=mydata)
sumLMClust=summary(resLMClust)
sumLMClust
