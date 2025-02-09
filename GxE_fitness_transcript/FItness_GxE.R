###This script will make basic plots (eg Fecundity comparison between env) and also test whether there is significant block effects 

library(lmerTest)
library(MuMIn)
library(forcats)
library(dplyr)
library(tidyr)

#INDICA
ind_salt_fitness_all = read.table("ricesalinity/data/Fitness/Indica_salt_fitness.txt",  header=TRUE, sep="\t", check.names = FALSE)
ind_salt_fitness_all = ind_salt_fitness_all %>% separate(`RNAseq Code`, c("RNAInd", "RNAEnv"), remove = FALSE)
ind_wet_fitness_all = read.table("ricesalinity/data/Fitness/Indica_normal_fitness.txt",  header=TRUE, sep="\t", check.names = FALSE)
ind_wet_fitness_all = ind_wet_fitness_all %>% separate(`RNAseq Code`, c("RNAInd", "RNAEnv"), remove = FALSE)


tmp1 = ind_wet_fitness_all[,c(7,14)]
names(tmp1) = c("IRGC", "Fitness_Wet")
tmp1$IRGC = sub(" ", "", tmp1$IRGC)
tmp1$IRGC = as.factor(tmp1$IRGC)
tmp1 = as.data.frame(tmp1 %>% group_by(IRGC)  %>% 
  summarize(Fecundity_Wet = mean(Fitness_Wet)))

tmp2 = ind_salt_fitness_all[,c(7,14)]
names(tmp2) = c("IRGC", "Fitness_Salt")
tmp2$IRGC = sub(" ", "", tmp2$IRGC)
tmp2$IRGC = as.factor(tmp2$IRGC)
tmp2 = as.data.frame(tmp2 %>% group_by(IRGC)  %>% 
                       summarize(Fecundity_Salt = mean(Fitness_Salt)))

plotData = merge(tmp1, tmp2, by="IRGC")
t.test(plotData$Fecundity_Wet, plotData$Fecundity_Salt, paired = T, alternative = "greater")
#(t = 6.0244, df = 129, p-value = 8.288e-09)
t.test(plotData$Fecundity_Wet, plotData$Fecundity_Salt, paired = T)
#(t = 6.0244, df = 129, p-value = 1.658e-08)
dim(plotData[which(plotData$Fecundity_Salt < plotData$Fecundity_Wet),])
prop.test(x=94, n = 130, p=0.5, correct=FALSE)
#p-value = 3.639e-07


####comparison of fecundity in the IND pop between wet and salt env
#merge the files based on the ind part of the RNASeq Code (represents the ind ID)
ind_compare_fecundity = merge(ind_salt_fitness_all, ind_wet_fitness_all, by="RNAInd")
names(ind_compare_fecundity)[14] = "Fecundity_Salt"; names(ind_compare_fecundity)[28] = "Fecundity_Wet"; 
plotData = ind_compare_fecundity[,c(1,14,28)]

plotData = plotData[order(plotData$Fecundity_Wet, decreasing = FALSE),]
plotData$RNAInd = as.factor(plotData$RNAInd)

#We have 383 paired individuals (not doing two way anova since have additional terms like blocks)
#I am updating block numbers since the blocks are distinct in wet and sal t but they are all numbered 1-6 instead of 1-12.
tmp = ind_salt_fitness_all$Block + 6
ind_salt_fitness_all$Block = tmp

tmp = ind_wet_fitness_all[,c(7,4,11,12,13,14)]
tmp2 = ind_salt_fitness_all[,c(7,4,11,12,13,14)]
ind_g.e = rbind(tmp, tmp2)
names(ind_g.e)[1] = "Geno_IRGCNr"; names(ind_g.e)[6] = "Fitness"; names(ind_g.e)[5] = "SubBlock";names(ind_g.e)[3] = "Replicate";
ind_g.e$Geno_IRGCNr = as.factor(ind_g.e$Geno_IRGCNr)
ind_g.e$RNAEnv = as.factor(ind_g.e$RNAEnv)
ind_g.e$Block = as.factor(ind_g.e$Block)
ind_g.e$SubBlock = as.factor(ind_g.e$SubBlock)
ind_g.e$Replicate = as.factor(ind_g.e$Replicate)

#since the blocks do not have both rep and/or env represented in the same block, including block in these models will compound the result
#Two way anova (with Env fixed and Geno random) with F-test 
GEind_aov2 <- aov(Fitness ~ Geno_IRGCNr * RNAEnv, data = ind_g.e) 
GEind_aov2
#Saving the sum of squares
SS_E = 6284012
SS_G = 48646270
SS_GE = 19014967
SS_R = 74858064

n=3; g=130; e=2
dfg=129; dfe=1; dfge=129; dfr=539

MS_G = SS_R + n*e*(SS_G)
MS_E = SS_R + n*(SS_GE) + n*g*(SS_E)
MS_GE = SS_R + n*(SS_GE)
MS_R = SS_R
F_G = MS_G/MS_R
F_GE = MS_GE/MS_R
F_E = MS_E/MS_GE
p_G = pf(F_G, dfg, dfr, lower.tail=F)   #3.050231e-39
p_GE = pf(F_GE, dfge, dfr, lower.tail=F)  #7.091438e-06
p_E = pf(F_E, dfe, dfge, lower.tail=F)  #2.031535e-05

