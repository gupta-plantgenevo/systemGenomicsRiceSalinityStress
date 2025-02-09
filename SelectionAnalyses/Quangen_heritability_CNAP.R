library(lme4)
library(dplyr)
library(tidyr)
library(MuMIn)
library(forcats)
library(ggplot2)
library(cowplot)
library(gridExtra)

###This script deals with estimating H2 and contional neutrality-antagonistic pleiotropy.

##Loading the transcript data and transposing the matrix to have gene exp as the columns instead of the current version where rows represent each transcript
ind_salt = read.table("ricesalinity/data/geneExp/Indica_salt_exp.txt", header=TRUE, sep="\t", check.names = FALSE)
geneID = ind_salt$gene_ID
rnaid = colnames(ind_salt)
ind_salt <- as.data.frame(t(ind_salt[,-1]))
colnames(ind_salt) = geneID
ind_salt$RNASeqCode = rnaid[-1]
ind_salt = ind_salt[,c(18142,1:18141)]

ind_wet = read.table("ricesalinity/data/geneExp/Indica_normal_exp.txt", header=TRUE, sep="\t", check.names = FALSE)
geneID = ind_wet$gene_ID
rnaid = colnames(ind_wet)
ind_wet <- as.data.frame(t(ind_wet[,-1]))
colnames(ind_wet) = geneID
ind_wet$RNASeqCode = rnaid[-1]
ind_wet = ind_wet[,c(18142,1:18141)]

identical(colnames(ind_wet), colnames(ind_salt))  ##Just to test that the column names are the same (and in the same order)
ind_all = rbind(ind_wet, ind_salt)            ##This has all the transcript counts for all the indica pop (along with a column of RNASeqCode)

##Now I will load file that has the Genotype and Env info along with RNASeqCode and then merge it with the above file (ind_all)
geno_info = read.table("ricesalinity/data/geneExp/INDICA-code_IRGC_Env.txt", header=TRUE, sep="\t")
names(geno_info)[1] = "RNASeqCode"
ind_quangenData = merge(geno_info, ind_all, by="RNASeqCode")
ind_quangenData$IRGC.Nr = as.factor(ind_quangenData$IRGC.Nr)
ind_quangenData$Environment = as.factor(ind_quangenData$Environment)


##Estimating H2
re = 3
varH2_an = read.table("ricesalinity/analysis/transcript_anova_GE/INDtrans_SS_MS_F_P.txt", header=TRUE, sep="\t")
Var_G = varH2_an$SS_G/varH2_an$df_g
Var_E = varH2_an$SS_E/varH2_an$df_e
Var_GE = varH2_an$SS_GE/varH2_an$df_ge
H2_an = (Var_G)/((Var_G)+((Var_GE)/2) +(Var_E/re))
varH2_an$H2 = H2_an
write.table(varH2_an, "ricesalinity/analysis/transcript_anova_GE/INDtrans_SS_MS_F_P_H2_wo0.5.txt", col.names = T, row.names = F, quote = F, sep="\t")



###Antagonistic Pleiotropy vs Conditional Neutrality
sel_coeff_wet = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_normal.w0fecun.nonStd.txt", header=TRUE, sep = "\t")
sel_coeff_salt = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_salt.w0fecun.nonStd.txt", header=TRUE, sep = "\t")
sel_coeff = merge(sel_coeff_wet, sel_coeff_salt, by="TrID")
#This now has 17962 transcripts
sel_coeff_S = sel_coeff[,c(1:2,4,10,12)]
names(sel_coeff_S)[2] = "S_wet"; names(sel_coeff_S)[4] = "S_salt"
names(sel_coeff_S)[3] = "P_wet"; names(sel_coeff_S)[5] = "P_salt"

#Conditional Neutrality
cn = sel_coeff_S[which((sel_coeff_S$P_wet < 0.025 & sel_coeff_S$P_salt > 0.05) | (sel_coeff_S$P_wet > 0.05 & sel_coeff_S$P_salt < 0.025)),]
psig = sel_coeff_S[which(sel_coeff_S$P_wet < 0.05 & sel_coeff_S$P_salt < 0.05),]
samedit_sel = psig[which((psig$S_wet < 0 & psig$S_salt < 0)  | (psig$S_wet > 0 & psig$S_salt > 0)),]
tmp = anti_join(psig, samedit_sel)
ap_wetben = tmp[which(tmp$S_wet  > 0),]
ap_saltben = tmp[which(tmp$S_salt > 0),]


#So we have 1807 CN, 68 AP (55 beneficial in wet and 13 beneficial in salt), 36 beneficial in both, 58 detrimental in both
cond = rep("Antagonistically Pleiotropic (Salt Beneficial)", times = nrow(ap_saltben))
cond = append(cond, rep("Antagonistically Pleiotropic (Wet Beneficial)", times = nrow(ap_wetben)))
t1 = rbind(ap_saltben, ap_wetben)
tmp = samedit_sel[which(samedit_sel$S_wet < 0),]
cond = append(cond, rep("Detrimental", times = nrow(tmp)))
t2 = rbind(t1, tmp)
tmp = samedit_sel[which(samedit_sel$S_wet > 0),]
cond = append(cond, rep("Beneficial", times = nrow(tmp)))
t1 = rbind(t2, tmp)
cond = append(cond, rep("Conditionally Neutral", times = nrow(cn)))
t2 = rbind(t1, cn)

trID = read.table("ricesalinity/data/geneExp/transcriptID_proxyID.txt", header=FALSE, sep="\t", check.names = FALSE)
tmp = subset(trID, !(V1 %in% t2$TrID))
names(tmp)[1] = "TrID"
t1 = merge(tmp, sel_coeff_S, by="TrID", all.x = TRUE)
t1 = t1[,c(1,3:6)]
cond = append(cond, rep("Neutral", times = nrow(t1)))
alldat = rbind(t2,t1)
alldat$Condition = cond
write.table(alldat, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/CNAP_all.infor.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")


n= c(1875, 1875)
x = c(1807, 68)
prop.test(x,n)
