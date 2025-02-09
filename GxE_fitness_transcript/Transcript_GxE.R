###This script deals with mixed ANOVA and thus idenitifies plastic genes and p-val for Differential Gene Expression
library(dplyr)
library(tidyr)

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

SS_GE = c()
SS_G = c()
SS_E = c()
SS_R = c()
for(i in 4:ncol(ind_quangenData)){
  tmp2 = summary(aov(ind_quangenData[,i] ~ Environment*IRGC.Nr, data = ind_quangenData))
  #Saving the sum of squares
  SS_E = append(SS_E, tmp2[[1]][1, 2])
  SS_G = append(SS_G, tmp2[[1]][2, 2])
  SS_GE = append(SS_GE, tmp2[[1]][3, 2])
  SS_R = append(SS_R, tmp2[[1]][4, 2])
}
tr_aov = data.frame(SS_G,SS_E,SS_GE,SS_R)
tr_aov$TrID = colnames(ind_quangenData)[-c(1:4)]
tr_aov = tr_aov[,c(5,1:4)]
write.table(tr_aov, "ricesalinity/analysis/transcript_anova_GE/IND_trans_SS.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


#For IND
sst = tr_aov
n=3; g=130; e=2
dfg=129; dfe=1; dfge=129; dfr=539

sst$MS_G = sst$SS_R + n*e*(sst$SS_G)
sst$MS_E = sst$SS_R + n*(sst$SS_GE) + n*g*(sst$SS_E)
sst$MS_GE = sst$SS_R + n*(sst$SS_GE)
sst$MS_R = sst$SS_R

sst$F_G = sst$MS_G/sst$MS_R
sst$F_GE = sst$MS_GE/sst$MS_R
sst$F_E = sst$MS_E/sst$MS_GE

sst$p_G = pf(sst$F_G, dfg, dfr, lower.tail=F)
sst$p_GE = pf(sst$F_GE, dfge, dfr, lower.tail=F)
sst$p_E = pf(sst$F_E, dfe, dfge, lower.tail=F)

sst$df_g = rep(dfg, times = nrow(sst))
sst$df_e = rep(dfe, times = nrow(sst))
sst$df_ge = rep(dfge, times = nrow(sst))
sst$df_r = rep(dfr, times = nrow(sst))

sst$fdr_g = p.adjust(sst$p_G, method="fdr")
sst$fdr_e = p.adjust(sst$p_E, method="fdr")
sst$fdr_ge = p.adjust(sst$p_GE, method="fdr")
sst = sst[,c(1,16,2,6,10,13,20,17,3,7,12,15,21,18,4,8,11,14,22,19,5,9)]

write.table(sst, "ricesalinity/analysis/transcript_anova_GE/INDtrans_SS_MS_F_P.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

