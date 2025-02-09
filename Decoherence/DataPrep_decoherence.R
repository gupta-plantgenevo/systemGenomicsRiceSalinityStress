###Prepping data for comparison of IND 
library(dplyr)
library(tidyr)
library(data.table)

geno1 = read.table("ricesalinity/data/geneExp/INDICA_Wet_GeneExp.transposed.txt", header = T, sep="\t", check.names = F)
geno2 = read.table("ricesalinity/data/geneExp/INDICA_Salt_GeneExp.transposed.txt", header = T, sep="\t", check.names = F)
sel_coeff_salt = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_salt.w0fecun.nonStd.txt", header = T, sep="\t")
sel_coeff_wet = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_wet.w0fecun.nonStd.txt", header = T, sep="\t")

t1 = sel_coeff_wet[which(abs(sel_coeff_wet$S) > 0.1),] #3953 transcripts
t2 = sel_coeff_salt[which(abs(sel_coeff_salt$S) > 0.1),] #4352 transcripts
t = unique(union(t1$TrID, t2$TrID))   #6804 transcripts

#Find all the transcripts that have over 50% missing data in either dataset
tmp = geno1
tmp[tmp == 0] <- NA 
cols_0count = as.data.frame(colSums(is.na(tmp[,-1])))
cols_0count$TrID = rownames(cols_0count)
rmtr1 = cols_0count[which(cols_0count$`colSums(is.na(tmp[, -1]))` > nrow(geno1)/2),]$TrID
#rmtr1 contains a list of transcripts that have over 50% NAs in wet and thus should be removed.

tmp = geno2
tmp[tmp == 0] <- NA 
cols_0count = as.data.frame(colSums(is.na(tmp[,-1])))
cols_0count$TrID = rownames(cols_0count)
rmtr2 = cols_0count[which(cols_0count$`colSums(is.na(tmp[, -1]))` > nrow(geno2)/2),]$TrID
#rmtr2 contains a list of transcripts that have over 50% NAs in salt and thus should be removed.

rmtr = unique(union(rmtr1, rmtr2)) #all tr with over 50% NAs in either
keeptr = setdiff(t, rmtr)   #all tr with |S| > 0.1 and less than 50% NAs
write.table(keeptr, "ricesalinity/data/module_pathway_diff/decoherence/ind_salt_wet/trwith_above01S_lessthan50percentNA.eitherEnv.list", col.names = TRUE, row.names = F, sep="\t", quote = F)

tmp1 = geno1[,(colnames(geno1) %in% keeptr)]
tmp2 = geno2[,(colnames(geno2) %in% keeptr)]
tmpexp = rbind(tmp1, tmp2)  #2051 transcripts 
#NOTE: these havent been filtered for 50% expression values since if we do that the number of transcripts reduce to less than 200
write.table(tmpexp, "ricesalinity/data/module_pathway_diff/decoherence/ind_salt_wet/trwith_above01S_lessthan50percentNA.eitherEnv.exp.txt", col.names = TRUE, row.names = F, sep="\t", quote = F)



#Pheno File (remains the same irrespective of the number of transcripts)
p1 = read.table("ricesalinity/data/Fitness/Indica_salt_fitness.txt", header = T, sep = "\t")
p2 = read.table("ricesalinity/data/Fitness/Indica_normal_fitness.txt", header = T, sep = "\t")
names(p1)[2] = "RNASeqCode"; names(p2)[2] = "RNASeqCode";
tmp = merge(expmat[,c(1,2053)], p1[,c(2,5)], by="RNASeqCode")
tmp2 = merge(expmat[,c(1,2053)], p2[,c(2,5)], by="RNASeqCode")
names(tmp2)[3] = "IRGCNr"
tmp3 = rbind(tmp,tmp2)
pheno = tmp3[,c(2,3)]
write.table(pheno, "ricesalinity/data/module_pathway_diff/decoherence/ind_salt_wet/pheno_forDeco.in.txt", col.names = TRUE, row.names = F, sep="\t", quote = F)

