###For running deco on terminal  
library(dplyr)
library(tidyr)
library(data.table)

tmpexp1 = read.table("ricesalinity/data/module_pathway_diff/decoherence/ind_salt_wet/trwith_above01S_lessthan50percentNA.eitherEnv.exp.txt", header = T, sep="\t", check.names = F)
tmp = as.data.frame(t(tmpexp))
tmpexpin = tmp

pheno = read.table("ricesalinity/data/module_pathway_diff/decoherence/ind_salt_wet/pheno_forDeco.in.txt", header = T, sep="\t", check.names = F)
pheno$IRGCNr = as.factor(pheno$IRGCNr)

input_matrix = tmpexpin
predictor_DF = pheno
predictor_column = 1 
covariates_column_start = 2
covariates_column_end = 2
class1_name = "Wet"
class2_name = "Salt"

# get all possible pairwise combinations of genes
genes=dim(input_matrix)[1]
e_all<-input_matrix
predictor<-predictor_DF[,predictor_column]
pairs<-as.data.frame(t(combn( 1:genes, 2)))
write.table(pairs, "ricesalinity/data/module_pathway_diff/decoherence/ind_salt_wet/trwith_above01S_lessthan50percentNA.pairs.txt", row.names = F, sep="\t", col.names = T, quote = F)

print('scaling data within the classes to be compared')
tmp<-matrix(nrow=dim(e_all)[1],ncol=dim(e_all)[2])
for (i in c(1:dim(tmp)[1])) {
  tmp[i,which(predictor==1)]<-scale( t(e_all[i,which(predictor==1)] ))
  tmp[i,which(predictor==0)]<-scale( t(e_all[i,which(predictor==0)] ))
} 
e_all_norm<-tmp


#Model with covrariates
print('testing for differences in correlation')
model<-model.matrix(~predictor_DF[,predictor_column] + as.matrix(predictor_DF[,covariates_column_start:covariates_column_end]) )

results <- as.data.frame(t(Reduce(cbind,lapply(1:dim(pairs)[1], function(i) {
  tmp1<-(e_all_norm[pairs$V1[i],])
  tmp2<-(e_all_norm[pairs$V2[i],])
  return( c( cor.test( e_all_norm[pairs$V1[i],which(predictor==0)] , e_all_norm[pairs$V2[i],which(predictor==0)] ,method='spearman', )$estimate, cor.test( e_all_norm[pairs$V1[i],which(predictor==1)] , e_all_norm[pairs$V2[i],which(predictor==1)] ,method='spearman')$estimate, summary(lm( scale(tmp1*tmp2) ~ model))$coefficients[2,4]))
}))))
names(results)<-c(paste('cor',class1_name,sep='_'),paste('cor',class2_name,sep='_'),'p_value')

results$tr1 = colnames(tmpexp)[pairs$V1]
results$tr2 = colnames(tmpexp)[pairs$V2]

library(qvalue)
require(sqldf)
results$qvalue<-qvalue(results$p_value)$qvalues

write.table(results, "ricesalinity/data/module_pathway_diff/decoherence/ind_salt_wet/deco_trwith_above01S_lessthan50percentNA.txt", row.names = F, col.names = T, sep="\t")



