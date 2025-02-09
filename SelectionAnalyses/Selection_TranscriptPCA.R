library(dplyr)
library(stringr)
library(factoextra)
library(lmerTest)

selestimate_fecun <- function(fecun){
  S = c()
  C = c()
  Pval_S = c()
  Serr = c()
  Pval_C = c()
  tval_S = c()
  tval_C = c()
  Cerr = c()
  ##This loop will fit a linear and quadratic model and then save the S and C values (regression coeff) and the respective P values in separate vectors
  ##Takes a while to run
  for(i in 3:ncol(fecun)){
    tr = fecun[,i]
    linear_model = lmer(Fecundity ~ tr + (1|Block), data = fecun)
    Serr = append(Serr, summary(linear_model)$coefficients[2,2])
    S = append(S, summary(linear_model)$coefficients[2,1])
    tval_S = append(tval_S, summary(linear_model)$coefficients[2,4])
    Pval_S = append(Pval_S, summary(linear_model)$coefficients[2,5])
    tr2 = tr^2
    quad_model = lmer(Fecundity ~ tr + tr2 + (1|Block), data = fecun)
    Cerr = append(Cerr, summary(quad_model)$coefficients[3,2])
    C = append(C, 2*(summary(quad_model)$coefficients[3,1]))
    tval_C = append(tval_C, (summary(quad_model)$coefficients[3,4]))
    Pval_C = append(Pval_C, summary(quad_model)$coefficients[3,5])
  }
  output = data.frame(S,Serr,tval_S,Pval_S,C,Cerr,tval_C,Pval_C)
  return(output)
}


##Reading the files:
fitness_all = read.table("ricesalinity/data/Fitness/Indica_normal_fitness.txt",  header=TRUE, sep="\t", check.names = FALSE)
trmat = read.table("ricesalinity/data/geneExp/INDICA_Wet_GeneExp.transposed.txt", header=TRUE, sep="\t", check.names = FALSE)
trID = read.table("ricesalinity/data/geneExp/transcriptID_proxyID.txt", header=FALSE, sep="\t", check.names = FALSE)

##The first step is removing individuals with 0 fecundity individuals:
fitness_filtered = fitness_all %>% filter(`Fecundity (Total Filled Grain Nr)` > 0 )
str(fitness_filtered)       ##removes 26 individuals with 0 fecundity; Has a total of 384 individuals
fitness_filtered = fitness_all
names(fitness_filtered)[2] = "RNASeqCode"; names(fitness_filtered)[12] = "Fecundity"

#Standardize the fecundity fitness
fitness_filtered$Fecundity = fitness_filtered$Fecundity/(mean(fitness_filtered$Fecundity))
fitness_filtered = fitness_filtered[,c(2,12,10)]

##We then need to remove individuals with 0 fecundity from the RNASeq DataSet as well
trmat_filt = trmat %>% filter(RNASeqCode %in% fitness_filtered$RNASeqCode)

####PCA Analysis Data Prep
pcin_wet = merge(fitness_filtered, trmat_filt, by=c("RNASeqCode"))
#This is a file which has ind with 0 fecundity removed and also  the fitness has been standardized.
#Now will do PCA on the transcripts with center and scale on 

###Trying the PC on ind_wet for multivariate selection analysis
pc_wet_in = prcomp(pcin_wet[,-c(1:3)])
#summary(pc_wet_in) will also give you the variance explained (but for all)
eig.val <- get_eigenvalue(pc_wet_in)
eig.val[1:15,]      #Just print for first 15 PCs, will just check to see which one is above 0.5% variance explained. In this case first 11 PCs.
write.table(eig.val[1:15,], "ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/IND_wet_variance.pca.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
eig.load = as.data.frame(pc_wet_in$rotation[,c(1:11)])
str(eig.load)
eig.load$TransciptID = trID$V1
eig.load = eig.load[,c(12,1:11)]
write.table(eig.load, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/IND_Wet_PCA_Loadings.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

eig.ind = as.data.frame(pc_wet_in$x[,c(1:11)])
eig.ind$Fecundity = pcin_wet$Fecundity
eig.ind$Block = pcin_wet$Block
eig.ind = eig.ind[,c(12:13,1:11)]
write.table(eig.ind, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/IND_Wet_PCA_Sel.input.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
eig.ind$RNASeqCode = pcin_wet$RNASeqCode
eig.ind = eig.ind[,c(14,3:13)]
write.table(eig.ind, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/IND_Wet_PCA.ind.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

#Running Selection
fecundity = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/IND_Wet_PCA_Sel.input.txt", header=TRUE, sep="\t", check.names = FALSE)
TrID = colnames(fecundity)
TrID = TrID[-c(1,2)]
fecundity$Block = as.factor(fecundity$Block)
#tmp = sel_coeff
sel_coeff = selestimate_fecun(fecundity)
sel_coeff$TrID = TrID
sel_coeff = sel_coeff[,c(9,1:8)]
#str(sel_coeff)
s_wet_correctedPval = p.adjust(sel_coeff$Pval_S, method = "bonferroni")
sel_coeff$Pval_S.corrected = s_wet_correctedPval
c_wet_correctedPval = p.adjust(sel_coeff$Pval_C, method = "bonferroni")
sel_coeff$Pval_C.corrected = c_wet_correctedPval
write.table(sel_coeff, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/PCA_SelCoeff_Indica_Wet_Fecundity.block.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

##Beta
linear_model = lmer(Fecundity ~ PC1 + PC2+ PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + (1|Block), data = fecundity)
summary(linear_model)$coefficients
wet_beta = data.frame(summary(linear_model)$coefficients)[-1,]
wet_beta$Pheno = rownames(wet_beta)
names(wet_beta) = c( "beta", "SE_b", "df_b", "tval_b", "Pval_b", "Pheno")
wet_beta = wet_beta[,c(6,1:5)]


#gamma
t1 = (fecundity$PC1)^2; t2 = (fecundity$PC2)^2; t3 = (fecundity$PC3)^2; t4 = (fecundity$PC4)^2; t5 = (fecundity$PC5)^2; t6 = (fecundity$PC6)^2
t7 = (fecundity$PC7)^2; t8 = (fecundity$PC8)^2; t9 = (fecundity$PC9)^2; t10 = (fecundity$PC10)^2; t11 = (fecundity$PC11)^2
quad_model = lmer(Fecundity ~ PC1 + t1+ PC2+t2+ PC3 +t3+ PC4 +t4+ PC5 +t5+ PC6 +t6+ PC7 +t7+ PC8 +t8+ PC9 +t9+ PC10 +t10+ PC11 +t11+ (1|Block), data = fecundity)
summary(quad_model)$coefficients
wet_gamma = data.frame(summary(quad_model)$coefficients)
t = wet_gamma[-c(1:2,4,6,8,10,12,14,16,18,20,22),]
t$Pheno = wet_beta$Pheno
wet_gamma = t[,c(6,1:5)]
names(wet_gamma) = c("Pheno", "gamma", "SE_g", "df_g", "tval_g", "Pval_g")

wet_beta_gamma =  merge(wet_beta, wet_gamma, by="Pheno")
write.table(wet_beta_gamma, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/PCA_beta_gamma_Indica_Wet_Fecundity.block.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")


#IND_Salt
fitness_all = read.table("ricesalinity/data/Fitness/Indica_salt_fitness.txt",  header=TRUE, sep="\t", check.names = FALSE)
trmat = read.table("ricesalinity/data/geneExp/INDICA_Salt_GeneExp.transposed.txt", header=TRUE, sep="\t", check.names = FALSE)

#fitness_filtered = fitness_all %>% filter(`Fecundity (Total Filled Grain Nr)` > 0 )
fitness_filtered = fitness_all
str(fitness_filtered)
names(fitness_filtered)[2] = "RNASeqCode"; names(fitness_filtered)[12] = "Fecundity"

#Standardize the fecundity fitness
fitness_filtered$Fecundity = fitness_filtered$Fecundity/(mean(fitness_filtered$Fecundity))
fitness_filtered = fitness_filtered[,c(2,12,10)]
trmat_filt = trmat %>% filter(RNASeqCode %in% fitness_filtered$RNASeqCode)

####PCA Analysis Data Prep
pcin_wet = merge(fitness_filtered, trmat_filt, by=c("RNASeqCode"))
#summary(pc_wet_in) will also give you the variance explained (but for all)
pc_wet_in = prcomp(pcin_wet[,-c(1:3)])
eig.val <- get_eigenvalue(pc_wet_in)
eig.val[1:15,]      #Just print for first 15 PCs, will just check to see which one is above 0.5% variance explained. In this case first 12 PCs.
write.table(eig.val[1:15,], "ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/IND_salt_variance.pca.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
eig.load = as.data.frame(pc_wet_in$rotation[,c(1:12)])
str(eig.load)
eig.load$TransciptID = trID
eig.load = eig.load[,c(13,1:12)]
write.table(eig.load, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/IND_Salt_PCA_Loadings.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

eig.ind = as.data.frame(pc_wet_in$x[,c(1:12)])
eig.ind$Fecundity = pcin_wet$Fecundity
eig.ind$Block = pcin_wet$Block
eig.ind = eig.ind[,c(13:14,1:12)]
write.table(eig.ind, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/IND_Salt_PCA_Sel.input.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
eig.ind$RNASeqCode = pcin_wet$RNASeqCode
eig.ind = eig.ind[,c(15,3:14)]
write.table(eig.ind, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/IND_Salt_PCA.ind.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

#Running Selection
fecundity = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/IND_Salt_PCA_Sel.input.txt", header=TRUE, sep="\t", check.names = FALSE)
TrID = colnames(fecundity)
TrID = TrID[-c(1,2)]
fecundity$Block = as.factor(fecundity$Block)
#tmp = sel_coeff
sel_coeff = selestimate_fecun(fecundity)
sel_coeff$TrID = TrID
sel_coeff = sel_coeff[,c(9,1:8)]
#str(sel_coeff)
s_wet_correctedPval = p.adjust(sel_coeff$Pval_S, method = "bonferroni")
sel_coeff$Pval_S.corrected = s_wet_correctedPval
c_wet_correctedPval = p.adjust(sel_coeff$Pval_C, method = "bonferroni")
sel_coeff$Pval_C.corrected = c_wet_correctedPval
write.table(sel_coeff, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/PCA_SelCoeff_Indica_Salt_Fecundity.block.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

#linear_model = lmer(Fecundity ~ PC1 + PC2+ PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + (1|Block), data = fecundity)
linear_model = lmer(Fecundity ~ PC1 + PC2+ PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + (1|Block), data = fecundity)
summary(linear_model)$coefficients
wet_beta = data.frame(summary(linear_model)$coefficients)[-1,]
wet_beta$Pheno = rownames(wet_beta)
names(wet_beta) = c( "beta", "SE_b", "df_b", "tval_b", "Pval_b", "Pheno")
wet_beta = wet_beta[,c(6,1:5)]

t1 = (fecundity$PC1)^2; t2 = (fecundity$PC2)^2; t3 = (fecundity$PC3)^2; t4 = (fecundity$PC4)^2; t5 = (fecundity$PC5)^2; t6 = (fecundity$PC6)^2; 
t7 = (fecundity$PC7)^2; t8 = (fecundity$PC8)^2; t9 = (fecundity$PC9)^2; t10 = (fecundity$PC10)^2; t11 = (fecundity$PC11)^2; t12 = (fecundity$PC12)^2;
quad_model = lmer(Fecundity ~ PC1 + t1+ PC2+t2+ PC3 +t3+ PC4 +t4+ PC5 +t5+ PC6 +t6+ PC7 +t7+ PC8 +t8+ PC9 +t9+ PC10 +t10+ PC11 +t11+ PC12 +t12+ (1|Block), data = fecundity)
summary(quad_model)$coefficients
wet_gamma = data.frame(summary(quad_model)$coefficients)
t = wet_gamma[-c(1:2,4,6,8,10,12,14,16,18,20,22,24),]
t$Pheno = wet_beta$Pheno
wet_gamma = t[,c(6,1:5)]
names(wet_gamma) = c("Pheno", "gamma", "SE_g", "df_g", "tval_g", "Pval_g")

wet_beta_gamma =  merge(wet_beta, wet_gamma, by="Pheno")
write.table(wet_beta_gamma, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/PCA/PCA_beta_gamma_Indica_Salt_Fecundity.block.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")


#Checking correlation
pc_ind1 = read.table("ricesalinity/data/SelectionAnalysis/InputFiles/IND_Wet_PCA_Sel.input.txt", header = T, sep = "\t")
load1 = read.table("ricesalinity/data/SelectionAnalysis/PCA/IND_Wet_PCA_Loadings.txt", header = T, sep="\t")
tail1 = read.table("ricesalinity/data/SelectionAnalysis/PCA/PCA_IND_wet_sigTail.txt", header = T, sep="\t")

cor(pc_ind1$Fecundity, pc_ind1$PC2)
tmp1 = load1[which(load1$TransciptID %in% tail1$PC2),]$PC2
cor(pc_ind1$Fecundity, tmp1)


#PCtails determining selectional direction
pctails_wet = read.table("ricesalinity/data/SelectionAnalysis/PCA/PCA_IND_wet_sigTail.txt", header = T, sep = "\t")
pctails_salt = read.table("ricesalinity/data/SelectionAnalysis/PCA/PCA_IND_salt_sigTail.2.txt", header = T, sep = "\t")
sel_coeff_wet = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_normal.w0fecun.nonStd.txt", header=TRUE, sep = "\t")
sel_coeff_salt = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_salt.w0fecun.nonStd.txt", header=TRUE, sep = "\t")

for (i in 1:ncol(pctails_wet)){
  a = (sel_coeff_wet[which(sel_coeff_wet$TrID %in% pctails_wet[,i]),]$S)
  x = length(a[which(a > 0)])/length(pctails_wet[,i])
  print(x)
}

for (i in 1:ncol(pctails_salt)){
  a = (sel_coeff_wet[which(sel_coeff_wet$TrID %in% pctails_salt[,i]),]$S)
  x = length(a[which(a > 0)])/length(pctails_salt[,i])
  print(x)
}

