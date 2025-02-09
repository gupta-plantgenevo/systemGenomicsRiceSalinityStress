####This script will prepare the transcripts for selection analysis -- keep 0 fitness ind, remove outliers, and estimate mean-standadized selection coeff

library(dplyr)
library(stringr)
library(data.table)
library(stats)
library(lmerTest)

###Functions that will be used by the script
#Function to find the outliers that are beyond 3SD from the normalized data
findOutlier <- function(data, cutoff = 3) {
  ## Calculate the sd
  sds <- apply(data, 2, sd, na.rm = TRUE)
  avg <- apply(data, 2, mean, na.rm =TRUE)
  ## Identify the cells with value greater than cutoff * sd (column wise)
  result <- mapply(function(d, m, s) {
    which(abs((d - m)/s) > cutoff)
  }, data, avg, sds)
  result
}

#Function to filter out the outliers that have been identified using the above function (will add NA in place of the outliers)
removeOutlier <- function(data, outliers) {
  result <- mapply(function(d, o) {
    res <- d
    res[o] <- NA
    return(res)
  }, data, outliers)
  return(as.data.frame(result))
}

###Doing it for INDICA -- Wet Fecundity:

##Reading the files:
ind_wet_fitness_all = read.table("ricesalinity/data/Fitness/Indica_normal_fitness.txt",  header=TRUE, sep="\t", check.names = FALSE)
ind_wet = read.table("ricesalinity/data/geneExp/Indica_normal_exp.txt", header=TRUE, sep="\t", check.names = FALSE)
trID = read.table("ricesalinity/data/geneExp/transcriptID_proxyID.txt", header=FALSE, sep="\t", check.names = FALSE)

#I will transpose the exp data to have the transcripts as columns:
geneID = ind_wet$gene_ID
rnaid = colnames(ind_wet)
ind_wet <- as.data.frame(t(ind_wet[,-1]))
colnames(ind_wet) = geneID
ind_wet$RNASeqCode = rnaid[-1]
ind_wet = ind_wet[,c(18142,1:18141)]
write.table(ind_wet, "ricesalinity/data/geneExp/INDICA_Wet_GeneExp.transposed.txt",row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
ind_wet = read.table("ricesalinity/data/geneExp/INDICA_Wet_GeneExp.transposed.txt", header = T, sep="\t", check.names = F)

m = apply(ind_wet[,-1],2,mean, na.rm=TRUE)
sd = apply(ind_wet[,-1],2,sd, na.rm=TRUE)
TrID = names(m)
meansd = data.frame(TrID, m, sd)
names(meansd) = c("TrID", "mean", "sd")
write.table(meansd, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/Indnormal_Transcript_mean_sd.txt", row.names = F, col.names = T, sep="\t", quote = F)

#Standardize the fecundity fitness (by mean fitness)
ind_wet_fitness_filtered = ind_wet_fitness_all
names(ind_wet_fitness_filtered)[2] = "RNASeqCode"; names(ind_wet_fitness_filtered)[12] = "Fecundity"
ind_wet_fitness_filtered$Fecundity = ind_wet_fitness_filtered$Fecundity/(mean(ind_wet_fitness_filtered$Fecundity))
ind_wet_fitness_filtered = ind_wet_fitness_filtered[,c(2,12,10)]

##Removing from each column the values that are outside of 3SD limits
ind_wet_sel = ind_wet[,-1]
ind_wet_sel[ind_wet_sel == 0] <- NA
outliers = findOutlier(ind_wet_sel)
tmp = data.frame(table(lengths(outliers)))
ind_wet_sel = removeOutlier(ind_wet_sel, outliers)

##Now we need to remove transcripts which do not have NA is more than 20 individuals (basically transcripts which had abundance in at least 20 individuals; 410-20 = 390)
cols_0count = as.data.frame(colSums(is.na(ind_wet_sel)))
cols_0count$TranscriptID = rownames(cols_0count)
removecols = cols_0count[which(cols_0count$`colSums(is.na(ind_wet_sel))` > 389),2]           ##remove 69 transcripts with less than 20 individuals with abundance
ind_wet_sel = ind_wet_sel[, !(names(ind_wet_sel) %in% removecols)]
##Left with 18072 transcripts (of 18141)
ind_wet_sel$RNASeqCode = ind_wet$RNASeqCode; 

##Now here I will merge the RNASeq data with fitness data
##NOTE: the extra steps regarding trID will just replace the transcriptIDs with their proxy since the python script acts up when using the complex naming system which is the original transcript ID
##I fixed the above comment, so I think I will just use the original trIDs to avoid further confusion
tmp = merge(ind_wet_fitness_filtered, ind_wet_sel, by=c("RNASeqCode"))
#str(tmp)
total = ncol(tmp)
tmp = tmp[,-1]
data_forSel_analysis_wet_IND = tmp
write.table(data_forSel_analysis_wet_IND, "ricesalinity/data/SelectionAnalysis/InputFiles/IND_normal.w0fecun.nonStd.Selinput.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

###Doing it for INDICA -- Saline Fecundity:
##Reading the files:
ind_salt_fitness_all = read.table("ricesalinity/data/Fitness/Indica_salt_fitness.txt",  header=TRUE, sep="\t", check.names = FALSE)
ind_salt = read.table("ricesalinity/data/geneExp/Indica_salt_exp.txt", header=TRUE, sep="\t", check.names = FALSE)
trID = read.table("ricesalinity/data/geneExp/transcriptID_proxyID.txt", header=FALSE, sep="\t", check.names = FALSE)

#I will transpose the exp data to have the transcripts as columns:
geneID = ind_salt$gene_ID
rnaid = colnames(ind_salt)
ind_salt <- as.data.frame(t(ind_salt[,-1]))
colnames(ind_salt) = geneID
ind_salt$RNASeqCode = rnaid[-1]
ind_salt = ind_salt[,c(18142,1:18141)]
write.table(ind_salt, "ricesalinity/data/geneExp/INDICA_Salt_GeneExp.transposed.txt",row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

m = apply(ind_salt[,-1],2,mean, na.rm=TRUE)
sd = apply(ind_salt[,-1],2,sd, na.rm=TRUE)
TrID = names(m)
tmp = data.frame(TrID, m, sd)
names(tmp) = c("TrID", "mean", "sd")
write.table(tmp, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/IndSalt_Transcript_mean_sd.txt", row.names = F, col.names = T, sep="\t", quote = F)

##The first step is removing individuals with 0 fecundity individuals:
#ind_salt_fitness_filtered = ind_salt_fitness_all %>% filter(`Fecundity (Total Filled Grain Nr)` > 0 )
#str(ind_salt_fitness_filtered)       ##removes 24 individuals with 0 fecundity of a total of 389 individuals; Now has a total of 365 individuals

#Standardize the fecundity fitness
ind_salt_fitness_filtered = ind_salt_fitness_all 
names(ind_salt_fitness_filtered)[2] = "RNASeqCode"; names(ind_salt_fitness_filtered)[12] = "Fecundity"
ind_salt_fitness_filtered$Fecundity = ind_salt_fitness_filtered$Fecundity/(mean(ind_salt_fitness_filtered$Fecundity))
ind_salt_fitness_filtered = ind_salt_fitness_filtered[,c(2,12,10)]

ind_salt_sel = ind_salt[,-1]
ind_salt_sel[ind_salt_sel == 0] <- NA     ##replace 0 values with NA
outliers = findOutlier(ind_salt_sel)
ind_salt_sel = removeOutlier(ind_salt_sel, outliers)

##Now we need to remove transcripts which do not have NA is more than 20 individuals (basically transcripts which had abundance in at least 20 individuals; 389-20 = 369)
cols_0count = as.data.frame(colSums(is.na(ind_salt_sel)))
cols_0count$TranscriptID = rownames(cols_0count)
removecols = cols_0count[which(cols_0count$`colSums(is.na(ind_salt_sel))` > 368),2]           ##remove 156 transcripts with less than 20 individuals with abundance
ind_salt_sel = ind_salt_sel[, !(names(ind_salt_sel) %in% removecols)]
##Left with 18021 transcripts (of 18141)
ind_salt_sel$RNASeqCode = ind_salt$RNASeqCode; 

##Now here I will merge the RNASeq data with fitness data
tmp = merge(ind_salt_fitness_filtered, ind_salt_sel, by=c("RNASeqCode"))
total = ncol(tmp)
tmp = tmp[,-1]
data_forSel_analysis_salt_IND = tmp
write.table(data_forSel_analysis_salt_IND, "ricesalinity/data/SelectionAnalysis/InputFiles/IND_salt.w0fecun.nonStd.Selinput.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")


data_forSel_analysis_wet_IND = read.table("ricesalinity/data/SelectionAnalysis/InputFiles/IND_wet.w0fecun.nonStd.Selinput.txt", header = T, sep = "\t")
data_forSel_analysis_salt_IND = read.table("ricesalinity/data/SelectionAnalysis/InputFiles/IND_salt.w0fecun.nonStd.Selinput.txt", header = T, sep = "\t")
tmp = intersect(colnames(data_forSel_analysis_wet_IND), colnames(data_forSel_analysis_salt_IND))
tmp = tmp[-2]
d1 = data_forSel_analysis_wet_IND[,(colnames(data_forSel_analysis_wet_IND) %in% tmp)]
d2 = data_forSel_analysis_salt_IND[,(colnames(data_forSel_analysis_salt_IND) %in% tmp)]

cov1 = cov2 = c()
for(i in c(2:ncol(d1))){
  cov1 = append(cov1, cov(d1$Fecundity, d1[,i], use="pairwise.complete.obs"))
  cov2 = append(cov2, cov(d2$Fecundity, d2[,i], use="pairwise.complete.obs"))
} 
t.test(cov1, cov2,paired = T)


###Just visualizing one of the transcripts with extreme difference
exp1 = read.table("ricesalinity/data/SelectionAnalysis/InputFiles/IND_salt_fecundity.input.txt", header = T, sep = "\t", check.names = F)
exp2 = read.table("ricesalinity/data/SelectionAnalysis/InputFiles/IND_salt_fecundity.w0fecun.input.txt", header = T, sep = "\t", check.names = F)
sum(is.na(df$col))
na_count1 <- data.frame(sapply(exp1[,-c(1:2)], function(y) sum(length(which(is.na(y))))))
na_count1$TrID = rownames(na_count1)
names(na_count1)[1] = "na_count1"
na_count2 <- data.frame(sapply(exp2[,-c(1:2)], function(y) sum(length(which(is.na(y))))))
na_count2$TrID = rownames(na_count2)
names(na_count2)[1] = "na_count2" 
tmp1 = merge(tmp, na_count1, by="TrID")
tmp = merge(tmp1, na_count2, by="TrID")

tr1 = exp1$`OS08T0500700-01`
linear_model1 = lmer(Fecundity ~ tr1 + (1|Block), data = exp1)
tr1 = seq(-3,3,by=0.1)
block = as.factor(sample(c(1:6), size=length(tr1), replace = TRUE))
predictedCounts1 = predict(linear_model1, newdata = data.frame(tr2 = tr2, Block=block))
fitcurve1 = data.frame(tr1, predictedCounts1)

tr2 = exp2$`OS08T0500700-01`
linear_model2 = lmer(Fecundity ~ tr2 + (1|Block), data = exp2)
tr2 = seq(-3,3,by=0.1)
block = as.factor(sample(c(1:6), size=length(tr2), replace = TRUE))
predictedCounts2 = predict(linear_model2, newdata = data.frame(tr2 = tr2, Block=block))
fitcurve2 = data.frame(tr2, predictedCounts2)

par(mfrow=c(1,2))
x1 = ggplot() + geom_point(data=exp1, aes(x=`OS08T0500700-01`, y=Fecundity), color="#F8766D", size=2, fill="#F8766D") +
  geom_smooth(data = fitcurve1, aes(x=tr1, y=predictedCounts1), color="red", size=2,method = "loess") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

x2 = ggplot() + geom_point(data=exp2, aes(x=`OS08T0500700-01`, y=Fecundity), color="lightblue3", size=2, fill="lightblue3") +
  geom_smooth(data = fitcurve2, aes(x=tr2, y=predictedCounts2), color="lightblue3", size=2,method = "loess") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggarrange(x1,x2)
