setwd("ricesalinity/data/eQTL/ind_data/fdr_trans1Mb/")
snpfile = read.csv("ricesalinity/data/eQTL/ind_data/snps.location.txt", header = T, sep="\t")
maf = read.table("../SNPid_maf.txt", header = F, sep="\t")
names(maf)[1] = "SNP_id"; names(maf)[2] = "maf"
snpfile = cbind(snpfile, maf[,2])
names(snpfile)[4] = "maf"

cis.salt.2017.lead = read.table("cis-eQTLs.salt.leadSNPs.txt", header = T, sep="\t")
trans.1mb.salt.2017.lead = read.table("trans-eQTLs.salt.1mb.leadSNPs.txt", header = T, sep="\t")
salt_gene = read.table("salt_cis_trans_compensating_reinfo.txt", header = T, sep="\t")

t1 = cis.salt.2017.lead[which(cis.salt.2017.lead$gene %in% salt_gene[which(salt_gene$Group == "Compensating"),]$TrID),]
t2 = trans.1mb.salt.2017.lead[which(trans.1mb.salt.2017.lead$gene %in% salt_gene[which(salt_gene$Group == "Compensating"),]$TrID),]
salt_comp = merge(t1[,c(1,5,6)], t2[,c(1,5,6)], by="gene")
names(salt_comp) = c("gene", "cis-snp", "cis-Chr", "trans-snp", "trans-Chr")
tmp = merge(salt_comp, maf, by.x='cis-snp', by.y='SNP_id'); names(tmp)[6] = 'maf_cis'
tmp = merge(tmp, maf, by.x='trans-snp', by.y='SNP_id'); names(tmp)[7] = 'maf_trans'
tmp1 = tmp[,c(3,2,4,6,1,5,7)]
salt_comp = tmp1
dim(salt_comp)  
dim(salt_comp[which(salt_comp$`cis-Chr` == salt_comp$`trans-Chr`),])  

#I might remove the same Chr for the ease of randomization, and just estimate ILD
tmp1 = salt_comp[which(salt_comp$`cis-Chr` != salt_comp$`trans-Chr`),]  
tmp1$`cis-snp` = ifelse(tmp1$`cis-Chr` < 10, paste("chr0", tmp1$`cis-snp`, sep=""), paste("chr", tmp1$`cis-snp`, sep=""))
tmp1$`trans-snp` = ifelse(tmp1$`trans-Chr` < 10, paste("chr0", tmp1$`trans-snp`, sep=""), paste("chr", tmp1$`trans-snp`, sep=""))
write.table(tmp1, "salt_compensating_diifChr_snps.txt", quote = F, row.names = F, col.names = T, sep = "\t")
salt_comp$`cis-snp` = ifelse(salt_comp$`cis-Chr` < 10, paste("chr0", salt_comp$`cis-snp`, sep=""), paste("chr", salt_comp$`cis-snp`, sep=""))
salt_comp$`trans-snp` = ifelse(salt_comp$`trans-Chr` < 10, paste("chr0", salt_comp$`trans-snp`, sep=""), paste("chr", salt_comp$`trans-snp`, sep=""))
write.table(salt_comp, "salt_compensating_all_snps.txt", quote = F, row.names = F, col.names = T, sep = "\t")


#estimating the min distance between the same chr SNPs
tmp1 = merge(t1[,c(1,5,6,7)], t2[,c(1,5,6,7)], by="gene")
names(tmp1) = c("gene", "cis-snp", "cis-Chr","cis-Pos", "trans-snp", "trans-Chr", "trans-Pos")
tmp1 = tmp1[which(tmp1$`cis-Chr` == tmp1$`trans-Chr`),]
tmp1$snpDist = abs(tmp1$`trans-Pos` - tmp1$`cis-Pos`)
min(tmp1$snpDist)  #967011

library(dplyr)
g1 = sample_n(snpfile, 200000)
g2 = sample_n(snpfile, 200000)
f1 = cbind(g1,g2)
names(f1) = c("snp_id1", "Chr1", "Pos1", "maf1", "snp_id2", "Chr2", "Pos2", "maf2")
f2 = f1[which(f1$Chr1 != f1$Chr2),]
t1 = salt_comp[which(salt_comp$`cis-Chr` != salt_comp$`trans-Chr`),]
#h1.1 = sample_n(f2, 1434)
h1.1 = f2[FALSE,]
for(i in c(1:nrow(t1))){
  tmp = f2[which((t1[i,]$maf_cis-0.05<f2$maf1 & f2$maf1<t1[i,]$maf_cis+0.05) & (t1[i,]$maf_trans-0.05 < f2$maf2 &  f2$maf2 < t1[i,]$maf_trans+0.05)),]
  h1.1[i,] = sample_n(tmp, 1)
}
f3 = f1[which(f1$Chr1 == f1$Chr2),]
f3$Dist = abs(f3$Pos2 - f3$Pos1)
t2 = salt_comp[which(salt_comp$`cis-Chr` == salt_comp$`trans-Chr`),]
h1.2 = f3[FALSE,]
for(i in c(1:nrow(t2))){
  tmp = f3[which(f3$maf1 > t2[i,]$maf_cis-0.05 & f3$maf1 < t2[i,]$maf_cis+0.05 & f3$maf2> t2[i,]$maf_trans-0.05 & f3$maf2< t2[i,]$maf_trans+0.05 & f3$Dist > 967011),]
  h1.2[i,] = sample_n(tmp, 1)
}
h1 = rbind(h1.1, h1.2[,-9])
h1$snp_id1 = ifelse(h1$Chr1 < 10, paste("chr0", h1$snp_id1, sep=""), paste("chr", h1$snp_id1, sep=""))
h1$snp_id2 = ifelse(h1$Chr2 < 10, paste("chr0", h1$snp_id2, sep=""), paste("chr", h1$snp_id2, sep=""))
write.table(h1, "salt-comp_randSNPpair_all.5.txt", row.names = F, col.names = F, quote = F, sep = "\t")




t1 = cis.salt.2017.lead[which(cis.salt.2017.lead$gene %in% salt_gene[which(salt_gene$Group == "Reinforcing"),]$TrID),]
t2 = trans.1mb.salt.2017.lead[which(trans.1mb.salt.2017.lead$gene %in% salt_gene[which(salt_gene$Group == "Reinforcing"),]$TrID),]
salt_re = merge(t1[,c(1,5,6)], t2[,c(1,5,6)], by="gene")
names(salt_re) = c("gene", "cis-snp", "cis-Chr", "trans-snp", "trans-Chr")
tmp = merge(salt_re, maf, by.x='cis-snp', by.y='SNP_id'); names(tmp)[6] = 'maf_cis'
tmp = merge(tmp, maf, by.x='trans-snp', by.y='SNP_id'); names(tmp)[7] = 'maf_trans'
tmp1 = tmp[,c(3,2,4,6,1,5,7)]
salt_re = tmp1
dim(salt_re)  #2107 SNP pairs
dim(salt_re[which(salt_re$`cis-Chr` == salt_re$`trans-Chr`),])  #394
#I might remove the same Chr for the ease of randomization, and just estimate ILD
tmp1 = salt_re[which(salt_re$`cis-Chr` != salt_re$`trans-Chr`),]  #1713
tmp1$`cis-snp` = ifelse(tmp1$`cis-Chr` < 10, paste("chr0", tmp1$`cis-snp`, sep=""), paste("chr", tmp1$`cis-snp`, sep=""))
tmp1$`trans-snp` = ifelse(tmp1$`trans-Chr` < 10, paste("chr0", tmp1$`trans-snp`, sep=""), paste("chr", tmp1$`trans-snp`, sep=""))
write.table(tmp1, "salt_reinforcing_diifChr_snps.txt", quote = F, row.names = F, col.names = T, sep = "\t")
salt_re$`cis-snp` = ifelse(salt_re$`cis-Chr` < 10, paste("chr0", salt_re$`cis-snp`, sep=""), paste("chr", salt_re$`cis-snp`, sep=""))
salt_re$`trans-snp` = ifelse(salt_re$`trans-Chr` < 10, paste("chr0", salt_re$`trans-snp`, sep=""), paste("chr", salt_re$`trans-snp`, sep=""))
write.table(salt_re, "salt_reinforcing_all_snps.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#estimating the min distance between the same chr SNPs
tmp1 = merge(t1[,c(1,5,6,7)], t2[,c(1,5,6,7)], by="gene")
names(tmp1) = c("gene", "cis-snp", "cis-Chr","cis-Pos", "trans-snp", "trans-Chr", "trans-Pos")
tmp1 = tmp1[which(tmp1$`cis-Chr` == tmp1$`trans-Chr`),]
tmp1$snpDist = abs(tmp1$`trans-Pos` - tmp1$`cis-Pos`)
min(tmp1$snpDist)  #909961

g1 = sample_n(snpfile, 200000)
g2 = sample_n(snpfile, 200000)
f1 = cbind(g1,g2)
names(f1) = c("snp_id1", "Chr1", "Pos1", "maf1", "snp_id2", "Chr2", "Pos2", "maf2")
f2 = f1[which(f1$Chr1 != f1$Chr2),]
t1 = salt_re[which(salt_re$`cis-Chr` != salt_re$`trans-Chr`),]
#h1.1 = sample_n(f2, 1434)
h1.1 = f2[FALSE,]
for(i in c(1:nrow(t1))){
  tmp = f2[which((t1[i,]$maf_cis-0.05<f2$maf1 & f2$maf1<t1[i,]$maf_cis+0.05) & (t1[i,]$maf_trans-0.05 < f2$maf2 &  f2$maf2 < t1[i,]$maf_trans+0.05)),]
  h1.1[i,] = sample_n(tmp, 1)
}
f3 = f1[which(f1$Chr1 == f1$Chr2),]
f3$Dist = abs(f3$Pos2 - f3$Pos1)
t2 = salt_re[which(salt_re$`cis-Chr` == salt_re$`trans-Chr`),]
h1.2 = f3[FALSE,]
for(i in c(1:nrow(t2))){
  tmp = f3[which(f3$maf1 > t2[i,]$maf_cis-0.05 & f3$maf1 < t2[i,]$maf_cis+0.05 & f3$maf2> t2[i,]$maf_trans-0.05 & f3$maf2< t2[i,]$maf_trans+0.05 & f3$Dist > 909961),]
  h1.2[i,] = sample_n(tmp, 1)
}
h1 = rbind(h1.1, h1.2[,-9])
h1$snp_id1 = ifelse(h1$Chr1 < 10, paste("chr0", h1$snp_id1, sep=""), paste("chr", h1$snp_id1, sep=""))
h1$snp_id2 = ifelse(h1$Chr2 < 10, paste("chr0", h1$snp_id2, sep=""), paste("chr", h1$snp_id2, sep=""))
write.table(h1, "salt-reinfo_randSNPpair_all.5.txt", row.names = F, col.names = F, quote = F, sep = "\t")
rm(g1, g2, f1, f2, f3, h1.1, h1.2, h1)


library(tidyr)
library(stringr)
###Need to edit the output file rom plink
h1 = read.table("filesforLD/all_leadSNP/salt-reinfo_randSNPpair_all.5.ld", header = F, sep = "\t")
tmp = h1 %>% separate(V3, c("V3", "V4"), sep="\n")
tmp1 = as.data.frame(tmp$V4); colnames(tmp1)[1] = "V4"
tmp1 = tmp1 %>% separate(V4, c("V1", "V2", "V3"), sep="\t")
tmp = rbind(tmp[,c(1:3)], tmp1)
tmp1 = tmp %>% separate(V3, c(NA, "R2", "D"), sep="= ")
tmp1$R2 = gsub("\\s+", "", str_trim(tmp1$R2))
tmp1$R2 = str_replace_all(tmp1$R2, "D", "")
h1 = tmp1
names(h1) = c("snp1", "snp2", "R2", "Dprime")
write.table(h1, "filesforLD/all_leadSNP/salt-reinfo_randSNPpair_all.5.ld", quote = F, sep = "\t", row.names = F, col.names = T)
##Done

h1 = read.table("filesforLD/all_leadSNP/salt_compensating_all_snps.ld", header = T, sep = "\t")
h1.1 = read.table("filesforLD/all_leadSNP/salt-comp_randSNPpair_all.1.ld", header = T, sep = "\t")
h1.2 = read.table("filesforLD/all_leadSNP/salt-comp_randSNPpair_all.2.ld", header = T, sep = "\t")
h1.3 = read.table("filesforLD/all_leadSNP/salt-comp_randSNPpair_all.3.ld", header = T, sep = "\t")
h1.4 = read.table("filesforLD/all_leadSNP/salt-comp_randSNPpair_all.4.ld", header = T, sep = "\t")
h1.5 = read.table("filesforLD/all_leadSNP/salt-comp_randSNPpair_all.5.ld", header = T, sep = "\t")

h2 = read.table("filesforLD/all_leadSNP/salt_reinforcing_all_snps.ld", header = T, sep = "\t")
h2.1 = read.table("filesforLD/all_leadSNP/salt-reinfo_randSNPpair_all.1.ld", header = T, sep = "\t")
h2.2 = read.table("filesforLD/all_leadSNP/salt-reinfo_randSNPpair_all.2.ld", header = T, sep = "\t")
h2.3 = read.table("filesforLD/all_leadSNP/salt-reinfo_randSNPpair_all.3.ld", header = T, sep = "\t")
h2.4 = read.table("filesforLD/all_leadSNP/salt-reinfo_randSNPpair_all.4.ld", header = T, sep = "\t")
h2.5 = read.table("filesforLD/all_leadSNP/salt-reinfo_randSNPpair_all.5.ld", header = T, sep = "\t")

summary(h1$R2)
summary(h1.1$R2); summary(h1.2$R2); summary(h1.3$R2); summary(h1.4$R2); summary(h1.5$R2);
boxplot(h1$R2, h1.1$R2, h1.2$R2, h1.3$R2, h1.4$R2, h1.5$R2)
summary(h2$R2)
summary(h2.1$R2); summary(h2.2$R2); summary(h2.3$R2); summary(h2.4$R2); summary(h2.5$R2);
boxplot(h2$R2, h2.1$R2, h2.2$R2, h2.3$R2, h2.4$R2, h2.5$R2)
boxplot(h1$R2, h1.1$R2, h1.2$R2, h1.3$R2, h1.4$R2, h1.5$R2, h2$R2, h2.1$R2, h2.2$R2, h2.3$R2, h2.4$R2, h2.5$R2)
##OK, so looks like the compensating and teh reinforcing SNPs both sets have higher r2 as compared to the 5 random background SNP datasets

#Might do like 1000 iterations or something to get a null dist and a p-value
#awk '{print $2}' indica123.SNP.final.bial.pruned.bim | shuf -n 10000 > snplistforBackgroundLD.txt
#plink --bfile indica123.SNP.final.bial.pruned --ld-snp-list snplistforBackgroundLD.txt --r2 inter-chr --ld-window-r2 0 --out backgroundLD_selectedSNPS.inter
#Using the backgroundLD_selectedSNPS.inter.ld file I then extracted random compensating and reinforcing LD with same number of pairs with equal numbers of variant pairs and a similar distribution of distances between the cis- and trans-variants (code as above)

rand = read.table("filesforLD/all_leadSNP/comp_rand_meanR2.txt", header = F, sep="\t")
compranddist = as.vector(rand$V1)
rand = read.table("filesforLD/all_leadSNP/rein_rand_meanR2.txt", header = F, sep="\t")
reinranddist = as.vector(rand$V1)

#nulldist = compranddist
#OGteststat =  mean(h1$R2)
nulldist = reinranddist
OGteststat =  mean(h2$R2)
m=1000
#hist(nulldist); abline(v=OGteststat, col="red")
lowtail = sum(nulldist <= OGteststat) + 1
uptail = sum(nulldist >= OGteststat) + 1
numerator = min(lowtail, uptail)
denom = m+1
pval = 2*numerator/denom
#two-sided pval compensating = 0.001998002
#two-sided pval reinforcing = 0.001998002



#for wet
cis.wet.2017.lead = read.table("cis-eQTLs.wet.leadSNPs.txt", header = T, sep="\t")
trans.1mb.wet.2017.lead = read.table("trans-eQTLs.wet.1mb.fdr001.txt", header = T, sep="\t")
wet_gene = read.table("wet_cis_trans_compensating_reinfo.txt", header = T, sep="\t")

t1 = cis.wet.2017.lead[which(cis.wet.2017.lead$gene %in% wet_gene[which(wet_gene$Group == "Compensating"),]$TrID),]
t2 = trans.1mb.wet.2017.lead[which(trans.1mb.wet.2017.lead$gene %in% wet_gene[which(wet_gene$Group == "Compensating"),]$TrID),]
wet_comp = merge(t1[,c(1,5,6)], t2[,c(1,5,6)], by="gene")
names(wet_comp) = c("gene", "cis-snp", "cis-Chr", "trans-snp", "trans-Chr")
dim(wet_comp)  #12434 SNP pairs
dim(wet_comp[which(wet_comp$`cis-Chr` == wet_comp$`trans-Chr`),])  #3995
#I might remove the same Chr for the ease of randomization, and just estimate ILD
tmp1 = wet_comp[which(wet_comp$`cis-Chr` != wet_comp$`trans-Chr`),]  #8439
write.table(tmp1, "wet_compensating_diifChr_snps.txt", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(wet_comp, "wet_compensating_all_snps.txt", quote = F, row.names = F, col.names = T, sep = "\t")


t1 = cis.wet.2017.lead[which(cis.wet.2017.lead$gene %in% wet_gene[which(wet_gene$Group == "Reinforcing"),]$TrID),]
t2 = trans.1mb.wet.2017.lead[which(trans.1mb.wet.2017.lead$gene %in% wet_gene[which(wet_gene$Group == "Reinforcing"),]$TrID),]
wet_re = merge(t1[,c(1,5,6)], t2[,c(1,5,6)], by="gene")
names(wet_re) = c("gene", "cis-snp", "cis-Chr", "trans-snp", "trans-Chr")
dim(wet_re)  #25178 SNP pairs
dim(wet_re[which(wet_re$`cis-Chr` == wet_re$`trans-Chr`),])  #12700
#I might remove the same Chr for the ease of randomization, and just estimate ILD
tmp2 = wet_re[which(wet_re$`cis-Chr` != wet_re$`trans-Chr`),]  #12478
write.table(tmp2, "wet_reinforcing_diifChr_snps.txt", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(wet_re, "wet_reinforcing_all_snps.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#estimating the min distance between the same chr SNPs
t1 = cis.wet.2017.lead[which(cis.wet.2017.lead$gene %in% wet_gene[which(wet_gene$Group == "Compensating"),]$TrID),]
t2 = trans.1mb.wet.2017.lead[which(trans.1mb.wet.2017.lead$gene %in% wet_gene[which(wet_gene$Group == "Compensating"),]$TrID),]
tmp1 = merge(t1[,c(1,5,6,7)], t2[,c(1,5,6,7)], by="gene")
names(tmp1) = c("gene", "cis-snp", "cis-Chr","cis-Pos", "trans-snp", "trans-Chr", "trans-Pos")
tmp1 = tmp1[which(tmp1$`cis-Chr` == tmp1$`trans-Chr`),]
tmp1$snpDist = abs(tmp1$`trans-Pos` - tmp1$`cis-Pos`)
min(tmp1$snpDist)  #909998

t1 = cis.wet.2017.lead[which(cis.wet.2017.lead$gene %in% wet_gene[which(wet_gene$Group == "Reinforcing"),]$TrID),]
t2 = trans.1mb.wet.2017.lead[which(trans.1mb.wet.2017.lead$gene %in% wet_gene[which(wet_gene$Group == "Reinforcing"),]$TrID),]
tmp1 = merge(t1[,c(1,5,6,7)], t2[,c(1,5,6,7)], by="gene")
names(tmp1) = c("gene", "cis-snp", "cis-Chr","cis-Pos", "trans-snp", "trans-Chr", "trans-Pos")
tmp1 = tmp1[which(tmp1$`cis-Chr` == tmp1$`trans-Chr`),]
tmp1$snpDist = abs(tmp1$`trans-Pos` - tmp1$`cis-Pos`)
min(tmp1$snpDist)  #904450

library(dplyr)
g1 = sample_n(snpfile, 70000)
g2 = sample_n(snpfile, 70000)
f1 = cbind(g1,g2)
names(f1) = c("snp_id1", "Chr1", "Pos1", "snp_id2", "Chr2", "Pos2")
f2 = f1[which(f1$Chr1 != f1$Chr2),]
f3 = f1[which(f1$Chr1 == f1$Chr2),]
f3$Dist = abs(f3$Pos2 - f3$Pos1)
h1.1 = sample_n(f2, 8439)
h1.2 = sample_n(f3[which(f3$Dist > 909998),], 3995)
h1 = rbind(h1.1, h1.2[,-7])
write.table(h1, "wet-comp_randSNPpair_all.txt", row.names = F, col.names = F, quote = F, sep = "\t")

g1 = sample_n(snpfile, 170000)
g2 = sample_n(snpfile, 170000)
f1 = cbind(g1,g2)
names(f1) = c("snp_id1", "Chr1", "Pos1", "snp_id2", "Chr2", "Pos2")
f2 = f1[which(f1$Chr1 != f1$Chr2),]
f3 = f1[which(f1$Chr1 == f1$Chr2),]
f3$Dist = abs(f3$Pos2 - f3$Pos1)
h1.1 = sample_n(f2, 12478)
h1.2 = sample_n(f3[which(f3$Dist > 904450),], 12700)
h1 = rbind(h1.1, h1.2[,-7])
write.table(h1, "wet-reinfo_randSNPpair_all.txt", row.names = F, col.names = F, quote = F, sep = "\t")
rm(g1,g2,f1,f2,f3,h1.1,h1.2,h1)
