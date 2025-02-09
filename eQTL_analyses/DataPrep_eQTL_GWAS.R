
library(tidyr)
library(stringr)

dat = read.table("ricesalinity/data/GWAS/Indica/ind123.SNP.bial.pruned.txt", header = T, sep="\t")
#prepping the genetic map/snp map files for eQTL and GWAS (with matrix-eQTL and GAPIT)
snpid = data.frame(dat[,c(1)])
names(snpid)[1] = "snp_id"
tmp = separate(data = snpid, col = snp_id, into = c("Chr", "Pos"), sep = "_", remove=FALSE)
snpid = tmp
write.table(snpid, "ricesalinity/data/eQTL/ind_data/snps.location.txt", row.names = F, col.names = T, quote = F, sep="\t")
snpid$Chr = str_replace(snpid$Chr, "chr0", "")
snpid$Chr = str_replace(snpid$Chr, "chr", "")
write.table(snpid, "ricesalinity/data/GWAS/Indica/snps.location.txt", row.names = F, col.names = T, quote = F, sep="\t")

#transposing the data for gapit in
newname=colnames(dat)
tmp = as.data.frame(transpose(dat[,-1]))
colnames(tmp) = snpid$snp_id
tmp$IRGC = newname[-1]
tmp = tmp[,c(246715, 1:246714)]
write.table(tmp, "ricesalinity/data/GWAS/Indica/ind123.SNP.bial.pruned.GAPITin.txt", row.names = F, col.names = T, quote = F, sep="\t")

#Pheno Data
pheno=  read.table("ricesalinity/data/Physio_otherMeasurements/wet_rawdata.txt", header = T, sep="\t")
orderind = read.table("ricesalinity/data/eQTL/ind_data/ordered.in.txt", header = FALSE, sep = "\t")
info1 = read.table("ricesalinity/data/Fitness/Indica_wet_fitness.txt", header = T, sep="\t")
info2 = read.table("ricesalinity/data/Fitness/Japonica_wet_fitness.txt", header = T, sep="\t")

#I am loading both Jao and Ind for indica since IRGC4541 is ind but has been recorded as jap in our datasets. 
info1$IRGC.Nr = str_replace(info1$IRGC.Nr, " ", "")
info2$IRGC.Nr = str_replace(info2$IRGC.Nr, " ", "")
tmp1 = info1[which(info1$IRGC.Nr %in% orderind$V1),]
tmp2 = info2[which(info2$IRGC.Nr %in% orderind$V1),]
names(tmp1)[12] = "Fecundity"; names(tmp2)[12] = "Fecundity"
tmp = rbind(tmp1[,c(5,8,12)],tmp2[,c(5,8,12)])
phen = merge(tmp, pheno, by="Plot")
#We have all the data now

#Next since we have multiple replicates per genotype, I will take the aggregate of these (mean) and then log-transform them to bring them close to normality for GWAS. 
tmp = aggregate(phen[, -c(1:2,4:5)], by = list(phen$IRGC.Nr), FUN = mean, na.rm = TRUE)
names(tmp)[1] = "IRGC"
write.table(tmp, "ricesalinity/data/GWAS/Indica/ind123_pheno_wet.mean.txt", row.names = F, col.names = T, sep = "\t", quote = F)
#Use the above command and not the below one as they treat NAs differently
#tmp = aggregate(. ~ IRGC.Nr, data = phen[,-c(1,4:5)], FUN = mean)
#apply(tmp[,-1],2,shapiro.test)
test =  as.data.frame(log10(tmp[,-1]))
phen = data.frame(tmp[,1])
phen = cbind(phen, test)
names(phen)[1] = "IRGC"
write.table(phen, "ricesalinity/data/GWAS/Indica/ind123_pheno_wet.norm.txt", row.names = F, col.names = T, quote = F, sep="\t")


pheno=  read.table("ricesalinity/data/Physio_otherMeasurements/salt_rawdata.txt", header = T, sep="\t")
orderind = read.table("ricesalinity/data/eQTL/ind_data/ordered.in.txt", header = FALSE, sep = "\t")
info1 = read.table("ricesalinity/data/Fitness/Indica_salt_fitness.txt", header = T, sep="\t")
info2 = read.table("ricesalinity/data/Fitness/Japonica_salt_fitness.txt", header = T, sep="\t")

#I am loading both Jap and Ind for indica since IRGC4541 is ind but has been recorded as jap in our datasets. 
info1$IRGC.Nr = str_replace(info1$IRGC.Nr, " ", "")
info2$IRGC.Nr = str_replace(info2$IRGC.Nr, " ", "")
tmp1 = info1[which(info1$IRGC.Nr %in% orderind$V1),]
tmp2 = info2[which(info2$IRGC.Nr %in% orderind$V1),]
names(tmp1)[12] = "Fecundity"; names(tmp2)[12] = "Fecundity"
tmp = rbind(tmp1[,c(5,8,12)],tmp2[,c(5,8,12)])
phen = merge(tmp, pheno, by="Plot")
phen$chl_b_T2 = as.numeric(phen$chl_b_T2)
#We have all the data now

#Next since we have multiple replicates per genotype, I will take the aggregate of these (mean) and then log-transform them to bring them close to normality for GWAS. 
tmp = aggregate(phen[, -c(1:2,4:5)], by = list(phen$IRGC.Nr), FUN = mean, na.rm = TRUE)
names(tmp)[1] = "IRGC"
write.table(tmp, "ricesalinity/data/GWAS/Indica/ind123_pheno_salt.mean.txt", row.names = F, col.names = T, sep = "\t", quote = F)
#Use the above command and not the below one as they treat NAs differently
#tmp = aggregate(. ~ IRGC.Nr, data = phen[,-c(1,4:5)], FUN = mean)
#apply(tmp[,-1],2,shapiro.test)
test =  as.data.frame(log10(tmp[,-1]))
phen = data.frame(tmp[,1])
phen = cbind(phen, test)
names(phen)[1] = "IRGC"

write.table(phen, "ricesalinity/data/GWAS/Indica/ind123_pheno_salt.norm.txt", row.names = F, col.names = T, quote = F, sep="\t")
#for IRGC88396, I changed the fecundity to NA

cova = read.csv("GAPIT.Genotype.PCA.csv", header = T, sep=",")
newcol = cova$taxa
covname = colnames(cova)
tmp = data.frame(t(cova[,-1]))
colnames(tmp) = newcol
rownames(tmp) = covname[-1]
write.table(tmp, "ricesalinity/data/eQTL/ind_data/pop_str_cov.ind.txt", row.names = T, col.names = T, sep = "\t", quote = F)

