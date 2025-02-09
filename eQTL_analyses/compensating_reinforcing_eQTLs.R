##Compensating and reinforcing eQTLs

setwd("ricesalinity/data/eQTL/ind_data/fdr_trans1Mb/")

cis.salt.2017 = read.table("cis-eQTLs.salt.fdr001.txt", header = T, sep="\t")
cis.wet.2017 = read.table("cis-eQTLs.wet.fdr001.txt", header = T, sep="\t")
trans.1mb.salt.2017 = read.table("trans-eQTLs.salt.1mb.fdr001.txt", header = T, sep="\t")
trans.1mb.wet.2017 = read.table("trans-eQTLs.wet.1mb.fdr001.txt", header = T, sep="\t")

##common eqtls between env
#ciscomm
dim(merge(cis.wet.2017, cis.salt.2017, by=c("gene","snps")))
dim(merge(trans.1mb.wet.2017, trans.1mb.salt.2017, by=c("gene","snps")))
dim(merge(cis.wet.2017.lead, cis.salt.2017.lead, by=c("gene","snps")))
dim(merge(trans.1mb.wet.2017.lead, trans.1mb.salt.2017.lead, by=c("gene","snps")))

summary(abs(cis.wet.2017$beta));
summary(abs(trans.1mb.wet.2017$beta))
t.test(abs(cis.wet.2017$beta), abs(trans.1mb.wet.2017$beta))

summary(abs(cis.salt.2017$beta));summary(abs(trans.1mb.salt.2017$beta))

ciscomm = merge(cis.wet.2017, cis.salt.2017, by=c("gene","snps"))
transcomm = merge(trans.1mb.wet.2017, trans.1mb.salt.2017, by=c("gene","snps"))
t.test((ciscomm$beta.x), (ciscomm$beta.y), paired = T)
t.test((transcomm$beta.x), (transcomm$beta.y), paired = T)

tmp = ciscomm[,c(2:4,1,5:8)]; colnames(tmp)= colnames(cis.wet.2017)
cisdiff_wet = setdiff(cis.wet.2017, tmp)
tmp = ciscomm[,c(2,9:10,1,11:14)]; colnames(tmp)= colnames(cis.salt.2017)
cisdiff_salt = setdiff(cis.salt.2017, tmp)

tmp = transcomm[,c(1,3:5,2,6:12)]; colnames(tmp)= colnames(trans.1mb.wet.2017)
transdiff_wet = setdiff(trans.1mb.wet.2017, tmp)
tmp = transcomm[,c(1,13:15,2,16:22)]; colnames(tmp)= colnames(trans.1mb.salt.2017)
transdiff_salt = setdiff(trans.1mb.salt.2017, tmp)

#Looking at the beta values for common and different eQTLs
median(abs(ciscomm$beta.x)); median(abs(cisdiff_wet$beta))
median(abs(ciscomm$beta.y)); median(abs(cisdiff_salt$beta))
median(abs(transcomm$beta.x)); median(abs(transdiff_wet$beta))
median(abs(transcomm$beta.y)); median(abs(transdiff_salt$beta))

#testing whether the proportions are significantly different (for fdr 0.001)
x = c(29623,25528)
n = c(59669, 137100)
prop.test(x,n)  #this is the test reported (this does not take leadSNPs)
x = c(30046, 111572)
prop.test(x,n)

cis.wet.2017.lead = read.table("cis-eQTLs.wet.leadSNPs.txt", header = T, sep="\t")
trans.1mb.wet.2017.lead = read.table("trans-eQTLs.wet.1mb.leadSNPs.txt", header = T, sep="\t")

#Testing for compensating vs stabilizing cis-trans
x = cis.wet.2017.lead
y = trans.1mb.wet.2017.lead
x$grp = rep("Cis.Wet", times=nrow(x))
y$grp = rep("trans.Wet", times=nrow(y))
tmp1 = x[which(x$gene %in% as.vector(intersect(x$gene, y$gene))),]
tmp1 = aggregate(beta ~ gene + grp, data = tmp1, FUN = mean, na.rm = TRUE)
tmp2 = y[which(y$gene %in% as.vector(intersect(x$gene, y$gene))),]
tmp2 = aggregate(beta ~ gene + grp, data = tmp2, FUN = mean, na.rm = TRUE)
wetcomm = merge(tmp1, tmp2, by="gene")
write.table(wetcomm, "wet_common_cistrans.txt", row.names = F, col.names = T, quote = T, sep="\t")

cis.salt.2017.lead = read.table("cis-eQTLs.salt.leadSNPs.txt", header = T, sep="\t")
trans.1mb.salt.2017.lead = read.table("trans-eQTLs.salt.1mb.leadSNPs.txt", header = T, sep="\t")

x = cis.salt.2017.lead
y = trans.1mb.salt.2017.lead
x$grp = rep("Cis.Salt", times=nrow(x))
y$grp = rep("trans.Salt", times=nrow(y))
tmp1 = x[which(x$gene %in% as.vector(intersect(x$gene, y$gene))),]
tmp1 = aggregate(beta ~ gene + grp, data = tmp1, FUN = mean, na.rm = TRUE)
tmp2 = y[which(y$gene %in% as.vector(intersect(x$gene, y$gene))),]
tmp2 = aggregate(beta ~ gene + grp, data = tmp2, FUN = mean, na.rm = TRUE)
saltcomm = merge(tmp1, tmp2, by="gene")
write.table(saltcomm, "salt_common_cistrans.txt", row.names = F, col.names = T, quote = T, sep="\t")

##Wetcomm is 434 genes and saltcomm has 388 genes with correlation about ~0.845.
#But another thing to consider is what combination of alleles are present in cis and trans 
#for this I will first for each gene, check whether each cis-trans pair is same direction, if yes I will estimate the proportion (00+22)/(00+22+02+20), if not then proportion (02+20)/(00+22+02+20).
#save the proportion of each cis-trans pair per gene and then take the mean per gene

snpfile = read.csv("ricesalinity/data/eQTL/ind_data/ind123.SNP.bial.pruned.txt", header = T, sep="\t")
cis.wet.2017.lead = read.table("cis-eQTLs.wet.leadSNPs.txt", header = T, sep="\t")
trans.1mb.wet.2017.lead = read.table("trans-eQTLs.wet.1mb.fdr001.txt", header = T, sep="\t")

cdfmain = cis.wet.2017.lead
tdfmain = trans.1mb.wet.2017.lead
a = wetdiffsign$gene

geneid = c()
meanpropsame = c()
meanpropdiff = c()

for (x in a){
  geneid = append(geneid, x)
  cdf = cdfmain[which(cdfmain$gene == x),]
  tdf = tdfmain[which(tdfmain$gene == x),]
  
  propsame = c()
  propdiff = c()
  for (i in 1:nrow(cdf)){
    tmp1 = snpfile[which(snpfile$snp_id == cdf$snps[i]),]
    for (j in 1:nrow(tdf)){
      tmp2 = snpfile[which(snpfile$snp_id == tdf$snps[j]),]
      tmp = rbind(tmp1, tmp2)
      c = colSums(tmp[,-1])
      p = ifelse(sign(cdf$beta[i]) == sign(tdf$beta[j]), (length(c[which(c==0 | c == 4)])/length(c)), (length(c[which(c==2)])/length(c)))
      propsame = append(propsame, p)
      p = ifelse(sign(cdf$beta[i]) == sign(tdf$beta[j]), (length(c[which(c==2)])/length(c)), (length(c[which(c==0 | c == 4)])/length(c)))
      propdiff = append(propdiff, p)
    }
  }
  meanpropsame = append(meanpropsame, mean(propsame))
  meanpropdiff = append(meanpropdiff, mean(propdiff))
}

df = data.frame(geneid, meanpropsame, meanpropdiff)

prop_wetsame = df
prop_wetdiff = df
names(prop_wetdiff) = c("geneid", "meanpropdiff", "meanpropsame")

prop_saltsame = df
prop_saltdiff = df
names(prop_saltdiff) = c("geneid", "meanpropdiff", "meanpropsame")

write.table(wetcomm, "wet_common_cistrans.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(saltcomm, "salt_common_cistrans.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(prop_wetsame, "wetReinfor_prop-revscomp.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(prop_saltsame, "saltReinfor_prop-revscomp.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(prop_wetdiff, "wetComp_prop-revscomp.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(prop_saltdiff, "saltComp_prop-revscomp.txt", row.names = F, col.names = T, sep="\t", quote = F)

wetsamesign = wetcomm[ (sign(wetcomm$beta.x) == sign(wetcomm$beta.y)), ]
wetdiffsign = wetcomm[ !(sign(wetcomm$beta.x) == sign(wetcomm$beta.y)), ]
saltsamesign = saltcomm[ (sign(saltcomm$beta.x) == sign(saltcomm$beta.y)), ]
saltdiffsign = saltcomm[ !(sign(saltcomm$beta.x) == sign(saltcomm$beta.y)), ]

wetsamesign = wetsamesign[which(wetsamesign$gene %in% prop_wetsame[which(prop_wetsame$meanpropsame >= 0.6),]$geneid),]
wetdiffsign = wetdiffsign[which(wetdiffsign$gene %in% prop_wetdiff[which(prop_wetdiff$meanpropdiff >= 0.6),]$geneid),]
saltsamesign = saltsamesign[which(saltsamesign$gene %in% prop_saltsame[which(prop_saltsame$meanpropsame >= 0.6),]$geneid),]
saltdiffsign = saltdiffsign[which(saltdiffsign$gene %in% prop_saltdiff[which(prop_saltdiff$meanpropdiff >= 0.6),]$geneid),]

wetsamesign$grp = rep("Same", times = nrow(wetsamesign))
wetdiffsign$grp = rep("Diff", times = nrow(wetdiffsign))
wetcomm = rbind(wetsamesign, wetdiffsign)
saltsamesign$grp = rep("Same", times = nrow(saltsamesign))
saltdiffsign$grp = rep("Diff", times = nrow(saltdiffsign))
saltcomm = rbind(saltsamesign, saltdiffsign)

write.table(wetsamesign, "wet_common_samesign_cistrans.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(wetdiffsign, "wet_common_diffsign_cistrans.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(saltsamesign, "salt_common_samesign_cistrans.txt", row.names = F, col.names = T, sep="\t", quote = F)
write.table(wetdiffsign, "salt_common_diffsign_cistrans.txt", row.names = F, col.names = T, sep="\t", quote = F)


##################################################################
####Getting Gene expression polymorphism
ind_wet = read.table("ricesalinity/data/geneExp/INDICA_Wet_GeneExp.transposed.txt", header = T, sep="\t")
ind_salt = read.table("ricesalinity/data/geneExp/INDICA_Salt_GeneExp.transposed.txt", header = T, sep="\t")

identical(colnames(ind_wet), colnames(ind_salt))  ##Just to test that the column names are the same (and in the same order)
ind_all = rbind(ind_wet, ind_salt)            ##This has all the transcript counts for all the indica pop (along with a column of RNASeqCode)

geno_info = read.table("ricesalinity/data/QuanGen/INDICA-code_IRGC_Env.txt", header=TRUE, sep="\t")
names(geno_info)[1] = "RNASeqCode"
ind_quangenData = merge(geno_info, ind_all, by="RNASeqCode")
ind_quangenData$IRGC.Nr = as.factor(ind_quangenData$IRGC.Nr)
ind_quangenData$Environment = as.factor(ind_quangenData$Environment)

ind_geno_wet = ind_quangenData[which(ind_quangenData$Environment == 'Wet'),]
ind_geno_salt = ind_quangenData[which(ind_quangenData$Environment == 'Salt'),]

polymor = ind_geno_salt %>% group_by(IRGC.Nr) %>% 
  summarise_at(vars(-c(1:2)), list(name=mean)) %>% 
  summarise_if(is.numeric, var) %>% 
  slice(1) %>% as.numeric()
exppoly_salt = data.frame(colnames(colnames(ind_salt[,-1]), polymor))
names(exppoly_salt) = c("geneID", "ExpPolymorphismSalt")                        
write.table(exppoly_salt, "exp_polymorphism_salt.txt", row.names = F, col.names = T, quote = F, sep = "\t")

polymor = ind_geno_wet %>% group_by(IRGC.Nr) %>% 
  summarise_at(vars(-c(1:2)), list(name=mean)) %>% 
  summarise_if(is.numeric, var) %>% 
  slice(1) %>% as.numeric()
exppoly_wet = data.frame(colnames(colnames(ind_wet[,-1]), polymor))
names(exppoly_wet) = c("geneID", "ExpPolymorphismwet")  
write.table(exppoly_wet, "exp_polymorphism_wet.txt", row.names = F, col.names = T, quote = F, sep = "\t")


###Looking at compensating and reinforcing and their expression polymorphisms
z = as.data.frame(exppoly_wet[which(exppoly_wet$geneID %in% wetsamesign$gene),]$ExpPolymorphismWet)
names(z)[1] = "ExpPol"
z$grp = rep("Reinforcing", times=nrow(z))
y = as.data.frame(exppoly_wet[which(exppoly_wet$geneID %in% wetdiffsign$gene),]$ExpPolymorphismWet)
names(y)[1] = "ExpPol"
y$grp = rep("Compensating", times=nrow(y))
x = rbind(z,y)
wilcox.test(z$ExpPol, y$ExpPol, alternative = "greater")
boxplot(wetcov$ExpPolymorphismWet, z, y)

z = as.data.frame(exppoly_salt[which(exppoly_salt$geneID %in% wetsamesign$gene),]$ExpPolymorphismSalt)
names(z)[1] = "ExpPol"
z$grp = rep("Reinforcing", times=nrow(z))
y = as.data.frame(exppoly_salt[which(exppoly_salt$geneID %in% wetdiffsign$gene),]$ExpPolymorphismSalt)
names(y)[1] = "ExpPol"
y$grp = rep("Compensating", times=nrow(y))
x = rbind(z,y)
wilcox.test(z$ExpPol, y$ExpPol, alternative = "greater")
boxplot(exppoly_wet$ExpPolymorphismSalt, z, y)

