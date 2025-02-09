setwd("ricesalinity/data/eQTL/ind_data/ems_paper_analysis/")
library(dplyr)
library(ggplot2)

maf = read.table("SNPid_maf.txt", header = F, sep="\t")
names(maf)[1] = "SNP_id"; names(maf)[2] = "maf"
cbind()

cis.salt.eqtl = read.table("../fdr_trans1Mb/cis-eQTLs.salt.leadSNPs.txt", header = T, sep="\t")
cis.wet.eqtl = read.table("../fdr_trans1Mb/cis-eQTLs.wet.leadSNPs.txt", header = T, sep="\t")
trans.salt.eqtl = read.table("../fdr_trans1Mb/trans-eQTLs.salt.1mb.leadSNPs.txt", header = T, sep="\t")
trans.wet.eqtl = read.table("../fdr_trans1Mb/trans-eQTLs.wet.1mb.leadSNPs.txt", header = T, sep="\t")
names(cis.salt.eqtl)[5] = "SNP_id"; names(cis.wet.eqtl)[5] = "SNP_id"; 
names(trans.salt.eqtl)[5] = "SNP_id"; names(trans.wet.eqtl)[5] = "SNP_id"; 

cis.wet = maf[which(maf$SNP_id %in% cis.wet.eqtl$SNP_id),]
cis.salt = maf[which(maf$SNP_id %in% cis.salt.eqtl$SNP_id),]
trans.wet = maf[which(maf$SNP_id %in% trans.wet.eqtl$SNP_id),]
trans.salt = maf[which(maf$SNP_id %in% trans.salt.eqtl$SNP_id),]
mafleft_wet = data.frame(setdiff(maf, rbind(cis.wet, trans.wet)))
mafleft_salt = data.frame(setdiff(maf, rbind(cis.salt, trans.salt)))
summary(trans.wet$maf); summary(cis.wet$maf); summary(mafleft_wet$maf)
boxplot(mafleft_wet$maf, cis.wet$maf, trans.wet$maf)

#Tesing for diff in maf (lower maf -- stronger negative selection)
t.test(cis.wet$maf, mafleft_wet$maf, alternative = "greater")
t.test(cis.salt$maf, mafleft_salt$maf, alternative = "greater")

t.test(cis.wet$maf, mafleft_wet$maf)
t.test(cis.salt$maf, mafleft_salt$maf)

t.test(mafleft_wet$maf, trans.wet[,2], alternative = "greater")     
t.test(mafleft_salt$maf, trans.salt[,2], alternative = "greater")   

mean(cis.wet[,2]); mean(trans.wet[,2]); mean(mafleft_wet$maf)    
t.test(cis.wet[,2], trans.wet[,2], alternative = "greater")       
mean(cis.salt[,2]); mean(trans.salt[,2]); mean(mafleft_salt$maf)    
t.test(cis.salt[,2], trans.salt[,2], alternative = "greater")        

cis.salt$group = rep("cis.salt", times=nrow(cis.salt))
trans.salt$group = rep("trans.salt", times=nrow(trans.salt))
cis.wet$group = rep("cis.wet", times=nrow(cis.wet))
trans.wet$group = rep("trans.wet", times=nrow(trans.wet))
mafleft_wet$group = rep("none.wet", times=nrow(mafleft_wet))


salt = rbind(cis.salt, trans.salt)
wet = rbind(cis.wet, trans.wet)
cis = rbind(cis.wet, cis.salt)
trans = rbind(trans.wet, trans.salt)

#Testing for the diff in effect sizes
mean(abs(cis.salt.eqtl$beta));mean(abs(trans.salt.eqtl$beta));    
mean(abs(cis.wet.eqtl$beta));mean(abs(trans.wet.eqtl$beta));      
wilcox.test(abs(cis.salt.eqtl$beta), abs(trans.salt.eqtl$beta))   
wilcox.test(abs(cis.wet.eqtl$beta), abs(trans.wet.eqtl$beta))     
#the effect sizes in trans are higher

#So just testing for just the upper quartile of beta values (to see whether the pattern hold even the upper quartile of beta values)
summary(abs(cis.salt.eqtl$beta))[5]
summary(abs(trans.salt.eqtl$beta))[5]
cis.salt.eqtl.3qu = cis.salt.eqtl[which(abs(cis.salt.eqtl$beta) >= 1.006836 ),]
trans.salt.eqtl.3qu = trans.salt.eqtl[which(abs(trans.salt.eqtl$beta) >= 1.281524  ),]
summary(abs(cis.salt.eqtl$beta)); summary(abs(trans.salt.eqtl.3qu$beta)); 
tmp1 = cis.salt[which(cis.salt$SNP_id %in% unique(cis.salt.eqtl.3qu$SNP_id)),]
tmp2 = trans.salt[which(trans.salt$SNP_id %in% unique(trans.salt.eqtl.3qu$SNP_id)),]
mean(tmp1[,2]); mean(tmp2[,2])    #0.1860632, 0.08169785
t.test(tmp1[,2], tmp2[,2], alternative = "greater")   #t = 17.87, df = 504.29, p-value < 2.2e-16

summary(abs(cis.wet.eqtl$beta))[5]
summary(abs(trans.wet.eqtl$beta))[5]
cis.wet.eqtl.3qu = cis.wet.eqtl[which(abs(cis.wet.eqtl$beta) >= 0.9672904 ),]
trans.wet.eqtl.3qu = trans.wet.eqtl[which(abs(trans.wet.eqtl$beta) >= 1.267311  ),]
summary(abs(cis.wet.eqtl$beta)); summary(abs(trans.wet.eqtl.3qu$beta)); 
tmp1 = cis.wet[which(cis.wet$SNP_id %in% unique(cis.wet.eqtl.3qu$SNP_id)),]
tmp2 = trans.wet[which(trans.wet$SNP_id %in% unique(trans.wet.eqtl.3qu$SNP_id)),]
mean(tmp1[,2]); mean(tmp2[,2])    #0.19007357, 0.08316512 
t.test(tmp1[,2], tmp2[,2], alternative = "greater")   #t = 19.226, df = 581.47, p-value < 2.2e-16
#So trans is under stronger seelction even when the uppper quartile of values are taken


#####The above shows that trans is under stonger selection than cis in both environments
cis.salt.eqtl = merge(cis.salt.eqtl, cis.salt, by="SNP_id")
cis.wet.eqtl = merge(cis.wet.eqtl, cis.wet, by="SNP_id")
trans.salt.eqtl = merge(trans.salt.eqtl, trans.salt, by="SNP_id")
trans.wet.eqtl = merge(trans.wet.eqtl, trans.wet, by="SNP_id")

###Now I will test for difference in selection within the same type of eQTL but between env
#same tests as above but for same type of eQTL in different env
t.test(cis.salt[,2], cis.wet[,2])   #no diff between cis maf between env
t.test(trans.salt[,2], trans.wet[,2], alternative = "less")   #t = -2.3749, df = 13587, p-value = 0.008783
t.test(abs(trans.salt.eqtl$beta), abs(trans.wet.eqtl$beta), alternative = "greater") #t = 2.1745, df = 21233, p-value = 0.01484
wilcox.test(abs(trans.salt.eqtl$beta), abs(trans.wet.eqtl$beta), alternative = "greater") #W = 58710892, p-value = 0.01315
#trans-eQTLs in the salt conditions are under stronger selection (lower maf) and so is the the effect size (bigger in salt conditions)

##the above tells us that no sig diff between cis, whereas both the effect size and -ve selection is stronger in salt than normal for trans

#So just testing for just the upper quartile of beta values (to see whether the pattern hold even the upper quartile of beta values)
summary(abs(trans.wet.eqtl$beta))[5]
summary(abs(trans.salt.eqtl$beta))[5]
trans.wet.eqtl.3qu = trans.wet.eqtl[which(abs(trans.wet.eqtl$beta) >= 1.267311 ),]
trans.salt.eqtl.3qu = trans.salt.eqtl[which(abs(trans.salt.eqtl$beta) >= 1.281524 ),]
summary(abs(trans.wet.eqtl.3qu$beta)); summary(abs(trans.salt.eqtl.3qu$beta)); 
tmp1 = trans.wet[which(trans.wet$SNP_id %in% unique(trans.wet.eqtl.3qu$SNP_id)),]
tmp2 = trans.salt[which(trans.salt$SNP_id %in% unique(trans.salt.eqtl.3qu$SNP_id)),]
mean(tmp1[,2]); mean(tmp2[,2])    #0.0839506, 0.08170797
t.test(tmp1[,2], tmp2[,2], alternative = "greater")   #t = 0.84515, df = 3684.8, p-value = 0.199
#It doesnt really hold up for higher values


##### Next, I wille explore the trans-eQTLs that impact one gene vs more than one gene
dim(trans.wet.eqtl)   #10054
dim(trans.salt.eqtl)  #11478
tmp1 = as.data.frame(table(trans.wet.eqtl$SNP_id))
tmp2 = as.data.frame(table(trans.salt.eqtl$SNP_id))

dim(tmp1)[1] - dim(tmp1[which(tmp1$Freq > 1),])[1]    #5306
dim(tmp1[which(tmp1$Freq > 1),])[1]     #1569
dim(tmp2)[1] - dim(tmp2[which(tmp2$Freq > 1),])[1]    #4986
dim(tmp2[which(tmp2$Freq > 1),])[1]      #1728
summary(tmp1$Freq); summary(tmp2$Freq)
summary(tmp1[which(tmp1$Freq >1),]$Freq); summary(tmp2[which(tmp2$Freq >1),]$Freq)

uniq.trans.wet.eqtl = trans.wet.eqtl[which(trans.wet.eqtl$SNP_id %in% (as.vector(tmp1[which(tmp1$Freq == 1),]$Var1))),]
dup.trans.wet.eqtl = trans.wet.eqtl[which(trans.wet.eqtl$SNP_id %in% (as.vector(tmp1[which(tmp1$Freq > 1),]$Var1))),]
uniq.trans.salt.eqtl = trans.salt.eqtl[which(trans.salt.eqtl$SNP_id %in% (as.vector(tmp2[which(tmp2$Freq == 1),]$Var1))),]
dup.trans.salt.eqtl = trans.salt.eqtl[which(trans.salt.eqtl$SNP_id %in% (as.vector(tmp2[which(tmp2$Freq > 1),]$Var1))),]
#tmp1 = trans.wet.eqtl[which(trans.wet.eqtl$SNP_id %in% (as.vector(tmp1[which(tmp1$Freq > 12),]$Var1))),]

#now-testing maf
uniq.trans.wet = trans.wet[which(trans.wet$SNP_id %in% unique(uniq.trans.wet.eqtl$SNP_id)),]
dup.trans.wet = trans.wet[which(trans.wet$SNP_id %in% unique(dup.trans.wet.eqtl$SNP_id)),]
uniq.trans.salt = trans.salt[which(trans.salt$SNP_id %in% unique(uniq.trans.salt.eqtl$SNP_id)),]
dup.trans.salt = trans.salt[which(trans.salt$SNP_id %in% unique(dup.trans.salt.eqtl$SNP_id)),]

t.test(uniq.trans.wet$maf, uniq.trans.salt$maf, alternative = "greater")  #t = 1.3959, df = 10256, p-value = 0.08139
t.test(dup.trans.wet$maf, dup.trans.salt$maf, alternative = "greater")  #t = 1.2486, df = 3201.8, p-value = 0.106
#no sig difference between maf for either single or dup trans eqtls between env

mean(uniq.trans.wet$maf); mean(dup.trans.wet$maf)
mean(uniq.trans.salt$maf); mean(dup.trans.salt$maf)


#effect size
mean(abs(uniq.trans.wet.eqtl$beta)); mean(abs(uniq.trans.salt.eqtl$beta))   #0.8738061, 0.8866009
mean(abs(dup.trans.wet.eqtl$beta)); mean(abs(dup.trans.salt.eqtl$beta))     #1.039581, 1.030649

t.test(abs(uniq.trans.wet.eqtl$beta), abs(uniq.trans.salt.eqtl$beta), alternative = "less")       #t = -1.8151, df = 15985, p-value = 0.03476
t.test(abs(dup.trans.wet.eqtl$beta), abs(dup.trans.salt.eqtl$beta), alternative = "greater")   #t = 1.1275, df = 14987, p-value = 0.1298
#no significant difference in effect size due to singletons

t.test(uniq.trans.wet$maf, dup.trans.wet$maf, alternative = "greater")  #t = 24.356, df = 6787.3, p-value < 2.2e-16
t.test(uniq.trans.salt$maf, dup.trans.salt$maf, alternative = "greater")  #t = 25.201, df = 7984.9, p-value < 2.2e-16

###############################################################################################################

###HOtspot detection
#below is binning the genome which is common ffor all conditions
library(Biostrings)
mygenome <- readDNAStringSet(list.files("ricesalinity/data/fromOnlineResources/other_probablyNOTused/","IRGSP-1.0_genomic.fna", full=TRUE))
chrSizes <- width(mygenome)
chrSizes = chrSizes[1:12]
names(chrSizes) = as.character(c(1:12))
#names(chrSizes) <- names(mygenome)
#print(chrSizes)
bins   <- tileGenome(chrSizes, tilewidth=100000, cut.last.tile.in.chrom=T)
test_bin = as.data.frame(bins)
test_bin$binno = c(1:nrow(test_bin))

library(tidyr)
#To change the condition replace the dataframe represented by a below:
a = trans.salt.eqtl
tmp1 = as.data.frame(table(a$SNP_id))
tmp = separate(data = tmp1, col = Var1, into = c("Chr", "Pos"), sep = "_", remove=FALSE)
tmp$Pos = as.numeric(tmp$Pos)

df = tmp
gr = GRanges(seqnames=df$Chr,
             IRanges(start=df$Pos,end=df$Pos),freq=df$Freq)
hits.df <- as.data.frame(findOverlaps(bins, gr))

binno = c()
nogenes = c()
for (i in (unique(hits.df$queryHits))){
  x = tmp[(hits.df[which(hits.df$queryHits == i),]$subjectHits),]$Var1
  y = length(unique(a[which(a$SNP_id %in% x),]$gene))
  binno=append(binno, i)
  nogenes = append(nogenes, y)
}

z = data.frame(binno, nogenes)
x = merge(test_bin,z, by="binno")
geneperbin_wet = x    #Change name of df here to save
geneperbin_salt = x    #Change name of df here to save

y1 = merge(test_bin,z, by="binno", all.x=T)
y2 = merge(test_bin,z, by="binno", all.x=T)
test = merge(y1,y2,by="binno")
test = test[!with(test,is.na(nogenes.x)& is.na(nogenes.y)),]
test = test[,c(1:5,7,13)]
names(test) = c("binno", "chr", "Start", "Stop", "Width", "Nogenes_Wet", "Nogenes_Salt")
geneperbin = test

write.table(geneperbin, "../fdr_trans1Mb/transGenePer100kb_fromlead1mbfile.txt", quote = F, row.names = F, col.names = T, sep="\t")

tmp1 = as.data.frame(table(trans.wet.eqtl[,5]))
tmp2 = as.data.frame(table(trans.salt.eqtl[,5]))

geneexp_wet = read.table("ricesalinity/data/GeneExp/INDICA_Wet_GeneExp.transposed.txt", header = T, sep="\t", check.names = F)
geneexp_salt = read.table("ricesalinity/data/GeneExp/INDICA_Salt_GeneExp.transposed.txt", header = T, sep="\t", check.names = F)

#R2R3 is expressed in both 
mean(geneexp_wet$`OS07T0627300-01`)     #2.44137
mean(geneexp_salt$`OS07T0627300-01`)    #2.24618

tmp1 = geneexp_wet[,(colnames(geneexp_wet) %in% trans.salt.eqtl[which(trans.salt.eqtl$snps == "7_25986894"),]$gene)]
tmp1$r2r3 = geneexp_wet$`OS07T0627300-01`
tmp2 = geneexp_salt[,(colnames(geneexp_salt) %in% trans.salt.eqtl[which(trans.salt.eqtl$snps == "7_25986894"),]$gene)]
tmp2$r2r3 = geneexp_salt$`OS07T0627300-01`

data_cor1 =  as.data.frame(cor(tmp1[ ,colnames(tmp1) != "r2r3"], tmp1$r2r3))  # Calculate correlations
data_cor2 =  as.data.frame(cor(tmp2[ ,colnames(tmp2) != "r2r3"], tmp2$r2r3))  # Calculate correlations
data_cor1$gene = rownames(data_cor1); data_cor2$gene = rownames(data_cor2)

datacor = merge(data_cor1, data_cor2, by="gene")
names(datacor) = c("gene", "corr.wet", "corr.salt")


####Getting the genes that are regulated by hotspots in salt (#genes >= 30)
df = geneperbin[which(geneperbin$Nogenes_Salt >=30),]
bins = GRanges(seqnames=df$chr,
             IRanges(start=df$Start,end=df$Stop), nogenes_salt=df$Nogenes_Salt)
test_bin = as.data.frame(bins)
test_bin$binno = c(1:nrow(test_bin))

a = trans.salt.eqtl
tmp1 = as.data.frame(table(a[,1]))
tmp = separate(data = tmp1, col = Var1, into = c("Chr", "Pos"), sep = "_", remove=FALSE)
tmp$Pos = as.numeric(tmp$Pos)
df = tmp
gr = GRanges(seqnames=df$Chr,
             IRanges(start=df$Pos,end=df$Pos),freq=df$Freq)
hits.df <- as.data.frame(findOverlaps(bins, gr))

binno = c()
geneid = list()
for (i in (unique(hits.df$queryHits))){
  x = tmp[(hits.df[which(hits.df$queryHits == i),]$subjectHits),]$Var1
  y = (unique(a[which(a$SNP_id %in% x),]$gene))
  binno=append(binno, i)
  geneid = append(geneid, list(y))
}

z = data.frame(binno)
z$geneid = geneid

test = merge(test_bin, z, by="binno")
test = test[,-6]
names(test) = c("binno", "chr", "Start", "Stop", "Width", "Nogenes_Salt", "regulatedGeneid")
hotspot_genes = test

library(tidyverse)
hotspot_genes %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ' ')) %>% 
  write.csv('../fdr_trans1Mb/saltTrans_hotspots.gene.txt', row.names = FALSE, sep="\t")


##Getting gene id of eQTL snps
#reading the annotation file first
anno = read.csv("ricesalinity/data/fromOnlineResources/other_probablyNOTused/IRGSP-1.0_representative/locus.gff", header = F, sep = "\t")
tmp = as.data.frame(separate(data = anno, col = V9, into = c("geneid", "genename", "Anno"), extra = "drop", sep = ";", remove=T))
tmp$geneid = sub("ID=","",tmp$geneid)
tmp$Anno = sub("Note=","",tmp$Anno)
tmp = tmp[,c(9,11)]
anno = tmp

library(GenomicRanges)
ref = read.table("ricesalinity/data/fromOnlineResources/other_probablyNOTused/gene_cord.bed", header = T, sep="\t")

df = (trans.salt.eqtl[,c(1,6,7)])
Freq = as.data.frame(table(df$SNP_id))
df = unique(df)
tmp = merge(df, Freq, by.x="SNP_id", by.y="Var1")
df = tmp
snps = GRanges(seqnames=df$Chr.y,
               IRanges(start=df$Pos,end=df$Pos), snpid=df$SNP_id)
reference <- GRanges(seqnames=ref$Chr,
                     IRanges(start=ref$Start,end=ref$Stop), geneid=ref$Gene)

olap <- findOverlaps(snps, reference)
tmp = rep("NA", times=nrow(df))
tmp[queryHits(olap)] = reference$geneid[subjectHits(olap)]
df$geneid = tmp

tmp = merge(df, anno, by="geneid", all.x=T)
tmp = tmp[,c(2:5,1,6)]
tmp = tmp[order(tmp$Chr.y, tmp$Pos),]

cis.salt.eqtlgenes = tmp
cis.wet.eqtlgenes = tmp
trans.salt.eqtlgenes = tmp
trans.wet.eqtlgenes = tmp

cis.salt.eqtlgenes$group = rep("cis.salt", times = nrow(cis.salt.eqtlgenes))
cis.wet.eqtlgenes$group = rep("cis.wet", times = nrow(cis.wet.eqtlgenes))
trans.salt.eqtlgenes$group = rep("trans.salt", times = nrow(trans.salt.eqtlgenes))
trans.wet.eqtlgenes$group = rep("trans.wet", times = nrow(trans.wet.eqtlgenes))

eqtlgenes = rbind(cis.salt.eqtlgenes, cis.wet.eqtlgenes, trans.salt.eqtlgenes, trans.wet.eqtlgenes)
eqtlgenes = eqtlgenes[,c(7,1:6)]
write.table(eqtlgenes, "../fdr_trans1Mb/eqtl_snpGenes.txt", row.names = F, col.names = T, sep="\t", quote = F)
