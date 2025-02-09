setwd("ricesalinity/data/eQTL/ind_data/fdr_trans1Mb/")

cis.salt.2017 = read.table("cis-eQTLs.salt.txt", header = T, sep="\t")
trans.salt.2017 = read.table("trans-eQTLs.salt.txt", header = T, sep="\t")
cis.wet.2017 = read.table("cis-eQTLs.wet.txt", header = T, sep="\t")
trans.wet.2017 = read.table("trans-eQTLs.wet.txt", header = T, sep="\t")

cis.wet.2017 = cis.wet.2017[which(cis.wet.2017$FDR < 0.001),]
cis.salt.2017 = cis.salt.2017[which(cis.salt.2017$FDR < 0.001),]
trans.wet.2017 = trans.wet.2017[which(trans.wet.2017$FDR < 0.001),]
trans.salt.2017 = trans.salt.2017[which(trans.salt.2017$FDR < 0.001),]

summary(data.frame(table(cis.wet.2017$gene))$Freq)
summary(data.frame(table(cis.salt.2017$gene))$Freq)
summary(data.frame(table(trans.wet.2017$gene))$Freq)
summary(data.frame(table(trans.salt.2017$gene))$Freq)

write.table(cis.wet.2017, "cis-eQTLs.wet.fdr001.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(cis.salt.2017, "cis-eQTLs.salt.fdr001.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(trans.1mb.wet.2017, "trans-eQTLs.wet.fdr001.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(trans.1mb.salt.2017, "trans-eQTLs.salt.fdr001.txt", quote = F, row.names = F, col.names = T, sep="\t")

#saltDffull and wetDffull estimated by matrixeQTL = 116
dfFull = 116
cis.wet.2017$R2 <- (cis.wet.2017$statistic^2) / (cis.wet.2017$statistic^2 + dfFull)
cis.salt.2017$R2 <- (cis.salt.2017$statistic^2) / (cis.salt.2017$statistic^2 + dfFull)
p1 <- ggplot() +
  geom_density( aes(x = abs(R2), y = ..density..), data = cis.wet.2017, fill="#69b3a2" ) +
  geom_density( aes(x = abs(R2), y = -..density..), data = cis.salt.2017, fill= "#404080") +
  xlab("Distribution of variance explained")

trans.wet.2017$R2 <- (trans.wet.2017$statistic^2) / (trans.wet.2017$statistic^2 + dfFull)
trans.salt.2017$R2 <- (trans.salt.2017$statistic^2) / (trans.salt.2017$statistic^2 + dfFull)
p2 <- ggplot() +
  geom_density( aes(x = abs(R2), y = ..density..), data = trans.wet.2017, fill="#69b3a2" ) +
  geom_density( aes(x = abs(R2), y = -..density..), data = trans.salt.2017, fill= "#404080") +
  xlab("Distribution of variance explained")


genepos = read.table("../gene_locations.txt", header = T, sep="\t")
cis.salt.2017 = merge(genepos, cis.salt.2017, by="gene")
trans.salt.2017 = merge(genepos, trans.salt.2017, by="gene")
cis.wet.2017 = merge(genepos, cis.wet.2017, by="gene")
trans.wet.2017 = merge(genepos, trans.wet.2017, by="gene")

cis.salt.2017$dist <- ifelse(cis.salt.2017$Pos < cis.salt.2017$s1, cis.salt.2017$s1-cis.salt.2017$Pos, cis.salt.2017$Pos-cis.salt.2017$s2)   
cis.wet.2017$dist <- ifelse(cis.wet.2017$Pos < cis.wet.2017$s1, cis.wet.2017$s1-cis.wet.2017$Pos, cis.wet.2017$Pos-cis.wet.2017$s2)   

tmp = trans.salt.2017[which(trans.salt.2017$Chr.y == trans.salt.2017$Chr.x),]
tmp$dist = ifelse(tmp$Pos < tmp$s1, tmp$s1-tmp$Pos, tmp$Pos-tmp$s2)
tmp2 = trans.salt.2017[which(!(trans.salt.2017$Chr.y == trans.salt.2017$Chr.x)),]
tmp2$dist = "NA"
trans.salt.2017 = rbind(tmp,tmp2)

tmp = trans.wet.2017[which(trans.wet.2017$Chr.y == trans.wet.2017$Chr.x),]
tmp$dist = ifelse(tmp$Pos < tmp$s1, tmp$s1-tmp$Pos, tmp$Pos-tmp$s2)
tmp2 = trans.wet.2017[which(!(trans.wet.2017$Chr.y == trans.wet.2017$Chr.x)),]
tmp2$dist = "NA"
trans.wet.2017 = rbind(tmp,tmp2)

cis.salt.2017$dist = abs(cis.salt.2017$dist); cis.wet.2017$dist = abs(cis.wet.2017$dist)
trans.salt.2017$dist = as.numeric(trans.salt.2017$dist); trans.wet.2017$dist = as.numeric(trans.wet.2017$dist)

#Defining trans as over 1Mb or on different chr and estimating whether the same gene effects are more compensating vs reinforcing
#I did multiple analysis and realized that the results from 1mb away are similar to diff chr (indicating the loss of linkage). So will take this
trans.1mb.wet.2017 = trans.wet.2017[which(trans.wet.2017$dist>1000000 | is.na(trans.wet.2017$dist)),]
trans.1mb.salt.2017 = trans.salt.2017[which(trans.salt.2017$dist>1000000 | is.na(trans.salt.2017$dist)),]
write.table(trans.1mb.wet.2017, "trans-eQTLs.wet.1mb.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(trans.1mb.salt.2017, "trans-eQTLs.salt.1mb.txt", quote = F, row.names = F, col.names = T, sep="\t")


cis.salt.2017 = read.table("cis-eQTLs.salt.fdr001.txt", header = T, sep="\t")
cis.wet.2017 = read.table("cis-eQTLs.wet.fdr001.txt", header = T, sep="\t")
trans.1mb.salt.2017 = read.table("trans-eQTLs.salt.1mb.fdr001.txt", header = T, sep="\t")
trans.1mb.wet.2017 = read.table("trans-eQTLs.wet.1mb.fdr001.txt", header = T, sep="\t")


##Now I will identify the lead SNPs and take that only to analyse the data:
library(GenomicRanges)
library(dplyr)
z = trans.1mb.wet.2017      #change name of df here
x = unique(z$gene)   
y = c()
a = c()
for(i in x){
  df = z[which(z$gene == i),]  
  df  = df[order(df$Chr.y, df$Pos),]
  gr = GRanges(seqnames=df$Chr.y,
               IRanges(start=df$Pos,end=df$Pos),pval=df$pvalue, geneid=df$gene)
  names(gr) = df$snps
  FLANK = 100000
  #flank your snps by 100kb and we merge all these regions together
  REGIONS <- reduce(flank(gr,FLANK,both=TRUE))
  # each SNP can only be matched to one merged region, so we just find overlap between region and snp and assign the snp to the region
  gr$region = subjectHits(findOverlaps(gr,REGIONS))
  # order by pvalue
  gr = gr[order(gr$pval),]
  # keep only the top snp in each region
  tmp = gr[!duplicated(gr$region)]
  y=append(y, tmp@ranges@NAMES)
  b = rep(i, times=length(tmp@ranges@NAMES))
  a = append(a, b)
}
zz=data.frame(a,y)
names(zz) = c("gene", "snps")
tmp = merge(zz, z, by=c("gene", "snps"))
tmp = tmp[c(1,3:5,2,6:ncol(tmp))]
trans.1mb.wet.lead = tmp   #change name of df here

write.table(cis.wet.2017.lead, "cis-eQTLs.wet.leadSNPs.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(cis.salt.2017.lead, "cis-eQTLs.salt.leadSNPs.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(trans.1mb.wet.2017.lead, "trans-eQTLs.wet.1mb.leadSNPs.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(trans.1mb.salt.2017.lead, "trans-eQTLs.salt.1mb.leadSNPs.txt", quote = F, row.names = F, col.names = T, sep="\t")
