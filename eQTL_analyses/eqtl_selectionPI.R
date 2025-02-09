###This script will test for the difference of pi between regions near cis and trans eqtl and the background genome-wide patterns.
library(GenomicRanges)
library(dplyr)

setwd("ricesalinity/data/eQTL/ind_data/")

pi_all = read.table("indica123.SNP.final.bial.pruned.sites.pi", header = T, sep="\t")
pi_all$CHROM = sub("chr0", "", pi_all$CHROM)
pi_all$CHROM = sub("chr", "", pi_all$CHROM)

pibase = GRanges(seqnames=pi_all$CHROM,
                 IRanges(start=pi_all$POS,end=pi_all$POS),pi=pi_all$PI)

#Reading the genome and dividing it into 100kb bins for the background pi (average over 100kb) -- doing only for normal conditions 
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

#Now estimating the mean of cis - normal pi in +/-50kb regions (100kb total)
df = read.table("fdr_trans1Mb/cis-eQTLs.wet.leadSNPs.txt", header = T, sep="\t")

gr = GRanges(seqnames=df$Chr.y,
             IRanges(start=df$Pos-50000,end=df$Pos+50000))
names(gr) = df$snps
tmp = as.data.frame(findOverlaps(gr,pibase))
cispi_wet = c()
for (i in unique(tmp$queryHits)){
  t = tmp[which(tmp$queryHits == i),]
  cispi_wet = append(cispi_wet, mean(pi_all[(t$subjectHits),]$PI))
}


remainnormal = pi_all[-(unique((tmp$subjectHits))),]
#The above has all snps that are not flanked by the cis lead snps
gr = GRanges(seqnames = remainnormal$CHROM, 
             IRanges(start=remainnormal$POS,end=remainnormal$POS),pi=remainnormal$PI)

tmp2= as.data.frame(findOverlaps(bins,pibase))
backgroundpi_normal = c()
for (i in unique(tmp2$queryHits)){
  t = tmp2[which(tmp2$queryHits == i),]
  backgroundpi_normal = append(backgroundpi_normal, mean(remainnormal[(t$subjectHits),]$PI))
}

#Estimating the diff between cis and background within 100kb regions
wilcox.test(cispi_wet, backgroundpi_normal, alternative = "greater")
t.test(cispi_wet, backgroundpi_normal, alternative = "greater")
t.test(cispi_wet, backgroundpi_normal)
mean(cispi_wet); mean(backgroundpi_normal,na.rm = T)

#Estimating the diff between cis and background with only SNPs (and noy 100kb region)
tmp = merge(pi_all, df[,c(6,7)], by.x=c("CHROM", "POS"), by.y=c("Chr.y", "Pos"))
tmp2 = anti_join(pi_all, tmp)
wilcox.test(tmp$PI, tmp2$PI, alternative = "greater")
mean(tmp$PI); mean(tmp2$PI, na.rm=T)


#Now estimating the mean of cis - salt pi in +/-50kb regions (100kb total)
df = read.table("fdr_trans1Mb/cis-eQTLs.salt.leadSNPs.txt", header = T, sep="\t")

gr = GRanges(seqnames=df$Chr.y,
             IRanges(start=df$Pos-50000,end=df$Pos+50000))
names(gr) = df$snps
tmp = as.data.frame(findOverlaps(gr,pibase))
cispi_salt = c()
for (i in unique(tmp$queryHits)){
  t = tmp[which(tmp$queryHits == i),]
  cispi_salt = append(cispi_salt, mean(pi_all[(t$subjectHits),]$PI))
}


remainsalt = pi_all[-(unique((tmp$subjectHits))),]
#The above has all snps that are not flanked by the cis lead snps
gr = GRanges(seqnames = remainsalt$CHROM, 
               IRanges(start=remainsalt$POS,end=remainsalt$POS),pi=remainsalt$PI)

tmp2= as.data.frame(findOverlaps(bins,pibase))
backgroundpi_salt = c()
for (i in unique(tmp2$queryHits)){
  t = tmp2[which(tmp2$queryHits == i),]
  backgroundpi_salt = append(backgroundpi_salt, mean(remainsalt[(t$subjectHits),]$PI))
}

wilcox.test(cispi_salt, backgroundpi_salt, alternative = "greater")
t.test(cispi_salt, backgroundpi_salt)
mean(cispi_salt); mean(backgroundpi_salt,na.rm = T)
mean(cispi_wet); mean(backgroundpi_normal,na.rm = T)

tmp = merge(pi_all, df[,c(6,7)], by.x=c("CHROM", "POS"), by.y=c("Chr.y", "Pos"))
tmp2 = anti_join(pi_all, tmp)
wilcox.test(tmp$PI, tmp2$PI, alternative = "greater")
mean(tmp$PI); mean(tmp2$PI, na.rm=T)


#######TRANS
#Now estimating the mean of trans - normal pi in +/-50kb regions (100kb total)
df = read.table("fdr_trans1Mb/trans-eQTLs.wet.1mb.leadSNPs.txt", header = T, sep="\t")

gr = GRanges(seqnames=df$Chr.y,
             IRanges(start=df$Pos-50000,end=df$Pos+50000))
names(gr) = df$snps
tmp = as.data.frame(findOverlaps(gr,pibase))
transpi_wet = c()
for (i in unique(tmp$queryHits)){
  t = tmp[which(tmp$queryHits == i),]
  transpi_wet = append(transpi_wet, mean(pi_all[(t$subjectHits),]$PI))
}


remainnormal = pi_all[-(unique((tmp$subjectHits))),]
#The above has all snps that are not flanked by the cis lead snps
gr = GRanges(seqnames = remainnormal$CHROM, 
             IRanges(start=remainnormal$POS,end=remainnormal$POS),pi=remainnormal$PI)

tmp2= as.data.frame(findOverlaps(bins,pibase))
backgroundpi_normal = c()
for (i in unique(tmp2$queryHits)){
  t = tmp2[which(tmp2$queryHits == i),]
  backgroundpi_normal = append(backgroundpi_normal, mean(remainnormal[(t$subjectHits),]$PI))
}

wilcox.test(transpi_wet, backgroundpi_normal, alternative = "less")
mean(transpi_wet); mean(backgroundpi_normal,na.rm = T)

tmp = merge(pi_all, df[,c(6,7)], by.x=c("CHROM", "POS"), by.y=c("Chr.y", "Pos"))
tmp2 = anti_join(pi_all, tmp)
wilcox.test(tmp$PI, tmp2$PI, alternative = "less")
mean(tmp$PI); mean(tmp2$PI, na.rm=T)


#Now estimating the mean of cis - salt pi in +/-50kb regions (100kb total)
df = read.table("fdr_trans1Mb/trans-eQTLs.salt.1mb.leadSNPs.txt", header = T, sep="\t")

gr = GRanges(seqnames=df$Chr.y,
             IRanges(start=df$Pos-50000,end=df$Pos+50000))
names(gr) = df$snps
tmp = as.data.frame(findOverlaps(gr,pibase))
transpi_salt = c()
for (i in unique(tmp$queryHits)){
  t = tmp[which(tmp$queryHits == i),]
  transpi_salt = append(transpi_salt, mean(pi_all[(t$subjectHits),]$PI))
}

remainsalt = pi_all[-(unique((tmp$subjectHits))),]
#The above has all snps that are not flanked by the cis lead snps
gr = GRanges(seqnames = remainsalt$CHROM, 
             IRanges(start=remainsalt$POS,end=remainsalt$POS),pi=remainsalt$PI)

tmp2= as.data.frame(findOverlaps(bins,pibase))
backgroundpi_salt = c()
for (i in unique(tmp2$queryHits)){
  t = tmp2[which(tmp2$queryHits == i),]
  backgroundpi_salt = append(backgroundpi_salt, mean(remainsalt[(t$subjectHits),]$PI))
}

wilcox.test(transpi_salt, backgroundpi_salt, alternative = "less")
mean(transpi_salt); mean(backgroundpi_salt,na.rm = T)
mean(transpi_wet); mean(backgroundpi_normal,na.rm = T)

wilcox.test(transpi_salt, transpi_wet, alternative = "less")

tmp = merge(pi_all, df[,c(6,7)], by.x=c("CHROM", "POS"), by.y=c("Chr.y", "Pos"))
tmp2 = anti_join(pi_all, tmp)
wilcox.test(tmp$PI, tmp2$PI, alternative = "less")
mean(tmp$PI); mean(tmp2$PI, na.rm=T)




