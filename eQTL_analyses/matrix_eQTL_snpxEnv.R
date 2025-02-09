library(MatrixEQTL)
setwd("ricesalinity/data/eQTL/ind_data/mat_eQTL_linearCross/")

#Prepping data for linearCross model which will have both env data in the same file
snps = read.table("../ind123.SNP.bial.pruned.txt", header = T, sep = "\t")
tmp = colnames(snps)[-1]
t1 = paste(tmp, "_w", sep = "")
t2 = paste(tmp, "_s", sep = "")
t = c(t1,t2)
tmp = snps[,c(1:124,2:124)]
colnames(tmp) = c(colnames(snps)[1], t)
snps = tmp
write.table(snps, "ind123_2.SNP.bial.pruned.txt", quote = F, row.names = F, col.names = T, sep = "\t")
rm(snps)

exp1 = read.table("../IND_wet_exp.eQTLin.txt", header = T, sep = "\t")
exp2 = read.table("../IND_salt_exp.eQTLin.txt", header = T, sep = "\t")
exp = merge(exp1, exp2, by="gene_id")
colnames(exp) = c(colnames(exp)[1], t)
write.table(exp, "IND_exp.eQTLin.txt", quote = F, row.names = F, col.names = T, sep = "\t")
rm(exp); rm(exp1); rm(exp2)

cov = read.table("../pop_str_cov.ind.txt", header = T, sep = "\t")
tmp = cov[,c(1:124,2:124)]
colnames(tmp) = c(colnames(cov)[1], t)
env = c("Env", rep("0", times = 123), rep("1", times = 123))
tmp = rbind(tmp, env)
cov = tmp
write.table(cov, "pop_str_cov_2.ind.txt", quote = F, row.names = F, col.names = T, sep = "\t")
rm(cov)


###getting data for running Matrix-eQTL
SNP_file_name = "ind123_2.SNP.bial.pruned.txt"
#snp file is in a numerical format (0/1/2) with ind as columns and snps as rows; first col is snp id followed by inds
expression_file_name = "IND_exp.eQTLin.txt"
#same as snp file but the num here are the expression counts and the first col is geneid
covariates_file_name = "pop_str_cov_2.ind.txt"
#covariates as rows

output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

pvOutputThreshold_cis = 1e-2;
pvOutputThreshold_tra = 1e-5;

errorCovariance = numeric();

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 4000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}


cisDist = 100000;
snpspos = read.table("snps.location.txt", header = TRUE, stringsAsFactors = FALSE);
genepos = read.table("gene_locations.txt", header = TRUE, stringsAsFactors = FALSE);
useModel = modelLINEAR_CROSS;

## Run the analysis
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:
## Plot the histogram of all p-values
plot(me)
##high inflation, meaning false-positive; will need to fix this

library(tidyr)
tmp = me$cis$eqtls
dim(tmp)  #25603
tmp = separate(data = me$cis$eqtls, col = snps, into = c("Chr", "Pos"), sep = "_", remove=FALSE)
write.table(tmp, "cis-eQTLs.txt", col.names = T, row.names = F, sep="\t", quote = F)
tmp = tmp[which(tmp$FDR <= 0.001),]
write.table(tmp, "cis-eQTLs.fdr001.txt", col.names = T, row.names = F, sep="\t", quote = F)
cis = tmp
dim(cis)  #92

tmp = me$trans$eqtls
dim(tmp)  #162260
tmp = separate(data = me$trans$eqtls, col = snps, into = c("Chr", "Pos"), sep = "_", remove=FALSE)
write.table(tmp, "trans-eQTLs.txt", col.names = T, row.names = F, sep="\t", quote = F)
tmp = tmp[which(tmp$FDR <= 0.001),]
write.table(tmp, "trans-eQTLs.fdr001.txt", col.names = T, row.names = F, sep="\t", quote = F)
trans = tmp
dim(trans)  #30323

cis = merge(genepos, cis, by="gene")
cis$Pos = as.integer(cis$Pos)
cis$dist <- ifelse(cis$Pos < cis$s1, cis$s1-cis$Pos, cis$Pos-cis$s2) 

trans = merge(genepos, trans, by="gene")
trans$Pos = as.integer(trans$Pos)
tmp = trans[which(trans$Chr.y == trans$Chr.x),]
tmp$dist = ifelse(tmp$Pos < tmp$s1, tmp$s1-tmp$Pos, tmp$Pos-tmp$s2)
tmp2 = trans[which(!(trans$Chr.y == trans$Chr.x)),]
tmp2$dist = "NA"
trans = rbind(tmp,tmp2)
trans.1mb = trans[which(trans$dist>1000000 | is.na(trans$dist)),]
dim(trans.1mb)  #29869

trans_old1 = read.table("../fdr_trans1Mb/trans-eQTLs.wet.1mb.fdr001.txt", header = T, sep="\t")
trans_old2 = read.table("../fdr_trans1Mb/trans-eQTLs.salt.1mb.fdr001.txt", header = T, sep="\t")
tmp = unique(union(setdiff(trans_old1[,c(1,5)], trans_old2[,c(1,5)]), setdiff(trans_old2[,c(1,5)], trans_old1[,c(1,5)])))
dim(tmp)  #111572
t1 = (merge(trans.1mb[,c(1,5)], tmp))  #882
write.table(t1, "transfdr001_1mb_envsp_intersectwithByEnvComparison.txt", quote = F, row.names = F, col.names = T, sep = "\t")
#total common gxe identified in the two -- 882


cis_old1 = read.table("../fdr_trans1Mb/cis-eQTLs.wet.fdr001.txt", header = T, sep="\t")
cis_old2 = read.table("../fdr_trans1Mb/cis-eQTLs.salt.fdr001.txt", header = T, sep="\t")
tmp = unique(union(setdiff(cis_old1[,c(4,1)], cis_old2[,c(4,1)]), setdiff(cis_old2[,c(4,1)], cis_old1[,c(4,1)])))
dim(tmp)  #30046
t2= (merge(cis[,c(1,5)], tmp))   #53
write.table(t2, "cisfdr001_envsp_intersectwithByEnvComparison.txt", quote = F, row.names = F, col.names = T, sep = "\t")

x = c(53,882)
n = c(30046, 111572)
prop.test(x, n, alternative = "two.sided")


