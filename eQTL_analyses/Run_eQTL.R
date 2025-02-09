##This code will run MatrixeQTL. Note: re-run the same code again with a change the expression file name to salt exp data; also change the output file names in write.table

library(MatrixEQTL)

setwd("ricesalinity/data/eQTL/ind_data/")
useModel = modelLINEAR;
SNP_file_name = "ind123.SNP.bial.pruned.txt"
#snp file is in a numerical format (0/1/2) with ind as columns and snps as rows; first col is snp id followed by inds
expression_file_name = "IND_wet_exp.eQTLin.txt"
#same as snp file but the num here are the expression counts and the first col is geneid
covariates_file_name = "pop_str_cov.ind.txt"
#covariates as rows

output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

pvOutputThreshold_cis = 1e-7;
pvOutputThreshold_tra = 1e-7;

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

#Setting cisDist to 100kb; in the downstream file I will remove any eQTL from 100kb to 1Mb; Less than 100kb will be defined cis and more than 1Mb on the same chr will be defined trans
cisDist = 100000;
snpspos = read.table("snps.location.txt", header = TRUE, stringsAsFactors = FALSE);
genepos = read.table("gene_locations.txt", header = TRUE, stringsAsFactors = FALSE);

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
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

## Plot the histogram of all p-values
plot(me)
##high inflation, meaning false-positive; will need to fix this

library(tidyr)
tmp = me$cis$eqtls
tmp = separate(data = me$cis$eqtls, col = snps, into = c("Chr", "Pos"), sep = "_", remove=FALSE)
write.table(tmp, "fdr_trans1Mb/cis-eQTLs.wet.txt", col.names = T, row.names = F, sep="\t", quote = F)

tmp = me$trans$eqtls
tmp = separate(data = me$trans$eqtls, col = snps, into = c("Chr", "Pos"), sep = "_", remove=FALSE)
write.table(tmp, "fdr_trans1Mb/trans-eQTLs.wet.txt", col.names = T, row.names = F, sep="\t", quote = F)

