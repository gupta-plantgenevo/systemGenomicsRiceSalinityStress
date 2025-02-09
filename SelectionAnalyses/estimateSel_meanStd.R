##This script will estimate the linear (S) and quadratic (C) selection coeff through linear regression (on Fecundity) taking Block as the random effect

library(stats)
library(lmerTest)
library(dplyr)

selestimate_fecun <- function(fecun){
  S = c()
  C = c()
  Pval_S = c()
  S_OK_or_singular = c()
  Pval_C = c()
  tval_S = c()
  tval_C = c()
  C_OK_or_singular = c()
  ##This loop will fit a linear and quadratic model and then save the S and C values (regression coeff) and the respective P values in separate vectors
  ##Takes a while to run
  for(i in 3:ncol(fecun)){
    tr = fecun[,i]
    linear_model = lmer(Fecundity ~ tr + (1|Block), data = fecun)
    S_OK_or_singular = append(S_OK_or_singular, ifelse(grepl("Singular",tail(summary(linear_model), n=1)),"Singular","OK"))
    S = append(S, summary(linear_model)$coefficients[2,1])
    tval_S = append(tval_S, summary(linear_model)$coefficients[2,4])
    Pval_S = append(Pval_S, summary(linear_model)$coefficients[2,5])
    tr2 = tr^2
    quad_model = lmer(Fecundity ~ tr + tr2 + (1|Block), data = fecun)
    C_OK_or_singular = append(C_OK_or_singular, ifelse(grepl("Singular",tail(summary(quad_model), n=1)),"Singular","OK"))
    #summary(quad_model)
    C = append(C, 2*(summary(quad_model)$coefficients[3,1]))
    tval_C = append(tval_C, (summary(quad_model)$coefficients[3,4]))
    Pval_C = append(Pval_C, summary(quad_model)$coefficients[3,5])
  }
  output = data.frame(S,tval_S,Pval_S,S_OK_or_singular,C,tval_C,Pval_C,C_OK_or_singular)
  return(output)
}


###Doing for Salt Env
ind_fecundity = read.table("ricesalinity/data/SelectionAnalysis/InputFiles/IND_salt.w0fecun.nonStd.Selinput.txt", header=TRUE, sep="\t", check.names = FALSE)
TrID = colnames(ind_fecundity)
TrID = TrID[-c(1,2)]
ind_fecundity$Block = as.factor(ind_fecundity$Block)
#tmp = sel_coeff
sel_coeff = selestimate_fecun(ind_fecundity)
sel_coeff$TrID = TrID
sel_coeff = sel_coeff[,c(9,1:8)]
#str(sel_coeff)
write.table(sel_coeff, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeffraw_IND_salt.w0fecun.nonStd.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

sel_coeff = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeffraw_IND_wet.w0fecun.nonStd.txt", header = T, sep = "\t")
trmean = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/IndWet_Transcript_mean_sd.txt", header = T, sep="\t")
str(sel_coeff); str(trmean)
tmp = merge(sel_coeff, trmean, by="TrID")
str(tmp)
tmp$S.meanStd = tmp$S*tmp$mean; tmp$C.meanStd = tmp$C*(tmp$mean*tmp$mean)
tmp$S.varStd = tmp$S*tmp$sd; tmp$C.varStd = tmp$C*(tmp$sd*tmp$sd)
tmp = tmp[,c(1:5,12,14,6:9,13,15)]
write.table(tmp, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_wet.w0fecun.nonStd.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")



##Normal Env
ind_fecundity = read.table("ricesalinity/data/SelectionAnalysis/InputFiles/IND_normal.w0fecun.nonStd.Selinput.txt", header=TRUE, sep="\t", check.names = FALSE)
TrID = colnames(ind_fecundity)
TrID = TrID[-c(1,2)]
ind_fecundity$Block = as.factor(ind_fecundity$Block)
#tmp = sel_coeff
sel_coeff = selestimate_fecun(ind_fecundity)
sel_coeff$TrID = TrID
sel_coeff = sel_coeff[,c(9,1:8)]
#str(sel_coeff)
write.table(sel_coeff, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeffraw_IND_normal.w0fecun.nonStd.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

sel_coeff = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeffraw_IND_normal.w0fecun.nonStd.txt", header = T, sep = "\t")
trmean = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/Indnormal_Transcript_mean_sd.txt", header = T, sep="\t")
str(sel_coeff); str(trmean)
tmp = merge(sel_coeff, trmean, by="TrID")
str(tmp)
tmp$S.meanStd = tmp$S*tmp$mean; tmp$C.meanStd = tmp$C*(tmp$mean*tmp$mean)
tmp$S.varStd = tmp$S*tmp$sd; tmp$C.varStd = tmp$C*(tmp$sd*tmp$sd)
tmp = tmp[,c(1:5,12,14,6:9,13,15)]
write.table(tmp, "ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_normal.w0fecun.nonStd.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")


