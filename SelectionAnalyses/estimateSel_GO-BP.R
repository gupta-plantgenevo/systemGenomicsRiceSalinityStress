library('stringr')
library(dplyr)
library(tidyr)
library(data.table)

quantile_cl <- function(y, q=0.5, conf.level = 0.95, na.rm=TRUE) {
  alpha <- 1 - conf.level
  if (na.rm) y <- y[!is.na(y)]
  n <- length(y)
  l <- qbinom(alpha/2, size=n, prob = q)
  u <- 1 + n - l
  ys <- sort.int(c(-Inf, y, Inf), partial = c(1 + l, 1 + u))
  data.frame(
    y = quantile(y, probs = q, na.rm=na.rm, type = 8),
    ymin = ys[1 + l],
    ymax = ys[1 + u]
  )
}

#Get the GO terms for third lvel
library(GO.db)
getAllBPChildren <- function(goids)
{
  ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names=FALSE))
  ans <- ans[!is.na(ans)]
}

level1_BP_terms <- getAllBPChildren("GO:0008150")     # 23 terms
level2_BP_terms <- getAllBPChildren(level1_BP_terms)  # 256 terms
level3_BP_terms <- getAllBPChildren(level2_BP_terms)  # 3059 terms
level4_BP_terms <- getAllBPChildren(level3_BP_terms)  # 9135 terms
#level5_BP_terms <- getAllBPChildren(level4_BP_terms)  # 15023 terms
write.table(level4_BP_terms, file = "ricesalinity/data/module_pathway_diff/GO_BP/GOid_BP.level4", quote = F, col.names = F, row.names = F)
write.table(level3_BP_terms, file = "ricesalinity/data/module_pathway_diff/GO_BP/GOid_BP.level3", quote = F, col.names = F, row.names = F)

z = fread("module_pathway_diff/GO_BP/go.osa.csv", sep="\t", header = T)
#x = read.table("module_pathway_diff/GO_BP/GOid_BP.all", header = F)
x = read.table("module_pathway_diff/GO_BP/GOid_BP.level4", header = F)
y = z[which(z$go %in% x$V1),]

y$TrID = y$gene_id
y$TrID <- str_replace_all(y$TrID, 'Os', 'OS')
y$TrID <- str_replace_all(y$TrID, 'g', 'T')
a = y[,c(10,1,3:8)]
write.table(a, "module_pathway_diff/GO_BP/go.osa.BP.txt", row.names = F, col.names = T, quote = F, sep = "\t")

tmp = y[,c(10,1,3)]
geno1 = fread("ricesalinity/data/geneExp/INDICA_Wet_GeneExp.transposed.txt", header = T, sep="\t", check.names = F)
z = as.data.frame(colnames(geno1)[-1])
#rm(geno1)
colnames(z)[1] = "TrID"
z = z %>% separate(TrID, into=c("Tr","no"), sep="-", remove = F)

x = unique(merge(tmp, z, by.x = "TrID", by.y = "Tr"))
length(unique(x$gene_id))    #9972 genes and transcripts
length(unique(x$go))
x = x[,c(4,2,3)]
names(x)[1] = "TrID"

sel_coeff_wet = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_normal.w0fecun.nonStd.txt", header=TRUE, sep = "\t")
sel_coeff_salt = read.table("ricesalinity/data/SelectionAnalysis/SelCoeff_output/SelCoeff_IND_salt.w0fecun.nonStd.txt", header=TRUE, sep = "\t")

z = x
tmp = merge(z, sel_coeff_wet[,c(1,2)], by="TrID", all.x = T)
names(tmp)[4] = "S.wet"
tmp$S.wet = abs(tmp$S.wet)
tmp2 = merge(tmp, sel_coeff_salt[,c(1,2)], by="TrID", all.x = T)
names(tmp2)[5] = "S.salt"
tmp2$S.salt = abs(tmp2$S.salt)
#write.table(tmp2, "module_pathway_diff/GO_BP/selTr_BP.GOid.level3.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(tmp2, "module_pathway_diff/GO_BP/selRaw.Tr_BP.GOid.level4.txt", row.names = F, col.names = T, sep = "\t", quote = F)
trGO.S = tmp2

y = as.data.frame(table(x$go))    #5584 Go terms
tmp = y[which(y$Freq > 20),]    #1490 GO terms have frequency of over 20
z = x[which(x$go %in% tmp$Var1),]
length(unique(z$TrID))      #Still 10742 transcripts

#z = x
tmp = merge(z, sel_coeff_wet[,c(1,2)], by="TrID", all.x = T)
names(tmp)[4] = "S.wet"
tmp$S.wet = abs(tmp$S.wet)
tmp2 = merge(tmp, sel_coeff_salt[,c(1,2)], by="TrID", all.x = T)
names(tmp2)[5] = "S.salt"
tmp2$S.salt = abs(tmp2$S.salt)
#write.table(tmp2, "module_pathway_diff/GO_BP/selTr_GOover20_withS.level3.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(tmp2, "module_pathway_diff/GO_BP/selRaw.Tr_GOover20_withS.level3.txt", row.names = F, col.names = T, sep = "\t", quote = F)

trGO.S = tmp2

ttestpval = c()
gon = c()
for (i in unique(tmp2$go)){
  o = tmp2[which(tmp2$go == i),]
  ttestpval = append(ttestpval, t.test(o$S.wet, o$S.salt)$p.value)
  gon = append(gon, i)
}
p = data.frame(gon,ttestpval)

trGO.S = read.table("module_pathway_diff/GO_BP/selTr_GOover20_withS.txt", header = T, sep="\t")
#Getting transcriptome-wide median (considering only the 10862 transcripts used here)
tmp2 = unique(trGO.S[,c(1,4)])
wetSmed = median(tmp2$S.wet, na.rm = T)
tmp2 = unique(trGO.S[,c(1,5)])
saltSmed = median(tmp2$S.salt, na.rm = T)

library(Rmisc)
#Group by GO, and getting median S by GO
tmp = as.data.frame(trGO.S %>% group_by(go) %>% 
                      summarise_if(is.numeric, list(~ sum(!is.na(.)), ~ median(., na.rm = TRUE), ~ mean(., na.rm = TRUE), ~ CI(na.omit(.), ci =0.95)[3])))
#the above gets the lower 95%CI around the mean
names(tmp)[8] = "L.ci.wet"
names(tmp)[9] = "L.ci.salt"
tmp = tmp[,c(1,2,4,6,8,3,5,7,9)]
GOmed_ci = tmp

temp = as.data.frame(trGO.S %>% group_by(go) %>% 
                      summarise_if(is.numeric, list(~ sum(!is.na(.)), ~ median(., na.rm = TRUE), ~ mean(., na.rm = TRUE), ~ quantile_cl(., conf.level =0.99, na.rm = T)$ymin)))
#the above gets the lower 95%CI around the median
names(temp)[8] = "L.ci.wet"
names(temp)[9] = "L.ci.salt"
temp = temp[,c(1,2,4,6,8,3,5,7,9)]
GOmed_ci = temp


wetSel = GOmed_ci[which(GOmed_ci$L.ci.wet>wetSmed),]
dim(wetSel)
saltSel = GOmed_ci[which(GOmed_ci$L.ci.salt>saltSmed),]  
dim(saltSel)

length(intersect(wetSel$go, saltSel$go)) 
length(setdiff(wetSel$go, saltSel$go))  
length(setdiff(saltSel$go, wetSel$go))  

a = GOmed_ci %>% mutate(wetCiPass = case_when(
    GOmed_ci$go %in% wetSel$go ~ "Yes",
    TRUE ~ "No"))

GOmed_ci = a %>% mutate(saltCiPass = case_when(
    GOmed_ci$go %in% saltSel$go ~ "Yes",
    TRUE ~ "No"))

a = GOmed_ci %>% mutate(WetSaltCiPass = case_when(
    GOmed_ci$go %in% (intersect(wetSel$go, saltSel$go)) ~ "Yes",
    TRUE ~ "No"))

GOmed_ci = a %>% mutate(WetonlyCiPass = case_when(
    GOmed_ci$go %in% (setdiff(wetSel$go, saltSel$go)) ~ "Yes",
    TRUE ~ "No"))

a = GOmed_ci %>% mutate(SaltonlyCiPass = case_when(
    GOmed_ci$go %in% (setdiff(saltSel$go, wetSel$go)) ~ "Yes",
    TRUE ~ "No"))

GOmed_ci = a
GOmed_ci= GOmed_ci[,c(1:5, 10,6:9,11:14)]

write.table(GOmed_ci, "module_pathway_diff/GO_BP/revisedSelEstimates/Sraw/golev4_median_CI0.99_sel.txt", row.names = F, col.names = T, sep="\t", quote = F)


#Getting rank-sum test stats and p-value 
#NOt sure whether I should do paired test or not -- I think we should (not mentioned in Neils' paper)
wilcox.test(tmp$S.wet_median, tmp$S.salt_median, alternative = "less", paired = T)
#to get the z-value, qnorm(pvalue)
qnorm(6.51e-05) 
#z = -3.82608; p-value = 6.51e-05

